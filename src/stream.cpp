/*#########################################################################

 Author: Thang V. Pham {t.pham@amsterdamumc.nl}

 All rights reserved.

 Streaming reader: builds the same tables as iq_filter (fast_read) from a
 sequence of data-frame batches, e.g. record batches of a parquet file read
 by the arrow package. Only one batch is held in memory at a time.

    state <- .Call("iq_stream_begin", options)
    .Call("iq_stream_add", state, batch)   # repeated, batch = data.frame / named list
    res <- .Call("iq_stream_end", state, verbose)   # same structure as fast_read(), or NULL

 Filters and missing values follow the text-file reader (iq_filter):
  - filter_double_less / greater keep x < value / x > value, and keep rows whose
    filter value is missing (NA, NaN, empty, or text that is not a number)
  - a quantitative value is missing if it is NA, NaN, empty, or equal to na_string

 Citation:
 
 Pham TV, Tran CTM, Henneman AA, Pham LHC, Le DG, Can AH, Bui PHL, Piersma SR, Jimenez CR. 
 Boosting the speed and accuracy of protein quantification algorithms in mass spectrometry-based proteomics.
 Journal of Proteome Research 2026 Feb 6; 25(2):1198-1203.

 Software version: 2.2.0

#########################################################################*/

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#define R_NO_REMAP

#include <R.h>
#include <Rinternals.h>

#include "utils.h"
#include "iq_writer.h"

using std::string;
using std::vector;
using std::runtime_error;

extern "C" SEXP iq_stream_begin(SEXP);
extern "C" SEXP iq_stream_add(SEXP, SEXP);
extern "C" SEXP iq_stream_end(SEXP, SEXP);
extern "C" SEXP iq_stream_write(SEXP, SEXP);

namespace {

// read-only view of one column of a batch
class column_view {
    SEXP x;
    SEXP levels;     // factor levels, or R_NilValue
    bool is_int64;   // bit64::integer64 stored in a double vector
    char buf[64];

public:
    column_view() : x(R_NilValue), levels(R_NilValue), is_int64(false) {
        buf[0] = 0;
    }

    explicit column_view(SEXP col) : x(col), levels(R_NilValue), is_int64(false) {
        buf[0] = 0;
        switch (TYPEOF(col)) {
            case STRSXP:
                break;
            case INTSXP:
                if (Rf_isFactor(col)) {
                    levels = Rf_getAttrib(col, R_LevelsSymbol);
                }
                break;
            case REALSXP:
                is_int64 = Rf_inherits(col, "integer64");
                break;
            case LGLSXP:
                break;
            default:
                throw runtime_error("unsupported column type");
        }
    }

    bool is_text() const {
        return TYPEOF(x) == STRSXP || levels != R_NilValue;
    }

    // text of a cell; NA becomes "" (an empty cell in a text file)
    const char* text(R_xlen_t i) {
        switch (TYPEOF(x)) {
            case STRSXP: {
                SEXP s = STRING_ELT(x, i);
                return s == NA_STRING ? "" : CHAR(s);
            }
            case INTSXP: {
                int v = INTEGER(x)[i];
                if (v == NA_INTEGER) {
                    return "";
                }
                if (levels != R_NilValue) {
                    return CHAR(STRING_ELT(levels, v - 1));
                }
                snprintf(buf, sizeof(buf), "%d", v);
                return buf;
            }
            case LGLSXP: {
                int v = LOGICAL(x)[i];
                return v == NA_LOGICAL ? "" : (v ? "TRUE" : "FALSE");
            }
            case REALSXP: {
                if (is_int64) {
                    int64_t v;
                    memcpy(&v, REAL(x) + i, sizeof(v));
                    if (v == INT64_MIN) {  // bit64 NA
                        return "";
                    }
                    snprintf(buf, sizeof(buf), "%lld", (long long)v);
                    return buf;
                }
                double v = REAL(x)[i];
                if (R_IsNA(v)) {
                    return "";
                }
                if (std::isnan(v)) {
                    return "NaN";
                }
                if (v == std::floor(v) && std::fabs(v) < 1e15) {
                    snprintf(buf, sizeof(buf), "%.0f", v);
                } else {
                    snprintf(buf, sizeof(buf), "%.15g", v);
                }
                return buf;
            }
        }
        return "";
    }

    // numeric value of a cell: NA -> NA_REAL, NaN -> NaN
    double number(R_xlen_t i) {
        switch (TYPEOF(x)) {
            case REALSXP: {
                if (is_int64) {
                    int64_t v;
                    memcpy(&v, REAL(x) + i, sizeof(v));
                    return v == INT64_MIN ? NA_REAL : (double)v;
                }
                return REAL(x)[i];
            }
            case INTSXP:
                if (levels == R_NilValue) {
                    int v = INTEGER(x)[i];
                    return v == NA_INTEGER ? NA_REAL : (double)v;
                }
                break;
            case LGLSXP: {
                int v = LOGICAL(x)[i];
                return v == NA_LOGICAL ? NA_REAL : (double)v;
            }
        }
        // text: the same rule as the text-file reader (utils.h)
        double v;
        return filter_number(text(i), &v) ? v : NA_REAL;
    }
};

struct stream_state {
    // options
    string col_sample;
    string col_primary_id;
    vector<string> col_secondary_ids;
    vector<string> col_quant;           // one column, or one column per fragment
    vector<string> col_quant_id;        // empty, or the same length as col_quant
    vector<string> col_annotations;

    vector<std::pair<string, string>> filter_string_equal;
    vector<std::pair<string, string>> filter_string_not_equal;
    vector<std::pair<string, double>> filter_double_less;
    vector<std::pair<string, double>> filter_double_greater;

    bool has_sep = false;
    char sep = ';';
    bool has_na_string = false;
    string na_string;
    bool na_string_is_number = false;
    double na_number = 0.0;

    // results
    vector<int> protein_index;
    vector<int> sample_index;
    vector<int> ion_index;
    vector<double> quant;

    vector<vector<string>> annotation;  // primary id, then annotation columns
    vector<string> samples;
    vector<string> ions;                // display name, parts joined with '_'
    vector<string> ion_prefixes;        // the secondary_id part of the display name, without the suffix

    std::unordered_map<string, int> map_proteins;
    std::unordered_map<string, int> map_samples;
    std::unordered_map<string, int> map_ions;   // key parts joined with '\001'

    size_t n_rows = 0;
    size_t n_rows_kept = 0;
    double n_values_written = 0;        // values passed on to an iq format writer
};

void finalize_state(SEXP ptr) {
    stream_state* st = static_cast<stream_state*>(R_ExternalPtrAddr(ptr));
    if (st) {
        delete st;
        R_ClearExternalPtr(ptr);
    }
}

stream_state* get_state(SEXP ptr) {
    if (TYPEOF(ptr) != EXTPTRSXP) {
        throw runtime_error("invalid reader state");
    }
    stream_state* st = static_cast<stream_state*>(R_ExternalPtrAddr(ptr));
    if (!st) {
        throw runtime_error("the reader has already been finished");
    }
    return st;
}

SEXP list_element(SEXP list, const char* name) {
    SEXP names = Rf_getAttrib(list, R_NamesSymbol);
    if (names == R_NilValue) {
        return R_NilValue;
    }
    for (R_xlen_t i = 0; i < Rf_xlength(list); i++) {
        if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
            return VECTOR_ELT(list, i);
        }
    }
    return R_NilValue;
}

vector<string> get_strings(SEXP options, const char* name, bool required) {
    SEXP x = list_element(options, name);
    vector<string> out;
    if (x == R_NilValue) {
        if (required) {
            throw runtime_error(string("missing option: ") + name);
        }
        return out;
    }
    if (TYPEOF(x) != STRSXP) {
        throw runtime_error(string("option must be character: ") + name);
    }
    for (R_xlen_t i = 0; i < Rf_xlength(x); i++) {
        if (STRING_ELT(x, i) == NA_STRING) {
            throw runtime_error(string("NA in option: ") + name);
        }
        out.emplace_back(CHAR(STRING_ELT(x, i)));
    }
    if (required && out.empty()) {
        throw runtime_error(string("empty option: ") + name);
    }
    return out;
}

// named character vector -> (column, value) pairs
vector<std::pair<string, string>> get_named_strings(SEXP options, const char* name) {
    SEXP x = list_element(options, name);
    vector<std::pair<string, string>> out;
    if (x == R_NilValue) {
        return out;
    }
    SEXP nm = Rf_getAttrib(x, R_NamesSymbol);
    if (TYPEOF(x) != STRSXP || nm == R_NilValue) {
        throw runtime_error(string("option must be a named character vector: ") + name);
    }
    for (R_xlen_t i = 0; i < Rf_xlength(x); i++) {
        if (STRING_ELT(x, i) == NA_STRING || STRING_ELT(nm, i) == NA_STRING || *CHAR(STRING_ELT(nm, i)) == 0) {
            throw runtime_error(string("missing name or value in option: ") + name);
        }
        out.emplace_back(CHAR(STRING_ELT(nm, i)), CHAR(STRING_ELT(x, i)));
    }
    return out;
}

// named numeric vector -> (column, value) pairs
vector<std::pair<string, double>> get_named_numbers(SEXP options, const char* name) {
    SEXP x = list_element(options, name);
    vector<std::pair<string, double>> out;
    if (x == R_NilValue) {
        return out;
    }
    SEXP nm = Rf_getAttrib(x, R_NamesSymbol);
    if (TYPEOF(x) != REALSXP || nm == R_NilValue) {
        throw runtime_error(string("option must be a named numeric vector: ") + name);
    }
    for (R_xlen_t i = 0; i < Rf_xlength(x); i++) {
        if (ISNAN(REAL(x)[i]) || STRING_ELT(nm, i) == NA_STRING || *CHAR(STRING_ELT(nm, i)) == 0) {
            throw runtime_error(string("missing name or value in option: ") + name);
        }
        out.emplace_back(CHAR(STRING_ELT(nm, i)), REAL(x)[i]);
    }
    return out;
}

stream_state* create_state(SEXP options) {
    if (TYPEOF(options) != VECSXP) {
        throw runtime_error("options must be a list");
    }

    std::unique_ptr<stream_state> st(new stream_state);

    st->col_sample = get_strings(options, "sample_id", true)[0];
    st->col_primary_id = get_strings(options, "primary_id", true)[0];
    st->col_secondary_ids = get_strings(options, "secondary_id", true);
    st->col_quant = get_strings(options, "intensity_col", true);
    st->col_quant_id = get_strings(options, "intensity_col_id", false);
    st->col_annotations = get_strings(options, "annotation_col", false);

    st->filter_string_equal = get_named_strings(options, "filter_string_equal");
    st->filter_string_not_equal = get_named_strings(options, "filter_string_not_equal");
    st->filter_double_less = get_named_numbers(options, "filter_double_less");
    st->filter_double_greater = get_named_numbers(options, "filter_double_greater");

    vector<string> sep = get_strings(options, "intensity_col_sep", false);
    if (!sep.empty()) {
        if (sep[0].size() != 1) {
            throw runtime_error("intensity_col_sep must be a single character");
        }
        st->has_sep = true;
        st->sep = sep[0][0];
    }

    vector<string> na = get_strings(options, "na_string", false);
    if (!na.empty()) {
        st->has_na_string = true;
        st->na_string = na[0];
        char* end = NULL;
        double v = strtod(na[0].c_str(), &end);
        if (!na[0].empty() && end && *end == 0) {
            st->na_string_is_number = true;
            st->na_number = v;
        }
    }

    if (st->col_quant.size() > 1 && st->has_sep) {
        throw runtime_error("intensity_col_sep can only be used with a single intensity column");
    }
    if (!st->col_quant_id.empty() && st->col_quant_id.size() != st->col_quant.size()) {
        throw runtime_error("intensity_col_id must have the same length as intensity_col");
    }
    if (st->col_quant.size() == 1 && !st->col_quant_id.empty() && !st->has_sep) {
        throw runtime_error("intensity_col_id with a single intensity column needs intensity_col_sep");
    }

    return st.release();
}

int find_column(SEXP names, const string& name) {
    for (R_xlen_t i = 0; i < Rf_xlength(names); i++) {
        if (STRING_ELT(names, i) != NA_STRING && name == CHAR(STRING_ELT(names, i))) {
            return (int)i;
        }
    }
    throw runtime_error("Cannot find column in the input: " + name);
}

column_view get_column(SEXP batch, SEXP names, const string& name, R_xlen_t nrow) {
    SEXP col = VECTOR_ELT(batch, find_column(names, name));
    if (Rf_xlength(col) != nrow) {
        throw runtime_error("columns of different lengths in the input batch");
    }
    try {
        return column_view(col);
    } catch (std::exception&) {
        throw runtime_error("Unsupported type of column: " + name);
    }
}

void split(const char* s, char sep, vector<string>* out) {
    out->clear();
    const char* start = s;
    for (const char* c = s;; c++) {
        if (*c == sep || *c == 0) {
            out->emplace_back(start, c - start);
            if (*c == 0) {
                break;
            }
            start = c + 1;
        }
    }
}

// one batch
void add_batch_rows(stream_state* st, SEXP batch);

void add_batch(stream_state* st, SEXP batch) {
    add_batch_rows(st, batch);
    if (st->quant.size() > (size_t)IQ_MAX_VALUES) {  // values kept in memory (passed on with a writer)
        throw runtime_error(too_many_values_message);
    }
}

void add_batch_rows(stream_state* st, SEXP batch) {
    if (TYPEOF(batch) != VECSXP) {
        throw runtime_error("a batch must be a data frame or a named list");
    }
    SEXP names = Rf_getAttrib(batch, R_NamesSymbol);
    if (names == R_NilValue) {
        throw runtime_error("a batch must have column names");
    }

    if (Rf_xlength(batch) == 0) {
        return;
    }
    R_xlen_t nrow = Rf_xlength(VECTOR_ELT(batch, 0));

    column_view sample = get_column(batch, names, st->col_sample, nrow);
    column_view primary = get_column(batch, names, st->col_primary_id, nrow);

    vector<column_view> secondary;
    for (auto& c : st->col_secondary_ids) {
        secondary.push_back(get_column(batch, names, c, nrow));
    }
    vector<column_view> quant_cols;
    for (auto& c : st->col_quant) {
        quant_cols.push_back(get_column(batch, names, c, nrow));
    }
    vector<column_view> quant_id_cols;
    for (auto& c : st->col_quant_id) {
        quant_id_cols.push_back(get_column(batch, names, c, nrow));
    }
    vector<column_view> annotation;
    for (auto& c : st->col_annotations) {
        annotation.push_back(get_column(batch, names, c, nrow));
    }

    vector<column_view> f_eq, f_neq, f_less, f_greater;
    for (auto& f : st->filter_string_equal) {
        f_eq.push_back(get_column(batch, names, f.first, nrow));
    }
    for (auto& f : st->filter_string_not_equal) {
        f_neq.push_back(get_column(batch, names, f.first, nrow));
    }
    for (auto& f : st->filter_double_less) {
        f_less.push_back(get_column(batch, names, f.first, nrow));
    }
    for (auto& f : st->filter_double_greater) {
        f_greater.push_back(get_column(batch, names, f.first, nrow));
    }

    // split mode: a single text column holding several values
    const bool split_mode = st->has_sep && quant_cols[0].is_text();

    vector<string> pieces, id_pieces;
    vector<double> values;          // per fragment: value, NAN if missing
    vector<string> suffix;          // per fragment: ion suffix
    string ion_prefix, key, display;
    char num[32];

    // is a text value missing? (same rule as the text-file reader)
    auto text_missing = [&](const char* t) {
        return *t == 0 || strcmp(t, "NaN") == 0 || (st->has_na_string && st->na_string == t);
    };

    for (R_xlen_t r = 0; r < nrow; r++) {
        st->n_rows++;

        // sample: every sample in the input is registered, also one whose rows are all
        // filtered out or missing (it is then a column of NA in the result)
        const char* sid = sample.text(r);
        int s;
        auto ss = st->map_samples.find(sid);
        if (ss != st->map_samples.end()) {
            s = ss->second;
        } else {
            s = (int)st->samples.size();
            st->samples.emplace_back(sid);
            st->map_samples.emplace(sid, s);
        }

        // filters
        bool ok = true;
        for (size_t i = 0; ok && i < f_eq.size(); i++) {
            ok = strcmp(f_eq[i].text(r), st->filter_string_equal[i].second.c_str()) == 0;
        }
        for (size_t i = 0; ok && i < f_neq.size(); i++) {
            ok = strcmp(f_neq[i].text(r), st->filter_string_not_equal[i].second.c_str()) != 0;
        }
        // a missing value (NA, NaN, empty or non-numeric text) keeps the row
        for (size_t i = 0; ok && i < f_less.size(); i++) {
            double v = f_less[i].number(r);
            if (!std::isnan(v) && v >= st->filter_double_less[i].second) {
                ok = false;
            }
        }
        for (size_t i = 0; ok && i < f_greater.size(); i++) {
            double v = f_greater[i].number(r);
            if (!std::isnan(v) && v <= st->filter_double_greater[i].second) {
                ok = false;
            }
        }
        if (!ok) {
            continue;
        }

        // quantitative values of this row
        values.clear();
        suffix.clear();
        if (split_mode) {
            split(quant_cols[0].text(r), st->sep, &pieces);
            if (!quant_id_cols.empty()) {
                split(quant_id_cols[0].text(r), st->sep, &id_pieces);
                if (id_pieces.size() != pieces.size()) {
                    throw runtime_error("Row " + std::to_string(st->n_rows) +
                                        ": the number of quant values and ids differ.");
                }
            }
            for (size_t k = 0; k < pieces.size(); k++) {
                const char* t = pieces[k].c_str();
                values.push_back(text_missing(t) ? NAN : atof(t));
                if (quant_id_cols.empty()) {
                    snprintf(num, sizeof(num), "%d", (int)(k + 1));
                    suffix.emplace_back(num);
                } else {
                    suffix.push_back(id_pieces[k]);
                }
            }
        } else {
            for (size_t k = 0; k < quant_cols.size(); k++) {
                double v;
                if (quant_cols[k].is_text()) {
                    const char* t = quant_cols[k].text(r);
                    v = text_missing(t) ? NAN : atof(t);
                } else {
                    v = quant_cols[k].number(r);
                    if (std::isnan(v) || (st->na_string_is_number && v == st->na_number)) {
                        v = NAN;
                    }
                }
                values.push_back(v);
                if (quant_id_cols.empty()) {
                    snprintf(num, sizeof(num), "%d", (int)(k + 1));
                    suffix.emplace_back(num);
                } else {
                    suffix.emplace_back(quant_id_cols[k].text(r));
                }
            }
        }

        int n_values = 0;
        for (double v : values) {
            if (!std::isnan(v)) {
                n_values++;
            }
        }
        if (n_values == 0) {
            continue;
        }
        st->n_rows_kept++;

        // protein
        const char* pid = primary.text(r);
        int p;
        auto pp = st->map_proteins.find(pid);
        if (pp != st->map_proteins.end()) {
            p = pp->second;
        } else {
            p = (int)st->annotation.size();
            vector<string> a;
            a.emplace_back(pid);
            for (auto& c : annotation) {
                a.emplace_back(c.text(r));
            }
            st->annotation.push_back(std::move(a));
            st->map_proteins.emplace(pid, p);
        }

        // ion prefix: secondary columns
        ion_prefix.clear();
        for (size_t j = 0; j < secondary.size(); j++) {
            if (j > 0) {
                ion_prefix += '\001';
            }
            ion_prefix += secondary[j].text(r);
        }

        for (size_t k = 0; k < values.size(); k++) {
            if (std::isnan(values[k])) {
                continue;
            }
            key = ion_prefix;
            key += '\001';
            key += suffix[k];

            int i;
            auto ii = st->map_ions.find(key);
            if (ii != st->map_ions.end()) {
                i = ii->second;
            } else {
                i = (int)st->ions.size();
                display = key;
                for (auto& c : display) {
                    if (c == '\001') {
                        c = '_';
                    }
                }
                st->ions.push_back(display);
                st->ion_prefixes.push_back(display.substr(0, ion_prefix.size()));  // '\001' -> '_' keeps the length
                st->map_ions.emplace(key, i);
            }

            st->protein_index.push_back(p + 1);  // 1-based for R
            st->sample_index.push_back(s + 1);
            st->ion_index.push_back(i + 1);
            st->quant.push_back(values[k]);
        }
    }
}

// Pass the values read so far to an iq format writer and free them (long_format_to_iq_format):
// only the protein, sample and ion lists stay in memory. Values without a protein or sample
// identifier are left out, as in process().
void write_values(stream_state* st, iq_writer* w) {
    for (size_t k = 0; k < st->quant.size(); k++) {
        int p = st->protein_index[k];
        int s = st->sample_index[k];
        if (st->annotation[p - 1][0].empty() || st->samples[s - 1].empty()) {
            w->n_unidentified++;
            continue;
        }
        w->add(p, s, st->ion_index[k], st->quant[k]);
    }
    st->n_values_written += (double)st->quant.size();
    vector<int>().swap(st->protein_index);
    vector<int>().swap(st->sample_index);
    vector<int>().swap(st->ion_index);
    vector<double>().swap(st->quant);
}

SEXP make_result(stream_state* st, bool verbose) {
    if (verbose) {
        Rprintf("\n# rows read                           = %.0f\n", (double)st->n_rows);
        Rprintf("# quantitative values after filtering = %.0f\n\n", (double)st->quant.size() + st->n_values_written);
        Rprintf("# samples  = %d\n", (int)st->samples.size());
        Rprintf("# proteins = %d\n", (int)st->annotation.size());
    }

    if (st->ions.empty()) {
        return R_NilValue;
    }

    size_t n = st->quant.size();

    SEXP p_index = PROTECT(Rf_allocMatrix(INTSXP, (int)n, 1));
    std::copy(st->protein_index.begin(), st->protein_index.end(), INTEGER(p_index));

    SEXP s_index = PROTECT(Rf_allocMatrix(INTSXP, (int)n, 1));
    std::copy(st->sample_index.begin(), st->sample_index.end(), INTEGER(s_index));

    SEXP i_index = PROTECT(Rf_allocMatrix(INTSXP, (int)n, 1));
    std::copy(st->ion_index.begin(), st->ion_index.end(), INTEGER(i_index));

    SEXP q = PROTECT(Rf_allocMatrix(REALSXP, (int)n, 1));
    std::copy(st->quant.begin(), st->quant.end(), REAL(q));

    SEXP ion_str = PROTECT(Rf_allocMatrix(STRSXP, (int)st->ions.size(), 1));
    for (size_t i = 0; i < st->ions.size(); i++) {
        SET_STRING_ELT(ion_str, i, Rf_mkCharCE(st->ions[i].c_str(), CE_UTF8));
    }

    SEXP ion_prefix_str = PROTECT(Rf_allocMatrix(STRSXP, (int)st->ion_prefixes.size(), 1));
    for (size_t i = 0; i < st->ion_prefixes.size(); i++) {
        SET_STRING_ELT(ion_prefix_str, i, Rf_mkCharCE(st->ion_prefixes[i].c_str(), CE_UTF8));
    }

    SEXP sample_str = PROTECT(Rf_allocMatrix(STRSXP, (int)st->samples.size(), 1));
    for (size_t i = 0; i < st->samples.size(); i++) {
        SET_STRING_ELT(sample_str, i, Rf_mkCharCE(st->samples[i].c_str(), CE_UTF8));
    }

    size_t np = st->annotation.size();
    size_t na = st->col_annotations.size() + 1;
    SEXP ann = PROTECT(Rf_allocMatrix(STRSXP, (int)np, (int)na));
    for (size_t i = 0; i < np; i++) {
        for (size_t j = 0; j < na; j++) {
            SET_STRING_ELT(ann, i + j * np, Rf_mkCharCE(st->annotation[i][j].c_str(), CE_UTF8));
        }
    }

    SEXP col_names = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t)na));
    SET_STRING_ELT(col_names, 0, Rf_mkCharCE(st->col_primary_id.c_str(), CE_UTF8));
    for (size_t j = 1; j < na; j++) {
        SET_STRING_ELT(col_names, j, Rf_mkCharCE(st->col_annotations[j - 1].c_str(), CE_UTF8));
    }
    SEXP dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dimnames, 1, col_names);
    Rf_setAttrib(ann, R_DimNamesSymbol, dimnames);

    SEXP quant_tab = PROTECT(Rf_allocVector(VECSXP, 4));
    SET_VECTOR_ELT(quant_tab, 0, p_index);
    SET_VECTOR_ELT(quant_tab, 1, s_index);
    SET_VECTOR_ELT(quant_tab, 2, i_index);
    SET_VECTOR_ELT(quant_tab, 3, q);

    SEXP quant_tab_names = PROTECT(Rf_allocVector(STRSXP, 4));
    SET_STRING_ELT(quant_tab_names, 0, Rf_mkChar("protein_list"));
    SET_STRING_ELT(quant_tab_names, 1, Rf_mkChar("sample_list"));
    SET_STRING_ELT(quant_tab_names, 2, Rf_mkChar("id"));
    SET_STRING_ELT(quant_tab_names, 3, Rf_mkChar("quant"));
    Rf_setAttrib(quant_tab, R_NamesSymbol, quant_tab_names);

    SEXP vec = PROTECT(Rf_allocVector(VECSXP, 5));
    SET_VECTOR_ELT(vec, 0, quant_tab);
    SET_VECTOR_ELT(vec, 1, ann);
    SET_VECTOR_ELT(vec, 2, sample_str);
    SET_VECTOR_ELT(vec, 3, ion_str);
    SET_VECTOR_ELT(vec, 4, ion_prefix_str);

    SEXP names = PROTECT(Rf_allocVector(STRSXP, 5));
    SET_STRING_ELT(names, 0, Rf_mkChar("quant_table"));
    SET_STRING_ELT(names, 1, Rf_mkChar("protein"));
    SET_STRING_ELT(names, 2, Rf_mkChar("sample"));
    SET_STRING_ELT(names, 3, Rf_mkChar("ion"));
    SET_STRING_ELT(names, 4, Rf_mkChar("ion_prefix"));
    Rf_setAttrib(vec, R_NamesSymbol, names);

    UNPROTECT(14);
    return vec;
}

// error text kept outside C++ scopes, so Rf_error (a longjmp) skips no destructor
char error_message[2048];

void set_error(const char* msg) {
    snprintf(error_message, sizeof(error_message), "%s", msg);
}

}  // namespace

// Entry points: all C++ objects live inside the try blocks; Rf_error is called
// only after they are destroyed.

SEXP iq_stream_begin(SEXP options) {
    stream_state* st = NULL;
    bool failed = false;
    try {
        st = create_state(options);
    } catch (std::exception& e) {
        set_error(e.what());
        failed = true;
    }
    if (failed) {
        Rf_error("%s", error_message);
    }
    SEXP ptr = PROTECT(R_MakeExternalPtr(st, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(ptr, finalize_state, TRUE);
    UNPROTECT(1);
    return ptr;
}

SEXP iq_stream_add(SEXP ptr, SEXP batch) {
    bool failed = false;
    try {
        add_batch(get_state(ptr), batch);
    } catch (std::exception& e) {
        set_error(e.what());
        failed = true;
    }
    if (failed) {
        Rf_error("%s", error_message);
    }
    return R_NilValue;
}

SEXP iq_stream_write(SEXP ptr, SEXP writer) {
    bool failed = false;
    try {
        write_values(get_state(ptr), get_writer(writer));
    } catch (std::exception& e) {
        set_error(e.what());
        failed = true;
    }
    if (failed) {
        Rf_error("%s", error_message);
    }
    return R_NilValue;
}

SEXP iq_stream_end(SEXP ptr, SEXP verbose) {
    stream_state* st = NULL;
    bool failed = false;
    try {
        st = get_state(ptr);
        if (st->quant.size() > (size_t)IQ_MAX_VALUES) {  // the R vectors of make_result have an int length
            throw runtime_error(too_many_values_message);
        }
    } catch (std::exception& e) {
        set_error(e.what());
        failed = true;
    }
    if (failed) {
        Rf_error("%s", error_message);
    }
    // the state stays owned by ptr (freed by its finalizer) until the result is built
    SEXP res = PROTECT(make_result(st, Rf_asLogical(verbose) == TRUE));
    finalize_state(ptr);
    UNPROTECT(1);
    return res;
}
