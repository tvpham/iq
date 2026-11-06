/*#########################################################################

 Author: Thang V. Pham {t.pham@amsterdamumc.nl}

 All rights reserved.

 Citation:

 Pham TV, Tran CTM, Henneman AA, Pham LHC, Le DG, Can AH, Bui PHL, 
 Piersma SR, Jimenez CR. Boosting the speed and accuracy of protein 
 quantification algorithms in mass spectrometry-based proteomics.
 Journal of Proteome Research 2026 Feb 6; 25(2):1198-1203.

 Software version: 2.2.0

#########################################################################*/

// R entry points of the iq format writer (iq_writer.h), used by long_format_to_iq_format.
// The readers add their values with iq_filter_write (text file, iq.cpp) and iq_stream_write
// (parquet file, stream.cpp). iq_read_sample_values reads the values of some samples back
// from the protein files, for the sample medians; iq_read_protein_batch reads a batch of
// protein files for process_iq_format.

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <vector>
#include <stdexcept>
#include <string>

#include "iq_writer.h"

extern "C" SEXP iq_writer_begin(SEXP);
extern "C" SEXP iq_writer_end(SEXP);
extern "C" SEXP iq_read_sample_values(SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP iq_read_protein_batch(SEXP, SEXP, SEXP, SEXP, SEXP);

namespace {

char writer_error[2048];

void finalize_writer(SEXP ptr) {
    iq_writer* w = static_cast<iq_writer*>(R_ExternalPtrAddr(ptr));
    if (w) {
        delete w;
        R_ClearExternalPtr(ptr);
    }
}

SEXP element(SEXP list, const char* name) {
    SEXP names = Rf_getAttrib(list, R_NamesSymbol);
    for (R_xlen_t i = 0; i < Rf_xlength(list); i++) {
        if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
            return VECTOR_ELT(list, i);
        }
    }
    return R_NilValue;
}

}  // namespace

// options: dir (ends with '/'), cutoff (NULL for none), buffer (megabytes), statistics
// (TRUE to keep the fragment statistics of each protein)
SEXP iq_writer_begin(SEXP options) {
    iq_writer* w = NULL;
    bool failed = false;
    try {
        SEXP dir = element(options, "dir");
        SEXP cutoff = element(options, "cutoff");
        SEXP buffer = element(options, "buffer");
        SEXP statistics = element(options, "statistics");
        if (TYPEOF(dir) != STRSXP || Rf_xlength(dir) != 1 ||
            TYPEOF(buffer) != REALSXP || Rf_xlength(buffer) != 1 || !(REAL(buffer)[0] > 0)) {
            throw std::runtime_error("iq_writer_begin: dir and buffer are needed");
        }
        bool has_cutoff = cutoff != R_NilValue;
        if (has_cutoff && (TYPEOF(cutoff) != REALSXP || Rf_xlength(cutoff) != 1)) {
            throw std::runtime_error("iq_writer_begin: cutoff must be one number or NULL");
        }
        w = new iq_writer(CHAR(STRING_ELT(dir, 0)), has_cutoff, has_cutoff ? REAL(cutoff)[0] : 0.0,
                          (size_t)(REAL(buffer)[0] * 1024.0 * 1024.0),
                          TYPEOF(statistics) == LGLSXP && Rf_xlength(statistics) == 1 && LOGICAL(statistics)[0] == TRUE);
    } catch (std::exception& e) {
        snprintf(writer_error, sizeof(writer_error), "%s", e.what());
        failed = true;
    }
    if (failed) {
        Rf_error("%s", writer_error);
    }
    SEXP ptr = PROTECT(R_MakeExternalPtr(w, Rf_install("iq_writer"), R_NilValue));
    R_RegisterCFinalizerEx(ptr, finalize_writer, TRUE);
    UNPROTECT(1);
    return ptr;
}

// Writes the remaining data and frees the writer. Returns a list: sample_count (values
// written per sample), protein_has_file, hist_bin and hist_count (histogram of all finite
// log2 values: bin b holds [b, b + 1) / hist_bins_per_unit), hist_bins_per_unit, log2_min,
// log2_max, n_written, n_not_positive, n_below_cutoff, n_unidentified, and with statistics
// protein_xic (per protein) and protein_ions (a list: the ions of each protein, first
// appearance order); NULL without.
SEXP iq_writer_end(SEXP ptr) {
    iq_writer* w = NULL;
    bool failed = false;
    try {
        w = get_writer(ptr);
        w->flush();
    } catch (std::exception& e) {
        snprintf(writer_error, sizeof(writer_error), "%s", e.what());
        failed = true;
    }
    if (failed) {
        Rf_error("%s", writer_error);
    }

    SEXP count = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)w->sample_count.size()));
    for (size_t s = 0; s < w->sample_count.size(); s++) {
        REAL(count)[s] = w->sample_count[s];
    }
    SEXP has_file = PROTECT(Rf_allocVector(LGLSXP, (R_xlen_t)w->protein_has_file.size()));
    for (size_t p = 0; p < w->protein_has_file.size(); p++) {
        LOGICAL(has_file)[p] = w->protein_has_file[p] ? TRUE : FALSE;
    }

    SEXP hist_bin = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)w->histogram.size()));
    SEXP hist_count = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)w->histogram.size()));
    R_xlen_t k = 0;
    for (const auto& h : w->histogram) {
        REAL(hist_bin)[k] = (double)h.first;
        REAL(hist_count)[k] = h.second;
        k++;
    }

    SEXP xic = R_NilValue;
    SEXP ions = R_NilValue;
    if (w->statistics) {
        size_t np = w->protein_xic.size();
        xic = Rf_allocVector(REALSXP, (R_xlen_t)np);
        PROTECT(xic);
        ions = Rf_allocVector(VECSXP, (R_xlen_t)np);
        PROTECT(ions);
        for (size_t p = 0; p < np; p++) {
            REAL(xic)[p] = (double)w->protein_xic[p];
            SEXP v = Rf_allocVector(INTSXP, (R_xlen_t)w->protein_ions[p].size());
            SET_VECTOR_ELT(ions, (R_xlen_t)p, v);
            std::copy(w->protein_ions[p].begin(), w->protein_ions[p].end(), INTEGER(v));
        }
    } else {
        PROTECT(xic);
        PROTECT(ions);
    }

    const int n = 13;
    const char* names[n] = {"sample_count", "protein_has_file", "hist_bin", "hist_count", "hist_bins_per_unit",
                            "log2_min", "log2_max", "n_written", "n_not_positive", "n_below_cutoff", "n_unidentified",
                            "protein_xic", "protein_ions"};
    SEXP res = PROTECT(Rf_allocVector(VECSXP, n));
    SET_VECTOR_ELT(res, 0, count);
    SET_VECTOR_ELT(res, 1, has_file);
    SET_VECTOR_ELT(res, 2, hist_bin);
    SET_VECTOR_ELT(res, 3, hist_count);
    SET_VECTOR_ELT(res, 4, Rf_ScalarReal(iq_writer::HIST_BINS_PER_UNIT));
    SET_VECTOR_ELT(res, 5, Rf_ScalarReal(w->log2_min));
    SET_VECTOR_ELT(res, 6, Rf_ScalarReal(w->log2_max));
    SET_VECTOR_ELT(res, 7, Rf_ScalarReal(w->n_written));
    SET_VECTOR_ELT(res, 8, Rf_ScalarReal(w->n_not_positive));
    SET_VECTOR_ELT(res, 9, Rf_ScalarReal(w->n_below_cutoff));
    SET_VECTOR_ELT(res, 10, Rf_ScalarReal(w->n_unidentified));
    SET_VECTOR_ELT(res, 11, xic);
    SET_VECTOR_ELT(res, 12, ions);
    SEXP res_names = PROTECT(Rf_allocVector(STRSXP, n));
    for (int j = 0; j < n; j++) {
        SET_STRING_ELT(res_names, j, Rf_mkChar(names[j]));
    }
    Rf_setAttrib(res, R_NamesSymbol, res_names);

    finalize_writer(ptr);
    UNPROTECT(8);
    return res;
}

// The log2 values of some samples, read from the protein files written by the writer, in one
// pass over the files: dir (ends with '/'), proteins (the protein numbers that have a file),
// samples (the sample numbers wanted), counts (the number of values written for each, as
// counted by the writer). Returns a list with a numeric vector for each sample in samples,
// filled in place (no second copy of the values).
SEXP iq_read_sample_values(SEXP dir, SEXP proteins, SEXP samples, SEXP counts) {
    if (TYPEOF(dir) != STRSXP || Rf_xlength(dir) != 1 || TYPEOF(proteins) != INTSXP ||
        TYPEOF(samples) != INTSXP || TYPEOF(counts) != REALSXP || Rf_xlength(counts) != Rf_xlength(samples)) {
        Rf_error("iq_read_sample_values: invalid arguments");
    }
    R_xlen_t n_samples = Rf_xlength(samples);
    int max_sample = 0;
    for (R_xlen_t j = 0; j < n_samples; j++) {
        if (INTEGER(samples)[j] < 1 || !(REAL(counts)[j] >= 0)) {
            Rf_error("iq_read_sample_values: invalid sample number or count");
        }
        if (INTEGER(samples)[j] > max_sample) {
            max_sample = INTEGER(samples)[j];
        }
    }

    // the result, allocated before any C++ object (Rf_error skips no destructor)
    SEXP res = PROTECT(Rf_allocVector(VECSXP, n_samples));
    for (R_xlen_t j = 0; j < n_samples; j++) {
        SET_VECTOR_ELT(res, j, Rf_allocVector(REALSXP, (R_xlen_t)REAL(counts)[j]));
    }

    bool failed = false;
    try {
        std::string d = CHAR(STRING_ELT(dir, 0));

        // sample number -> position in samples, -1 if not wanted
        std::vector<R_xlen_t> position(max_sample + 1, -1);
        std::vector<double*> out(n_samples);
        std::vector<R_xlen_t> filled(n_samples, 0);
        for (R_xlen_t j = 0; j < n_samples; j++) {
            position[INTEGER(samples)[j]] = j;
            out[j] = REAL(VECTOR_ELT(res, j));
        }

        char line[256];
        for (R_xlen_t k = 0; k < Rf_xlength(proteins); k++) {
            std::string name = d + std::to_string(INTEGER(proteins)[k]) + ".tsv";
            FILE* f = fopen(name.c_str(), "r");
            if (f == NULL) {
                throw std::runtime_error("Cannot read " + name);
            }
            while (fgets(line, sizeof(line), f)) {  // sample <tab> ion <tab> intensity
                char* end;
                long s = strtol(line, &end, 10);
                if (s < 1 || s > max_sample || position[s] < 0) {
                    continue;
                }
                R_xlen_t j = position[s];
                if (filled[j] >= Rf_xlength(VECTOR_ELT(res, j))) {
                    fclose(f);
                    throw std::runtime_error("More values than written for a sample in " + name);
                }
                char* ion_end;
                strtol(end, &ion_end, 10);
                out[j][filled[j]++] = std::log2(strtod(ion_end, NULL));
            }
            fclose(f);
        }
        for (R_xlen_t j = 0; j < n_samples; j++) {
            if (filled[j] != Rf_xlength(VECTOR_ELT(res, j))) {
                throw std::runtime_error("Fewer values than written for a sample in " + d);
            }
        }
    } catch (std::exception& e) {
        snprintf(writer_error, sizeof(writer_error), "%s", e.what());
        failed = true;
    }
    if (failed) {
        Rf_error("%s", writer_error);
    }

    UNPROTECT(1);
    return res;
}

namespace {

// one line of a protein file: sample <tab> ion <tab> intensity
bool parse_value_line(const char* line, long* s, long* i, double* v) {
    char* end;
    *s = strtol(line, &end, 10);
    if (end == line) {
        return false;
    }
    const char* p = end;
    *i = strtol(p, &end, 10);
    if (end == p) {
        return false;
    }
    p = end;
    *v = strtod(p, &end);
    return end != p;
}

}  // namespace

// A batch of protein files of the iq format, for iq_process (process_iq_format): starting at
// position start (1-based) of proteins (the protein numbers with a file), the files are read
// until the batch holds at least max_values values (at least one protein). The values are
// returned as iq_process takes them: protein_index (1, 2, ... within the batch), ion_index,
// sample_index and quant = log2(intensity) + shift[sample] (the normalization of the sample);
// next_start is the position of the first protein not read. The files are read twice, first to count
// the values, so that the R vectors are allocated once, at their size.
SEXP iq_read_protein_batch(SEXP dir, SEXP proteins, SEXP start, SEXP max_values, SEXP shift) {
    if (TYPEOF(dir) != STRSXP || Rf_xlength(dir) != 1 || TYPEOF(proteins) != INTSXP ||
        TYPEOF(start) != INTSXP || Rf_xlength(start) != 1 || TYPEOF(max_values) != REALSXP ||
        Rf_xlength(max_values) != 1 || TYPEOF(shift) != REALSXP) {
        Rf_error("iq_read_protein_batch: invalid arguments");
    }
    R_xlen_t n = Rf_xlength(proteins);
    R_xlen_t first = (R_xlen_t)INTEGER(start)[0] - 1;
    if (first < 0 || first >= n) {
        Rf_error("iq_read_protein_batch: start out of range");
    }
    long n_samples = (long)Rf_xlength(shift);

    // pass 1: the proteins of the batch and their numbers of values
    bool failed = false;
    double total = 0;
    R_xlen_t last = first;  // one past the last protein of the batch
    try {
        std::string d = CHAR(STRING_ELT(dir, 0));
        char line[256];
        while (last < n && (last == first || total < REAL(max_values)[0])) {
            std::string name = d + std::to_string(INTEGER(proteins)[last]) + ".tsv";
            FILE* f = fopen(name.c_str(), "r");
            if (f == NULL) {
                throw std::runtime_error("Cannot read " + name);
            }
            long s, i;
            double v;
            while (fgets(line, sizeof(line), f)) {
                if (parse_value_line(line, &s, &i, &v)) {
                    total++;
                }
            }
            fclose(f);
            last++;
        }
    } catch (std::exception& e) {
        snprintf(writer_error, sizeof(writer_error), "%s", e.what());
        failed = true;
    }
    if (failed) {
        Rf_error("%s", writer_error);
    }

    R_xlen_t m = (R_xlen_t)total;
    SEXP p_index = PROTECT(Rf_allocVector(INTSXP, m));
    SEXP i_index = PROTECT(Rf_allocVector(INTSXP, m));
    SEXP s_index = PROTECT(Rf_allocVector(INTSXP, m));
    SEXP quant = PROTECT(Rf_allocVector(REALSXP, m));

    // pass 2: the values
    try {
        std::string d = CHAR(STRING_ELT(dir, 0));
        char line[256];
        R_xlen_t k = 0;
        for (R_xlen_t j = first; j < last; j++) {
            std::string name = d + std::to_string(INTEGER(proteins)[j]) + ".tsv";
            FILE* f = fopen(name.c_str(), "r");
            if (f == NULL) {
                throw std::runtime_error("Cannot read " + name);
            }
            long s, i;
            double v;
            while (fgets(line, sizeof(line), f)) {
                if (!parse_value_line(line, &s, &i, &v)) {
                    continue;
                }
                if (k >= m || s < 1 || s > n_samples || i < 1 || i > INT_MAX) {
                    fclose(f);
                    throw std::runtime_error("Invalid or changed protein file " + name);
                }
                INTEGER(p_index)[k] = (int)(j - first + 1);
                INTEGER(i_index)[k] = (int)i;
                INTEGER(s_index)[k] = (int)s;
                REAL(quant)[k] = std::log2(v) + REAL(shift)[s - 1];
                k++;
            }
            fclose(f);
        }
        if (k != m) {
            throw std::runtime_error("The protein files changed while they were read in " + d);
        }
    } catch (std::exception& e) {
        snprintf(writer_error, sizeof(writer_error), "%s", e.what());
        failed = true;
    }
    if (failed) {
        Rf_error("%s", writer_error);
    }

    const char* names[] = {"protein_index", "ion_index", "sample_index", "quant", "next_start"};
    SEXP res = PROTECT(Rf_allocVector(VECSXP, 5));
    SET_VECTOR_ELT(res, 0, p_index);
    SET_VECTOR_ELT(res, 1, i_index);
    SET_VECTOR_ELT(res, 2, s_index);
    SET_VECTOR_ELT(res, 3, quant);
    SET_VECTOR_ELT(res, 4, Rf_ScalarInteger((int)(last + 1)));
    SEXP res_names = PROTECT(Rf_allocVector(STRSXP, 5));
    for (int j = 0; j < 5; j++) {
        SET_STRING_ELT(res_names, j, Rf_mkChar(names[j]));
    }
    Rf_setAttrib(res, R_NamesSymbol, res_names);
    UNPROTECT(6);
    return res;
}
