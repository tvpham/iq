/*#########################################################################

 Author: Thang V. Pham {t.pham@amsterdamumc.nl}

 All rights reserved.

 Citations:

 Pham TV, Tran CTM, Henneman AA, Pham LHC, Le DG, Can AH, Bui PHL, 
 Piersma SR, Jimenez CR. Boosting the speed and accuracy of protein 
 quantification algorithms in mass spectrometry-based proteomics.
 Journal of Proteome Research 2026 Feb 6; 25(2):1198-1203.
 
 Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative
 protein abundances from ion quantification in DIA-MS-based proteomics,
 Bioinformatics 2020 Apr 15;36(8):2611-2613.

 Software version: 2.2.0

#########################################################################*/

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <utility>
#include <exception>
#include <functional>
#include <stdexcept>

// The original MaxLFQ (Eigen) is compiled only with IQ_USE_EIGEN, see src/Makevars. Without it,
// "maxlfq", "weighted_maxlfq" and fast_MaxLFQ() use the solver of "maxlfq_bit".
#ifdef IQ_USE_EIGEN
#include <Eigen/Dense>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#define R_NO_REMAP

#include <R.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>

#include "utils.h"
#ifdef IQ_USE_EIGEN
#include "maxlfq.h"
#endif
#include "maxlfq_bit.h"
#include "coordinate_descent.h"
#include "iq_writer.h"

extern "C" SEXP iq_filter(SEXP);
extern "C" SEXP iq_filter_write(SEXP, SEXP);
extern "C" SEXP iq_MaxLFQ(SEXP);
extern "C" SEXP iq_quant(SEXP);
extern "C" SEXP iq_process(SEXP);
extern "C" SEXP iq_rescale(SEXP);
extern "C" SEXP iq_eigen(void);

// streaming reader, src/stream.cpp
extern "C" SEXP iq_stream_begin(SEXP);
extern "C" SEXP iq_stream_add(SEXP, SEXP);
extern "C" SEXP iq_stream_end(SEXP, SEXP);
extern "C" SEXP iq_stream_write(SEXP, SEXP);

// iq format writer, src/writer.cpp
extern "C" SEXP iq_writer_begin(SEXP);
extern "C" SEXP iq_writer_end(SEXP);
extern "C" SEXP iq_read_sample_values(SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP iq_read_protein_batch(SEXP, SEXP, SEXP, SEXP, SEXP);

using namespace std;

//------------------------------ FILTER -------------------------------------

// header
#define BUFFER_SIZE 4096

#define FREAD_BUFSIZE 2097152

class utils {
public:

    // count the number of unique values in a large INT array
    static int count_unique(int *vec, size_t nrow) {
        int *m = std::max_element(vec, vec + nrow);
        auto tmp = new vector<bool>(*m + 1, false);  // include 0, for a reasonable size m
        for (size_t i = 0; i < nrow; i++) {
            (*tmp)[vec[i]] = true;
        }
        int total = 0;
        for (auto b : (*tmp)) {
            if (b) {
                total++;
            }
        }
        delete tmp;
        return total;
    }

    // map the column name to index value
    static void map_header(const char **column_names, size_t len, size_t *index, char *line, size_t *tab_pos, size_t n_tab_pos) {
        for (size_t i = 0; i < len; i++) {
            size_t j = 0;
            while (j < n_tab_pos && strcmp(column_names[i], line + tab_pos[j]) != 0) {
                j++;
            }

            if (j < n_tab_pos) {
                index[i] = j;
            } else {
                throw runtime_error((string("Cannot find column in the header : ") + column_names[i]).c_str());
            }
        }
    }

    // concatenate strings with '_' in between
    static string concatenate(const vector<string> &v) {
        string tmp = v[0];
        for (size_t i = 1; i < v.size(); i++) {
            tmp += ("_" + v[i]);
        }
        return tmp;
    }

    // Read a tab-separated line from buffer
    static int getline_tab(char **buffer, size_t *n, FILE *stream, size_t **positions, size_t *npos) {
        if (*buffer == NULL) {
            *buffer = (char *)malloc(BUFFER_SIZE);
            if (*buffer == NULL) {
                throw runtime_error("Cannot allocate memory");
            }
            *n = BUFFER_SIZE;

            size_t *new_positions = (size_t *)realloc(*positions, (BUFFER_SIZE + 1) * sizeof(size_t));  // if NULL, realloc = malloc
            if (new_positions == NULL) {
                free(*buffer);
                free(*positions);
                *buffer = NULL;
                *positions = NULL;
                throw runtime_error("Cannot allocate memory");
            }
            *positions = new_positions;
        }

        *npos = 0;

        size_t remaining_buffer_size = *n;

        char *ptr = *buffer;

        (*positions)[*npos] = 0;
        (*positions)[*npos + 1] = (*positions)[*npos];
        (*npos)++;

        int end_of_file = -1;  // nothing left

        while (true) {
            if (!fgets(ptr, (int)remaining_buffer_size, stream) || ferror(stream) || feof(stream)) {
                return end_of_file;
            }

            end_of_file = 0;  // there is something left

            // search for new line
            while (*ptr) {
                (*positions)[*npos]++;
                if (*ptr == '\t') {
                    *ptr = 0;
                    (*positions)[*npos + 1] = (*positions)[*npos];
                    (*npos)++;
                }
                if (*ptr == '\n' || *ptr == '\r') {
                    *ptr = 0;
                    return 0;  // all is good
                }
                ptr++;
                remaining_buffer_size--;
            }

            if (remaining_buffer_size > 1) {
                // stop because of end of file
                return 0;
            }

            // OK, more data available extend buffer
            // (grow into a temporary: a failed realloc returns NULL but keeps the old block)
            char *new_buffer = (char *)realloc(*buffer, *n * 2);
            if (new_buffer == NULL) {
                free(*buffer);
                free(*positions);
                *buffer = NULL;
                *positions = NULL;
                throw runtime_error("Cannot allocate memory");
            }
            *buffer = new_buffer;
            remaining_buffer_size = *n + 1;
            ptr = *buffer + (*n - 1);
            *n *= 2;

            size_t *new_positions = (size_t *)realloc(*positions, (*n + 1) * sizeof(size_t));  // if NULL, realloc = malloc
            if (new_positions == NULL) {
                free(*buffer);
                free(*positions);
                *buffer = NULL;
                *positions = NULL;
                throw runtime_error("Cannot allocate memory");
            }
            *positions = new_positions;
        }
    }

    static void split(char* c, const char* sep, vector<char*>* out) {
        out->clear();
        out->push_back(c);
        while (*c) {
            if (*c == *sep) {
                *c++ = 0;
                out->push_back(c);
            }
            else {
                c++;
            }
        }
    }

};

class string_vector_view {
   public:
    vector<const char *> str_vec;

    string_vector_view(const char *str) {
        str_vec = {str};
    }

    string_vector_view(const string_vector_view& vv) {
        str_vec = vv.str_vec;
    }

    bool operator==(const string_vector_view &rhs) const {
        if (str_vec.size() != rhs.str_vec.size()) {
            return false;
        }
        for (size_t i = 0; i < str_vec.size(); i++) {
            if (strcmp(str_vec[i], rhs.str_vec[i]) != 0) {
                return false;
            }
        }
        return true;
    }
};


// a simple hash function for unordered_map keys, K&R version 2
struct hash_fn {
    size_t operator()(const string_vector_view &svv) const {
        size_t hash_val = 0;
        for (auto s : svv.str_vec) {
            for (const char *p = s; *p; p++) {
                hash_val = *p + 31 * hash_val;
            }
        }
        return hash_val;
    }
};

class double_buffering_file_t {
    FILE *f;

    size_t n = FREAD_BUFSIZE;

    int current;

    size_t right_pos;
    size_t left_pos;

public:
    char *buf_0;
    char *buf_1;
    bool eof;

    // error checking in the calling function
    double_buffering_file_t(FILE *f) : f(f) {
        eof = false;

        buf_0 = (char *)malloc(n);

        buf_1 = (char *)malloc(n);

        current = 0;
        right_pos = n;
        left_pos = n;  // nothing left over from the other buffer
    }

    ~double_buffering_file_t() {
        if (buf_0) {
            free(buf_0);
        }
        if (buf_1) {
            free(buf_1);
        }
    }

    // reading multiple lines
    int fgetss(vector<char *> &str_vec) {
        size_t pos;  // current position

        // check if the last buffer was full. If so, go back and increase buffer.
        if (left_pos == 0) {
            current = 1 - current;  // go back

            // grow into a temporary: a failed realloc returns NULL but keeps the old block,
            // which the destructor frees
            char *new_buf = (char *)realloc(buf_0, n * 2);
            if (new_buf == NULL) {
                return -1;
            }
            buf_0 = new_buf;

            new_buf = (char *)realloc(buf_1, n * 2);
            if (new_buf == NULL) {
                return -1;
            }
            buf_1 = new_buf;

            pos = right_pos;  // not 'n' because we did not read to the full of buffer
            n *= 2;
        } else {
            // copy leftover part
            pos = right_pos - left_pos;
            if (current == 0) {
                (void)memcpy((void *)buf_0, (void *)(buf_1 + left_pos), pos);
            } else {
                (void)memcpy((void *)buf_1, (void *)(buf_0 + left_pos), pos);
            }
        }

        // fill data from pos
        unsigned char *active_buffer;
        size_t last;

        if (current == 0) {
            last = fread(buf_0 + pos, 1, n - pos - 1, f);  // leave 1 for terminnating 0 for eof
            active_buffer = (unsigned char *)buf_0;
        } else {
            last = fread(buf_1 + pos, 1, n - pos - 1, f);
            active_buffer = (unsigned char *)buf_1;
        }
        eof = (last < (n - pos - 1)) || ferror(f) || feof(f);

        right_pos = pos + last;
        left_pos = 0;

        str_vec.clear();

        unsigned char *t;

        while ((left_pos < right_pos) && (t = (unsigned char *)memchr((void *)(active_buffer + left_pos), '\n', right_pos - left_pos))) {
            *t = 0;
            str_vec.push_back((char *)active_buffer + left_pos);
            left_pos = t - active_buffer + 1;
        }

        if (eof && (left_pos < right_pos)) {
            active_buffer[right_pos] = 0;
            str_vec.push_back((char *)active_buffer + left_pos);
            left_pos = right_pos;
        }
        current = 1 - current;  // switch buffer

        return 0;
    }
};


void process(const vector<string> &argv,
             vector<int> *protein_index,
             vector<int> *sample_index,
             vector<int> *ion_index,
             vector<double> *quant,
             vector<vector<string>*> *annotation,
             vector<string> *annotation_colnames,
             vector<string*> *samples,
             vector<vector<string>*> *ions,
             const std::function<void()> &on_chunk = std::function<void()>()) {
    // on_chunk: if set, called after each chunk of lines, with the values of that chunk in
    // protein_index, sample_index, ion_index and quant (iq format writer)
    //--- example input parameters
    /*
    const char* col_sample = "R.Condition";
    const char* col_primary_id = "PG.ProteinGroups";
    vector<const char*> col_secondary_ids({ "EG.ModifiedSequence", "FG.Charge", "F.FrgIon", "F.Charge" });
    const char* col_quant = "F.PeakArea";
    vector<const char*> col_annotations({ "PG.Genes", "PG.ProteinNames" });

    vector<pair<const char*, const char*>> filter_string_equal({ pair<const char*, const char*>("F.ExcludedFromQuantification", "False"),
                                                                  pair<const char*, const char*>("F.FrgLossType", "noloss") });

    vector<pair<const char*, double>> filter_double_less({ pair<const char*, double>("PG.Qvalue", 0.01),
                                                           pair<const char*, double>("EG.Qvalue", 0.01) });
    */

    const char *col_sample = "";
    const char *col_primary_id = "";
    vector<const char *> col_secondary_ids;
    vector<const char *> col_quants;      // one column (split by intensity_col_sep), or one per fragment
    vector<const char *> col_annotations;

    vector<pair<const char *, const char *>> filter_string_equal;
    vector<pair<const char*, const char*>> filter_string_not_equal;

    vector<pair<const char *, double>> filter_double_less;
    vector<pair<const char*, double>> filter_double_greater;

    const char* intensity_col_sep = 0;
    vector<const char *> col_quant_ids;   // as col_quants: one column (split), or one per intensity column
    const char* na_string = 0;

    if (argv.size() < 1) {
        /*
        std::cout << "\nUsage: [OPTION] input_file_name\n\n";
        std::cout << "Convert input file to a concise representation for iq.\n\n";
        std::cout << "  --header\n\tOnly show header of the input file\n";
        std::cout << "  --no-annotation\n\tNo annotation column\n";
        std::cout << "  --annotation COL(s) --|input_file_name\n\tFollowing entries are annotation columns until the next switch or input_file_name\n\t[PG.Genes PG.ProteinNames]\n";
        std::cout << "  --sample COL\n\tSample column [R.Condition]\n";
        std::cout << "  --quant COL\n\tQuantitative column [F.PeakArea]\n";
        std::cout << "  --primary COL\n\tProtein column [PG.ProteinGroups]\n";
        std::cout << "  --secondary COL(s) --|input_file_name\n\tFollowing entries are ion columns until the next switch or input_file_name\n\t[EG.ModifiedSequence FG.Charge F.FrgIon F.Charge]\n";
        std::cout << "  --clear-all-filters\n\tClear all filterings\n";
        std::cout << "  --filter-string-equal COL VALUE\n\tAdd a string filtering (COL == VALUE entries are kept)\n\t[F.ExcludedFromQuantification False] & [F.FrgLossType noloss]\n";
        std::cout << "  --filter-double-less COL VALUE\n\tAdd a double filtering (NaN COL and COL < VALUE entries are kept)\n\t[PG.Qvalue 0.01] & [EG.Qvalue 0.01]\n";
        std::cout << "  --output-dir VALUE\n\tDirectory for output files\n";
        */
        /*
        1.9.4
        --intensity_col_sep [;]
        --intensity_col_id [Fragment.Info]
        --na_string [0]
        */
        /*
        1.9.7
        --filter-string-not-equal
        --filter-double-greater
        */
        throw runtime_error("Not enough arguments");
    }

    bool show_header = false;
    bool verbose = true;   // --quiet: no command, columns, thread count or counts
    string output_dir = "";

    size_t p = 0;  // parameter
    while (p < (argv.size() - 1)) {
        if (argv[p].compare("--no-annotation") == 0) {
            col_annotations.clear();
            p++;
        } else if (argv[p].compare("--header") == 0) {
            show_header = true;
            p++;
        } else if (argv[p].compare("--quiet") == 0) {
            verbose = false;
            p++;
        } else if (argv[p].compare("--sample") == 0) {
            col_sample = argv[++p].c_str();
            p++;
        } else if (argv[p].compare("--annotation") == 0) {
            col_annotations.clear();
            p++;
            while (p < (argv.size() - 1)) {
                if (strncmp(argv[p].c_str(), "--", 2) == 0) {
                    break;
                }
                col_annotations.push_back(argv[p++].c_str());
            }
        } else if (argv[p].compare("--primary") == 0) {
            col_primary_id = argv[++p].c_str();
            p++;
        } else if (argv[p].compare("--quant") == 0) {
            col_quants.clear();
            p++;
            while (p < (argv.size() - 1)) {
                if (strncmp(argv[p].c_str(), "--", 2) == 0) {
                    break;
                }
                col_quants.push_back(argv[p++].c_str());
            }
        } else if (argv[p].compare("--secondary") == 0) {
            col_secondary_ids.clear();
            p++;
            while (p < (argv.size() - 1)) {
                if (strncmp(argv[p].c_str(), "--", 2) == 0) {
                    break;
                }
                col_secondary_ids.push_back(argv[p++].c_str());
            }
        } else if (argv[p].compare("--clear-all-filters") == 0) {
            filter_string_equal.clear();
            filter_string_not_equal.clear();
            filter_double_less.clear();
            filter_double_greater.clear();
            p++;
        } else if (argv[p].compare("--filter-string-equal") == 0) {
            if (p + 2 < (argv.size() - 1)) {
                filter_string_equal.push_back(make_pair(argv[p + 1].c_str(), argv[p + 2].c_str()));
                p += 3;
            } else {
                throw runtime_error("--filter-string-equal needs two parameters");
            }
        }
        else if (argv[p].compare("--filter-string-not-equal") == 0) {
            if (p + 2 < (argv.size() - 1)) {
                filter_string_not_equal.push_back(make_pair(argv[p + 1].c_str(), argv[p + 2].c_str()));
                p += 3;
            }
            else {
                throw runtime_error("--filter-string-not-equal needs two parameters");
            }
        } else if (argv[p].compare("--filter-double-less") == 0) {
            if (p + 2 < (argv.size() - 1)) {
                filter_double_less.push_back(make_pair(argv[p + 1].c_str(), atof(argv[p + 2].c_str())));
                p += 3;
            } else {
                throw runtime_error("--filter-double-less needs two parameters");
            }
        }
        else if (argv[p].compare("--filter-double-greater") == 0) {
            if (p + 2 < (argv.size() - 1)) {
                filter_double_greater.push_back(make_pair(argv[p + 1].c_str(), atof(argv[p + 2].c_str())));
                p += 3;
            }
            else {
                throw runtime_error("--filter-double-greater needs two parameters");
            }
        } else if (argv[p].compare("--output-dir") == 0) {
            if (p + 1 < (argv.size() - 1)) {
                output_dir = argv[p + 1];
                if (output_dir[output_dir.size() - 1] != '/') {
                    output_dir = output_dir + "/";
                }
                p += 2;
            } else {
                throw runtime_error("--output-dir needs one parameter");
            }
        } else if (argv[p].compare("--intensity_col_sep") == 0) {
            intensity_col_sep = argv[++p].c_str();
            p++;
        } else if (argv[p].compare("--intensity_col_id") == 0) {
            col_quant_ids.clear();
            p++;
            while (p < (argv.size() - 1)) {
                if (strncmp(argv[p].c_str(), "--", 2) == 0) {
                    break;
                }
                col_quant_ids.push_back(argv[p++].c_str());
            }
        } else if (argv[p].compare("--na_string") == 0) {
            na_string = argv[++p].c_str();
            p++;
        } else {
            throw runtime_error((string("Do not know what to do with ") + argv[p]).c_str());
        }
    }

    if (col_quants.empty()) {
        throw runtime_error("No intensity column");
    }
    if (col_secondary_ids.empty()) {  // the ions are made from these columns (col_secondary_ids_ind[0] below)
        throw runtime_error("No secondary_id column");
    }
    if (col_quants.size() > 1 && intensity_col_sep) {
        throw runtime_error("intensity_col_sep can only be used with a single intensity column");
    }
    if (!col_quant_ids.empty() && col_quant_ids.size() != col_quants.size()) {
        throw runtime_error("intensity_col_id must have the same length as intensity_col");
    }
    if (col_quants.size() == 1 && !col_quant_ids.empty() && !intensity_col_sep) {
        throw runtime_error("intensity_col_id with a single intensity column needs intensity_col_sep");
    }
    const bool has_quant_ids = !col_quant_ids.empty();

    //--- headers

    const char *filename = argv[(argv.size() - 1)].c_str();

    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        throw runtime_error((string("Check file name and path. Cannot open: ") + filename).c_str());
    }

    char *line = NULL;
    size_t len = 0;

    size_t *tab_pos = NULL;
    size_t n_tab_pos = 0;

    int header_status;
    try {
        header_status = utils::getline_tab(&line, &len, fp, &tab_pos, &n_tab_pos);
    } catch (...) {
        fclose(fp);  // getline_tab has already freed its buffers
        throw;
    }

    if (header_status == -1) {
        fclose(fp);
        if (line) {
            free(line);
        }
        if (tab_pos) {
            free(tab_pos);
        }
        throw runtime_error((string("Cannot parse header of ") + filename).c_str());
    }

    // files the reader cannot read: say why, instead of "Cannot find column in the header"
    {
        string problem;
        const unsigned char *b = (const unsigned char *)line;
        if (b[0] == 0x1f && b[1] == 0x8b) {
            problem = " is gzip-compressed; decompress it first, or read it into R and pass the data frame";
        } else {
            for (size_t i = 0; i < n_tab_pos; i++) {
                if (line[tab_pos[i]] == '"') {
                    problem = " has quoted fields (e.g. written by write.table with quote = TRUE); write it with "
                              "quote = FALSE, or read it into R and pass the data frame";
                    break;
                }
            }
        }
        if (!problem.empty()) {
            fclose(fp);
            free(line);
            if (tab_pos) {
                free(tab_pos);
            }
            throw runtime_error(string(filename) + problem);
        }
    }

    // a UTF-8 byte-order mark before the first column name (e.g. "CSV UTF-8" from Excel)
    if (n_tab_pos > 0 && strncmp(line + tab_pos[0], "\xEF\xBB\xBF", 3) == 0) {
        tab_pos[0] += 3;
    }

    if (show_header) {
        for (size_t i = 0; i < n_tab_pos; i++) {
            Rprintf("%s\n", line + tab_pos[i]);
        }
        //--- finish reading
        fclose(fp);
        if (line) {
            free(line);
        }
        if (tab_pos) {
            free(tab_pos);
        }
        //exit(0);
        return;
    }

    // display command
    if (verbose) {
        Rprintf("\nCommand: ");
        for (size_t i = 0; i < argv.size(); i++) {
            Rprintf("%s ", argv[i].c_str());
        }
        Rprintf("\n\nSample column:\n    %s\n", col_sample);
        Rprintf("Protein column:\n    %s\n", col_primary_id);
        Rprintf("Ion column(s):\n   ");
        for (auto str : col_secondary_ids) {
            Rprintf(" %s", str);
        }
        Rprintf("\nQuant column(s):\n   ");
        for (auto str : col_quants) {
            Rprintf(" %s", str);
        }
        Rprintf("\n");
        if (!col_annotations.empty()) {
            Rprintf("Annotation column(s):\n   ");
            for (auto str : col_annotations) {
                Rprintf(" %s", str);
            }
            Rprintf("\n");
        }

        if (!filter_string_equal.empty()) {
            Rprintf("String equal filter(s):\n");
            for (auto str : filter_string_equal) {
                Rprintf("    %s == %s\n", str.first, str.second);
            }
        }

        if (!filter_string_not_equal.empty()) {
            Rprintf("String not equal filter(s):\n");
            for (auto str : filter_string_not_equal) {
                Rprintf("    %s != %s\n", str.first, str.second);
            }
        }

        if (!filter_double_less.empty()) {
            Rprintf("Double less filter(s):\n");
            for (auto str : filter_double_less) {
                Rprintf("    %s < %f\n", str.first, str.second);
            }
        }

        if (!filter_double_greater.empty()) {
            Rprintf("Double greater filter(s):\n");
            for (auto str : filter_double_greater) {
                Rprintf("    %s > %f\n", str.first, str.second);
            }
        }

        if (intensity_col_sep) {
            Rprintf("Quant column separator:\n    %s\n", intensity_col_sep);
        }

        if (na_string) {
            Rprintf("Quant column NA string (in addition to NaN):\n    %s\n", na_string);
        }

        if (has_quant_ids) {
            Rprintf("Quant id column(s):\n   ");
            for (auto str : col_quant_ids) {
                Rprintf(" %s", str);
            }
            Rprintf("\n");
        }

    }

    //--- check header
    size_t n_cols = n_tab_pos;

    // close the file and free the header buffers if a column is not found
    auto map_header = [&](const char **column_names, size_t len, size_t *index) {
        try {
            utils::map_header(column_names, len, index, line, tab_pos, n_tab_pos);
        } catch (...) {
            fclose(fp);
            free(line);
            free(tab_pos);
            throw;
        }
    };

    size_t col_sample_ind;
    map_header(&col_sample, 1, &col_sample_ind);

    size_t col_primary_id_ind;
    map_header(&col_primary_id, 1, &col_primary_id_ind);

    vector<size_t> col_secondary_ids_ind(col_secondary_ids.size());
    map_header(col_secondary_ids.data(), col_secondary_ids.size(), col_secondary_ids_ind.data());

    vector<size_t> col_quants_ind(col_quants.size());
    map_header(col_quants.data(), col_quants.size(), col_quants_ind.data());

    vector<size_t> col_annotations_ind(col_annotations.size());
    map_header(col_annotations.data(), col_annotations.size(), col_annotations_ind.data());

    vector<size_t> col_filter_string_equal(filter_string_equal.size());
    for (size_t i = 0; i < filter_string_equal.size(); i++) {
        map_header(&(filter_string_equal[i].first), 1, &(col_filter_string_equal[i]));
    }

    vector<size_t> col_filter_string_not_equal(filter_string_not_equal.size());
    for (size_t i = 0; i < filter_string_not_equal.size(); i++) {
        map_header(&(filter_string_not_equal[i].first), 1, &(col_filter_string_not_equal[i]));
    }

    vector<size_t> col_filter_double_less(filter_double_less.size());
    for (size_t i = 0; i < filter_double_less.size(); i++) {
        map_header(&(filter_double_less[i].first), 1, &(col_filter_double_less[i]));
    }

    vector<size_t> col_filter_double_greater(filter_double_greater.size());
    for (size_t i = 0; i < filter_double_greater.size(); i++) {
        map_header(&(filter_double_greater[i].first), 1, &(col_filter_double_greater[i]));
    }

    annotation_colnames->push_back(string(col_primary_id));
    for (auto n : col_annotations) {
        annotation_colnames->push_back(string(n));
    }

    vector<size_t> col_quant_ids_ind(col_quant_ids.size());
    map_header(col_quant_ids.data(), col_quant_ids.size(), col_quant_ids_ind.data());


    if (line) {
        free(line);
    }
    if (tab_pos) {
        free(tab_pos);
    }

    //--- start building

    auto map_proteins = new unordered_map<string_vector_view, size_t, hash_fn>;
    auto map_ions = new unordered_map<string_vector_view, size_t, hash_fn>;
    auto map_samples = new unordered_map<string_vector_view, size_t, hash_fn>;

    size_t n_ion_ind = col_secondary_ids_ind.size();

    size_t line_no = 1;
    size_t n_after_filtered = 0;
    size_t print_threshold = 20;

    auto _line = new vector<char *>;
    auto _line_new = new vector<char *>;
    auto _tab = new vector<vector<char *>>;

    auto quant_vec = new vector< vector<char *> >; // vector of quant values
    auto quant_id_vec = new vector< vector<char *> >; // vector of quant values
    auto quant_vec_is_na = new vector< vector<bool> >; // vector of NA values
    auto quant_vec_is_na_n = new vector< int >; // vector of NA values
    auto row_valid = new vector< char >; // 1: the line has all columns (not blank); char, not bool: set in parallel

    vector<string*> numbers;

    auto do_protein = [&]() {
        for (size_t n = 0; n < _line->size(); n++) {

            if ((*quant_vec_is_na_n)[n] == 0) {
                continue;
            }

            auto pp = map_proteins->find(string_vector_view((*_tab)[n][col_primary_id_ind]));
            size_t p;

            if (pp != map_proteins->end()) {
                p = pp->second;
            } else {
                p = annotation->size();

                annotation->emplace_back(new vector<string>);

                (*annotation)[p]->emplace_back((*_tab)[n][col_primary_id_ind]);
                for (auto i : col_annotations_ind) {
                    (*annotation)[p]->emplace_back((*_tab)[n][i]);
                }

                map_proteins->emplace(((*annotation)[p]->at(0)).c_str(), p);
            }

            for (int i = 0; i < (*quant_vec_is_na_n)[n]; i++) {
                protein_index->push_back(p + 1);
            }
        }
    };

    auto do_sample = [&]() {

        for (size_t n = 0; n < _line->size(); n++) {

            // every sample in the input is registered, also one whose lines are all filtered
            // out or missing (it is then a column of NA in the result); not a blank line
            if (!(*row_valid)[n]) {
                continue;
            }

            auto ss = map_samples->find(string_vector_view((*_tab)[n][col_sample_ind]));

            size_t s;

            if (ss != map_samples->end()) {
                s = ss->second;
            } else {

                s = samples->size();

                samples->emplace_back(new string((*_tab)[n][col_sample_ind]));

                map_samples->emplace((*samples)[s]->c_str(), s);
            }

            if ((*quant_vec_is_na_n)[n] == 0) {
                continue;
            }

            for (int i = 0; i < (*quant_vec_is_na_n)[n]; i++) {
                sample_index->push_back(s + 1); // 1-based vector for R
            }

            for (size_t i = 0; i < (*quant_vec)[n].size(); i++) {
                if (((*quant_vec_is_na)[n])[i]) {
                    quant->push_back(atof(((*quant_vec)[n])[i]));
                }
            }
        }
    };

    auto do_ion = [&]() {
        for (size_t n = 0; n < _line->size(); n++) {
            if ((*quant_vec_is_na_n)[n] == 0) {
                continue;
            }

            for (size_t k = 0; k < (*quant_vec)[n].size(); k++) {
                if (((*quant_vec_is_na)[n])[k]) {

                    string_vector_view iv((*_tab)[n][col_secondary_ids_ind[0]]);
                    for (size_t j = 1; j < n_ion_ind; j++) {
                        iv.str_vec.push_back((*_tab)[n][col_secondary_ids_ind[j]]);
                    }

                    if (has_quant_ids) {
                        iv.str_vec.push_back(((*quant_id_vec)[n])[k]);
                    }
                    else {
                        iv.str_vec.push_back(numbers[k]->c_str());
                    }

                    auto ii = map_ions->find(iv);
                    size_t i;

                    if (ii != map_ions->end()) {
                        i = ii->second;
                    }
                    else {
                        i = ions->size();

                        ions->emplace_back(new vector<string>);

                        for (size_t j = 0; j < n_ion_ind; j++) {
                            (*ions)[i]->emplace_back((*_tab)[n][col_secondary_ids_ind[j]]);
                        }
                        if (has_quant_ids) {
                            (*ions)[i]->emplace_back(((*quant_id_vec)[n])[k]);
                        }
                        else {
                            (*ions)[i]->emplace_back(numbers[k]->c_str());
                        }

                        for (size_t j = 0; j <= n_ion_ind; j++) { // 1 extra at the end. Needed because it now point to a
                            iv.str_vec[j] = (*ions)[i]->at(j).c_str();
                        }

                        map_ions->emplace(iv, i);
                    }

                    ion_index->push_back(i + 1);
                }
            }
        }
    };


    #ifdef _OPENMP
        omp_set_dynamic(0);
        int no_threads = 4;
        omp_set_num_threads(no_threads);
        if (verbose) {
            Rprintf("\nUsing %d threads ...\n", no_threads);
        }
    #else
        if (verbose) {
            Rprintf("Using a single CPU core ...\n");
        }
    #endif

    auto double_buffers = new double_buffering_file_t(fp);

    auto cleanup = [&]() {

        for (auto& a : numbers) {
            if (a) {
                delete a;
            }
        }

        delete map_proteins;
        delete map_ions;
        delete map_samples;

        delete _line;
        delete _line_new;
        delete _tab;
        fclose(fp);
        delete double_buffers;

        delete quant_vec;
        delete quant_id_vec;
        delete quant_vec_is_na;
        delete quant_vec_is_na_n;
        delete row_valid;
    };

    if (double_buffers->buf_0 == NULL || double_buffers->buf_1 == NULL) {
        cleanup();
        throw runtime_error("Cannot allocate memory.");
    }

    bool memory_ok = true;

    while (true) {

        swap(_line, _line_new);

        if (!_line->empty()) {

            if (_tab->size() < _line->size()) {
                _tab->resize(_line->size(), vector<char *>(n_cols, 0));
            }

            if (quant_vec->size() < _line->size()) {
                quant_vec->resize(_line->size());
            }

            if (quant_id_vec->size() < _line->size()) {
                quant_id_vec->resize(_line->size());
            }

            if (quant_vec_is_na->size() < _line->size()) {
                quant_vec_is_na->resize(_line->size());
            }

            if (quant_vec_is_na_n->size() < _line->size()) {
                quant_vec_is_na_n->resize(_line->size());
            }

            if (row_valid->size() < _line->size()) {
                row_valid->resize(_line->size());
            }

            size_t line_mismatched = SIZE_MAX;  // first bad line number, SIZE_MAX if none

            #pragma omp parallel for schedule(static) reduction(min:line_mismatched)
            for (size_t n = 0; n < _line->size(); n++) {

                (*quant_vec)[n].clear();
                (*quant_id_vec)[n].clear();
                (*quant_vec_is_na)[n].clear();
                (*quant_vec_is_na_n)[n] = 0;
                (*row_valid)[n] = 0;

                (*_tab)[n][0] = (*_line)[n];
                size_t pos = 0;
                char *c = (*_line)[n];
                while (*c) {
                    if (*c == '\r') {
                        *c++ = 0;
                    } else if (*c == '\t') {
                        *c++ = 0;
                        pos++;
                        if (pos < n_cols) {
                            (*_tab)[n][pos] = c;
                        }
                    } else {
                        c++;
                    }
                }

                // blank line, ignored
                if (pos == 0 && *(*_line)[n] == 0) {
                    continue;
                }

                // stop before touching columns the line does not have
                if ((pos + 1) != n_cols) {
                    line_mismatched = min(line_mismatched, line_no + n + 1);
                    continue;
                }

                (*row_valid)[n] = 1;  // its sample is registered, even if no value is kept

                size_t i = 0;
                while (i < filter_string_equal.size() && strcmp((*_tab)[n][col_filter_string_equal[i]], filter_string_equal[i].second) == 0) {
                    i++;
                }

                if (i >= filter_string_equal.size()) {

                    i = 0;
                    while (i < filter_string_not_equal.size() && strcmp((*_tab)[n][col_filter_string_not_equal[i]], filter_string_not_equal[i].second) != 0) {
                        i++;
                    }

                    if (i >= filter_string_not_equal.size()) {

                        bool ok = true;
                        i = 0;
                        while (i < filter_double_less.size() && ok) {
                            double v;
                            if (filter_number((*_tab)[n][col_filter_double_less[i]], &v) && v >= filter_double_less[i].second) {
                                ok = false;
                            }
                            i++;
                        }

                        if (ok) {

                            i = 0;
                            while (i < filter_double_greater.size() && ok) {
                                double v;
                                if (filter_number((*_tab)[n][col_filter_double_greater[i]], &v) && v <= filter_double_greater[i].second) {
                                    ok = false;
                                }
                                i++;
                            }

                            if (ok) {

                                if (intensity_col_sep) {
                                    // a single column holding several values
                                    utils::split((*_tab)[n][col_quants_ind[0]], intensity_col_sep, &((*quant_vec)[n]));
                                    if (has_quant_ids) {
                                        utils::split((*_tab)[n][col_quant_ids_ind[0]], intensity_col_sep, &((*quant_id_vec)[n]));
                                        if ((*quant_vec)[n].size() != (*quant_id_vec)[n].size()) {
                                            line_mismatched = min(line_mismatched, line_no + n + 1);
                                        }
                                    }
                                }
                                else {
                                    // one value per column
                                    for (size_t k = 0; k < col_quants_ind.size(); k++) {
                                        (*quant_vec)[n].push_back((*_tab)[n][col_quants_ind[k]]);
                                        if (has_quant_ids) {
                                            (*quant_id_vec)[n].push_back((*_tab)[n][col_quant_ids_ind[k]]);
                                        }
                                    }
                                }

                                for (auto& s : (*quant_vec)[n]) {

                                    if (*s != 0 && strcmp(s, "NaN") != 0 && (!na_string || strcmp(s, na_string) != 0)) {
                                        (*quant_vec_is_na)[n].push_back(true);
                                        (*quant_vec_is_na_n)[n]++;
                                    }
                                    else {
                                        (*quant_vec_is_na)[n].push_back(false);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            for (size_t n = 0; n < _line->size(); n++) {
                if (numbers.size() < (*quant_vec)[n].size()) {
                    for (auto k = numbers.size(); k < (*quant_vec)[n].size(); k++) {
                        numbers.push_back(new string(to_string(k + 1)));
                    }
                }
            }

            if (line_mismatched != SIZE_MAX) { // no exception in parallel regions
                cleanup();
                throw runtime_error(("Line " + to_string(line_mismatched) + ": the number of columns does not match header or quant values and ids.\n").c_str());
            }

            line_no += _line->size();
            for (size_t n = 0; n < _line->size(); n++) {
                n_after_filtered += (*quant_vec_is_na_n)[n];
            }
        }

        if (double_buffers->eof) {  // wrapping up the last round
            do_ion();
            do_sample();
            do_protein();
            if (on_chunk) {
                try {
                    on_chunk();
                } catch (...) {
                    cleanup();
                    throw;
                }
            }
            break;
        }

        #pragma omp parallel sections
        {
            #pragma omp section
            {
                if (!_line->empty()) {
                    do_ion();
                }
            }
            #pragma omp section
            {
                if (!_line->empty()) {
                    do_sample();
                }
            }
            #pragma omp section
            {
                if (!_line->empty()) {
                    do_protein();
                }
            }

            #pragma omp section
            {
                if (double_buffers->fgetss(*_line_new) != 0) {
                    memory_ok = false;
                }
            }
        }

        if (!memory_ok) {
            cleanup();
            throw runtime_error("Cannot allocate memory.");
        }

        if (on_chunk && !_line->empty()) {
            try {
                on_chunk();
            } catch (...) {
                cleanup();
                throw;
            }
        }

        if (quant->size() > (size_t)IQ_MAX_VALUES) {  // values kept in memory (none with on_chunk)
            cleanup();
            throw runtime_error(too_many_values_message);
        }

        if (verbose && samples->size() >= print_threshold) {
            Rprintf("%d samples read\n", (int)samples->size());
            print_threshold = samples->size() + 20;
        }
    }
    if (verbose) {
        Rprintf("\n# lines read (excluding headers)      = %.0f\n", (double)(line_no-1));
        Rprintf("# quantitative values after filtering = %.0f\n\n", (double)n_after_filtered);

        Rprintf("# samples  = %d\n", (int)samples->size());
        Rprintf("# proteins = %d\n", (int)annotation->size());
    }

    cleanup();
}

// writer: NULL to return all values (fast_read), or an iq format writer that receives the
// values chunk by chunk (long_format_to_iq_format); the result then has no values, only the
// lists of proteins, samples and ions.
// The error of the text reader, kept outside the C++ scopes: iq_filter_impl raises it with
// Rf_error (a longjmp, which skips destructors) only after iq_filter_body has returned.
static char filter_error[2048];

static SEXP iq_filter_body(SEXP cmd, iq_writer* writer, bool* failed) {

    const char *c = CHAR(STRING_ELT(cmd, 0));

    // the arguments, separated by \001 (.fast_read_command in R); an empty argument, e.g. an
    // empty filter value, is kept: two separators in a row give an empty argument between them
    vector<string> argv;

    const char sep = 1;

    const char *p = c;
    if (*p) {
        while (true) {
            const char *start = p;
            while (*p && *p != sep) {
                p++;
            }
            argv.emplace_back(start, p - start);
            if (*p == 0) {
                break;
            }
            p++;  // past the separator
        }
    }

    auto protein_index = new vector<int>;
    auto sample_index = new vector<int>;
    auto ion_index = new vector<int>;
    auto quant = new vector<double>;

    auto annotation = new vector<vector<string>*>;
    auto col_annotation = new vector<string>;
    auto samples = new vector<string*>;
    auto ions = new vector<vector<string>*>;

    // values without a protein or sample identifier are left out, as in process()
    std::function<void()> on_chunk;
    if (writer) {
        on_chunk = [&]() {
            for (size_t k = 0; k < quant->size(); k++) {
                int p = (*protein_index)[k];
                int s = (*sample_index)[k];
                if ((*annotation)[p - 1]->at(0).empty() || (*samples)[s - 1]->empty()) {
                    writer->n_unidentified++;
                    continue;
                }
                writer->add(p, s, (*ion_index)[k], (*quant)[k]);
            }
            protein_index->clear();
            sample_index->clear();
            ion_index->clear();
            quant->clear();
        };
    }

    try {
        process(argv, protein_index, sample_index, ion_index, quant, annotation, col_annotation, samples, ions, on_chunk);
    } catch (exception & e) {
        snprintf(filter_error, sizeof(filter_error), "%s", e.what());
        *failed = true;

        for (auto& a : (*ions)) {
            if (a) {
                delete a;
            }
        }
        delete ions;

        for (auto& a : (*samples)) {
            if (a) {
                delete a;
            }
        }
        delete samples;

        for (auto& a : (*annotation)) {
            if (a) {
                delete a;
            }
        }
        delete annotation;

        delete col_annotation;
        delete quant;
        delete ion_index;
        delete sample_index;
        delete protein_index;
        return (R_NilValue);
    }

    // the R vectors below have an int length
    const bool too_many = protein_index->size() > (size_t)IQ_MAX_VALUES;
    if (too_many) {
        snprintf(filter_error, sizeof(filter_error), "%s", too_many_values_message);
        *failed = true;
    } else if (ions->size() == 0) {
        snprintf(filter_error, sizeof(filter_error), "No quantitative values left after filtering.");
        *failed = true;
    }

    if (*failed) {

        for (auto& a : (*ions)) {
            if (a) {
                delete a;
            }
        }
        delete ions;

        for (auto& a : (*samples)) {
            if (a) {
                delete a;
            }
        }
        delete samples;

        for (auto& a : (*annotation)) {
            if (a) {
                delete a;
            }
        }
        delete annotation;

        delete col_annotation;
        delete quant;
        delete ion_index;
        delete sample_index;
        delete protein_index;
        return (R_NilValue);
    }

    SEXP p_index = PROTECT(Rf_allocMatrix(INTSXP, protein_index->size(), 1));
    copy(protein_index->begin(), protein_index->end(), INTEGER(p_index));

    SEXP s_index = PROTECT(Rf_allocMatrix(INTSXP, sample_index->size(), 1));
    copy(sample_index->begin(), sample_index->end(), INTEGER(s_index));

    SEXP i_index = PROTECT(Rf_allocMatrix(INTSXP, ion_index->size(), 1));
    copy(ion_index->begin(), ion_index->end(), INTEGER(i_index));

    SEXP q = PROTECT(Rf_allocMatrix(REALSXP, quant->size(), 1));
    copy(quant->begin(), quant->end(), REAL(q));

    SEXP ion_str = PROTECT(Rf_allocMatrix(STRSXP, ions->size(), 1));
    for (size_t i = 0; i < ions->size(); i++) {
        // concatenate
        SET_STRING_ELT(ion_str, i, Rf_mkChar(utils::concatenate(*(ions->at(i))).c_str()));
    }

    // the secondary_id part of each ion, without the suffix of the intensity value (1, 2, ... or its id)
    SEXP ion_prefix_str = PROTECT(Rf_allocMatrix(STRSXP, ions->size(), 1));
    for (size_t i = 0; i < ions->size(); i++) {
        const vector<string>& parts = *(ions->at(i));
        string prefix = parts[0];
        for (size_t j = 1; j + 1 < parts.size(); j++) {
            prefix += "_" + parts[j];
        }
        SET_STRING_ELT(ion_prefix_str, i, Rf_mkChar(prefix.c_str()));
    }

    SEXP sample_str = PROTECT(Rf_allocMatrix(STRSXP, samples->size(), 1));
    for (size_t i = 0; i < samples->size(); i++) {
        SET_STRING_ELT(sample_str, i, Rf_mkChar(samples->at(i)->c_str()));
    }

    SEXP ann = PROTECT(Rf_allocMatrix(STRSXP, annotation->size(), annotation->at(0)->size()));
    for (size_t i = 0; i < annotation->size(); i++) {
        for (size_t j = 0; j < annotation->at(i)->size(); j++) {
            SET_STRING_ELT(ann, i + j * annotation->size(), Rf_mkChar(annotation->at(i)->at(j).c_str()));
        }
    }

    // name for annotation tables
    SEXP col_names = PROTECT(Rf_allocVector(STRSXP, col_annotation->size()));
    for (size_t i = 0; i < col_annotation->size(); i++) {
        SET_STRING_ELT(col_names, i, Rf_mkChar(col_annotation->at(i).c_str()));
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

    for (auto& a : (*ions)) {
        if (a) {
            delete a;
        }
    }
    delete ions;

    for (auto& a : (*samples)) {
        if (a) {
            delete a;
        }
    }
    delete samples;

    for (auto& a : (*annotation)) {
        if (a) {
            delete a;
        }
    }
    delete annotation;

    delete col_annotation;
    delete quant;
    delete ion_index;
    delete sample_index;
    delete protein_index;

    return vec;
}

static SEXP iq_filter_impl(SEXP cmd, iq_writer* writer) {
    bool failed = false;
    SEXP res = iq_filter_body(cmd, writer, &failed);
    if (failed) {
        Rf_error("%s", filter_error);
    }
    return res;
}

SEXP iq_filter(SEXP cmd) {
    return iq_filter_impl(cmd, NULL);
}

SEXP iq_filter_write(SEXP cmd, SEXP writer) {
    iq_writer* w = NULL;
    static char message[1024];
    bool failed = false;
    try {
        w = get_writer(writer);
    } catch (exception& e) {
        snprintf(message, sizeof(message), "%s", e.what());
        failed = true;
    }
    if (failed) {
        Rf_error("%s", message);
    }
    return iq_filter_impl(cmd, w);
}

//------------------------------ MaxLFQ -------------------------------------

#ifdef IQ_USE_EIGEN

using Eigen::FullPivHouseholderQR;
using Eigen::HouseholderQR;
using Eigen::MatrixXd;
using Eigen::VectorXd;

class ion_table {

    static int full_connection;
    static FullPivHouseholderQR<MatrixXd> full_qr;

public:
    static void init(int N) {
        MatrixXd AtA = MatrixXd::Zero(N + 1, N + 1);

        for (int j = 0; j < (N - 1); j++) {
            for (int k = j + 1; k < N; k++) {
                AtA(j, k) = -1.0;
                AtA(k, j) = -1.0;
                AtA(j, j) += 1.0;
                AtA(k, k) += 1.0;
            }
        }

        AtA *= 2.0;

        for (int j = 0; j < N; j++) {
            AtA(j, N) = 1.0;
            AtA(N, j) = 1.0;
        }
        AtA(N, N) = 0.0;

        full_qr = FullPivHouseholderQR<MatrixXd>(AtA);

        full_connection = N * (N - 1) / 2;
    }

    unordered_map<int, vector<double>> map;
    int ncol;

    ion_table(int ncol) : ncol(ncol) {
    }

    // add an entry (log2) to a map, increase size if necessary; more values of the same ion and
    // sample: log2 of the summed intensities, as in iq_process
    void add_to_table(int ion, int col, double quant) {
        auto r = map.emplace(ion, vector<double>());
        if (r.second) {
            //n++;
            r.first->second.resize(ncol, NAN);
        }
        double& cell = r.first->second[col];
        cell = std::isnan(cell) ? quant : log2(pow(2.0, cell) + pow(2.0, quant));
    }

    void spread(int *g, int i, int val) {
        g[i] = val;

        for (const auto &m : map) {
            if (!isnan(m.second[i])) {
                for (int j = 0; j < ncol; j++) {
                    if (g[j] < 0 && !isnan(m.second[j])) {
                        spread(g, j, val);
                    }
                }
            }
        }
    }

    inline double get_median(double *median_buffer, int median_size) {
        sort(median_buffer, median_buffer + median_size);
        int mid = median_size / 2;
        return median_size % 2 == 0 ? (median_buffer[mid] + median_buffer[mid - 1]) / 2.0 : median_buffer[mid];
    }

    void maxLFQ(double *buffer, int *g) {
        if (map.size() < 1) {
            fill_n(buffer, ncol, NAN);
            fill_n(g, ncol, 0);
            return;
        }

        if (map.size() == 1) {
            for (int i = 0; i < ncol; i++) {
                buffer[i] = map.begin()->second[i];
            }
            fill_n(g, ncol, 0);
            return;
        }

        double *median_buffer = new double[map.size()];
        int median_size = 0;

        fill_n(g, ncol, -1);

        int val = 0;
        for (int i = 0; i < ncol; i++) {
            if (g[i] < 0) {
                spread(g, i, val++);
            }
        }

        fill_n(buffer, ncol, NAN);

        for (int i = 0; i < val; i++) {
            vector<int> ind;
            for (int j = 0; j < ncol; j++) {
                if (g[j] == i) {
                    ind.push_back(j);
                }
            }

            if (ind.size() == 1) {
                median_size = 0;

                for (const auto &m : map) {
                    if (!isnan(m.second[ind[0]])) {
                        median_buffer[median_size++] = m.second[ind[0]];
                    }
                }
                if (median_size == 0) {
                    buffer[ind[0]] = NAN;
                } else {
                    buffer[ind[0]] = get_median(median_buffer, median_size);
                }
            } else {

                int N = ind.size();

                MatrixXd AtA = MatrixXd::Zero(N + 1, N + 1);
                VectorXd Atb = VectorXd::Zero(N + 1);

                int n_connection = 0;

                for (int j = 0; j < (N - 1); j++) {
                    for (int k = j + 1; k < N; k++) {
                        median_size = 0;
                        for (const auto &m : map) {
                            if (!isnan(m.second[ind[j]]) && !isnan(m.second[ind[k]])) {
                                median_buffer[median_size++] = m.second[ind[k]] - m.second[ind[j]];
                            }
                        }
                        if (median_size > 0) {
                            n_connection++;

                            double r_i_j = get_median(median_buffer, median_size);
                            AtA(j, k) = -1.0;
                            AtA(k, j) = -1.0;

                            AtA(j, j) += 1.0;
                            AtA(k, k) += 1.0;

                            Atb(j) -= r_i_j;
                            Atb(k) += r_i_j;
                        }
                    }
                }

                AtA *= 2.0;
                Atb *= 2.0;

                for (int j = 0; j < N; j++) {
                    AtA(j, N) = 1.0;
                    AtA(N, j) = 1.0;
                }
                AtA(N, N) = 0.0;

                // mean data
                double sum = 0.0;
                int count = 0;
                for (const auto &m : map) {
                    for (const auto &s : ind) {
                        if (!isnan(m.second[s])) {
                            sum += m.second[s];
                            count++;
                        }
                    }
                }
                Atb(N) = sum * (double)N / (double)count;

                if (n_connection == full_connection) {
                    VectorXd x = full_qr.solve(Atb);
                    for (int j = 0; j < N; j++) {
                        buffer[ind[j]] = x(j);
                    }
                } else {
                    FullPivHouseholderQR<MatrixXd> qr(AtA);
                    VectorXd x = qr.solve(Atb);
                    for (int j = 0; j < N; j++) {
                        buffer[ind[j]] = x(j);
                    }
                }
            }
        }

        delete[] median_buffer;
    }
};

int ion_table::full_connection;
FullPivHouseholderQR<MatrixXd> ion_table::full_qr;

#else  // IQ_USE_EIGEN

// The ions of one protein (fast_MaxLFQ), quantified with the solver of maxlfq_bit.
class ion_table {

    unordered_map<int, int> rows;  // ion -> row of the matrix
    vector<int> row_of, col_of;    // each value added: its row and column
    vector<double> quant_of;       // and the value
    int ncol;

public:
    static void init(int) {
    }

    ion_table(int ncol) : ncol(ncol) {
    }

    // an ion with a value (log2) in a sample; more values of the same ion and sample are summed
    // in linear space when the matrix is built (maxLFQ below)
    void add_to_table(int ion, int col, double quant) {
        row_of.push_back(rows.emplace(ion, (int)rows.size()).first->second);
        col_of.push_back(col);
        quant_of.push_back(quant);
    }

    // estimates (NaN without data) and the connected component of each sample (g)
    void maxLFQ(double *buffer, int *g) {
        int M = (int)rows.size();
        if (M < 1) {
            fill_n(buffer, ncol, NAN);
            fill_n(g, ncol, 0);
            return;
        }

        // ions x samples, column-major, NaN = missing; an ion with more than one value in a
        // sample gets log2 of the summed intensities, as in iq_process
        vector<double> X((size_t)M * ncol, NAN);
        vector<int> count((size_t)M * ncol, 0);
        for (size_t k = 0; k < quant_of.size(); k++) {
            size_t cell = (size_t)col_of[k] * M + row_of[k];
            if (count[cell] == 0) {
                X[cell] = quant_of[k];
            } else if (count[cell] == 1) {
                X[cell] = pow(2.0, X[cell]) + pow(2.0, quant_of[k]);  // linear sum
            } else {
                X[cell] += pow(2.0, quant_of[k]);
            }
            count[cell]++;
        }
        for (size_t c = 0; c < X.size(); c++) {
            if (count[c] > 1) {
                X[c] = log2(X[c]);
            }
        }

        maxlfq_bit(X.data(), M, ncol, buffer, 1, 0);

        if (M == 1) {
            fill_n(g, ncol, 0);
        } else {
            int n_components;
            connected_components(X.data(), M, ncol, g, n_components);
        }
    }
};

#endif  // IQ_USE_EIGEN

SEXP get_list_element(SEXP list, const char *str) {
    SEXP names = Rf_getAttrib(list, R_NamesSymbol);

    for (int i = 0; i < Rf_length(list); i++) {
        if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
            return (VECTOR_ELT(list, i));
        }
    }

    throw runtime_error((string("Cannot find list element: ") + str).c_str());
}


static void tp_user(void *dummy) {
    R_CheckUserInterrupt();
}

// call the above in a top-level context
int tp_check() {
    return (R_ToplevelExec(tp_user, NULL) == FALSE);
}

/******************
iq format 
*******************/

double dist(double *a, int step, double *b, int n) {
    
    double sum = 0;
    int t = 0;
    
    double *ap = a;
    for (int i = 0; i < n; i++) {
        if (!isnan(*ap) && !isnan(b[i])) {
            t++;
            sum += (*ap - b[i]);
        }
        ap+= step;
    }
    
    double m = sum / t;
    
    double d2 = 0;
    
    ap = a;
    for (int i = 0; i < n; i++) {
        if (!isnan(*ap) && !isnan(b[i])) {
            double d =  (*ap - b[i] - m);
            d2 += d*d;
        }
        ap+= step;
    }
    
    if (t == 0) {
        return NAN;
    }
    
    if (t == 1) {
        //std::cout << "d2 = " << d2 << endl;
        return NAN;
    }
    
    /* NOTE: sqrt is correct. Eq. 4 in the paper misses a SQRT*/

    return sqrt((double)n / t * (d2));
}

void missing_value_weight(double* X, int* M, int* N,
                          double* w) {
    
    for (int i = 0; i < *M; i++) {
        
        double* pX = X + i;
        int cc = *N;
        
        for(int j = 0; j < *N; j++) {
            if (isnan(*pX)) {
                cc--;
            }
            pX += *M;
        }
        
        w[i] = (double)cc / (double)*N;
        
        /*
         if (cc == 0) {
         w[i] = 0;
         }
         else {
         w[i] = (double)cc / (double)*N;
         }
         */
    }
}

void huber_weight(double* X, int* M, int* N,
                  double *b,
                  double* w,                  //w pre-allocated
                  double k = 1.4826,
                  double c = 1.345) {                
    
    // absolute deviation
    double* abs_deviation = new double[*M];
    int mm = 0;
    
    for (int i = 0; i < *M; i++) {
        w[i] = dist(X+i, *M, b, *N); // residual
        
        if (!isnan(w[i])) {
            abs_deviation[mm++] = w[i];
        }
    }
    
    // mad of residuals
    if (mm == 0) {
        for (int i = 0; i < *M; i++) {
            w[i] = 1.0;
        }
        delete [] abs_deviation;
        return;
    }
    
    double m = quick_median(abs_deviation, mm);
    
    int mm2 = 0;
    for (int i = 0; i < *M; i++) {
        if (!isnan(w[i])) {
            abs_deviation[mm2++] = abs(w[i] - m);
        }
    }
    
    // mad
    m = quick_median(abs_deviation, mm2) * k;
    
    if (m < 1e-10) {
        for (int i = 0; i < *M; i++) {
            w[i] = 1.0;
        }
        delete [] abs_deviation;
        return;
    }
    
    //cout << "m=" << m << endl;
    
    delete [] abs_deviation;
    
    // huber psi
    for (int i = 0; i < *M; i++) {
        if (!isnan(w[i])) {
            w[i] = abs(w[i] / m); // rescale
            
            if (w[i] <= c) {
                w[i] = 1;
            }
            else {
                w[i] = c / w[i];
            }
        }
        else {
            w[i] = 1;
        }
    }
}

SEXP iq_MaxLFQ(SEXP list) {

    int stop_sig = 0;

    int* proteins;
    int* ions;
    int* samples;
    double* quants;
    bool verbose;
    R_xlen_t n_p = 0, n_i = 0, n_s = 0, n_q = 0;

    try {
        proteins = INTEGER(get_list_element(list, "protein_index"));
        ions = INTEGER(get_list_element(list, "ion_index"));
        samples = INTEGER(get_list_element(list, "sample_index"));
        quants = REAL(get_list_element(list, "quant"));
        verbose = Rf_asLogical(get_list_element(list, "verbose")) == TRUE;
        n_p = Rf_xlength(get_list_element(list, "protein_index"));
        n_i = Rf_xlength(get_list_element(list, "ion_index"));
        n_s = Rf_xlength(get_list_element(list, "sample_index"));
        n_q = Rf_xlength(get_list_element(list, "quant"));
    }
    catch (exception & e) {
        Rprintf("%s\n", e.what());
        return (R_NilValue);
    }

    // also checked by fast_MaxLFQ in R: no values, vectors of different lengths, or an index
    // below 1 (or NA) would read or write outside the arrays below
    if (n_p == 0) {
        Rf_error("iq_MaxLFQ: no values");
    }
    if (n_i != n_p || n_s != n_p || n_q != n_p) {
        Rf_error("iq_MaxLFQ: protein_index, ion_index, sample_index and quant must have the same length");
    }
    for (R_xlen_t r = 0; r < n_p; r++) {
        if (proteins[r] < 1 || samples[r] < 1 || ions[r] == NA_INTEGER) {
            Rf_error("iq_MaxLFQ: protein_index and sample_index must be 1 or more, ion_index not NA");
        }
    }

    size_t nrow = (size_t)n_p;

    int n_proteins = utils::count_unique(proteins, nrow);
    int n_samples = utils::count_unique(samples, nrow);

    SEXP table = PROTECT(Rf_allocMatrix(REALSXP, n_proteins, n_samples));
    double* buffer = REAL(table);

    auto group_annotation = new vector<string>(n_proteins, "");
    int* row_names = new int[n_proteins];
    int* col_names = new int[n_samples];

    if (verbose) {
        Rprintf("nrow = %d, # proteins = %d, # samples = %d\n", (int)nrow, (int)n_proteins, (int)n_samples);
    }

    //auto protein_index = new vector< vector<int> >(n_proteins);
    auto protein_index = new vector<vector<int>>(*(std::max_element(proteins, proteins + nrow)));  // allowing for missing proteins

    for (auto i : *protein_index) {
        i.clear();
    }

    vector<int> sample_set;  // sample index, in case a subset of samples is quantified

    for (size_t i = 0; i < nrow; i++) {
        (*protein_index)[proteins[i] - 1].push_back(i);

        if (samples[i] >= (int)sample_set.size()) {
            sample_set.resize(samples[i] + 1, -1);
            sample_set[samples[i]] = 0;
        }
        else {
            sample_set[samples[i]] = 0;
        }
    }

    int cc = 0;
    for (int i = 0; i < (int)sample_set.size(); i++) {
        if (sample_set[i] > -1) {
            col_names[cc] = i;
            sample_set[i] = cc++;
        }
    }

    ion_table::init(n_samples);  // QR for full matrix

    int nr = 0;
    vector<int> map_back((*protein_index).size(), -1);
    for (size_t i = 0; i < (*protein_index).size(); i++) {
        if (!(*protein_index)[i].empty()) {
            map_back[i] = nr;
            row_names[nr++] = i + 1;
        }
    }

    size_t thres_display = 0;

    #ifdef _OPENMP

        omp_set_dynamic(0);

        int no_threads = omp_get_num_procs() - 1;

        if (no_threads < 1) {
            no_threads = 1;
        }

        omp_set_num_threads(no_threads);

        // the solver of a protein runs single-threaded inside the loop over proteins
        int saved_levels = omp_get_max_active_levels();
        omp_set_max_active_levels(1);

        if (verbose) {
            Rprintf("Using %d threads...\n", no_threads);
        }

    #else

        if (verbose) {
            Rprintf("Using a single CPU core...\n");
        }

    #endif

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < (*protein_index).size(); i++) {
        if (stop_sig) {
            continue;
        }

        int thread_id = 0;

        #ifdef _OPENMP
            thread_id = omp_get_thread_num();
        #endif

        if (!(*protein_index)[i].empty()) {

            ion_table *tab = new ion_table(n_samples);

            double *w = new double[n_samples];
            int *gr = new int[n_samples];

            for (auto j : (*protein_index)[i]) {
                tab->add_to_table(ions[j], sample_set[samples[j]], quants[j]);
            }

            tab->maxLFQ(w, gr);

            //--- annotation
            bool single_compoment = true;
            int component = -1;

            for (int j = 0; j < n_samples; j++) {
                if (!isnan(w[j])) {
                    if (component == -1) {
                        component = gr[j];
                    } else {
                        if (component != gr[j]) {
                            single_compoment = false;
                            break;
                        }
                    }
                }
            }

            if (component == -1) {
                (*group_annotation)[map_back[i]] = "NA";
            } else {
                if (!single_compoment) {  // all NOT in the same component, otherwise the default empty string is good
                    for (int j = 0; j < n_samples; j++) {
                        if (j > 0) {
                            (*group_annotation)[map_back[i]].push_back(';');
                        }
                        if (isnan(w[j])) {
                            (*group_annotation)[map_back[i]].append("NA");
                        } else {
                            (*group_annotation)[map_back[i]].append(to_string(gr[j] + 1));
                        }
                    }
                }
            }

            //--- filling buffer
            for (int j = 0; j < n_samples; j++) {  // byrow = TRUE
                size_t k = j * n_proteins + map_back[i];
                if (isnan(w[j])) {
                    buffer[k] = NA_REAL;
                } else {
                    buffer[k] = w[j];
                }
            }

            delete[] w;
            delete[] gr;
            delete tab;
        }

        if (thread_id == 0) {
            if (verbose && i > thres_display) {
                Rprintf("%d%%\n", (int)(i * 100 / (*protein_index).size()));
                R_FlushConsole();
                thres_display = i + (*protein_index).size() / 20;
            }

            if (tp_check()) {  // user interrupted ...
                stop_sig = 1;
                #pragma omp flush(stop_sig)
            }
        }
    }

    #ifdef _OPENMP
        omp_set_max_active_levels(saved_levels);
    #endif

    if (stop_sig) {
        Rprintf("Cancelled.\n");

        UNPROTECT(1);

        delete protein_index;
        delete group_annotation;
        delete[] col_names;
        delete[] row_names;

        return (R_NilValue);
    }

    // estimate names
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, n_proteins));
    for (int i = 0; i < n_proteins; i++) {
        SET_STRING_ELT(r_names, i, Rf_mkChar(to_string(row_names[i]).c_str()));
    }

    SEXP c_names = PROTECT(Rf_allocVector(STRSXP, n_samples));
    for (int i = 0; i < n_samples; i++) {
        SET_STRING_ELT(c_names, i, Rf_mkChar(to_string(col_names[i]).c_str()));
    }

    //SEXP dimnames = getAttrib(ann, R_DimNamesSymbol);
    SEXP dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dimnames, 0, r_names);
    SET_VECTOR_ELT(dimnames, 1, c_names);
    Rf_setAttrib(table, R_DimNamesSymbol, dimnames);

    // annotation
    SEXP ann = PROTECT(Rf_allocVector(STRSXP, n_proteins));
    for (int i = 0; i < n_proteins; i++) {
        SET_STRING_ELT(ann, i, Rf_mkChar(group_annotation->at(i).c_str()));
    }

    // return list
    SEXP vec = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(vec, 0, table);
    SET_VECTOR_ELT(vec, 1, ann);

    SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, Rf_mkChar("estimate"));
    SET_STRING_ELT(names, 1, Rf_mkChar("annotation"));

    Rf_setAttrib(vec, R_NamesSymbol, names);

    if (verbose) {
        Rprintf("Completed.\n");
    }

    UNPROTECT(7);

    delete protein_index;
    delete group_annotation;
    delete[] col_names;
    delete[] row_names;

    return (vec);
}

//------------------------------ quantification of one protein ---------------------

// p1 / p2 from R: an empty vector (a method without the parameter) gives a default of 0
// instead of R's pointer for an empty vector, which must not be dereferenced
static double p1_default = 0.0;
static int p2_default = 0;

static double* p1_pointer(SEXP p1) {
    return Rf_xlength(p1) > 0 ? REAL(p1) : &p1_default;
}

static int* p2_pointer(SEXP p2) {
    return Rf_xlength(p2) > 0 ? INTEGER(p2) : &p2_default;
}

// The method of iq_quant (process_matrix) for one protein.
// Direct summaries of a protein, per sample, on the scale of the data (X: M ions x N samples,
// column-major, log2, NaN = missing; x: N values, NaN for a sample without data):
//   method 7 "sum": log2 of the summed intensities
//   method 8 "mean": the mean log2 intensity
//   method 9 "topN": the mean log2 intensity of the top_n ions with the highest mean log2
//     intensity over the samples (equal means: in input order), as topN() in R
static void summarize_matrix(const double* X, int M, int N, int method, int top_n, double* x) {

    vector<int> rows(M);
    for (int i = 0; i < M; i++) {
        rows[i] = i;
    }

    if (method == 9) {
        vector<double> row_mean(M);
        vector<double> v;
        for (int i = 0; i < M; i++) {
            v.clear();
            for (int j = 0; j < N; j++) {
                if (!isnan(X[i + (size_t)j * M])) {
                    v.push_back(X[i + (size_t)j * M]);
                }
            }
            row_mean[i] = mean_value(v.data(), v.size());  // NaN for an ion without values
        }
        std::stable_sort(rows.begin(), rows.end(), [&](int a, int b) {
            if (isnan(row_mean[b])) {
                return !isnan(row_mean[a]);  // ions without values last
            }
            return !isnan(row_mean[a]) && row_mean[a] > row_mean[b];
        });
        rows.resize(top_n < M ? (top_n > 0 ? top_n : 0) : M);
    }

    vector<double> v;
    for (int j = 0; j < N; j++) {
        v.clear();
        for (int i : rows) {
            if (!isnan(X[i + (size_t)j * M])) {
                v.push_back(X[i + (size_t)j * M]);
            }
        }
        x[j] = (method == 7) ? log2_sum(v.data(), v.size()) : mean_value(v.data(), v.size());
    }
}

// X: M ions x N samples, column-major, NaN = missing.
// x: N sample estimates; s: M ion offsets (NaN unless median polish or rlm).
// After the method, x is shifted so that its mean equals the mean of the data (not for the
// direct summaries, methods 7 to 9).
// p1 / p2 point to one value each (a default of 0 for an empty vector from R, see p1_pointer).
// Returns false for an unknown method (x and s are then NaN).
// Used in a parallel loop by iq_process with no_threads = 1 (nested regions inactive).
static bool quantify_matrix(double* X, int M, int N, int method, double* p1, int* p2,
                            double k, int min_M, int no_threads, double* x, double* s) {

    fill_n(s, M, NAN);  // only median polish and rlm compute ion offsets

    if (method >= 7 && method <= 9) {
        summarize_matrix(X, M, N, method, *p2, x);
        return true;
    }

    // mean data
    int no_t = 0;
    double sum = 0.0;
    for (int i = 0; i < M * N; i++) {
        if (!isnan(X[i])) {
            no_t++;
            sum += X[i];
        }
    }

    if (no_t == 0) {  // all NA
        fill_n(x, N, NAN);
        return true;
    }

    bool known = true;

    switch (method) {

        case 0:  // maxlfq
        {
            #ifdef IQ_USE_EIGEN
            std::unordered_map<int, std::vector<double>> map;
            for (int i = 0; i < M; i++) {
                map[i].resize(N);
                for (int j = 0; j < N; j++) {
                    map[i][j] = X[j * M + i];
                }
            }
            vector<int> gr(N);
            maxlfq(map, x, gr.data());
            #else
            maxlfq_bit(X, M, N, x, no_threads, *p2);  // p2: memory level, as "maxlfq_bit"
            #endif
            break;
        }
        case 1:  // maxlfq_bit
        {
            maxlfq_bit(X, M, N, x, no_threads, *p2);
            break;
        }
        case 2:  // weighted maxlfq
        {
            #ifdef IQ_USE_EIGEN
            std::unordered_map<int, std::vector<double>> map;
            for (int i = 0; i < M; i++) {
                map[i].resize(N);
                for (int j = 0; j < N; j++) {
                    map[i][j] = X[j * M + i];
                }
            }
            vector<int> gr(N);
            maxlfq(map, x, gr.data());

            if (M >= min_M) {
                vector<double> w(M);
                huber_weight(X, &M, &N, x, w.data(), 1.4826, k);
                weighted_maxlfq(map, x, gr.data(), w.data());
            }
            #else
            maxlfq_bit(X, M, N, x, no_threads, *p2);  // p2: memory level, as "maxlfq_bit"

            if (M >= min_M) {
                vector<double> w(M);
                huber_weight(X, &M, &N, x, w.data(), 1.4826, k);
                weighted_maxlfq_bit(X, M, N, x, w.data(), no_threads, *p2);
            }
            #endif
            break;
        }
        case 3:  // median polish
        {
            median_polish(X, M, N, s, x, no_threads, *p1, *p2);
            break;
        }
        case 4:  // weighted median polish
        {
            if (M == 1) {
                for (int j = 0; j < N; j++) {
                    x[j] = X[j];
                }
            } else {
                median_polish(X, M, N, s, x, no_threads, *p1, *p2);

                if (M >= min_M) {
                    vector<double> w(M);
                    huber_weight(X, &M, &N, x, w.data(), 1.4826, k);
                    weighted_median_polish(X, M, N, s, x, w.data(), *p1, *p2);
                }
            }
            break;
        }
        case 5:  // huber
        {
            huber(X, M, N, s, x, no_threads, *p1, *p2);
            break;
        }
        case 6:  // weighted huber
        {
            huber(X, M, N, s, x, no_threads, *p1, *p2);

            if (M >= min_M) {
                vector<double> w(M);
                huber_weight(X, &M, &N, x, w.data(), 1.4826, k);
                weighted_huber(X, M, N, s, x, w.data(), no_threads, *p1, *p2);
            }
            break;
        }
        default:
            fill_n(x, N, NAN);
            fill_n(s, M, NAN);
            known = false;
            break;
    }

    int no_x = 0;
    double sum_x = 0.0;
    for (int j = 0; j < N; j++) {
        if (!isnan(x[j])) {
            no_x++;
            sum_x += x[j];
        }
    }

    if (no_x > 0) {
        // mean data
        double mean_data = sum * (double)no_x / (double)no_t;

        // scaling
        double delta = (mean_data - sum_x) / no_x;
        for (int j = 0; j < N; j++) {
            if (!isnan(x[j])) {
                x[j] += delta;
            }
        }
    }

    return known;
}

// number of threads from n_threads: > 0 at most that many, <= 0 all cores + n_threads
static int thread_count(int n_threads) {
    int no_threads = 1;
    #ifdef _OPENMP
    int n_procs = omp_get_num_procs();
    if (n_threads > 0) {
        no_threads = n_procs >= n_threads ? n_threads : n_procs;
    } else {
        no_threads = n_procs + n_threads;
    }
    if (no_threads < 1) {
        no_threads = 1;
    }
    #else
    (void)n_threads;
    #endif
    return no_threads;
}

extern "C" {
    
SEXP iq_quant(SEXP list) {
        
    double* X;
    int* M;
    int* N;
    int* method;
    double* p1;
    int* p2;
    double* k;
    int* min_M;
    int* threads;
        
    try {
        X = REAL(get_list_element(list, "X"));
        M = INTEGER(get_list_element(list, "M"));
        N = INTEGER(get_list_element(list, "N"));
        method = INTEGER(get_list_element(list, "method"));
        p1 = p1_pointer(get_list_element(list, "p1"));
        p2 = p2_pointer(get_list_element(list, "p2"));
        k = REAL(get_list_element(list, "k"));
        min_M = INTEGER(get_list_element(list, "min_M"));
        threads = INTEGER(get_list_element(list, "n_threads"));
    }
    catch (exception & e) {
        Rprintf("%s\n", e.what());
        return (R_NilValue);
    }
        
    SEXP x = PROTECT(Rf_allocVector(REALSXP, *N));
    double* _x = REAL(x);
    SEXP s = PROTECT(Rf_allocVector(REALSXP, *M));
    double* _s = REAL(s);

    int no_threads = thread_count(*threads);

    #ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(no_threads);
    #endif

    if (!quantify_matrix(X, *M, *N, *method, p1, p2, *k, *min_M, no_threads, _x, _s)) {
        Rprintf("Unknown method.\n");
    }

    for (int i = 0; i < *M; i++) {
        if (isnan(_s[i])) {
            _s[i] = NA_REAL;
        }
    }
    for (int i = 0; i < *N; i++) {
        if (isnan(_x[i])) {
            _x[i] = NA_REAL;
        }
    }
        
    // return list
    SEXP vec = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(vec, 0, x);
    SET_VECTOR_ELT(vec, 1, s);
        
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, Rf_mkChar("x"));
    SET_STRING_ELT(names, 1, Rf_mkChar("s"));
        
    Rf_setAttrib(vec, R_NamesSymbol, names);
        
    UNPROTECT(4);
        
    return (vec);
        
}

// All proteins at once, in parallel over proteins (like iq_MaxLFQ), each with the
// method of iq_quant (process_matrix). The matrix of a protein is built as in
// process_iq_format: ions (rows) in input order, all samples (columns); an ion with
// more than one value in a sample gets log2 of the summed intensities; non-finite
// values are missing.
//
// list: protein_index, ion_index, sample_index (1-based), quant (log2),
//       n_proteins, n_samples, method, p1, p2, k, min_M, n_threads,
//       rescale_method (0 none, 1-4 as in rescale_components)
// Returns an n_proteins x n_samples matrix (NA rows for proteins without data),
// or NULL if interrupted.
SEXP iq_process(SEXP list) {

    char error_message[256] = "";

    int *proteins = nullptr, *ions = nullptr, *samples = nullptr;
    double* quants = nullptr;
    R_xlen_t nrow = 0;
    int n_proteins = 0, n_samples = 0, method = 0, min_M = 0, n_threads = 0, rescale_method = 0;
    double* p1 = nullptr;
    int* p2 = nullptr;
    double k = 0;
    int verbose = 1;  // 0: no thread count, progress or "Completed." (process_iq_format, per batch)

    try {
        SEXP p = get_list_element(list, "protein_index");
        SEXP i = get_list_element(list, "ion_index");
        SEXP s = get_list_element(list, "sample_index");
        SEXP q = get_list_element(list, "quant");
        nrow = Rf_xlength(p);
        if (TYPEOF(p) != INTSXP || TYPEOF(i) != INTSXP || TYPEOF(s) != INTSXP || TYPEOF(q) != REALSXP ||
            Rf_xlength(i) != nrow || Rf_xlength(s) != nrow || Rf_xlength(q) != nrow) {
            throw runtime_error("protein_index, ion_index, sample_index (integer) and quant (double) must have the same length");
        }
        proteins = INTEGER(p);
        ions = INTEGER(i);
        samples = INTEGER(s);
        quants = REAL(q);

        n_proteins = INTEGER(get_list_element(list, "n_proteins"))[0];
        n_samples = INTEGER(get_list_element(list, "n_samples"))[0];
        method = INTEGER(get_list_element(list, "method"))[0];
        p1 = p1_pointer(get_list_element(list, "p1"));
        p2 = p2_pointer(get_list_element(list, "p2"));
        k = REAL(get_list_element(list, "k"))[0];
        min_M = INTEGER(get_list_element(list, "min_M"))[0];
        n_threads = INTEGER(get_list_element(list, "n_threads"))[0];
        rescale_method = INTEGER(get_list_element(list, "rescale_method"))[0];
        verbose = Rf_asLogical(get_list_element(list, "verbose")) == TRUE;

        if (method < 0 || method > 9) {
            throw runtime_error("unknown method");
        }
        for (R_xlen_t r = 0; r < nrow; r++) {
            if (proteins[r] < 1 || proteins[r] > n_proteins || samples[r] < 1 || samples[r] > n_samples ||
                ions[r] == NA_INTEGER) {
                throw runtime_error("protein_index, sample_index or ion_index out of range");
            }
        }
    } catch (exception& e) {
        snprintf(error_message, sizeof(error_message), "%s", e.what());
    }
    if (error_message[0]) {
        Rf_error("%s", error_message);
    }

    SEXP table = PROTECT(Rf_allocMatrix(REALSXP, n_proteins, n_samples));
    double* buffer = REAL(table);
    fill_n(buffer, (size_t)n_proteins * n_samples, NA_REAL);

    // rows of each protein, in input order
    vector<vector<R_xlen_t>> protein_rows(n_proteins);
    for (R_xlen_t r = 0; r < nrow; r++) {
        protein_rows[proteins[r] - 1].push_back(r);
    }

    int no_threads = thread_count(n_threads);

    #ifdef _OPENMP
    // OpenMP settings of the session, restored at the end
    int saved_dynamic = omp_get_dynamic();
    int saved_threads = omp_get_max_threads();
    int saved_levels = omp_get_max_active_levels();

    omp_set_dynamic(0);
    omp_set_num_threads(no_threads);
    omp_set_max_active_levels(1);  // the methods run single-threaded inside the loop
    if (verbose) {
        Rprintf("Using %d threads...\n", no_threads);
    }
    #endif

    int stop_sig = 0;
    int failed = 0;
    size_t thres_display = 0;

    #pragma omp parallel for schedule(dynamic)
    for (int p = 0; p < n_proteins; p++) {

        if (stop_sig || failed || protein_rows[p].empty()) {
            continue;
        }

        try {
            const vector<R_xlen_t>& rows = protein_rows[p];

            // ion -> row, in order of appearance
            std::unordered_map<int, int> ion_row;
            for (auto r : rows) {
                ion_row.emplace(ions[r], (int)ion_row.size());
            }
            int M = (int)ion_row.size();
            int N = n_samples;

            vector<double> X((size_t)M * N, NAN);
            vector<int> count((size_t)M * N, 0);

            for (auto r : rows) {
                size_t cell = ion_row[ions[r]] + (size_t)(samples[r] - 1) * M;
                double v = quants[r];
                if (count[cell] == 0) {
                    X[cell] = v;
                } else if (count[cell] == 1) {
                    X[cell] = pow(2.0, X[cell]) + pow(2.0, v);  // linear sum
                } else {
                    X[cell] += pow(2.0, v);
                }
                count[cell]++;
            }
            for (size_t c = 0; c < X.size(); c++) {
                if (count[c] > 1) {
                    X[c] = log2(X[c]);
                }
                if (!std::isfinite(X[c])) {
                    X[c] = NAN;
                }
            }

            vector<double> x(N), s(M);
            if (M == 1) {  // as process_matrix: the values of the only ion
                std::copy(X.begin(), X.end(), x.begin());
            } else {
                quantify_matrix(X.data(), M, N, method, p1, p2, k, min_M, 1, x.data(), s.data());
            }

            // as rescale() after process_matrix in process_iq_format (X is not changed by the methods)
            rescale_components(X.data(), M, N, x.data(), rescale_method);

            for (int j = 0; j < N; j++) {
                buffer[p + (size_t)j * n_proteins] = isnan(x[j]) ? NA_REAL : x[j];
            }
        } catch (...) {
            #pragma omp atomic write
            failed = 1;
        }

        #ifdef _OPENMP
        if (omp_get_thread_num() == 0) {
        #endif
            if (verbose && (size_t)p > thres_display) {
                Rprintf("%d%%\n", (int)((size_t)p * 100 / n_proteins));
                R_FlushConsole();
                thres_display = p + n_proteins / 20;
            }
            if (tp_check()) {  // user interrupted ...
                #pragma omp atomic write
                stop_sig = 1;
            }
        #ifdef _OPENMP
        }
        #endif
    }

    #ifdef _OPENMP
    omp_set_dynamic(saved_dynamic);
    omp_set_num_threads(saved_threads);
    omp_set_max_active_levels(saved_levels);
    #endif

    vector<vector<R_xlen_t>>().swap(protein_rows);  // free before Rf_error, which skips destructors

    if (failed) {
        UNPROTECT(1);
        Rf_error("Quantification failed (out of memory?).");
    }

    if (stop_sig) {
        Rprintf("Cancelled.\n");
        UNPROTECT(1);
        return R_NilValue;
    }

    if (verbose) {
        Rprintf("Completed.\n");
    }
    UNPROTECT(1);
    return table;
}

    
    
// R's rescale(): list(X, M, N, x, method), method 1-4 as in rescale_components
SEXP iq_rescale(SEXP list) {

    SEXP X = get_list_element(list, "X");
    SEXP x = get_list_element(list, "x");
    int M = INTEGER(get_list_element(list, "M"))[0];
    int N = INTEGER(get_list_element(list, "N"))[0];
    int method = INTEGER(get_list_element(list, "method"))[0];

    if (TYPEOF(X) != REALSXP || TYPEOF(x) != REALSXP || M < 1 || N < 1 ||
        Rf_xlength(X) != (R_xlen_t)M * N || Rf_xlength(x) != N) {
        Rf_error("rescale: x must have one value per column of X");
    }

    SEXP out = PROTECT(Rf_allocVector(REALSXP, N));
    double* _out = REAL(out);
    std::copy(REAL(x), REAL(x) + N, _out);

    rescale_components(REAL(X), M, N, _out, method);

    for (int j = 0; j < N; j++) {
        if (isnan(_out[j])) {
            _out[j] = NA_REAL;
        }
    }

    UNPROTECT(1);
    return out;
}

SEXP iq_connected_components(SEXP list) {
        
    double* X;
    int* M;
    int* N;
        
    try {
        X = REAL(get_list_element(list, "X"));
        M = INTEGER(get_list_element(list, "M"));
        N = INTEGER(get_list_element(list, "N"));
    }
    catch (exception & e) {
        Rprintf("%s\n", e.what());
        return (R_NilValue);
    }
        
    SEXP components = PROTECT(Rf_allocVector(INTSXP, *N));
    int* _components = INTEGER(components);
    int no_component;
        
    connected_components(X, *M, *N, _components, no_component);
        
    UNPROTECT(1);
        
    return (components);
}
    
SEXP iq_huber_weight(SEXP list) {
        
    double* X;
    int* M;
    int* N;
        
    double* x;
    
    double* c;
    double* k;
        
    try {
        X = REAL(get_list_element(list, "X"));
        M = INTEGER(get_list_element(list, "M"));
        N = INTEGER(get_list_element(list, "N"));
        x = REAL(get_list_element(list, "x"));
        c = REAL(get_list_element(list, "c"));
        k = REAL(get_list_element(list, "k"));
    }
        
    catch (exception & e) {
        Rprintf("%s\n", e.what());
        return (R_NilValue);
    }
        
    SEXP w = PROTECT(Rf_allocVector(REALSXP, *M));
    double* _w = REAL(w);
        
    huber_weight(X, M, N, x, _w, *c, *k); // not a mistake, c and k are swapped
        
    for (int i = 0; i < *M; i++) {
        if (isnan(_w[i])) {
            _w[i] = 0;
        }
    }
        
    UNPROTECT(1);
        
    return (w);
}
}


// TRUE if built with the original MaxLFQ (Eigen, IQ_USE_EIGEN), see src/Makevars
SEXP iq_eigen(void) {
    #ifdef IQ_USE_EIGEN
    return Rf_ScalarLogical(TRUE);
    #else
    return Rf_ScalarLogical(FALSE);
    #endif
}

static const R_CallMethodDef callMethods[]  = {
    {"iq_filter", (DL_FUNC) &iq_filter, 1},
    {"iq_filter_write", (DL_FUNC) &iq_filter_write, 2},
    {"iq_MaxLFQ", (DL_FUNC) &iq_MaxLFQ, 1},
    {"iq_quant", (DL_FUNC) &iq_quant, 1},
    {"iq_process", (DL_FUNC) &iq_process, 1},
    {"iq_rescale", (DL_FUNC) &iq_rescale, 1},
    {"iq_eigen", (DL_FUNC) &iq_eigen, 0},
    {"iq_connected_components", (DL_FUNC) &iq_connected_components, 1},
    {"iq_huber_weight", (DL_FUNC) &iq_huber_weight, 1},
    {"iq_stream_begin", (DL_FUNC) &iq_stream_begin, 1},
    {"iq_stream_add", (DL_FUNC) &iq_stream_add, 2},
    {"iq_stream_end", (DL_FUNC) &iq_stream_end, 2},
    {"iq_stream_write", (DL_FUNC) &iq_stream_write, 2},
    {"iq_writer_begin", (DL_FUNC) &iq_writer_begin, 1},
    {"iq_writer_end", (DL_FUNC) &iq_writer_end, 1},
    {"iq_read_sample_values", (DL_FUNC) &iq_read_sample_values, 4},
    {"iq_read_protein_batch", (DL_FUNC) &iq_read_protein_batch, 5},
    {NULL, NULL, 0}
};

// called by R when the library is loaded; the name must be R_init_<package name>
extern "C" void R_init_iq(DllInfo *info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
