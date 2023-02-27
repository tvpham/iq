/*
#########################################################################

Author: Thang V. Pham, t.pham@amsterdamumc.nl

All rights reserved.

Citation:

Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative
protein abundances from ion quantification in DIA-MS-based proteomics,
Bioinformatics 2020 Apr 15;36(8):2611-2613.

Software version: 1.9.9

#########################################################################
*/

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <utility>
#include <exception>
#include <stdexcept>

#include <RcppEigen.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define R_NO_REMAP

#include <R.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>

extern "C" SEXP iq_filter(SEXP);
extern "C" SEXP iq_MaxLFQ(SEXP);

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

            *positions = (size_t *)realloc(*positions, (BUFFER_SIZE + 1) * sizeof(size_t));  // if NULL, realloc = malloc
            if (*positions == NULL) {
                if (*buffer) {
                    free(*buffer);
                }
                throw runtime_error("Cannot allocate memory");
            }
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
            *buffer = (char *)realloc(*buffer, *n * 2);
            if (*buffer == NULL) {
                if (*positions) {
                    free(*positions);
                }
                throw runtime_error("Cannot allocate memory");
            }
            remaining_buffer_size = *n + 1;
            ptr = *buffer + (*n - 1);
            *n *= 2;

            *positions = (size_t *)realloc(*positions, (*n + 1) * sizeof(size_t));  // if NULL, realloc = malloc
            if (*positions == NULL) {
                if (*buffer) {
                    free(*buffer);
                }
                throw runtime_error("Cannot allocate memory");
            }
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
            buf_0 = (char *)realloc(buf_0, n * 2);

            if (buf_0 == NULL) {
                free(buf_1);
                return -1;
            }

            buf_1 = (char *)realloc(buf_1, n * 2);

            if (buf_1 == NULL) {
                free(buf_0);
                return -1;
            }

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
             vector<vector<string>*> *ions) {
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
    const char *col_quant = "";
    vector<const char *> col_annotations;

    vector<pair<const char *, const char *>> filter_string_equal;
    vector<pair<const char*, const char*>> filter_string_not_equal;

    vector<pair<const char *, double>> filter_double_less;
    vector<pair<const char*, double>> filter_double_greater;

    const char* intensity_col_sep = 0;
    const char* intensity_col_id = 0;
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
    string output_dir = "";

    size_t p = 0;  // parameter
    while (p < (argv.size() - 1)) {
        if (argv[p].compare("--no-annotation") == 0) {
            col_annotations.clear();
            p++;
        } else if (argv[p].compare("--header") == 0) {
            show_header = true;
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
            col_quant = argv[++p].c_str();
            p++;
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
            intensity_col_id = argv[++p].c_str();
            p++;
        } else if (argv[p].compare("--na_string") == 0) {
            na_string = argv[++p].c_str();
            p++;
        } else {
            throw runtime_error((string("Do not know what to do with ") + argv[p]).c_str());
        }
    }

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

    if (utils::getline_tab(&line, &len, fp, &tab_pos, &n_tab_pos) == -1) {
        fclose(fp);
        if (line) {
            free(line);
        }
        if (tab_pos) {
            free(tab_pos);
        }
        throw runtime_error((string("Cannot parse header of ") + filename).c_str());
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
    Rprintf("\nQuant column:\n    %s\n", col_quant);
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

    if (intensity_col_id) {
        Rprintf("Quant id column:\n    %s\n", intensity_col_id);
    }

    //--- check header
    size_t n_cols = n_tab_pos;

    size_t col_sample_ind;
    utils::map_header(&col_sample, 1, &col_sample_ind, line, tab_pos, n_tab_pos);

    size_t col_primary_id_ind;
    utils::map_header(&col_primary_id, 1, &col_primary_id_ind, line, tab_pos, n_tab_pos);

    vector<size_t> col_secondary_ids_ind(col_secondary_ids.size());
    utils::map_header(col_secondary_ids.data(), col_secondary_ids.size(), col_secondary_ids_ind.data(), line, tab_pos, n_tab_pos);

    size_t col_quant_ind;
    utils::map_header(&col_quant, 1, &col_quant_ind, line, tab_pos, n_tab_pos);

    vector<size_t> col_annotations_ind(col_annotations.size());
    utils::map_header(col_annotations.data(), col_annotations.size(), col_annotations_ind.data(), line, tab_pos, n_tab_pos);

    vector<size_t> col_filter_string_equal(filter_string_equal.size());
    for (size_t i = 0; i < filter_string_equal.size(); i++) {
        utils::map_header(&(filter_string_equal[i].first), 1, &(col_filter_string_equal[i]), line, tab_pos, n_tab_pos);
    }

    vector<size_t> col_filter_string_not_equal(filter_string_not_equal.size());
    for (size_t i = 0; i < filter_string_not_equal.size(); i++) {
        utils::map_header(&(filter_string_not_equal[i].first), 1, &(col_filter_string_not_equal[i]), line, tab_pos, n_tab_pos);
    }

    vector<size_t> col_filter_double_less(filter_double_less.size());
    for (size_t i = 0; i < filter_double_less.size(); i++) {
        utils::map_header(&(filter_double_less[i].first), 1, &(col_filter_double_less[i]), line, tab_pos, n_tab_pos);
    }

    vector<size_t> col_filter_double_greater(filter_double_greater.size());
    for (size_t i = 0; i < filter_double_greater.size(); i++) {
        utils::map_header(&(filter_double_greater[i].first), 1, &(col_filter_double_greater[i]), line, tab_pos, n_tab_pos);
    }

    annotation_colnames->push_back(string(col_primary_id));
    for (auto n : col_annotations) {
        annotation_colnames->push_back(string(n));
    }

    size_t intensity_col_id_ind;
    if (intensity_col_id) {
        utils::map_header(&intensity_col_id, 1, &intensity_col_id_ind, line, tab_pos, n_tab_pos);
    }


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

            if ((*quant_vec_is_na_n)[n] == 0) {
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

                    if (intensity_col_id) {
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
                        if (intensity_col_id) {
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
        Rprintf("\nUsing %d threads ...\n", no_threads);
    #else
        Rprintf("Using a single CPU core ...\n");
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

            int line_mismatched = 0;

            #pragma omp parallel for schedule(static)
            for (size_t n = 0; n < _line->size(); n++) {

                (*_tab)[n][0] = (*_line)[n];
                size_t pos = 0;
                char *c = (*_line)[n];
                while (*c) {
                    if (*c == '\r') {
                        *c++ = 0;
                    } else if (*c == '\t') {
                        *c = 0;
                        pos++;
                        if (pos == n_cols) {
                            line_mismatched = line_no + n + 1;
                        }
                        (*_tab)[n][pos] = ++c;
                    } else {
                        c++;
                    }
                }

                if ((pos + 1) != n_cols) {
                    line_mismatched = line_no + n + 1;
                }

                (*quant_vec)[n].clear();
                (*quant_id_vec)[n].clear();
                (*quant_vec_is_na)[n].clear();
                (*quant_vec_is_na_n)[n] = 0;

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
                            if (strcmp((*_tab)[n][col_filter_double_less[i]], "NaN") != 0 && atof((*_tab)[n][col_filter_double_less[i]]) >= filter_double_less[i].second) { // gcc atof does not work for NaN
                                ok = false;
                            }
                            i++;
                        }

                        if (ok) {

                            i = 0;
                            while (i < filter_double_greater.size() && ok) {
                                if (strcmp((*_tab)[n][col_filter_double_greater[i]], "NaN") != 0 && atof((*_tab)[n][col_filter_double_greater[i]]) <= filter_double_greater[i].second) { // gcc atof does not work for NaN
                                    ok = false;
                                }
                                i++;
                            }

                            if (ok) {

                                if (intensity_col_sep) {
                                    utils::split((*_tab)[n][col_quant_ind], intensity_col_sep, &((*quant_vec)[n]));
                                    if (intensity_col_id) {
                                        utils::split((*_tab)[n][intensity_col_id_ind], intensity_col_sep, &((*quant_id_vec)[n]));
                                        if ((*quant_vec)[n].size() != (*quant_id_vec)[n].size()) {
                                            line_mismatched = line_no + n + 1;
                                        }
                                    }
                                }
                                else {
                                    (*quant_vec)[n].push_back((*_tab)[n][col_quant_ind]);
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

            if (line_mismatched > 0) { // no exception in parallel regions
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

        if (samples->size() >= print_threshold) {
            Rprintf("%d samples read\n", samples->size());
            print_threshold = samples->size() + 20;
        }
    }
    Rprintf("\n# lines read (excluding headers)      = %d\n", line_no-1);
    Rprintf("# quantitative values after filtering = %d\n\n", n_after_filtered);

    Rprintf("# samples  = %d\n", samples->size());
    Rprintf("# proteins = %d\n", annotation->size());

    cleanup();
}

// [[Rcpp::export]]
SEXP iq_filter(SEXP cmd) {

    const char *c = CHAR(STRING_ELT(cmd, 0));

    vector<string> argv;

    const char *p = c;

    const char sep = 1;

    while (*p == sep) {
        p++;
    }

    while (*p) {
        string s = "";
        while (*p && *p != sep) {
            s = s + *p++;
        }
        argv.push_back(s);

        while (*p == sep) {
            p++;
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

    try {
        process(argv, protein_index, sample_index, ion_index, quant, annotation, col_annotation, samples, ions);
    } catch (exception & e) {
        Rprintf("%s\n", e.what());

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

    if (ions->size() == 0) {

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

    SEXP vec = PROTECT(Rf_allocVector(VECSXP, 4));
    SET_VECTOR_ELT(vec, 0, quant_tab);
    SET_VECTOR_ELT(vec, 1, ann);
    SET_VECTOR_ELT(vec, 2, sample_str);
    SET_VECTOR_ELT(vec, 3, ion_str);

    SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
    SET_STRING_ELT(names, 0, Rf_mkChar("quant_table"));
    SET_STRING_ELT(names, 1, Rf_mkChar("protein"));
    SET_STRING_ELT(names, 2, Rf_mkChar("sample"));
    SET_STRING_ELT(names, 3, Rf_mkChar("ion"));

    Rf_setAttrib(vec, R_NamesSymbol, names);

    UNPROTECT(13);

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

//------------------------------ MaxLFQ -------------------------------------

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

    // add an entry to a map, increase size if necessary
    void add_to_table(int ion, int col, double quant) {
        auto r = map.emplace(ion, vector<double>());
        if (r.second) {
            //n++;
            r.first->second.resize(ncol, NAN);
        }
        r.first->second[col] = quant;
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

// [[Rcpp::export]]
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

// [[Rcpp::export]]
SEXP iq_MaxLFQ(SEXP list) {

    int stop_sig = 0;

    int* proteins;
    int* ions;
    int* samples;
    double* quants;

    try {
        proteins = INTEGER(get_list_element(list, "protein_index"));
        ions = INTEGER(get_list_element(list, "ion_index"));
        samples = INTEGER(get_list_element(list, "sample_index"));
        quants = REAL(get_list_element(list, "quant"));
    }
    catch (exception & e) {
        Rprintf("%s\n", e.what());
        return (R_NilValue);
    }

    size_t nrow = Rf_xlength(get_list_element(list, "protein_index"));

    int n_proteins = utils::count_unique(proteins, nrow);
    int n_samples = utils::count_unique(samples, nrow);

    SEXP table = PROTECT(Rf_allocMatrix(REALSXP, n_proteins, n_samples));
    double* buffer = REAL(table);

    auto group_annotation = new vector<string>(n_proteins, "");
    int* row_names = new int[n_proteins];
    int* col_names = new int[n_samples];

    Rprintf("nrow = %d, # proteins = %d, # samples = %d\n", nrow, n_proteins, n_samples);

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

        Rprintf("Using %d threads...\n", no_threads);

    #else

        Rprintf("Using a single CPU core...\n");

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
            if (i > thres_display) {
                Rprintf("%d%%\n", i * 100 / (*protein_index).size());
                R_FlushConsole();
                thres_display = i + (*protein_index).size() / 20;
            }

            if (tp_check()) {  // user interrupted ...
                stop_sig = 1;
                #pragma omp flush(stop_sig)
            }
        }
    }

    if (stop_sig) {
        Rprintf("Canceled.\n");

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

    Rprintf("Completed.\n");

    UNPROTECT(7);

    delete protein_index;
    delete group_annotation;
    delete[] col_names;
    delete[] row_names;

    return (vec);
}

static const R_CallMethodDef callMethods[]  = {
    {"iq_filter", (DL_FUNC) &iq_filter, 1},
    {"iq_MaxLFQ", (DL_FUNC) &iq_MaxLFQ, 1},
    {NULL, NULL, 0}
};

void R_init_myLib(DllInfo *info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}
