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

// Writer of the iq format, filled incrementally by the readers (long_format_to_iq_format):
// the values of each chunk of the input are appended to per-protein buffers, written to the
// files <dir><protein>.tsv ("sample <tab> ion <tab> intensity", 1-based indices) when the
// buffered data reach a budget, so that the input is never held in memory as a whole.
//
// Each value goes through the per-value steps of fast_preprocess: a value whose log2 is not
// finite (zero or negative intensity) is left out, and only a value whose log2 is above the
// cut-off (if any) is written to its protein file, with its original intensity. The writer
// counts the values written for each sample (the sample medians for the normalization are
// computed afterwards from the protein files, iq_read_sample_values in writer.cpp) and keeps
// a histogram of all finite log2 values, also those not above the cut-off, for the QC plot.

#ifndef __IQ_WRITER_H
#define __IQ_WRITER_H

#include <cmath>
#include <cstdio>
#include <map>
#include <stdexcept>
#include <unordered_set>
#include <string>
#include <vector>

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>

class iq_writer {

public:
    std::string dir;          // ends with '/'
    bool has_cutoff;
    double cutoff;
    size_t budget;            // bytes buffered before a flush

    std::vector<std::string> protein_buffer;      // per protein: lines of its file
    std::vector<char> protein_has_file;
    std::vector<double> sample_count;             // per sample: values written
    size_t buffered = 0;

    // histogram of all finite log2 values, bins of width 1 / HIST_BINS_PER_UNIT
    static constexpr double HIST_BINS_PER_UNIT = 64.0;
    std::map<long long, double> histogram;

    double log2_min = INFINITY;
    double log2_max = -INFINITY;
    double n_written = 0;        // values written to protein files
    double n_not_positive = 0;   // zero, negative or missing intensity
    double n_below_cutoff = 0;
    double n_unidentified = 0;   // without a protein or sample identifier

    // fragment statistics (statistics = true), as in process(): for each protein, the summed
    // intensities of all values read (before preprocessing; summed in long double, as R's sum)
    // and its ions in order of first appearance
    bool statistics;
    std::vector<long double> protein_xic;
    std::vector<std::vector<int>> protein_ions;
    std::unordered_set<unsigned long long> protein_ion_seen;   // (protein, ion)

    iq_writer(const std::string& dir, bool has_cutoff, double cutoff, size_t budget, bool statistics)
        : dir(dir), has_cutoff(has_cutoff), cutoff(cutoff), budget(budget), statistics(statistics) {
    }

    std::string protein_file(size_t p) const {    // p: 0-based
        return dir + std::to_string(p + 1) + ".tsv";
    }

    // one value; p, s, i: 1-based indices of the reader; v: intensity
    void add(int p, int s, int i, double v) {
        if (statistics) {  // all values read, also those left out below
            size_t pi = (size_t)p - 1;
            if (protein_xic.size() <= pi) {
                protein_xic.resize(pi + 1, 0.0L);
                protein_ions.resize(pi + 1);
            }
            protein_xic[pi] += v;
            if (protein_ion_seen.insert(((unsigned long long)(unsigned)p << 32) | (unsigned)i).second) {
                protein_ions[pi].push_back(i);
            }
        }

        double l = std::log2(v);
        if (!std::isfinite(l)) {
            n_not_positive++;
            return;
        }

        histogram[(long long)std::floor(l * HIST_BINS_PER_UNIT)]++;

        if (l < log2_min) {
            log2_min = l;
        }
        if (l > log2_max) {
            log2_max = l;
        }

        if (has_cutoff && !(l > cutoff)) {
            n_below_cutoff++;
        } else {
            size_t pi = (size_t)p - 1;
            if (protein_buffer.size() <= pi) {
                protein_buffer.resize(pi + 1);
                protein_has_file.resize(pi + 1, 0);
            }
            char line[96];
            int n = snprintf(line, sizeof(line), "%d\t%d\t%.17g\n", s, i, v);
            protein_buffer[pi].append(line, (size_t)n);
            buffered += (size_t)n;
            n_written++;

            size_t si = (size_t)s - 1;
            if (sample_count.size() <= si) {
                sample_count.resize(si + 1, 0.0);
            }
            sample_count[si]++;
        }

        if (buffered >= budget) {
            flush();
        }
    }

    // append all buffered data to the files and free the buffers
    void flush() {
        for (size_t p = 0; p < protein_buffer.size(); p++) {
            std::string& b = protein_buffer[p];
            if (b.empty()) {
                continue;
            }
            std::string name = protein_file(p);
            FILE* f = fopen(name.c_str(), "ab");
            if (f == NULL) {
                throw std::runtime_error("Cannot write " + name);
            }
            size_t written = fwrite(b.data(), 1, b.size(), f);
            int closed = fclose(f);
            if (written != b.size() || closed != 0) {
                throw std::runtime_error("Cannot write " + name);
            }
            protein_has_file[p] = 1;
            std::string().swap(b);
        }

        buffered = 0;
    }
};

// the writer of an external pointer made by iq_writer_begin (writer.cpp)
inline iq_writer* get_writer(SEXP ptr) {
    if (TYPEOF(ptr) != EXTPTRSXP || R_ExternalPtrTag(ptr) != Rf_install("iq_writer")) {
        throw std::runtime_error("invalid iq format writer");
    }
    iq_writer* w = static_cast<iq_writer*>(R_ExternalPtrAddr(ptr));
    if (!w) {
        throw std::runtime_error("the iq format writer has already been finished");
    }
    return w;
}

#endif
