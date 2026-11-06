/*#########################################################################
 
 Authors: Thang V. Pham {t.pham@amsterdamumc.nl}, Pham Huy Chau Long, 
          Tran Thi Minh Chau
 
 All rights reserved.
 
 Citation:
 
 Pham TV, Tran CTM, Henneman AA, Pham LHC, Le DG, Can AH, Bui PHL, 
 Piersma SR, Jimenez CR. Boosting the speed and accuracy of protein 
 quantification algorithms in mass spectrometry-based proteomics.
 Journal of Proteome Research 2026 Feb 6; 25(2):1198-1203.
 
 Software version: 2.2.0
 
#########################################################################*/

#ifndef __UTILS_H
#define __UTILS_H

#include <iostream>
#include <vector>
#include <algorithm> 
#include <numeric>   
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>

struct weighted_sample_t {

    double value;
    double weight;

    bool operator<(const weighted_sample_t& other) const {
        if (value == other.value) {
            return weight < other.weight;
        }
        else {
            return value < other.value;
        }
    }
};


inline double quick_median(double *median_buffer, size_t n) {
    
    if (n == 0) {
        return NAN;
    }

    const size_t mid = n / 2;
    
    std::nth_element(median_buffer, median_buffer + mid, median_buffer + n);
    
    if (n % 2 == 1) {
        return median_buffer[mid];
    } else {
        double mid_val = median_buffer[mid];
        std::nth_element(median_buffer, median_buffer + mid - 1, median_buffer + n);
        return (median_buffer[mid - 1] + mid_val) / 2.0;
    }
}


// to check: equal values, different weight
inline double weighted_median(std::vector<weighted_sample_t>& points, size_t n = 0) {
    
    if (n ==0) {
        n = points.size();
    }
    
    sort(points.begin(), points.begin() + n);
    
    // normalize w
    double sum_w = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum_w += points[i].weight;
        points[i].weight *= 2;  // maybe more stable than sum_w/2 for integer weight
    }
    
    // half the weight on each side is a tie up to rounding: weights that are equal in exact
    // arithmetic (e.g. the two ions linking two samples, whose residuals mirror each other)
    // can differ in the last digits, and an exact comparison would then pick either value
    const double tolerance = 1e-12 * sum_w;

    double left_sum = 0;

    for (size_t pos = 0; pos < n; pos++) {
        left_sum += points[pos].weight;
        if (left_sum >= sum_w - tolerance) {
            // half the weight on each side: average the two middle values
            if (left_sum <= sum_w + tolerance && pos < (n-1)) {
                return (points[pos].value + points[pos+1].value) / 2.0;
            }
            return points[pos].value;
        }
    }
    
    return points[n - 1].value;
}


inline void connected_components(double* X, int M, int N, int *components, int& no_component) {
    
    std::fill_n(components, N, -1);
    
    int* next = new int[N];
    int* previous = new int[N];
    
    int start = 0;
    for (int i = 0; i < N; i++) {
        next[i] = i + 1;
        previous[i] = i-1;
    }
    next[N-1] = -1;
    
    no_component = 0;
    for (int i = 0; i < N; i++) {
        if (components[i] < 0) {
            
            components[i] = no_component;
            
            previous[i] = -1;
            start = next[i];
            
            std::vector<int> stack;
            stack.reserve(N); 
            stack.push_back(i);
            
            while (!stack.empty()) {
                int current = stack.back();
                stack.pop_back();
                
                // For each peptide measurement (row) involving current sample
                for (int r = 0; r < M; r++) {
                    if (!std::isnan(X[r + current * M])) {
                        // Find all other samples measured in this peptide
                        int k = start;
                        //for (int k = 0; k < N; k++) {
                        while (k >= 0) {
                            if (components[k] == -1 && !std::isnan(X[r + k * M])) {
                                // Found a new connected sample
                                components[k] = no_component;
                                
                                //remove k
                                if (previous[k] >= 0) {
                                    next[previous[k]] = next[k];
                                    if (next[k] >=0) {
                                        previous[next[k]] = previous[k];
                                    }
                                }
                                else {
                                    if (next[k] >=0) {
                                        previous[next[k]] = previous[k];
                                    }
                                    //start = next[k];
                                }
                                
                                stack.push_back(k);
                                
                            }
                            k = next[k];
                        }
                    }
                }
            }
            
            no_component++;
        }
    }
    
    delete[] next;
    delete[] previous;
}


inline void rescale(double* X, int M, int N, double *x, double *x_out) {
    
    std::copy_n(x, N, x_out);
    
    if (M == 1) {
        return;
    }
    
    int* components = new int[N];
    int n_components;
    
    connected_components(X, M, N, components, n_components);
    
    for (int c = 0; c < n_components; c++) {
        
        // mean data
        double sum = 0.0;
        int count = 0;
        double sum_x = 0.0;
        int count_x = 0;
        
        int col;
        for (int j = 0; j < N; j++) {
            if (components[j] == c) {
                
                if (!std::isnan(x[j])) {
                    sum_x += x[j];
                    count_x++;
                    col = j;
                }
                
                double* jp = X + j * M;
                for (int i = 0; i < M; i++) {
                    if (!std::isnan(*jp)) {
                        sum += *jp;
                        count++;
                    }
                    jp++;
                }
            }
        }
        
        if (count > 0) {
            if (count_x > 1) {
                double m = sum / (double)count;
                double m_x = sum_x / (double)count_x;
                for (int j = 0; j < N; j++) {
                    if (components[j] == c) {
                        x_out[j] = x[j] - m_x + m;
                    }
                }
            }
            else {
                
                double* buffer = new double[M];
                int buffer_size = 0;
                double* jp = X + col * M;
                for (int i = 0; i < M; i++) {
                    if (!std::isnan(*jp)) {
                        buffer[buffer_size++] = *jp;
                    }
                    jp++;
                }
                
                x_out[col] = quick_median(buffer, buffer_size);
                
                delete [] buffer;
            }
        }
        else {
            for (int j = 0; j < N; j++) {
                if (components[j] == c) {
                    x_out[j] = NAN;
                }
            }
        }
    }

    delete[] components;
}


// The readers return their values as R vectors with an int length (matrices of one column):
// at most INT_MAX values in memory (process(), fast_read*); the iq format has no such limit.
#ifndef IQ_MAX_VALUES
#define IQ_MAX_VALUES INT_MAX  // smaller only in a test build (-DIQ_MAX_VALUES=...)
#endif
static const char* const too_many_values_message =
    "Too many values to hold in memory (more than 2147483647); use long_format_to_iq_format and process_iq_format for data of this size";

// the mean of n values; NaN if there are none
inline double mean_value(const double* v, size_t n) {
    if (n == 0) {
        return NAN;
    }
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum += v[i];
    }
    return sum / n;
}

// log2(sum(2^v)), as m + log2(sum(2^(v - m))) in R; NaN if there are no values
inline double log2_sum(const double* v, size_t n) {
    if (n == 0) {
        return NAN;
    }
    double m = *std::max_element(v, v + n);
    long double s = 0.0;
    for (size_t i = 0; i < n; i++) {
        s += pow(2.0, v[i] - m);
    }
    return m + log2((double)s);
}

// The rescaling of the quantification x (N samples) of one protein with data X
// (M ions x N samples, column-major, NaN = missing), in place. Each connected
// component of samples is rescaled separately; a component whose estimates are
// all missing is left as it is.
// method: 1 "median_mean", 2 "to_median", 3 "to_mean", 4 "to_sum"; other values: no rescaling
//   median_mean, to_median, to_mean: a single-sample component gets the median
//     (mean for to_mean) of its data; a larger one is shifted so that the mean (median
//     for to_median) of its estimates equals that of its data
//   to_sum: a component is shifted so that its estimates in linear scale (2^x) add up to the
//     summed intensities of its data (2^X)
// Used by iq_process and by rescale() in R (iq_rescale).
inline void rescale_components(double* X, int M, int N, double* x, int method) {

    if (method < 1 || method > 4) {
        return;
    }

    std::vector<int> components(N);
    int n_components;
    connected_components(X, M, N, components.data(), n_components);

    std::vector<double> data;       // the data of the component, column by column
    std::vector<double> estimates;  // its non-missing estimates

    for (int c = 0; c < n_components; c++) {

        data.clear();
        estimates.clear();
        int size = 0;

        for (int j = 0; j < N; j++) {
            if (components[j] == c) {
                size++;
                if (!std::isnan(x[j])) {
                    estimates.push_back(x[j]);
                }
                for (int i = 0; i < M; i++) {
                    if (!std::isnan(X[i + (size_t)j * M])) {
                        data.push_back(X[i + (size_t)j * M]);
                    }
                }
            }
        }

        if (estimates.empty()) {
            continue;
        }

        if (method == 4) {  // sum
            double a = log2_sum(data.data(), data.size());
            double b = log2_sum(estimates.data(), estimates.size());
            for (int j = 0; j < N; j++) {
                if (components[j] == c) {
                    x[j] = (x[j] + a) - b;
                }
            }
        } else if (size == 1) {
            double value = (method == 3) ? mean_value(data.data(), data.size())
                                         : quick_median(data.data(), data.size());
            for (int j = 0; j < N; j++) {
                if (components[j] == c) {
                    x[j] = value;
                }
            }
        } else {
            double m, m_x;
            if (method == 2) {
                m = quick_median(data.data(), data.size());
                m_x = quick_median(estimates.data(), estimates.size());
            } else {
                m = mean_value(data.data(), data.size());
                m_x = mean_value(estimates.data(), estimates.size());
            }
            for (int j = 0; j < N; j++) {
                if (components[j] == c) {
                    x[j] = (x[j] - m_x) + m;
                }
            }
        }
    }
}

// The value of a filter column (filter_double_less / filter_double_greater) given as text.
// Returns false if the value is missing: empty, NaN, or text that is not a number
// (e.g. "NA"); such rows are kept. Used by both readers (iq.cpp and stream.cpp).
inline bool filter_number(const char* s, double* v) {
    if (strcmp(s, "NaN") == 0) {  // some atof/strtod do not parse NaN
        return false;
    }
    char* end = NULL;
    double x = strtod(s, &end);
    if (end == s || std::isnan(x)) {
        return false;
    }
    *v = x;
    return true;
}

#endif // __UTILS_H