/*#########################################################################
 
 Author: Thang V. Pham, t.pham@amsterdamumc.nl, Pham Huy Chau Long, 
         Tran Thi Minh Chau
 
 All rights reserved.
 
 Citation:
 
 Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative
 protein abundances from ion quantification in DIA-MS-based proteomics,
 Bioinformatics 2020 Apr 15;36(8):2611-2613.
 
 Software version: 2.0
 
#########################################################################*/

#ifndef __UTILS_H
#define __UTILS_H

#include <iostream>
#include <vector>
#include <algorithm> 
#include <numeric>   

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


double quick_median(double *median_buffer, size_t n) {
    
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
    
    double left_sum = 0;

    for (size_t pos = 0; pos < n; pos++) {
        left_sum += points[pos].weight;
        if (left_sum >= sum_w) {
            if (pos < (n-1)) {
                return (left_sum * points[pos].value + (2*sum_w-left_sum)* points[pos+1].value)/((2*sum_w));
            }
            return points[pos].value;
        }
    }
    
    return points[n - 1].value;
}


void connected_components(double* X, int M, int N, int *components, int& no_component) {
    
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


void rescale(double* X, int M, int N, double *x, double *x_out) {
    
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
}


#endif // __UTILS_H