/*#########################################################################
 
 Author: Thang V. Pham, t.pham@amsterdamumc.nl
 
 All rights reserved.
 
 Citation:
 
 Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative
 protein abundances from ion quantification in DIA-MS-based proteomics,
 Bioinformatics 2020 Apr 15;36(8):2611-2613.
 
 Software version: 2.0
 
#########################################################################*/

#ifndef __CDW_H
#define __CDW_H

#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "utils.h"

using namespace std;

#define EPSILON 1e-10

void median_polish(double *X, const int M, const int N, double *s, double *x, 
                   int no_threads = 1, 
                   double eps = 1e-4,
                   int max_iteration = 1000) {
    // X column-based matrix with M rows, N columns
    double *big_median_buffer = new double[(M > N ? M : N) * no_threads];

    fill_n(s, M, NAN);
    fill_n(x, N, NAN);

    double t = NAN;

    double* a = new double[M*N];
    
    double* ap = a;
    double* xp = X;
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < M; i++) {
            *ap = *xp;
            if (!isnan(*ap)) {
                s[i] = 0;
                x[j] = 0;
                t = 0;
            }
            ap++;
            xp++;
        }
    }
    
    if (isnan(t)) {
        return;
    }
    
    double old_sum = 0;

    for (int iteration = 0; iteration < max_iteration; iteration++) {
        
        // row median
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < M; i++) {
            
            if (isnan(s[i])) continue;
            
            int thread_id = 0;
            #ifdef _OPENMP
            thread_id = omp_get_thread_num();
            #endif
            
            double* median_buffer = big_median_buffer + N * thread_id;
            size_t num_of_coords = 0;
    
            double* ap = a + i;        
            for (int j = 0; j < N; j++) {
                if (!isnan(*ap)) {
                    median_buffer[num_of_coords++] = *ap;
                }
                ap += M;
            }

            double m = quick_median(median_buffer, num_of_coords);

            ap = a + i;
            for (int j = 0; j < N; j++) {
                if (!isnan(*ap)) {
                    *ap -= m;
                }
                ap += M;
            }
            
            s[i] += m;
            
        }

        // update t
        size_t num_of_coords = 0;
        for (int j = 0; j < N; j++) {
            if (!isnan(x[j])) {
                big_median_buffer[num_of_coords++] = x[j];
            }
        }

        double m = quick_median(big_median_buffer, num_of_coords);

        for (int j = 0; j < N; j++) {
            if (!isnan(x[j])) {
                x[j] -= m;
            }
        }

        t += m;

        // col median
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < N; j++) {
            if (isnan(x[j])) continue;

            int thread_id = 0;
            #ifdef _OPENMP
            thread_id = omp_get_thread_num();
            #endif
            
            double* median_buffer = big_median_buffer + M * thread_id;
            size_t num_of_coords = 0;
            
            double* ap = a + j * M;
            for (int i = 0; i < M; i++) {
                if (!isnan(*ap)) {
                    median_buffer[num_of_coords++] = *ap;
                }
                ap++;
            }

            double m = quick_median(median_buffer, num_of_coords);

            ap = a + j * M;
            for (int i = 0; i < M; i++) {
                if (!isnan(*ap)) {
                    *ap -= m;
                }
                ap++;
            }

            x[j] += m;
        }
        
        // update t
        num_of_coords = 0;
        for (int i = 0; i < M; i++) {
            if (!isnan(s[i])) {
                big_median_buffer[num_of_coords++] = s[i];
            }
        }

        m = quick_median(big_median_buffer, num_of_coords);

        for (int i = 0; i < M; i++) {
            if (!isnan(s[i])) {
                s[i] -= m;
            }
        }

        t += m;

        // calculate residual
        double new_sum = 0;
        double* ap = a;
        for (int i = 0; i < M*N; i++) {
            if (!isnan(*ap)) {
                new_sum += abs(*ap);
            }
            ap++;
        }


        bool converged = (new_sum == 0) || (abs(new_sum - old_sum) < eps * new_sum);
        if (converged) {
            break;
        }

        old_sum = new_sum;
    }

    delete[] a;
    delete[] big_median_buffer;
}

void weighted_median_polish(double *X, const int M, const int N, 
                              double *s, double *x, double *w, 
                              double eps = 1e-4,
                              int max_iteration = 1000) {
    // X column-based matrix with M rows, N columns

    double *median_buffer = new double[M > N ? M : N];

    std::vector<weighted_sample_t> points(M);

    fill_n(s, M, NAN);
    fill_n(x, N, NAN);

    double t = NAN;

    double* a = new double[M*N];
    
    double* ap = a;
    double* xp = X;
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < M; i++) {
            *ap = *xp;
            if (!isnan(*ap)) {
                s[i] = 0;
                x[j] = 0;
                t = 0;
            }
            ap++;
            xp++;
        }
    }
    

    if (isnan(t)) {
        return;
    }

    double old_sum = 0;

    for (int iteration = 0; iteration < max_iteration; iteration++) {

        // row median
        for (int i = 0; i < M; i++) {

            if (isnan(s[i])) continue;

            size_t num_of_coords = 0;
            double* ap = a + i;        
            for (int j = 0; j < N; j++) {
                if (!isnan(*ap)) {
                    median_buffer[num_of_coords++] = *ap;
                }
                ap += M;
            }

            double m = quick_median(median_buffer, num_of_coords);

            ap = a + i;
            for (int j = 0; j < N; j++) {
                if (!isnan(*ap)) {
                    *ap -= m;
                }
                ap += M;
            }
            
            s[i] += m;
        }

        // update t
        size_t num_of_coords = 0;
        for (int j = 0; j < N; j++) {
            if (!isnan(x[j])) {
                median_buffer[num_of_coords++] = x[j];
            }
        }

        double m = quick_median(median_buffer, num_of_coords);

        for (int j = 0; j < N; j++) {
            if (!isnan(x[j])) {
                x[j] -= m;
            }
        }

        t += m;

        // col median
        for (int j = 0; j < N; j++) {
            if (isnan(x[j])) continue;

            size_t median_size = 0;
            double* ap = a + j * M;
            for (int i = 0; i < M; i++) {
                if (!isnan(*ap)) {
                    points[median_size].value = *ap;
                    points[median_size].weight = w[i];
                    median_size++;
                }
                ap++;
            }

            double m = weighted_median(points, median_size);


            ap = a + j * M;
            for (int i = 0; i < M; i++) {
                if (!isnan(*ap)) {
                    *ap -= m;
                }
                ap++;
            }
            

            x[j] += m;
            
        }

        // update t
        num_of_coords = 0;
        for (int i = 0; i < M; i++) {
            if (!isnan(s[i])) {
                median_buffer[num_of_coords++] = s[i];
            }
        }

        m = quick_median(median_buffer, num_of_coords);

        for (int i = 0; i < M; i++) {
            if (!isnan(s[i])) {
                s[i] -= m;
            }
        }

        t += m;
        

        // calculate residual
        double new_sum = 0;
        double* ap = a;
        for (int i = 0; i < M*N; i++) {
            if (!isnan(*ap)) {
                new_sum += abs(*ap);
            }
            ap++;
        }
        
        bool converged = (new_sum == 0) || (abs(new_sum - old_sum) < eps * new_sum);
        if (converged) {
            break;
        }
        
        old_sum = new_sum;
    }

    delete[] a;
    delete[] median_buffer;
    
    return;
}


void weighted_least_squares(double *X, const int M, const int N, double *W, double *s, double *x, 
                            int no_threads = 1,
                            int max_iteration = 1000) {
    // X column-based matrix with M rows, N columns
    
    double t = NAN;
    
    double* a = new double[M * N];
    
    double* ap = a;
    double* xp = X;
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < M; i++) {
            *ap = *xp;
            if (!isnan(*ap)) {
                *ap = *ap - s[i] - x[j];
                t = 0;
            }
            ap++;
            xp++;
        }
    }
    
    if (isnan(t)) {
        return;
    }
    
    double old_sum = 0;
    
    for (int iteration = 0; iteration < max_iteration; iteration++) {
        
        // row
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < M; i++) {
            if (isnan(s[i])) continue;
            
            double sum_w = 0;
            double sum_aw = 0;
            double* wp = W + i;
        
            double* ap = a + i;    
            for (int j = 0; j < N; j++) {
                if (!isnan(*ap)) {
                    sum_w += *wp;
                    sum_aw += *wp * (*ap);
                }
                wp += M;
                ap += M;
            }
            
            double delta = sum_aw / sum_w;
            
            ap = a + i;
            for (int j = 0; j < N; j++) {
                if (!isnan(*ap)) {
                    *ap -= delta;
                }
                ap += M;
            }
            s[i] += delta;
        }
        
        // col
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < N; j++) {
            
            if (isnan(x[j])) continue;
           
            double sum_w = 0;
            double sum_aw = 0;
            double* wp = W + j * M;
            double* ap = a + j * M;
            
            for (int i = 0; i < M; i++) {
                if (!isnan(*ap)) {
                    sum_w += *wp;
                    sum_aw += *wp * (*ap);
                }
                wp++;
                ap++;
            }
            
            double delta = sum_aw / sum_w;
            
            ap = a + j * M;
            for (int i = 0; i < M; i++) {
                if (!isnan(*ap)) {
                    *ap -= delta;
                }
                ap++;
            }
            
            x[j] += delta;
        }
        
        // calculate residual
        double new_sum = 0;
        double* ap = a;
        for (int i = 0; i < M*N; i++) {
            if (!isnan(*ap)) {
                new_sum += abs(*ap);
            }
            ap++;
        }
        
        bool converged = (new_sum == 0) || (abs(new_sum - old_sum) < EPSILON * new_sum);
        if (converged) {
            break;
        }
        old_sum = new_sum;
    }
    delete[] a;
}


double mad(double *X, const int M, const int N, double *W, double *s, double *x) {
    
    int cc = 0;
    
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            if (!isnan(X[j * M + i])) {
                W[cc++] = X[j * M + i] - s[i] - x[j];
            }
        }
    }
    
    double center = quick_median(W, cc);
    
    for (int i = 0; i < cc; i++) {
        W[i] = abs(W[i] - center);
    }
    
    return (1.4826 * quick_median(W, cc));
    
}


void huber(double *X, const int M, const int N, double *s, double *x, 
           int no_threads = 1,
           const double k = 1.345, int max_iteration = 1000) {
    // X column-based matrix with M rows, N columns
    
    fill_n(s, M, NAN);
    fill_n(x, N, NAN);
    
    double *W = new double[M*N];
    fill_n(W, M*N, 1.0);
    
    double t = NAN;
    
    //vector<vector<double>> a(M, vector<double>(N, 0.0));
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            if (!isnan(X[j * M + i])) {
                s[i] = 0;
                x[j] = 0;
                //a[i][j] = 1.0;
                t = 0;
            }
        }
    }
    
    if (isnan(t)) {
        return;
    }
    
    weighted_least_squares(X, M, N, W, s, x, no_threads, max_iteration);
    
    double sigma = mad(X, M, N, W, s, x);
    if (sigma < EPSILON) {
        return;
    }
    
    double* old_s = new double[M];
    double* old_x = new double[N];
    
    for (int iteration = 0; iteration < max_iteration; iteration++) {
        std::copy(s, s + M, old_s);
        std::copy(x, x + N, old_x);
        
        //standardize residuals & apply Huber weight
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                if (!isnan(X[j * M + i])) {
                    W[j * M + i] = (X[j * M + i] - s[i] - x[j]) / sigma;
                    
                    if (abs(W[j * M + i]) > k) {
                        W[j * M + i] = k / abs(W[j * M + i]);
                    }
                    else {
                        W[j * M + i] = 1;
                    }
                }
                else {
                    W[j * M + i] = NAN;
                }
            }
        }
        
        // solve
        weighted_least_squares(X, M, N, W, s, x, no_threads, max_iteration);
        
        sigma = mad(X, M, N, W, s, x);
        
        if (sigma < EPSILON) {
            break;
        }
        
        bool stop = true;
        for (int i = 0; i < M; i++) {
            if (abs(s[i] - old_s[i]) > EPSILON) {
                stop = false;
            }
        }
        for (int j = 0; j < N; j++) {
            if (abs(x[j] - old_x[j]) > EPSILON) {
                stop = false;
            }
        }
        if (stop) {
            break;
        }
    }
    
    delete [] old_x;
    delete [] old_s;
    delete [] W;
    
}

void weighted_huber(double *X, const int M, const int N, double *s, double *x, 
                    double* w,
                    int no_threads = 1,
                    const double k = 1.345, int max_iteration = 1000) {
    // X column-based matrix with M rows, N columns
    
    fill_n(s, M, NAN);
    fill_n(x, N, NAN);
    
    double *W = new double[M*N];
    fill_n(W, M*N, 1.0);
    
    
    double t = NAN;
    
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            if (!isnan(X[j * M + i])) {
                s[i] = 0;
                x[j] = 0;
                t = 0;
            }
        }
    }
    
    if (isnan(t)) {
        return;
    }
    
    weighted_least_squares(X, M, N, W, s, x, no_threads, max_iteration);
    
    double sigma = mad(X, M, N, W, s, x);
    if (sigma < EPSILON) {
        return;
    }
    
    double* old_s = new double[M];
    double* old_x = new double[N];
    
    for (int iteration = 0; iteration < max_iteration; iteration++) {
        std::copy(s, s + M, old_s);
        std::copy(x, x + N, old_x);
        
        //standardize residuals & apply Huber weight
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                if (!isnan(X[j * M + i])) {
                    W[j * M + i] = (X[j * M + i] - s[i] - x[j]) / sigma;
                    
                    if (abs(W[j * M + i]) > k) {
                        W[j * M + i] = w[i] * k / abs(W[j * M + i]);
                    }
                    else {
                        W[j * M + i] = w[i];
                    }
                }
                else {
                    W[j * M + i] = NAN;
                }
            }
        }
        
        // solve
        weighted_least_squares(X, M, N, W, s, x, no_threads, max_iteration);
        
        sigma = mad(X, M, N, W, s, x);
        
        if (sigma < EPSILON) {
            break;
        }
        
        bool stop = true;
        for (int i = 0; i < M; i++) {
            if (abs(s[i] - old_s[i]) > EPSILON) {
                stop = false;
            }
        }
        for (int j = 0; j < N; j++) {
            if (abs(x[j] - old_x[j]) > EPSILON) {
                stop = false;
            }
        }
        if (stop) {
            break;
        }
    }
    
    delete [] old_x;
    delete [] old_s;
    delete [] W;
    
}


#endif  // __CDW_H
