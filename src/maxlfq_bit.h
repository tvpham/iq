/*#########################################################################
 
 Author: Thang V. Pham, t.pham@amsterdamumc.nl, Pham Huy Chau Long
 
 All rights reserved.
 
 Citation:
 
 Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative
 protein abundances from ion quantification in DIA-MS-based proteomics,
 Bioinformatics 2020 Apr 15;36(8):2611-2613.
 
 Software version: 2.0
 
#########################################################################*/

#ifndef __MAXLFQ1_H
#define __MAXLFQ1_H

#include <algorithm>

#include "utils.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static constexpr double EPS = 1e-10;  // Convergence threshold for iterative solver


void Ax(double* A, int M, int N, double* x, double* out) {

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < (size_t)N; i++) {
        
        int nn = 0;      
        double vv = 0.0; 

        for (size_t j = 0; j < (size_t)N; j++) {
            if (j != i) {

                bool has_ratio = false;

                double* ip = A + M * i;
                double* jp = A + M * j;
                for (size_t r = 0; r < (size_t)M; r++) {
                    if (!std::isnan(*ip) && !std::isnan(*jp)) {
                        has_ratio = true;
                        break;
                    }
                    ip++;
                    jp++;
                }

                if (has_ratio) {
                    nn++;
                    vv -= x[j]; // Accumulate negative contributions from connected samples
                }
            }
        }

        out[i] = 2.0 * (vv + (double)nn * x[i]) + x[N];
        
    }

    out[N] = std::accumulate(x, x + N, 0.0);
}


void Ax(double* AtA, int n1, double* x, double* out) {
    
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < (size_t)n1; i++) {
        double* r = AtA + i*n1;
        out[i] = std::inner_product(r, r + n1, x, 0.0);
    }
}


void Ax(unsigned char** AtA, int* d_AtA, int n, double* x, double* out) {
    
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < (size_t)n; i++) {
        
        out[i] = x[i] * (double)d_AtA[i];
        
        unsigned char b = 0;
        unsigned char *pv = AtA[i];
        unsigned char v = *pv;
        
        for (size_t j = 0; j < (size_t)n; j++) {
            if (v & 1) {
                out[i] = out[i] - x[j];
            }
            if (b == 7) {
                pv++;
                v = *pv;
                b = 0;
            }
            else {
                b++;
                v = v >> 1;
            }
        }
        out[i] = 2.0 * out[i] + x[n];
    }
    
    double sum = 0.0;
    
    #pragma omp parallel for reduction(+:sum)
    for (size_t i = 0; i < (size_t)n; i++) {
        sum += x[i];
    }
    out[n] = sum;
    
}


void solve(double* X, int M, int N,
           double* b,
           int n1,
           double* x,
           double* r,
           double* p,
           double* Ap,
           void* AtA,
           int* d_AtA,
           int max_iter = -1) {

    if (AtA == NULL) {
        Ax(X, M, N, x, r);
    }
    else {
        if (d_AtA == NULL) {
            Ax((double*)AtA, N+1, x, r);
        }
        else {
            Ax((unsigned char **)AtA, d_AtA, N, x, r);
        }
    }

    for (size_t i = 0; i < (size_t)n1; ++i) {
        r[i] = b[i] - r[i];
    }

    std::copy_n(r, n1, p);

    double rsold = std::inner_product(r, r + n1, r, 0.0);

    int n_rounds = max_iter < 0 ? n1 : max_iter;

    for (int iter = 0; iter < n_rounds && rsold > EPS; iter++) {

        if (AtA == NULL) {
            Ax(X, M, N, p, Ap);  // Matrix-vector multiply: Ap = A*p
        }
        else {
            if (d_AtA == NULL) {
                Ax((double*)AtA, N+1, p, Ap);
            }
            else {
                Ax((unsigned char **)AtA, d_AtA, N, p, Ap);
            }
        }
        
        double pAp = std::inner_product(p, p + n1, Ap, 0.0);

        if (std::abs(pAp) < EPS) {
            break;
        }

        double alpha = rsold / pAp;  

        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < (size_t)n1; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }
        
        double rsnew = std::inner_product(r, r + n1, r, 0.0);

        // Convergence check
        if (rsnew < EPS) break;

        double beta = rsnew / rsold; 

        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < (size_t)n1; i++) {
            p[i] = r[i] + beta * p[i];
        }
        rsold = rsnew;
    }
}


void maxlfq_bit(double* X, int M, int N, double* x, int no_threads = 1, int memory_level = 0) {

    if (M == 1) {
        for (int i = 0; i < N; i++) {
            x[i] = X[i];
        }
        return;
    }

    double* Atb = new double[N + 1];
    std::fill_n(Atb, (size_t)N+1, 0.0);

    // Start with zero initial guess
    double* x0 = new double[(size_t)N + 1];
    std::fill_n(x0, N+1, 0.0);

    double* r = new double[(size_t)N + 1];
    double* p = new double[(size_t)N + 1];
    double* Ap = new double[(size_t)N + 1];

    double *big_median_buffer = new double[M * no_threads];

    // default if memory_level < 0
    void* AtA = NULL;
    int* d_AtA = NULL;

    if (memory_level > 0) {
        
        AtA = malloc((size_t)(N + 1) * (size_t)(N + 1) * (size_t)sizeof(double));
        
        if (AtA == NULL) {
            //std::cerr << "Not enough memory, lower memory setting activated.\n";
        }
        else {
            std::fill_n((double *)AtA, (size_t)(N + 1) * (size_t)(N + 1), 0.0);
        }
    }
    
    if (memory_level == 0 || (memory_level > 0 && AtA == NULL)) {
        
        AtA = (unsigned char **)malloc((size_t)N  * (size_t)sizeof(unsigned char *));
        d_AtA = (int *)malloc((size_t)N  * (size_t)sizeof(int));

        if (AtA == NULL || d_AtA == NULL) {
            //std::cerr << "Not enough memory, lowest memory setting activated.\n";
            if (AtA != NULL) {
                free(AtA);
                AtA = NULL;
            }
            if (d_AtA != NULL) {
                free(d_AtA);
                d_AtA = NULL;
            }
        }
        else {
            std::fill_n(d_AtA, N, 0);
            
            for (size_t i = 0; i < (size_t)N; i++) {
                ((unsigned char **)AtA)[i] = (unsigned char *)malloc((size_t)(N / 8 + 1)  * sizeof(unsigned char));
                if (((unsigned char **)AtA)[i] == NULL) {
                    //std::cerr << "Not enough memory, lowest memory setting activated.\n";
                    for (size_t j = 0; j < i; j++) {
                        free(((unsigned char **)AtA)[j]);
                    }
                    free(AtA);
                    AtA = NULL;
                    
                    if (d_AtA != NULL) {
                        free(d_AtA);
                        d_AtA = NULL;
                    }
                    break;
                }
                std::fill_n(((unsigned char **)AtA)[i], (size_t)(N / 8 + 1), 0);
            }
        }
    }
    
    
    unsigned char bit_mask[8] = {1, 2, 4, 8, 16, 32, 64, 128};
    

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < (size_t)N - 1; i++) {
        
        int thread_id = 0;
        #ifdef _OPENMP
            thread_id = omp_get_thread_num();
        #endif

        double *median_buffer = big_median_buffer + M * thread_id;

        size_t i1 = i >> 3;
        size_t i2 = i & 7;
        
        for (size_t j = i + 1; j < (size_t)N; j++) {

            size_t median_size = 0;

            double* ip = X + M * i;
            double* jp = X + M * j;
            for (size_t r = 0; r < (size_t)M; r++) {
                if (!std::isnan(*ip) && !std::isnan(*jp)) {
                    median_buffer[median_size++] = *jp - *ip;
                }
                ip++;
                jp++;
            }

            if (median_size > 0) {
                double r_i_j;
                
                r_i_j = quick_median(median_buffer, median_size);

                #pragma omp atomic
                Atb[i] -= r_i_j;

                #pragma omp atomic
                Atb[j] += r_i_j;

                if (AtA != NULL) {
                    
                    if (d_AtA != NULL) {
                        size_t j1 = j >> 3;
                        size_t j2 = j & 7;
                        
                        #pragma omp atomic
                        ((unsigned char **)AtA)[i][j1] = ((unsigned char **)AtA)[i][j1] | bit_mask[j2];
                    
                        #pragma omp atomic
                        ((unsigned char **)AtA)[j][i1] = ((unsigned char **)AtA)[j][i1] | bit_mask[i2];
                    
                        #pragma omp atomic
                        d_AtA[i]++;
                    
                        #pragma omp atomic
                        d_AtA[j]++;
                    }
                    else {
                        ((double*)AtA)[i * (size_t)(N+1) + j] = -1;
                        
                        ((double*)AtA)[j * (size_t)(N+1) + i] = -1;
                        
                        #pragma omp atomic
                        ((double*)AtA)[i * (size_t)(N+1) + i]++;
                        
                        #pragma omp atomic
                        ((double*)AtA)[j * (size_t)(N+1) + j]++;
                    }
                }
            }
        }
    }

    delete [] big_median_buffer;

    if (d_AtA == NULL && AtA != NULL) {
        for (size_t i = 0; i < (size_t)N; i++) {
            for (size_t j = 0; j < (size_t)N; j++) {
                ((double*)AtA)[i * (size_t)(N+1) + j] *= 2;
            }

            ((double*)AtA)[(size_t)N * (size_t)(N+1) + i] = 1;
            ((double*)AtA)[(size_t)i * (size_t)(N+1) + (size_t)N] = 1;
        }
    }
    
    // Set up augmented right-hand side: [2*Atb; global_constraint]
    for (size_t j = 0; j < (size_t)N; j++) {
        Atb[j] *= 2.0;
    }

    solve(X, M, N, Atb, N+1, x0, r, p, Ap, AtA, d_AtA);

    rescale(X, M, N, x0, x);

    if (d_AtA != NULL) {
        for (size_t i = 0; i < (size_t)N; i++) {
            free(((unsigned char **)AtA)[i]);
        }
        free(AtA);
        free(d_AtA);
    }
    else {
        if (AtA != NULL) {
            free(AtA);
        }
    }

    delete [] Ap;
    delete [] p;
    delete [] r;
    delete [] x0;
    delete [] Atb;
}

#endif
