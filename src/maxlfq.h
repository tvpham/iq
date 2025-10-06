/*#########################################################################
 
 Author: Thang V. Pham, t.pham@amsterdamumc.nl
 
 All rights reserved.
 
 Citation:
 
 Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative
 protein abundances from ion quantification in DIA-MS-based proteomics,
 Bioinformatics 2020 Apr 15;36(8):2611-2613.
 
 Software version: 2.0
 
#########################################################################*/

#ifndef __MAXLFQ_H
#define __MAXLFQ_H

#include <RcppEigen.h>

#include <unordered_map>
#include <vector>

#include "utils.h"


using Eigen::FullPivHouseholderQR;
using Eigen::HouseholderQR;
using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;

int full_connection;
FullPivHouseholderQR<MatrixXd> full_qr;


void init_maxlfq(int N) {
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

/* MaxLFQ */
void spread(const unordered_map<int, vector<double>> &map, int *g, int i, int val) {
    int ncol = map.begin()->second.size();
    g[i] = val;

    for (const auto &m : map) {
        if (!isnan(m.second[i])) {
            for (int j = 0; j < ncol; j++) {
                if (g[j] < 0 && !isnan(m.second[j])) {
                    spread(map, g, j, val);
                }
            }
        }
    }
}

void maxlfq(const unordered_map<int, vector<double>> &map, double *buffer, int *g) {
    int ncol = map.begin()->second.size();

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
            spread(map, g, i, val++);
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
                buffer[ind[0]] = quick_median(median_buffer, median_size);
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

                        double r_i_j = quick_median(median_buffer, median_size);
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


void weighted_maxlfq(unordered_map<int, vector<double>> &map, double *buffer, int *g, double *w) {
    
    int ncol = map.begin()->second.size();
    
    if (map.size() == 1) {
        for (int i = 0; i < ncol; i++) {
            buffer[i] = map.begin()->second[i];
        }
        fill_n(g, ncol, 0);
        return;
    }
    
    int median_size = 0;
    std::vector<weighted_sample_t> points(map.size());
    
    fill_n(g, ncol, -1);
    
    int val = 0;
    for (int i = 0; i < ncol; i++) {
        if (g[i] < 0) {
            spread(map, g, i, val++);
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
            
            for (size_t j = 0; j < map.size(); j++) {
                if (!isnan(map[j][ind[0]])) {
                    points[median_size].value = map[j][ind[0]];
                    points[median_size].weight = w[j];
                    median_size++;
                }
            }
            if (median_size == 0) {
                buffer[ind[0]] = NAN;
            } else {
                buffer[ind[0]] = weighted_median(points, median_size);
            }
        } else {
            int N = ind.size();
            
            MatrixXd AtA = MatrixXd::Zero(N + 1, N + 1);
            VectorXd Atb = VectorXd::Zero(N + 1);
            
            int n_connection = 0;
            
            for (int j = 0; j < (N - 1); j++) {
                for (int k = j + 1; k < N; k++) {
                    median_size = 0;
                    for (size_t ell = 0; ell < map.size(); ell++) {
                        if (!isnan(map[ell][ind[j]]) && !isnan(map[ell][ind[k]])) {
                            points[median_size].value = map[ell][ind[k]] - map[ell][ind[j]];
                            points[median_size].weight = w[ell];
                            median_size++;
                        }
                    }
                    if (median_size > 0) {
                        n_connection++;
                        
                        double r_i_j = weighted_median(points, median_size);
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
}

#endif  // __MAXLFQ_H