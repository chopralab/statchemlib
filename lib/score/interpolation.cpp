/* This is interpolation.cpp and is part of StatChemLIB
 * Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

#include "statchem/score/interpolation.hpp"
#include <assert.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_vector.h>
#include <functional>
#include "statchem/fileio/inout.hpp"
#include "statchem/helper/help.hpp"
using namespace std;

namespace statchem {

namespace score {
namespace Interpolation {
vector<double> derivative(const vector<double>& y, const double step) {
    assert(y.size() > 0);
    vector<double> deriva;

    deriva.push_back((y[1] - y[0]) / step);  // first

    for (size_t i = 1; i < y.size() - 1; ++i) {
        deriva.push_back((y[i + 1] - y[i - 1]) /
                         (2 * step));  // symmetric difference quotient
    }

    deriva.push_back(0);  // last

    return deriva;
}

/**
 * Interpolate through every point
 *
 */

vector<double> interpolate(const vector<double>& dataX,
                           const vector<double>& dataY, const double step) {
    assert(dataX.size() > 0 && dataX.size() == dataY.size());
    vector<double> pot;

    gsl_interp_accel* acc = gsl_interp_accel_alloc();
    gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, dataX.size());
    gsl_spline_init(spline, &dataX[0], &dataY[0], dataX.size());

    for (double xi = dataX.front(); xi <= dataX.back(); xi += step) {
        dbgmsg("x = " << xi << " yeval = " << gsl_spline_eval(spline, xi, acc));
        pot.push_back(gsl_spline_eval(spline, xi, acc));
    }

    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);

    return pot;
}

struct BSplineFitWorkspace {
    gsl_bspline_workspace* bw;
    gsl_vector* B;
    gsl_vector* c;
    gsl_vector* w;
    gsl_vector* x;
    gsl_vector* y;
    gsl_matrix* X;
    gsl_matrix* cov;
    gsl_multifit_linear_workspace* mw;
    ~BSplineFitWorkspace() {
        gsl_bspline_free(bw);
        gsl_vector_free(B);
        gsl_vector_free(x);
        gsl_vector_free(y);
        gsl_matrix_free(X);
        gsl_vector_free(c);
        gsl_vector_free(w);
        gsl_matrix_free(cov);
        gsl_multifit_linear_free(mw);
    }
};

/**
 * B-spline interpolation to get smooth curve
 * ------------------------------------------
 * k : number of data points to fit
 * ncoeff : number of fit coefficients
 */

BSplineFit::BSplineFit(const vector<double>& dataX, const vector<double>& dataY,
                       const size_t k, const size_t ncoeffs)
    : n(dataX.size()),
      ncoeffs2((ncoeffs < n) ? ncoeffs : n),
      nbreak(ncoeffs2 + 2 - k),  // nbreak = ncoeffs2 + 2 - k
      bsplinewsp(new BSplineFitWorkspace) {
    assert(dataX.size() > 0 && dataX.size() == dataY.size());

    bsplinewsp->bw = gsl_bspline_alloc(
        k, nbreak);  // allocate a cubic bspline workspace (k = 4)
    bsplinewsp->B = gsl_vector_alloc(ncoeffs2);
    bsplinewsp->c = gsl_vector_alloc(ncoeffs2);
    bsplinewsp->w = gsl_vector_alloc(n);
    bsplinewsp->x = gsl_vector_alloc(n);
    bsplinewsp->y = gsl_vector_alloc(n);
    bsplinewsp->X = gsl_matrix_alloc(n, ncoeffs2);
    bsplinewsp->cov = gsl_matrix_alloc(ncoeffs2, ncoeffs2);
    bsplinewsp->mw = gsl_multifit_linear_alloc(n, ncoeffs2);

    double chisq;

    for (size_t i = 0; i < dataX.size(); ++i) {
        const double xi = dataX[i];
        const double yi = dataY[i];
        const double sigma = 0.1 * yi;

        gsl_vector_set(bsplinewsp->x, i, xi);
        gsl_vector_set(bsplinewsp->y, i, yi);
        gsl_vector_set(bsplinewsp->w, i, 1.0 / (sigma * sigma));
    }

    /* use uniform breakpoints on [0, 15] */
    gsl_bspline_knots_uniform(dataX.front(), dataX.back(), bsplinewsp->bw);

    /* construct the fit matrix X */
    for (size_t i = 0; i < n; ++i) {
        double xi = gsl_vector_get(bsplinewsp->x, i);

        /* compute B_j(xi) for all j */
        gsl_bspline_eval(xi, bsplinewsp->B, bsplinewsp->bw);

        /* fill in row i of X */
        for (size_t j = 0; j < ncoeffs2; ++j) {
            double Bj = gsl_vector_get(bsplinewsp->B, j);
            gsl_matrix_set(bsplinewsp->X, i, j, Bj);
        }
    }

    /* do the fit */
    gsl_multifit_wlinear(bsplinewsp->X, bsplinewsp->w, bsplinewsp->y,
                         bsplinewsp->c, bsplinewsp->cov, &chisq,
                         bsplinewsp->mw);
#ifndef STATCHEM_DEBUG_MESSAGES
    gsl_stats_wtss(bsplinewsp->w->data, 1, bsplinewsp->y->data, 1,
                   bsplinewsp->y->size);
#else
    double tss = gsl_stats_wtss(bsplinewsp->w->data, 1, bsplinewsp->y->data, 1,
                                bsplinewsp->y->size);
    dbgmsg("chisq/dof = " << chisq / (n - ncoeffs)
                          << ", Rsq = " << 1.0 - chisq / tss << "\n");
#endif
}

vector<double> BSplineFit::interpolate_bspline(const double start,
                                               const double end,
                                               const double step) {
    vector<double> pot;

    /* output the smoothed curve */
    for (double xi = start; xi <= end; xi += step) {
        double yi, yerr;
        gsl_bspline_eval(xi, bsplinewsp->B, bsplinewsp->bw);
        gsl_multifit_linear_est(bsplinewsp->B, bsplinewsp->c, bsplinewsp->cov,
                                &yi, &yerr);

        pot.push_back(yi);
    }

    return pot;
}

BSplineFit::~BSplineFit() { delete bsplinewsp; }
}  // namespace Interpolation
}  // namespace score
}  // namespace statchem
