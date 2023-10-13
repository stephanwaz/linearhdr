/**
 *
 * simplified hdr generation assuming raw linear response, fixes possible bugs in old code that lost absolute
 * calibration
 * @author Stephen Wasilewski stephen.wasilewski@epfl.ch
 *
 * Forked from: version pfstools 2.2.0:
 * @brief Robertson02 algorithm for automatic self-calibration.
 *
 * 
 * This file is derived from a part of PFS CALIBRATION package.
 * ---------------------------------------------------------------------- 
 * Copyright (C) 2004 Grzegorz Krawczyk
 * 
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * ---------------------------------------------------------------------- 
 * 
 * @author Grzegorz Krawczyk, <gkrawczyk@users.sourceforge.net>
 *         Ivo Ihrke        , <ihrke@mmci.uni-saarland.de>
 *
 * $Id: robertson02.cpp,v 1.11 2011/02/25 13:45:14 ihrke Exp $
 */

/*
 * Copyright (c) 2023 Stephen Wasilewski, EPFL
 *  =======================================================================
 *  This program is free software: you can redistribute it and/or
 *  modify it under the terms of theGNU Lesser General Public License
 *  as published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <https://www.gnu.org/licenses/>.
 *  =======================================================================
 */


#include <config.h>

#include <iostream>
#include <vector>

#include <cmath>
#include <linearhdr.h>

using namespace std;


#define PROG_NAME "linearhdr"

extern bool verbose; /* verbose should be declared for each standalone code */

inline float max(float a, float b) {
    return (a > b) ? a : b;
}

int grid_color(int i, int j, int g0, int r0) {
    // is green
    if (j % 2 == ((g0 + i) % 2))
        return 1;
    // is red row
    if (i % 2 == r0)
        return 0;
    // is blue row
    return 2;
}

float get_exposure_compensationX(const Exposure &ex) {
    return  100.0f * ex.aperture * ex.aperture / ( ex.exposure_time * ex.iso );
}

float get_weight(const double x, const double b = 25, const double m = 5){
    // w[i] = eG.exposure_time;
    // y\ =\ \frac{1}{1+e^{-b\cdot\min\left(1-x,x\right)+m}}
    return x / (1+ std::exp(-b * min(1-x, x) + m));
}

int linear_Response(pfs::Array2D *out[],
                   const ExposureList *imgs[],
                   const float opt_saturation_offset,
                   const float opt_black_offset,
                   float deghosting_value,
                   const float scale,
                   const float vlambda[3],
                   const float rgb_corr[3][3],
                   const float oor_high,
                   float oor_low,
                   bool isbayer,
                   const bool demosaic){

    isbayer = isbayer or demosaic;
    int g0 = first_non_zero(out[1]);
    int r0 = first_non_zero_row(out[0]);

    // number of exposures
    int N = imgs[0]->size();

    // frame size
    int width = out[0]->getCols();
    int height = out[0]->getRows();

    int NP = width * height;

    // number of saturated pixels
    int saturated_pixels = 0;
    int under_pixels = 0;

    float mmax[3] = {1e-30, 1e-30, 1e-30};
    float mmin[3] = {1e30, 1e30, 1e30};

    float aavg[3] = {0.0, 0.0, 0.0};
    float ravg[3] = {0.0, 0.0, 0.0};

    int skipped_deghost = 0;
    float lummax = 1e-30;
    float lummin = 1e30;
    float lum;
    int k;
    int ccf = 0;
    int cinc = 1;
    if (isbayer)
        cinc = 3;

    int mi = 0;
    int mj = 0;

    // loop over each pixel
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            k = j + i * width;

            float X[3][N];
            float w[N];
            bool saturated_exp[N];
            bool under_exp[N];
            bool skip_exp[N];
            fill(w, w + N, 1.0);
            fill(saturated_exp, saturated_exp + N, false);
            fill(under_exp, under_exp + N, false);
            fill(skip_exp, skip_exp + N, false);

            // track reference frame for deghosting
            int ref = 0;
            bool foundref = false;

            // only operate on relevant color when isbayer (modifies function of all cc loops)
            if (isbayer) {
                ccf = grid_color(i, j, g0, r0);
            }

            // First compute scene values, weights and reference frames for deghosting
            for (int n = 0; n < N; n++) {
                for (int cc = ccf; cc < 3; cc += cinc) {
                    X[cc][n] = (*(*imgs[cc])[n].yi)(k) * get_exposure_compensationX((*imgs[cc])[n]) * scale;
                    w[n] = min(w[n], get_weight((*(*imgs[cc])[n].yi)(k)));
                    saturated_exp[n] = ((*(*imgs[cc])[n].yi)(k) >= 1 - opt_saturation_offset) || saturated_exp[n];
                    under_exp[n] = ((*(*imgs[cc])[n].yi)(k) <= opt_black_offset) || under_exp[n];
                }
                if (!foundref && (!saturated_exp[n] && !under_exp[n])) {
                    // use first (slowest) exposure (least noise as reference frame for deghosting)
                    ref = n;
                    foundref = true;
                }
            }

            // flag frames to skip for deghosting
            if (deghosting_value > 0)
                for (int n = 0; n < N; n++)
                    if (n != ref) {
                        float deviation_from_ref = 0;
                        for (int cc = ccf; cc < 3; cc += cinc) {
                            //use absolute deviation
                            if (deghosting_value >= 1)
                                deviation_from_ref = max(deviation_from_ref, fabs((float) (X[cc][n] - X[cc][ref])));
                            else
                                deviation_from_ref = max(deviation_from_ref, fabs((float) (1 - X[cc][n] / X[cc][ref])));
                        }
                        skip_exp[n] = deviation_from_ref > deghosting_value;
                        if (skip_exp[n])
                            skipped_deghost++;
                    }

            float div = 0.0f;
            bool all_under = true;
            bool all_over = true;

            // loop over exposures to compute combined value
            for (int n = 0; n < N; n++) {
                all_under &= under_exp[n];
                all_over &= saturated_exp[n];
                if (!(saturated_exp[n] || under_exp[n] || skip_exp[n])) {
                    div += w[n];
                    for (int cc = ccf; cc < 3; cc += cinc) {
                        (*out[cc])(k) += X[cc][n] * w[n];
                    }
                }
            }

            // flag oob values, update output and store min/max value
            for (int cc = ccf; cc < 3; cc += cinc) {
                if (all_under) {
                    (*out[cc])(k) = -2;
                    under_pixels++;
                } else if (all_over) {
                    (*out[cc])(k) = -1;
                    saturated_pixels++;
                } else {
                    (*out[cc])(k) = (*out[cc])(k) / div;
                    if (isbayer) {
                        mmax[cc] = (mmax[cc] > (*out[cc])(k)) ? mmax[cc] : (*out[cc])(k);
                        mmin[cc] = (mmin[cc] < (*out[cc])(k)) ? mmin[cc] : (*out[cc])(k);
                    }
                }
            }

            // store min/max luminance and running color averages
            if (!isbayer && (!all_over && !all_under)) {
                lum = (*out[0])(k) * vlambda[0] + (*out[1])(k) * vlambda[1] + (*out[2])(k) * vlambda[2];
                lummin = (lum < lummin) ? lum : lummin;
                lummax = (lum > lummax) ? lum : lummax;
                for (int cc = ccf; cc < 3; cc += cinc) {
                    aavg[cc] += (*out[cc])(k) / NP;
                    ravg[cc] += ((*out[cc])(k) / lum) / NP;
                    mmin[cc] = (lum < lummin) ? (*out[cc])(k) : mmin[cc];
                    mmax[cc] = (lum < lummax) ? (*out[cc])(k) : mmax[cc];
                }
            }
        }
    }

    if (oor_high >= 0)
        lummax = oor_high;
    if (oor_low >= 0)
        lummin = oor_low;

    // correct min/max values to luminance (avoids magenta highlights for out of bounds values)
    if (!isbayer) {
        for (int cc = ccf; cc < 3; cc += cinc) {
            aavg[cc] = aavg[cc] * NP / (NP - under_pixels - saturated_pixels);
            ravg[cc] = ravg[cc] * NP / (NP - under_pixels - saturated_pixels);
        }
        int alum = aavg[0] * vlambda[0] + aavg[1] * vlambda[1] + aavg[2] * vlambda[2];
        int rlum = ravg[0] * vlambda[0] + ravg[1] * vlambda[1] + ravg[2] * vlambda[2];


        for (int cc = ccf; cc < 3; cc += cinc) {
            mmin[cc] = ravg[cc] * lummin / rlum;
            mmax[cc] = aavg[cc] * lummax / alum;
        }

    }

    // Fill in out of range values
    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++) {
            k = j + i * width;
            // only operate on relevant color when isbayer
            if (isbayer)
                ccf = grid_color(i, j, g0, r0);
            for (int cc = ccf; cc < 3; cc += cinc) {
                (*out[cc])(k) = ((*out[cc])(k) == -1) ? mmax[cc] : (*out[cc])(k);
                (*out[cc])(k) = ((*out[cc])(k) == -2) ? mmin[cc] : (*out[cc])(k);
            }
    }
    // demosaic after merge
    if (demosaic){
        dht_interpolate(out[0], out[1], out[2]);
    }
    // apply color transformation
    float x,y,z;
    for (int j = 0; j < NP; j++) {
        x = (*out[0])(j);
        y = (*out[1])(j);
        z = (*out[2])(j);
        (*out[0])(j) = x * rgb_corr[0][0] + y * rgb_corr[0][1] + z * rgb_corr[0][2];
        (*out[1])(j) = x * rgb_corr[1][0] + y * rgb_corr[1][1] + z * rgb_corr[1][2];
        (*out[2])(j) = x * rgb_corr[2][0] + y * rgb_corr[2][1] + z * rgb_corr[2][2];
        (*out[0])(j) = max((*out[0])(j), 0);
        (*out[1])(j) = max((*out[1])(j), 0);
        (*out[2])(j) = max((*out[2])(j), 0);
    }

    VERBOSE_STR << "Maximum Value: " << mmax[0] << ", " << mmax[1] << ", " << mmax[2] << std::endl;
    VERBOSE_STR << "Exposure pixels skipped due to deghosting: " <<
                (float) skipped_deghost * 100.f / (float) (NP * N) << "%" << std::endl;

    if (under_pixels > 0) {
        float perc = ceilf(100.0f * under_pixels / NP);
        VERBOSE_STR << "under-exposed pixels found in " << perc << "% of the image!" << endl;
    }
    return saturated_pixels;
}
