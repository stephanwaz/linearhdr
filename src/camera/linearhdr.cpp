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
 * This file is a part of PFS CALIBRATION package.
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


#include <config.h>

#include <iostream>
#include <vector>

#include <cmath>

#include <responses.h>


using namespace std;


#define PROG_NAME "linearhdr"

extern bool verbose; /* verbose should be declared for each standalone code */

inline float max3(float a[3]) {
    float max = (a[0] > a[1]) ? a[0] : a[1];
    return (a[2] > max) ? a[2] : max;
}

inline float max(float a, float b) {
    return (a > b) ? a : b;
}


float get_exposure_compensationX(const Exposure &ex) {
    return  100.0f * ex.aperture * ex.aperture / ( ex.exposure_time * ex.iso );
}

float get_weight(const double x, const double b = 25, const double m = 5){
    // w[i] = eG.exposure_time;
    // y\ =\ \frac{1}{1+e^{-b\cdot\min\left(1-x,x\right)+m}}
    return 1 / (1+ std::exp(-b * min(1-x, x) + m));
}

int linear_Response(pfs::Array2D *out[],
                   const ExposureList *imgs[],
                   const float opt_saturation_offset,
                   const float opt_black_offset,
                   float deghosting_value,
                   const float scale,
                   const int lead_channel = -1,
                   const float oor_high = 1e-30,
                   const float oor_low = 1e30){

    // number of exposures
    int N = imgs[0]->size();

    // frame size
    int width = out[0]->getCols();
    int height = out[0]->getRows();

    // number of saturated pixels
    int saturated_pixels = 0;
    int under_pixels = 0;

    float mmax[3] = {oor_high, oor_high, oor_high};
    float mmin[3] = {oor_low, oor_low, oor_low};

    // For each pixel
    int skipped_deghost = 0;
    for (int j = 0; j < width * height; j++) {
        bool saturated[3] = {false, false, false};
        bool under[3] = {false, false, false};

        float X[3][N];
        float w[3][N];
        bool saturated_exp[3][N];
        bool under_exp[3][N];
        bool skip_exp[N];
        fill(skip_exp, skip_exp + N, false);

        // track reference frame for deghosting
        int ref = 0;
        bool foundref = false;

        // First compute scene radiances
        for (int i = 0; i < N; i++) {
            for (int cc = 0; cc < 3; cc++) {
                w[cc][i] = 0;
                const Exposure &expo = (*imgs[cc])[i];
                if (lead_channel < 0 || cc == lead_channel) {
                    X[cc][i] = (*expo.yi)(j) * get_exposure_compensationX(expo) * scale;
                    w[cc][i] += get_weight((*expo.yi)(j));

                    saturated_exp[cc][i] = (*expo.yi)(j) >= 1 - opt_saturation_offset;
                    under_exp[cc][i] = (*expo.yi)(j) <= opt_black_offset;
                    if (saturated_exp[cc][i] || under_exp[cc][i])
                        w[cc][i] = 0;
                    else
                        ref = i;
                } else {
                    X[cc][i] = (*expo.yi)(j);
                    saturated_exp[cc][i] = false;
                    under_exp[cc][i] = false;
                }

            }
            // use first (slowest) exposure (least noise as reference frame for deghosting)
            if (!foundref && !saturated_exp[0][i] && !saturated_exp[1][i] && !saturated_exp[2][i] &&
                    !under_exp[0][i] && !under_exp[1][i] && !under_exp[2][i]){
                ref = i;
                foundref = true;
            }

        }

        if (deghosting_value > 0) {
            // Compute the number of pixels to be skipped due to deghosting
            for (int i = 0; i < N; i++)
                if (i != ref) {

                    float deviation_from_ref = 0;
                    for (int cc = 0; cc < 3; cc++) {
                        //deghost only on luminance unless RGB
                        if (lead_channel < 0 || cc == lead_channel) {
                            //use absolute deviation
                            if (deghosting_value >= 1)
                                deviation_from_ref = max(deviation_from_ref,
                                                         fabs((float) (X[cc][i] - X[cc][ref])));
                            else
                                deviation_from_ref = max(deviation_from_ref,
                                                         fabs((float) (1 - X[cc][i] / X[cc][ref])));
                        }
                    }
                    skip_exp[i] = deviation_from_ref > deghosting_value;
                    if (skip_exp[i])
                        skipped_deghost++;
                }
        }

        for (int cc = 0; cc < 3; cc++) {

            float sum = 0.0f;
            float div = 0.0f;
            bool all_under = true;
            bool all_over = true;
            int pc = lead_channel;
            if (lead_channel < 0)
                pc = cc;
            // For each exposure
            for (int i = 0; i < N; i++) {
                all_under &= under_exp[pc][i];
                all_over &= saturated_exp[pc][i];
                if (!(saturated_exp[pc][i] || under_exp[pc][i] || (deghosting_value != -1 && skip_exp[i]))) {
                    sum += X[cc][i] * w[pc][i];
                    div += w[pc][i];
                }
            }
            if (all_under) {
                (*out[cc])(j) = -2;
                under[cc] = true;
            } else if (all_over) {
                (*out[cc])(j) = -1;
                saturated[cc] = true;
            } else if (div >= 1e-15) {
                (*out[cc])(j) = sum / div;
                mmax[cc] = (mmax[pc] > (*out[pc])(j)) ? mmax[cc] : (*out[cc])(j);
                mmin[cc] = (mmin[pc] < (*out[pc])(j)) ? mmin[cc] : (*out[cc])(j);
            }
        }

        if (under[0] || under[1] || under[2]) {
            under_pixels++;
        }
        if (saturated[0] || saturated[1] || saturated[2]) {
            saturated_pixels++;
        }

    }
    // Fill in nan values NOTE: removed normalization here
    for (int j = 0; j < width * height; j++)
        for (int cc = 0; cc < 3; cc++) {
            if ((*out[cc])(j) == -1)
                (*out[cc])(j) = mmax[cc];
            else if ((*out[cc])(j) == -2)
                (*out[cc])(j) = mmin[cc];
        }

    VERBOSE_STR << "Maximum Value: " << mmax[0] << ", " << mmax[1] << ", " << mmax[2] << std::endl;
    VERBOSE_STR << "Exposure pixels skipped due to deghosting: " <<
                (float) skipped_deghost * 100.f / (float) (width * height * N) << "%" << std::endl;

    if (under_pixels > 0) {
        float perc = ceilf(100.0f * under_pixels / (width * height));
        VERBOSE_STR << "under-exposed pixels found in " << perc << "% of the image!" << endl;
    }
    return saturated_pixels;
}
