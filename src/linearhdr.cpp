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
    return 1 / (1+ std::exp(-b * min(1-x, x) + m));
}

int linear_Response(pfs::Array2D *out[],
                   const ExposureList *imgs[],
                   const float opt_saturation_offset,
                   const float opt_black_offset,
                   float deghosting_value,
                   const float scale,
                   const float rgb_corr[3][3],
                   const float oor_high,
                   float oor_low,
                   bool isbayer,
                   const bool demosaic){

    isbayer = isbayer or demosaic;
    // number of exposures
    int N = imgs[0]->size();

    // frame size
    int width = out[0]->getCols();
    int height = out[0]->getRows();

    // number of saturated pixels
    int saturated_pixels = 0;
    int under_pixels = 0;

    float mmax[3] = {1e-30, 1e-30, 1e-30};
    float mmin[3] = {1e30, 1e30, 1e30};

    // For each pixel
    int skipped_deghost = 0;
    for (int j = 0; j < width * height; j++) {
        bool saturated[3] = {false, false, false};
        bool under[3] = {false, false, false};

        float X[3][N];
        float w[3][N];
        bool saturated_exp[N];
        bool under_exp[N];
        bool skip_exp[N];
        fill(saturated_exp, saturated_exp + N, false);
        fill(under_exp, under_exp + N, false);
        fill(skip_exp, skip_exp + N, false);

        // track reference frame for deghosting
        int ref = 0;
        bool foundref = false;

        // First compute scene radiances
        for (int i = 0; i < N; i++) {
            for (int cc = 0; cc < 3; cc++) {
                const Exposure &expo = (*imgs[cc])[i];
                X[cc][i] = (*expo.yi)(j) * get_exposure_compensationX(expo) * scale;
                w[cc][i] = get_weight((*expo.yi)(j));
                saturated_exp[i] = (*expo.yi)(j) >= 1 - opt_saturation_offset || saturated_exp[i];
                // when isbayer, this value tracks not under-exposed
                if (isbayer)
                    under_exp[i] = (*expo.yi)(j) > opt_black_offset || under_exp[i];
                else
                    under_exp[i] = (*expo.yi)(j) <= opt_black_offset || under_exp[i];
            }
            under_exp[i] = under_exp[i] != isbayer;
            if (saturated_exp[i] || under_exp[i]){
                w[0][i] = 0;
                w[1][i] = 0;
                w[2][i] = 0;
            } else if (!foundref){
                // use first (slowest) exposure (least noise as reference frame for deghosting)
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
                        //use absolute deviation
                        if (deghosting_value >= 1)
                            deviation_from_ref = max(deviation_from_ref,
                                                     fabs((float) (X[cc][i] - X[cc][ref])));
                        else
                            deviation_from_ref = max(deviation_from_ref,
                                                     fabs((float) (1 - X[cc][i] / X[cc][ref])));
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

            // For each exposure
            for (int i = 0; i < N; i++) {
                all_under &= under_exp[i];
                all_over &= saturated_exp[i];
                if (!(saturated_exp[i] || under_exp[i] || (deghosting_value != -1 && skip_exp[i]))) {
                    sum += X[cc][i] * w[cc][i];
                    div += w[cc][i];
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
                mmax[cc] = (mmax[cc] > (*out[cc])(j)) ? mmax[cc] : (*out[cc])(j);
                if (!demosaic || (*out[cc])(j) > 0)
                    mmin[cc] = (mmin[cc] < (*out[cc])(j)) ? mmin[cc] : (*out[cc])(j);
            }
        }

        if (under[0] || under[1] || under[2]) {
            under_pixels++;
        }
        if (saturated[0] || saturated[1] || saturated[2]) {
            saturated_pixels++;
        }

    }
    if (oor_high >= 0)
        mmax[0] = mmax[1] = mmax[2] = oor_high;
    if (oor_low >= 0)
        mmin[0] = mmin[1] = mmin[2] = oor_low;
    // Fill in nan values NOTE: removed normalization here
    float x,y,z;

    int g0 = first_non_zero(out[1]);
    int r0 = first_non_zero_row(out[0]);

    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++)
            for (int cc = 0; cc < 3; cc++) {
                if ((*out[cc])(j, i) == -1)
                    (*out[cc])(j, i) = mmax[cc];
                else if ((*out[cc])(j, i) == -2) {
                    if (!demosaic || grid_color(i, j, g0, r0) == cc)
                        (*out[cc])(j, i) = mmin[cc];
                }
            }
    // demosaic after merge
    if (demosaic){
        dht_interpolate(out[0], out[1], out[2]);
    }
    // apply color transformation
    for (int j = 0; j < width * height; j++) {
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
                (float) skipped_deghost * 100.f / (float) (width * height * N) << "%" << std::endl;

    if (under_pixels > 0) {
        float perc = ceilf(100.0f * under_pixels / (width * height));
        VERBOSE_STR << "under-exposed pixels found in " << perc << "% of the image!" << endl;
    }
    return saturated_pixels;
}
