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

// based on the first green column (0 or 1) and first red row (0 or 1) return color at pixel coordinate
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

// basic photographic reciprocity (assumes calibrated aperture/exposure/iso
float get_exposure_compensationX(const Exposure &ex) {
    return  100.0f * ex.aperture * ex.aperture / ( ex.exposure_time * ex.iso );
}

// Poisson Photon Noise Estimator(PPNE) with rounded peak to reduce
float get_weight(const Exposure &ex, const double x, const double s){
    // old hat function:
    // y\ =\ \frac{1}{1+e^{-b\cdot\min\left(1-x,x\right)+m}}
//    return ex.exposure_time / (1+ std::exp(-25 * min(1-x, x) + 5));
    // new PPNE with easing at top
    // y\ =\ \frac{\max\left(x-t,0\right)}{1+3e^{-\frac{\left(1-s-x\right)}{.02}}}
    return ex.exposure_time  / (1+ 3*std::exp(-(1-x-s)/0.02));
}

void apply_color_transform(int j, pfs::Array2D *out[], const float rgb_corr[3][3]) {
    float x, y, z;
    x = (*out[0])(j);
    y = (*out[1])(j);
    z = (*out[2])(j);
    (*out[0])(j) = max(0, x * rgb_corr[0][0] + y * rgb_corr[0][1] + z * rgb_corr[0][2]);
    (*out[1])(j) = max(0, x * rgb_corr[1][0] + y * rgb_corr[1][1] + z * rgb_corr[1][2]);
    (*out[2])(j) = max(0, x * rgb_corr[2][0] + y * rgb_corr[2][1] + z * rgb_corr[2][2]);
}

std::tuple<long, long> linear_response(pfs::Array2D *out[],
                                     const ExposureList *imgs[],
                                     const float opt_saturation_offset,
                                     const float opt_black_offset,
                                     const float scale,
                                     const float vlambda[3],
                                     const float rgb_corr[3][3],
                                     const float oor_high,
                                     float oor_low,
                                     bool isbayer,
                                     const bool demosaic){

    // governs primary merging, do we check all color channels, or merge channels independently
    isbayer = isbayer or demosaic;
    // establish bayer grid coordinates
    int g0 = first_non_zero((*imgs[1])[0].yi);
    int r0 = first_non_zero_row((*imgs[0])[0].yi);

    if (isbayer)
        VERBOSE_STR << "1st green column: " << g0 << " 1st red row: " << r0 << std::endl;

    // number of exposures
    int N = imgs[0]->size();

    // frame size
    int width = out[0]->getCols();
    int height = out[0]->getRows();

    int NP = width * height;

    // track out of bounds pixels
    int saturated_pixels = 0;
    int under_pixels = 0;

    // track maximum in range values (in target color space)
    float mmax[3] = {0,0,0};
    float mmin[3] = {1e30, 1e30, 1e30};
    float lummax = 0;
    float lummin = 1e30;
    float lum;
    bool saturated_exp; // track sat/under
    bool under_exp;

    // we will calculate the absolute and relative color averages to better match out of bounds values
    // in target color space)
    double aavg[3] = {0.0, 0.0, 0.0};
    double ravg[3] = {0.0, 0.0, 0.0};


    int k; // used to ravel i,j pixel coordinates
    int ccf = 0; // first color index in loop, adjusted at each pixel when isbayer
    int cinc = isbayer ? 3: 1; // whether color loop will increment through all colors or only first

    // loop over image
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            k = j + i * width;
            for (int cc = 0; cc < 3; cc++) {
                (*out[cc])(k) = 0;
            }

            float X[3][N]; // stores the value for each frame
            float w[N]; // stores the weight for each frame
            fill(w, w + N, 1.0);

            // only operate on relevant color when isbayer (modifies function of all cc loops)
            if (isbayer) {
                ccf = grid_color(i, j, g0, r0);
            }

            float div = 0.0f; // track sum of weights for averaging
            bool all_under = true;
            bool all_over = true;

            // loop over individual frames, calculating weights
            for (int n = 0; n < N; n++) {
                saturated_exp = false;
                under_exp = false;
                for (int cc = ccf; cc < 3; cc += cinc) {
                    X[cc][n] = (*(*imgs[cc])[n].yi)(k) * get_exposure_compensationX((*imgs[cc])[n]) * scale;
                    // weight by worst color channel
                    w[n] = min(w[n], get_weight((*imgs[cc])[n], (*(*imgs[cc])[n].yi)(k), opt_saturation_offset));
                    saturated_exp = ((*(*imgs[cc])[n].yi)(k) >= 1 - opt_saturation_offset) || saturated_exp;
                    under_exp = ((*(*imgs[cc])[n].yi)(k) <= opt_black_offset) || under_exp;
                }
                all_under &= under_exp;
                all_over &= saturated_exp;
                if (!(saturated_exp || under_exp)) {
                    div += w[n];
                    for (int cc = ccf; cc < 3; cc += cinc) {
                        (*out[cc])(k) += X[cc][n] * w[n];
                    }
                }
            }
            // just in case (can end up with a zero weight when finding worst color channel)
            // but this should mean it is under-exposed
            if (div <= 0 && !all_over)
                all_under = true;

            // after looping over all frames, now we can
            // flag oob values, correctly weight output and store min/max value
            for (int cc = ccf; cc < 3; cc += cinc) {
                if (all_under) {
                    (*out[cc])(k) = -2;
                    under_pixels++;
                } else if (all_over) {
                    (*out[cc])(k) = -1;
                    saturated_pixels++;
                } else {
                    (*out[cc])(k) = (*out[cc])(k) / div;
                    // in this case directly store each channel max
                    if (isbayer) {
                        mmax[cc] = (mmax[cc] > (*out[cc])(k)) ? mmax[cc] : (*out[cc])(k);
                        mmin[cc] = (mmin[cc] < (*out[cc])(k)) ? mmin[cc] : (*out[cc])(k);
                    }
                }
            }

            // when operating on a demosaiced images, we want to balance max values to luminance
            // (based on the target color space)
            // store min/max luminance and running color averages
            if (!isbayer && !(all_over || all_under)) {
                // apply color transformation
                apply_color_transform(k, out, rgb_corr);

                lum = (*out[0])(k) * vlambda[0] + (*out[1])(k) * vlambda[1] + (*out[2])(k) * vlambda[2];
                for (int cc = 0; cc < 3; cc += 1) {
                    aavg[cc] += (*out[cc])(k);
                    ravg[cc] += (*out[cc])(k) / lum;
                    mmin[cc] = (lum < lummin) ? (*out[cc])(k) : mmin[cc];
                    mmax[cc] = (lum > lummax) ? (*out[cc])(k) : mmax[cc];
                }
                lummin = (lum < lummin) ? lum : lummin;
                lummax = (lum > lummax) ? lum : lummax;
            }
        }
    } // end of first loop over pixels

    // apply manual overrides
    if (oor_high >= 0)
        lummax = oor_high;
    if (oor_low >= 0)
        lummin = oor_low;

    // just in case no values were in range
    if (lummin > 1e29) {
        lummin = 0;
    }
    for (int cc = 0; cc < 3; cc++) {
        if (mmin[cc] > 1e29)
            mmin[cc] = 0;
    }

    // correct min/max values to luminance (tends to avoid magenta highlights for out of bounds values)
    if (!isbayer) {
        // these get triple counted
        under_pixels = under_pixels / 3;
        saturated_pixels = saturated_pixels / 3;
        for (int cc = 0; cc < 3; cc += 1) {
            aavg[cc] = aavg[cc] / (NP - under_pixels - saturated_pixels);
            ravg[cc] = ravg[cc] / (NP - under_pixels - saturated_pixels);
        }
        float alum = aavg[0] * vlambda[0] + aavg[1] * vlambda[1] + aavg[2] * vlambda[2];
        float rlum = ravg[0] * vlambda[0] + ravg[1] * vlambda[1] + ravg[2] * vlambda[2];

        VERBOSE_STR << "Average Value: " << aavg[0] << ", " << aavg[1] << ", " << aavg[2] << ", " << alum << std::endl;
        VERBOSE_STR << "Relative Average Value: " << ravg[0] << ", " << ravg[1] << ", " << ravg[2] << ", " << rlum << std::endl;

        for (int cc = ccf; cc < 3; cc += cinc) {
            if (rlum > 0)
                mmin[cc] = ravg[cc] * lummin / rlum;
            if (alum > 0)
                mmax[cc] = aavg[cc] * lummax / alum;
        }

    } else {
        // only correct bayer input if setting to user override
        for (int cc = 0; cc < 3; cc++) {
            if (oor_high >= 0)
                mmax[cc] = oor_high;
            if (oor_low >= 0)
                mmin[cc] = oor_low;
        }
    }

    // fill in out of range pixels
    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++) {
            k = j + i * width;
            // only operate on relevant color when isbayer need to set ccf again in new loop
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
        // in this case data not yet mapped to output colorspace, do this now
        for (int j = 0; j < NP; j++)
            apply_color_transform(j, out, rgb_corr);
    }

    VERBOSE_STR << "Maximum Value: " << mmax[0] << ", " << mmax[1] << ", " << mmax[2] << ", " << lummax << std::endl;
    VERBOSE_STR << "Minimum Value: " << mmin[0] << ", " << mmin[1] << ", " << mmin[2] << ", " << lummin << std::endl;

    return {saturated_pixels, under_pixels};
}
