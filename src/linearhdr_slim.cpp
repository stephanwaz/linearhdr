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
#include <cmath>
#include <linearhdr.h>

using namespace std;


#define PROG_NAME "linearhdr_slim"

extern bool verbose; /* verbose should be declared for each standalone code */

inline float max(float a, float b) {
    return (a > b) ? a : b;
}

inline float min(float a, float b) {
    return (a < b) ? a : b;
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
float get_exposure_compensation(const Exposure &ex) {
    return  ( ex.exposure_time * ex.iso ) / (100.0f * ex.aperture * ex.aperture);
}

float get_weight(const double x, const double s){
    return 1.0  / (1+ std::exp(-(0.88-x-s)/0.012));
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

// merge a mosaic image considering adjacent pixels
// helps to detect sensor blooming
void merge(pfs::Array2D *out[],
           const ExposureList *imgs[],
           int i, int j, int r0, int g0,
           const float opt_black_offset,
           const float opt_saturation_offset,
           const float scale,
           float& mmax,
           float& mmin) {
    bool saturated_exp;
    bool under_exp;
    bool all_under = true;
    bool all_over = true;
    int N = imgs[0]->size();
    int width = out[0]->getCols();
    int height = out[0]->getRows();
    int k = j + i * width;
    float w = 0.0; // stores the weight for each frame
    float div = 0.0;
    float saturation;
    float under;
    int cc = grid_color(i, j, g0, r0);
    const int bloommap[4][9] = {{2, 1, 2, 1, 0, 1, 2, 1, 2},
                                {1, 0, 1, 2, 1, 2, 1, 0, 1},
                                {0, 1, 0, 1, 2, 1, 0, 1, 0},
                                {1, 2, 1, 0, 1, 0, 1, 2, 1}};
    int bm = cc;
    if (cc == 1 && ((i-r0) % 2) == 0)
        bm = 3;
    int bloomk[9];
    fill(bloomk, bloomk + 9, -1.0);
    for (int bi = -1; bi < 2; bi++)
        if (i+bi >= 0 && i+bi < height)
            for (int bj = -1; bj < 2; bj++) {
                if (j+bj >= 0 && j+bj < width)
                    bloomk[(bi + 1) * 3 + bj + 1] = (j + bj) + (i + bi) * width;
            }
    // loop over individual frames, calculating weights
    for (int n = 0; n < N; n++) {
        saturation = 0.0;
        under = 1.0;
        for (int b = 0; b < 9; b++) {
            if (bloomk[b] > -1) {
                saturation = max(saturation, (*(*imgs[bloommap[bm][b]])[n].yi)(bloomk[b]));
                under = min(under, (*(*imgs[bloommap[bm][b]])[n].yi)(bloomk[b]));
            }
        }
        saturated_exp = saturation >= 1 - opt_saturation_offset;
        under_exp = under <= opt_black_offset;
        all_under &= under_exp;
        all_over &= saturated_exp;
        if (!(saturated_exp || under_exp)) {
            // use worst weight from surrounding pixels
            w = min(get_weight(saturation, opt_saturation_offset),
                       get_weight(under, opt_saturation_offset));
            div += w  * get_exposure_compensation((*imgs[cc])[n]);
            (*out[cc])(k) += w * (*(*imgs[cc])[n].yi)(k);
        }
    }
    // after looping over all frames, now we can
    // flag oob values, correctly weight output and store min/max value
    if (all_under) {
        (*out[cc])(k) = -2;
    } else if (all_over) {
        (*out[cc])(k) = -1;
    } else {
        (*out[cc])(k) = scale * (*out[cc])(k) / div;
        // in this case directly store each channel max
        mmax = (mmax > (*out[cc])(k)) ? mmax : (*out[cc])(k);
        mmin = (mmin < (*out[cc])(k)) ? mmin : (*out[cc])(k);
    }
}

std::tuple<long, long> linear_response_slim(pfs::Array2D *out[],
                                       const ExposureList *imgs[],
                                       const float opt_saturation_offset,
                                       const float opt_black_offset,
                                       const float scale,
                                       const float rgb_corr[3][3],
                                       const float oor_high,
                                       float oor_low,
                                       const bool demosaic){

    // establish bayer grid coordinates
    int g0 = first_non_zero((*imgs[1])[0].yi);
    int r0 = first_non_zero_row((*imgs[0])[0].yi);

    VERBOSE_STR << "1st green column: " << g0 << " 1st red row: " << r0 << std::endl;

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


    int k; // used to ravel i,j pixel coordinates
    int ccf; // first color index in loop, adjusted at each pixel when isbayer

    // loop over image
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            k = j + i * width;
            for (int cc = 0; cc < 3; cc++)
                (*out[cc])(k) = 0;


            ccf = grid_color(i, j, g0, r0);
            merge(out, imgs, i, j, r0, g0, opt_black_offset, opt_saturation_offset, scale, mmax[ccf], mmin[ccf]);

        }
    } // end of first loop over pixels
    // just in case no values were in range
    for (float & m : mmin) {
        if (m > 1e29)
            m = 0;
    }

    // only correct bayer input if setting to user override
    for (int cc = 0; cc < 3; cc++) {
        if (oor_high >= 0)
            mmax[cc] = oor_high;
        if (oor_low >= 0)
            mmin[cc] = oor_low;
    }

    // fill in out of range pixels
    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++) {
            k = j + i * width;
            ccf = grid_color(i, j, g0, r0);
            if ((*out[ccf])(k) == -2) {
                under_pixels++;
                (*out[ccf])(k) = mmin[ccf];
            } else if ((*out[ccf])(k) == -1) {
                saturated_pixels++;
                (*out[ccf])(k) = mmax[ccf];
            }
        }

    // demosaic after merge
    if (demosaic){
        dht_interpolate(out[0], out[1], out[2]);
        // in this case data not yet mapped to output colorspace, do this now
        for (int j = 0; j < NP; j++)
            apply_color_transform(j, out, rgb_corr);
    }

    return {saturated_pixels, under_pixels};
}
