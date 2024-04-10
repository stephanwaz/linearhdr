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

float get_weight(const double x, const double s, const double b){
    return std::exp(-b/x - 0.1/(1-s-x)) + 1e-9;
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
char merge(pfs::Array2D *out[],
           const ExposureList *imgs[],
           int i, int j, int r0, int g0,
           const float opt_black_offset,
           const float opt_saturation_offset,
           const float scale,
           const float wsp[3]) {
    int N = imgs[0]->size();
    int width = out[0]->getCols();
    int height = out[0]->getRows();
    int k = j + i * width;
    float w; // stores the weight for each frame
    float div = 0.0;
    for (int cc = 0; cc < 3; cc++) //make sure output is clear
        (*out[cc])(k) = 0;
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
    bool all_under = true;
    bool all_over = true;
    bool saturated_exp;
    bool under_exp;
    float saturation_b; // highest raw value among surrounding pixels of different color
    float saturation_a; // highest raw value among surrounding pixels of same color
    float under; // lowest raw value among surrounding pixels
    float high = 0.0; // max value (including oor)
    float low = 1e30; // min value (including oor)
    float ec;
    float pvalue;
    for (int n = 0; n < N; n++) {
        ec = get_exposure_compensation((*imgs[cc])[n]);
        saturation_a = saturation_b = 0.0;
        under = 1.0;
        // first check surrounding pixels for saturation
        for (int b = 0; b < 9; b++) {
            if (bloomk[b] > -1) {
                if (bloommap[bm][b] == cc)
                    saturation_a = max(saturation_a, (*(*imgs[cc])[n].yi)(bloomk[b]));
                else
                    saturation_b = max(saturation_b, (*(*imgs[bloommap[bm][b]])[n].yi)(bloomk[b]));
                under = min(under, (*(*imgs[bloommap[bm][b]])[n].yi)(bloomk[b]));
            }
        }
        saturated_exp = saturation_b >= 1 - opt_saturation_offset;
        under_exp = under <= opt_black_offset;

        pvalue = (*(*imgs[cc])[n].yi)(k);
        high = max(high, pvalue / ec);
        low = min(low, pvalue / ec);
        // if other channels are out of range, cap to white point of destination to avoid (likely) magenta highlights
        if (saturated_exp)
            pvalue = min(pvalue, wsp[cc]);

        saturated_exp |= saturation_a >= 1 - opt_saturation_offset;
        all_under &= under_exp;
        all_over &= saturated_exp;

        if (!(saturated_exp || under_exp)) {
            // use worst weight from surrounding pixels
            w = min(get_weight(max(saturation_a, saturation_b), opt_saturation_offset, opt_black_offset),
                       get_weight(under, opt_saturation_offset, opt_black_offset));
            div += w  * ec;
            (*out[cc])(k) += w * pvalue;
        }
    }
    // after looping over all frames, now we can
    // correct oor values and correctly weight output for in range
    if (all_under) {
        (*out[cc])(k) = low;
    } else if (all_over) {
        (*out[cc])(k) = min(high, wsp[cc] / ec);
    } else {
        (*out[cc])(k) = scale * (*out[cc])(k) / div;
    }
    return all_under + all_over * 2;
}

std::tuple<long, long> linear_response_slim(pfs::Array2D *out[],
                                       const ExposureList *imgs[],
                                       const float opt_saturation_offset,
                                       const float opt_black_offset,
                                       const float scale,
                                       const float rgb_corr[3][3],
                                       const float wsp[3],
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

    int k; // used to ravel i,j pixel coordinates
    char *oor = (char *) malloc(NP  * sizeof(char)); // track out of range pixels for masking
    // loop over image
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            k = j + i * width;
            oor[k] = merge(out, imgs, i, j, r0, g0, opt_black_offset, opt_saturation_offset, scale, wsp);
            // track out of range
            saturated_pixels += oor[k] == 2;
            under_pixels += oor[k] == 1;
        }
    }
    // demosaic after merge
    if (demosaic){
        dht_interpolate(out);
        // in this case data not yet mapped to output colorspace, do this now
        for (int j = 0; j < NP; j++) {
            apply_color_transform(j, out, rgb_corr);
            // mask out of range
            for (int cc = 0; cc < 3; cc++) {
                if (oor[j] == 2 && oor_high >= 0)
                    (*out[cc])(j) = oor_high;
                if (oor[j] == 1 && oor_low >= 0)
                    (*out[cc])(j) = oor_low;
            }

        }
    } else {
        int cc;
        // fill in out of range pixels
        for (int i = 0; i < height; i++)
            for (int j = 0; j < width; j++) {
                k = j + i * width;
                cc = grid_color(i, j, g0, r0);
                if (oor[j] == 2 && oor_high >= 0)
                    (*out[cc])(k) = oor_high;
                if (oor[j] == 1 && oor_low >= 0)
                    (*out[cc])(k) = oor_low;
            }
    }
    return std::make_tuple(saturated_pixels, under_pixels);
}
