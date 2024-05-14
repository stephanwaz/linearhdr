/**
 *
 * simplified hdr generation assuming raw linear response
 * @author Stephen Wasilewski stephen.wasilewski@epfl.ch
 *
 * Forked from: version pfstools 2.2.0:
 * @brief Robertson02 algorithm for automatic self-calibration.
 *
 * 
 * This file is derived from a part of PFS CALIBRATION package (but at this
 * point has mostly been entirely rewritten, but theseus' ship and all):
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


#define PROG_NAME "linearhdr"

extern bool verbose; /* verbose should be declared for each standalone code */

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
    return std::exp(-0.01/x - 0.01/(1-s-x)) + 1e-6;
}

void apply_color_transform(float irgb[3], float orgb[3], const float rgb_corr[3][3]){
    float x, y, z;
    x = irgb[0];
    y = irgb[1];
    z = irgb[2];
    orgb[0] = max(0, x * rgb_corr[0][0] + y * rgb_corr[0][1] + z * rgb_corr[0][2]);
    orgb[1] = max(0, x * rgb_corr[1][0] + y * rgb_corr[1][1] + z * rgb_corr[1][2]);
    orgb[2] = max(0, x * rgb_corr[2][0] + y * rgb_corr[2][1] + z * rgb_corr[2][2]);
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

float get_efc(int pixh, int height, float etime, const float efc[3][4]) {
    // choose right efc coefficients
    float efcc = 1.0;
    int efci;
    bool match;
    for (efci = 0; efci < 3; efci++) {
        match = etime < efc[0][efci];
        if (match)
            break;
    }
    if (match) {
        float pheight = height - pixh;
        efcc = 1 / (1 + (efc[efci][1]*pheight*pheight + efc[efci][2]*pheight + efc[efci][3]) / etime);
    }
    return efcc;
}

// merge a mosaic image considering adjacent pixels
// helps to detect sensor blooming
unsigned char merge(pfs::Array2D *out[],
                    const ExposureList *imgs[],
                    int i, int j, int r0, int g0,
                    const float blk_off,
                    const float sat_off,
                    const float scale,
                    const float wsp[3],
                    const float efc[3][4],
                    const bool best = false,
                    const bool onep = false) {
    int N = imgs[0]->size();
    int width = out[0]->getCols();
    int height = out[0]->getRows();
    int k = j + i * width;
    float w; // stores the weight for each frame
    float div = 0.0;
    for (int cc = 0; cc < 3; cc++) //make sure output is clear
        (*out[cc])(k) = 0;
    int cc = grid_color(i, j, g0, r0);
    int fp0 = 0;
    int fp = 5;
    int bloomk[5][5];
    fill(bloomk[0], bloomk[0] + 25, j + i * width);
    int ccb;
    int bm = cc;
    if (cc == 1 && ((i-r0) % 2) == 0) // g2
        bm = 3;
    // 1D indices for image lookups
    for (int bi = -2; bi < 3; bi++)
        if (i+bi >= 0 && i+bi < height)
            for (int bj = -2; bj < 3; bj++) {
                if (j+bj >= 0 && j+bj < width)
                    bloomk[bi + 2][bj + 2] = (j + bj) + (i + bi) * width;
            }
    // the surrounding color channels for the 4 positions in the bayer CFA (r, g1, b, g2)
    const int bloommap[4][5][5] = {{{0, 1, 0, 1, 0},
                                    {1, 2, 1, 2, 1},
                                    {0, 1, 0, 1, 0},
                                    {1, 2, 1, 2, 1},
                                    {0, 1, 0, 1, 0}},
                                 {{1, 2, 1, 2, 1},
                                  {0, 1, 0, 1, 0},
                                  {1, 2, 1, 2, 1},
                                  {0, 1, 0, 1, 0},
                                  {1, 2, 1, 2, 1}},
                                 {{2, 1, 2, 1, 2},
                                  {1, 0, 1, 0, 1},
                                  {2, 1, 2, 1, 2},
                                  {1, 0, 1, 0, 1},
                                  {2, 1, 2, 1, 2}},
                                 {{1, 0, 1, 0, 1},
                                  {2, 1, 2, 1, 2},
                                  {1, 0, 1, 0, 1},
                                  {2, 1, 2, 1, 2},
                                  {1, 0, 1, 0, 1}}};
    if (onep) { // to use a smaller footprint
        fp = 4;
        fp0 = 1;
    }

    // loop over individual frames, calculating weights
    bool all_under = true;
    bool all_over = true;
    bool saturated_exp;
    bool under_exp;
    float saturation_b; // highest raw value among surrounding pixels of different color
    float saturation_a; // highest raw value among surrounding pixels of same color
    float pvalue; // current pixel saturation
    float under; // lowest raw value among surrounding pixels
    float high = 0.0; // max value (including oor)
    float low = 1e30; // min value (including oor)
    float ec;
    float bw = -1;
    float bv = 0.0;
    for (int n = 0; n < N; n++) {
        ec = get_exposure_compensation((*imgs[cc])[n]);
        saturation_b = saturation_a = 0.0;
        under = 1.0;
        // first check surrounding pixels for saturation
        for (int b = fp0; b < fp; b++)
            for (int bj = fp0; bj < fp; bj++) {
                if (bloomk[b][bj] > -1) {
                    ccb = bloommap[bm][b][bj];
                    if (ccb != cc) {
                        saturation_b = max(saturation_b, (*(*imgs[ccb])[n].yi)(bloomk[b][bj]));

                    } else {
                        saturation_a = max(saturation_a, (*(*imgs[ccb])[n].yi)(bloomk[b][bj]));
                    }
                    under = min(under, (*(*imgs[ccb])[n].yi)(bloomk[b][bj]));
                }
        }
        pvalue = (*(*imgs[cc])[n].yi)(k);
        // saturation of other colors
        saturated_exp = saturation_b >= 1 - sat_off;
        under_exp = under < blk_off;

        high = max(high, min(pvalue, wsp[cc]) / ec);
        low = min(low, pvalue / ec);

        // if other channels are out of range, cap to white point of destination to avoid colored highlights
        if (saturated_exp)
            pvalue = min(pvalue, wsp[cc]);

        // saturation of same color
        saturated_exp |= saturation_a >= 1 - sat_off;

        all_under &= under_exp;
        all_over &= saturated_exp;

        if (!(saturated_exp || under_exp)) {
            // use worst weight from surrounding pixels to include some IPC
            w = min(get_weight(max(saturation_a, saturation_b), sat_off),
                    get_weight(under, sat_off));
            div += w  * ec;
            if (efc[0][0] > 0)
                pvalue *= get_efc(i, height, (*imgs[cc])[n].exposure_time, efc);
            (*out[cc])(k) += w * pvalue;
            // nullify poisson noise weighting (exchange with two lines above)
            // div += w;
            //(*out[cc])(k) += w * irgb[cc] / ec;
            if (w * ec > bw) {
                bw = w * ec;
                bv = pvalue / ec;
            }
        }
    }
    // after looping over all frames, now we can
    // correct oor values and correctly weight output for in range
    if (all_under) {
        (*out[cc])(k) = low;
    } else if (all_over) {
        (*out[cc])(k) = high;
    } else if (best) {
        (*out[cc])(k) = scale * bv;
    } else {
        (*out[cc])(k) = scale * (*out[cc])(k) / div;
    }
    return all_under + all_over * 2;
}

// merge an RGB demosaiced sequence, weights by worst channel
unsigned char merge_rgb(pfs::Array2D *out[],
                        const ExposureList *imgs[],
                        int i, int j,
                        const float blk_off,
                        const float sat_off,
                        const float scale,
                        const float wsp[3],
                        const float rgb_corr[3][3],
                        const float efc[3][4],
                        const bool best = false) {
    int N = imgs[0]->size();
    int width = out[0]->getCols();
    int height = out[0]->getRows();
    int k = j + i * width;
    float w; // stores the weight for each frame
    float div = 0.0;
    for (int cc = 0; cc < 3; cc++) //make sure output is clear
        (*out[cc])(k) = 0;
    // loop over individual frames, calculating weights
    bool all_under = true;
    bool all_over = true;
    bool saturated_exp;
    bool under_exp;
    float saturation; // highest raw value among 3 colors
    float under; // lowest raw value among 3 colors
    float high[3] = {0.0, 0.0, 0.0}; // max value (including oor)
    float low[3] = {1e30, 1e30, 1e30}; // min value (including oor)
    float ec;
    float bw = -1;
    float bv[3] = {0.0, 0.0, 0.0};
    float irgb[3];
    float orgb[3];
    for (int n = 0; n < N; n++) {
        ec = get_exposure_compensation((*imgs[0])[n]);
        irgb[0] = (*(*imgs[0])[n].yi)(k);
        irgb[1] = (*(*imgs[1])[n].yi)(k);
        irgb[2] = (*(*imgs[2])[n].yi)(k);
        saturation = max(irgb[0], irgb[1], irgb[2]);
        under = min(irgb[0], irgb[1], irgb[2]);;

        saturated_exp = saturation >= 1 - sat_off;
        under_exp = under < blk_off;
        // need to also check that destination space is not out of range otherwise
        // magenta highlights will creep in
        apply_color_transform(irgb, orgb, rgb_corr);
        saturated_exp |= max(orgb[0], orgb[1], orgb[2]) >= 1 - sat_off;
        under_exp |= min(orgb[0], orgb[1], orgb[2]) < blk_off;

        all_under &= under_exp;
        all_over &= saturated_exp;

        // use worst weight from all channels
        w = min(get_weight(saturation, sat_off),
                get_weight(under, sat_off));
        if (!(saturated_exp || under_exp)) {
            // no noise weighting: div += w
            div += w * ec;
            if (w * ec > bw) {
                bw = w * ec;
                bv[0] = (*(*imgs[0])[n].yi)(k) / ec;
                bv[1] = (*(*imgs[1])[n].yi)(k) / ec;
                bv[2] = (*(*imgs[2])[n].yi)(k) / ec;
            }
        }
        float efcc = 1.0;
        if (efc[0][0] > 0)
            efcc = get_efc(i, height, (*imgs[0])[n].exposure_time, efc);

        for (int cc = 0; cc < 3; cc++) {
            irgb[cc] *= efcc;
            high[cc] = max(high[cc], min(irgb[cc], wsp[cc])  / ec);
            low[cc] = min(low[cc], irgb[cc] / ec);
            // if other channels are out of range, cap to white point of destination to avoid (likely) magenta highlights
            if (saturated_exp && irgb[cc] < 1 - sat_off)
                irgb[cc] = min(irgb[cc], wsp[cc]);
            if (!(saturated_exp || under_exp)) {
                // no noise weighting: (*out[cc])(k) += w * pvalue[cc] / ec;
                (*out[cc])(k) += w * irgb[cc];
            }
        }
    }
    // if all frames are either over or under and none are in range (shouldn't happen with good exposure sequence
    // in camera raw space unless very peaky spectrum
    if (div == 0)
        all_over = true;

    // after looping over all frames, now we can
    // correct oor values and correctly weight output for in range
    for (int cc = 0; cc < 3; cc++) {
        if (all_under) {
            (*out[cc])(k) = low[cc];
        } else if (all_over) {
            (*out[cc])(k) = high[cc];
        } else if (best) {
            (*out[cc])(k) = scale * bv[cc];
        } else {
            (*out[cc])(k) = scale * (*out[cc])(k) / div;
        }
    }
    return all_under + all_over * 2;
}

// merge a single pixel channel (can be bayer or not)
// does not track saturation or under flow, should only be used experimentally
unsigned char merge_independent(pfs::Array2D *out[],
                    const ExposureList *imgs[],
                    int i, int j, int cc,
                    const float blk_off,
                    const float sat_off,
                    const float scale,
                    const float wsp,
                    const float efc[3][4],
                    const bool best = false) {
    int N = imgs[0]->size();
    int width = out[0]->getCols();
    int height = out[0]->getRows();
    int k = j + i * width;
    float w; // stores the weight for each frame
    float div = 0.0;
    (*out[cc])(k) = 0; // reset pixel
    // loop over individual frames, calculating weights
    bool all_under = true;
    bool all_over = true;
    bool saturated_exp;
    bool under_exp;
    float high = 0.0; // max value (including oor)
    float low = 1e30; // min value (including oor)
    float ec;
    float pvalue;
    float bw = -1;
    float bv = 0.0;
    for (int n = 0; n < N; n++) {
        ec = get_exposure_compensation((*imgs[cc])[n]);
        pvalue = (*(*imgs[cc])[n].yi)(k);
        saturated_exp = pvalue >= 1 - sat_off;
        under_exp = pvalue < blk_off;

        high = max(high, min(pvalue, wsp) / ec);
        low = min(low, pvalue / ec);

        all_under &= under_exp;
        all_over &= saturated_exp;

        if (!(saturated_exp || under_exp)) {
            w = get_weight(pvalue, sat_off);
            div += w  * ec;
            if (efc[0][0] > 0)
                pvalue *= get_efc(i, height, (*imgs[cc])[n].exposure_time, efc);
            (*out[cc])(k) += w * pvalue;
            if (w * ec > bw) {
                bw = w * ec;
                bv = pvalue / ec;
            }
        }
    }
    // after looping over all frames, now we can
    // correct oor values and correctly weight output for in range
    if (all_under) {
        (*out[cc])(k) = low;
    } else if (all_over) {
        (*out[cc])(k) = high;
    } else if (best) {
        (*out[cc])(k) = scale * bv;
    } else {
        (*out[cc])(k) = scale * (*out[cc])(k) / div;
    }
    return all_under + all_over * 2;
}

std::tuple<long, long> linear_response(pfs::Array2D *out[],
                                       const ExposureList *imgs[],
                                       const float sat_off,
                                       const float blk_off,
                                       const float scale,
                                       const float rgb_corr[3][3],
                                       const float wsp[3],
                                       const float efc[3][4],
                                       const float oor_high,
                                       float oor_low,
                                       bool isbayer,
                                       const bool demosaic,
                                       bool mergeeach,
                                       bool usebest,
                                       const bool median){

    // governs primary merging, do we merge all color channels, or according to CFA
    isbayer = isbayer or demosaic;

    // frame size
    int width = out[0]->getCols();
    int height = out[0]->getRows();

    int NP = width * height;

    // track out of bounds pixels
    int saturated_pixels = 0;
    int under_pixels = 0;

    int k; // used to ravel i,j pixel coordinates
    auto *oor = (unsigned char *) malloc(NP  * sizeof(unsigned char)); // track out of range pixels for masking

    // establish bayer grid coordinates
    int g0 = first_non_zero((*imgs[1])[0].yi);
    int r0 = first_non_zero_row((*imgs[0])[0].yi);

    // three different merging algorithms, vary by what information is used to test and weight pixels
    // (how many surrounding (CFA) or coincident (RGB) to include)
    if (mergeeach) {
        int ccf = 0; // first color index in loop, adjusted at each pixel when isbayer
        int cinc = isbayer ? 3: 1; // whether color loop will increment through all colors or only first
        for (int i = 0; i < height; i++)
            for (int j = 0; j < width; j++) {
                k = j + i * width;
                if (isbayer)
                    ccf = grid_color(i, j, g0, r0);
                for (int cc = ccf; cc < 3; cc += cinc) {
                    oor[k] = merge_independent(out, imgs, i, j, cc, blk_off, sat_off, scale, wsp[cc], efc, usebest);
                    // track out of range
                    saturated_pixels += oor[k] == 2;
                    under_pixels += oor[k] == 1;
                }
            }
    } else {
        if (isbayer) {
            VERBOSE_STR << "1st green column: " << g0 << " 1st red row: " << r0 << std::endl;
            // loop over bayer image
            for (int i = 0; i < height; i++)
                for (int j = 0; j < width; j++) {
                    k = j + i * width;
                    oor[k] = merge(out, imgs, i, j, r0, g0, blk_off, sat_off, scale, wsp, efc, usebest);
                    // track out of range
                    saturated_pixels += oor[k] == 2;
                    under_pixels += oor[k] == 1;
                }
        } else {
            //loop over rgb image
            for (int i = 0; i < height; i++)
                for (int j = 0; j < width; j++) {
                    k = j + i * width;
                    oor[k] = merge_rgb(out, imgs, i, j, blk_off, sat_off, scale, wsp, rgb_corr, efc, usebest);
                    // track out of range
                    saturated_pixels += oor[k] == 2;
                    under_pixels += oor[k] == 1;
                }
        }
    }
    if (isbayer && !demosaic) {
        int cc;
        // fill in out of range pixels
        for (int i = 0; i < height; i++)
            for (int j = 0; j < width; j++) {
                k = j + i * width;
                cc = grid_color(i, j, g0, r0);
                if (oor[k] == 2 && oor_high >= 0)
                    (*out[cc])(k) = oor_high;
                if (oor[k] == 1 && oor_low >= 0)
                    (*out[cc])(k) = oor_low;
            }
    } else {
        if (demosaic) // demosaic after merge
            dht_interpolate(out, median);
        // in this case data not yet mapped to output colorspace, do this now
        for (int j = 0; j < NP; j++) {
            apply_color_transform(j, out, rgb_corr);
            // mask out of range
            for (int cc = 0; cc < 3; cc++) {
                if (oor[k] == 2 && oor_high >= 0)
                    (*out[cc])(j) = oor_high;
                if (oor[k] == 1 && oor_low >= 0)
                    (*out[cc])(j) = oor_low;
            }

        }
    }
    return std::make_tuple(saturated_pixels, under_pixels);
}
