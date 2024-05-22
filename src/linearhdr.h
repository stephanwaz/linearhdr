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
 *
 * $Id: robertson02.h,v 1.4 2011/02/15 15:46:27 ihrke Exp $
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

#ifndef _linearhdr_h_
#define _linearhdr_h_

#include <stdio.h>

#include <vector>
#include <array2d.h>
#include <pfs.h>
#include <tuple>

/**
 * @brief Container for images taken with different exposures
 */
class Exposure
{
public:
    float exposure_time;          // exposure time
    float iso;                    // sensor gain iso value
    float aperture;
    pfs::Array2D* yi;		// exposed pixel value (camera output)
};


/**
 * @brief Container for a list of exposures
 */
typedef std::vector<Exposure> ExposureList;

struct FrameInfo {float etime; float iso; float aperture; float factor;};

std::tuple<long, long> linear_response(pfs::Array2D *out[],
                                    const ExposureList *imgs[],
                                    const float sat_off,
                                    const float blk_off,
                                    const float scale,
                                    const float rgb_corr[3][3],
                                    const float wsp[3], // white saturation point in target color space of raw values
                                    const float efc[3][4],
                                    const float oor_high,
                                    float oor_low,
                                    bool isbayer,
                                    const bool demosaic,
                                    bool mergeeach = false,
                                    bool usebest = false,
                                    const bool median = true,
                                    const bool usemax = false, // only used if isbayer and not mergeeach
                                    const bool usemin = false, // only used if isbayer and not mergeeach
                                    const bool onep = false); // only used if isbayer and not mergeeach


void dht_interpolate(pfs::Array2D *imgdata[], bool median = true);

int first_non_zero_row(pfs::Array2D *X);

int first_non_zero(pfs::Array2D *X);

int grid_color(int i, int j, int g0, int r0);

void apply_color_transform(int j, pfs::Array2D *out[], const float rgb_corr[3][3]);

void apply_color_transform(float irgb[3], float orgb[3], float rgb_corr[3][3]);

float get_efc(int pixh, int height, float etime, const float efc[3][4]);

inline float max(float a, float b) {
    return (a > b) ? a : b;
}

inline float min(float a, float b) {
    return (a > b) ? b : a;
}

inline float max(float a, float b, float c) {
    return max(max(a, b), c);
}

inline float min(float a, float b, float c) {
    return min(min(a, b), c);
}

inline float max(float a, float b, float c, float d) {
    return max(max(max(a, b), c), d);
}

#endif /* #ifndef _linearhdr_h_ */
