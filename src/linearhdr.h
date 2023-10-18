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

/**
 * @brief Container for images taken with different exposures
 */
class Exposure
{
public:
    float ti;			// exposure value (eg time) - including iso and apperture
    float exposure_time;          // exposure time
    float iso;                    // sensor gain iso value
    float aperture;
    pfs::Array2D* yi;		// exposed pixel value (camera output)

    // to be able to sort exposures using stl
    static bool msort( Exposure const & a, Exposure const & b)
    { return a.ti < b.ti; }
};


/**
 * @brief Container for a list of exposures
 */
typedef std::vector<Exposure> ExposureList;


int linear_response(pfs::Array2D *rgb_out[],
                    const ExposureList *imgs[],
                    const float opt_saturation_offset,
                    const float opt_black_offset,
                    float deghosting_value,
                    const float scale,
                    const float vlambda[3],
                    const float rgb_corr[3][3],
                    const float oor_high = -1,
                    const float oor_low = -1,
                    bool isbayer = false,
                    const bool demosaic = false);

void dht_interpolate(pfs::Array2D *Xj, pfs::Array2D *Yj, pfs::Array2D *Zj);

int first_non_zero_row(pfs::Array2D *X);

int first_non_zero(pfs::Array2D *X);

int grid_color(int i, int j, int g0, int r0);

#endif /* #ifndef _linearhdr_h_ */
