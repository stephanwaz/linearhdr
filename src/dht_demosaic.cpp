/* -*- C++ -*-
 * Copyright 2023 Stephen Wasilewski
 * Modified from:
 *
 * File: dht_demosaic.cpp
 * Copyright 2013 Anton Petrusevich
 * Created: Tue Apr  9, 2013
 *
 * This code is licensed under one of two licenses as you choose:
 *
 * 1. GNU LESSER GENERAL PUBLIC LICENSE version 2.1
 *    (See file LICENSE.LGPL provided in LibRaw distribution archive for
 * details).
 *
 * 2. COMMON DEVELOPMENT AND DISTRIBUTION LICENSE (CDDL) Version 1.0
 *    (See file LICENSE.CDDL provided in LibRaw distribution archive for
 * details).
 *
 */

/*
 * Copyright (c) 2023 Stephen Wasilewski, EPFL
 * primary modification is to make compatible with pfsarray data type rather than libraw
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

/*
 * функция вычисляет яркостную дистанцию.
 * если две яркости отличаются в два раза, например 10 и 20, то они имеют такой
 * же вес при принятии решения об интерполировании, как и 100 и 200 --
 * фотографическая яркость между ними 1 стоп. дистанция всегда >=1
 */

/* Deepl translation:
 * * the function calculates the luminance distance.
 * If two luminances differ by a factor of two, such as 10 and 20, they have the same weight in the interpolation decision as 100 and 200.
 * the same weight in the interpolation decision as 100 and 200 --
 * photographic brightness between them is 1 stop. distance is always >=1
 *
 */

//#include <internal/dmp_include.h>
#include <cmath>  //for sqrt()
#include <pfs.h> //use pfs::array instead in libraw.imgdata
#include <iostream>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

// always return ratio as >= 1
static inline float calc_dist(float c1, float c2) {
    return c1 > c2 ? c1 / c2 : c2 / c1;
}

//because floating point input, algorithm needs a non-zero base value
//this is added to the raw input but then subtracted at final copy
//edit: better to pass non-zero values when that is the desired effect
//this results in NAN when one channel in denominator is 0, wiping out any pixels
//take information from a zero pixel.
static const float FLOOR = 0.0;

struct DHT {
    int nr_height, nr_width, iwidth, iheight; //SW
    int r0, g0;  //SW to locate bayer grid after potential cropping
    static const int nr_topmargin = 4, nr_leftmargin = 4;

    float (*nraw)[3];
    float channel_maximum[3]; //SW change to float
    float channel_minimum[3];
//  LibRaw &libraw; //SW don't use

    // these are the integer codes for the different interpolation schemes
    enum {
        HVSH = 1, // make odd when there is a big difference
        HOR = 2, // add 2 for horizontal
        VER = 4, // add 4 for vertical
        HORSH = HOR | HVSH,
        VERSH = VER | HVSH,
        DIASH = 8, // add 8 when there is a big difference
        LURD = 16, // left-up right-down
        RULD = 32, // right-up left-down
        LURDSH = LURD | DIASH,
        RULDSH = RULD | DIASH,
    };

    static inline float Tg(void) throw() { return 256.0f; } // big difference in horizontal/vertical gradients

    static inline float T(void) throw() { return 1.4f; } // big difference in diagonal gradients

    char *ndir; // array for interpolation codes

    //SW renamed to shorten expressions
    //for uniformity would be better to change nraw to 2Darray
    inline int nro(int row, int col) throw() {
        return (row * nr_width + col);
    }

    // pick horizontal or vertical at non-green pixels
    int get_hv_grb(int x, int y, int kc) {
        // this makes weights green from whichever color (red blue) is at interpolation point
        //ratio of green:color above and below
        float hv1 = 2 * nraw[nro(y - 1, x)][1] /
                (nraw[nro(y - 2, x)][kc] + nraw[nro(y, x)][kc]);
        float hv2 = 2 * nraw[nro(y + 1, x)][1] /
                (nraw[nro(y + 2, x)][kc] + nraw[nro(y, x)][kc]);
        //relative difference in greenness up/down * relative color difference center:up/down
        float kv = calc_dist(hv1, hv2) *
                calc_dist(nraw[nro(y, x)][kc] * nraw[nro(y, x)][kc],
                          (nraw[nro(y - 2, x)][kc] * nraw[nro(y + 2, x)][kc]));
        // fourth power
        kv *= kv;
        kv *= kv;
        kv *= kv;
        //kv * relative green difference between farther and closer greens above and below
        float dv = kv * calc_dist(nraw[nro(y - 3, x)][1] * nraw[nro(y + 3, x)][1],
                                  nraw[nro(y - 1, x)][1] * nraw[nro(y + 1, x)][1]);

        // ratio of green:color left and right
        float hh1 = 2 * nraw[nro(y, x - 1)][1] /
                (nraw[nro(y, x - 2)][kc] + nraw[nro(y, x)][kc]);
        float hh2 = 2 * nraw[nro(y, x + 1)][1] /
                (nraw[nro(y, x + 2)][kc] + nraw[nro(y, x)][kc]);
        //relative green difference * relative color difference left and right
        float kh = calc_dist(hh1, hh2) *
                calc_dist(nraw[nro(y, x)][kc] * nraw[nro(y, x)][kc],
                          (nraw[nro(y, x - 2)][kc] * nraw[nro(y, x + 2)][kc]));
        // fourth power
        kh *= kh;
        kh *= kh;
        kh *= kh;
        //kh * relative green difference between farther and closer greens left and right
        float dh = kh * calc_dist(nraw[nro(y, x - 3)][1] * nraw[nro(y, x + 3)][1],
                                  nraw[nro(y, x - 1)][1] * nraw[nro(y, x + 1)][1]);
        float e = calc_dist(dh, dv);
        // assign to code depending on strength of difference
        // if horizontal is much larger than vertical: 5
        // if H is larger than V: 4
        // if V in much larger than H: 3
        // if V is larger than H: 2
        char d = dh < dv ? (e > Tg() ? HORSH : HOR) : (e > Tg() ? VERSH : VER);
        return d;
    }

    // pick horizontal/vertical at green pixels
    int get_hv_rbg(int x, int y, int hc) {
        // this makes weights for red/blue when green is at the interpolation point
        // ratio other_color:green above and below
        float hv1 = 2 * nraw[nro(y - 1, x)][hc ^ 2] /
                (nraw[nro(y - 2, x)][1] + nraw[nro(y, x)][1]);
        float hv2 = 2 * nraw[nro(y + 1, x)][hc ^ 2] /
                (nraw[nro(y + 2, x)][1] + nraw[nro(y, x)][1]);
        // relative green difference above and below
        float kv = calc_dist(hv1, hv2) * calc_dist(nraw[nro(y, x)][1] * nraw[nro(y, x)][1],
                                                   (nraw[nro(y - 2, x)][1] * nraw[nro(y + 2, x)][1]));

        // fourth power
        kv *= kv;
        kv *= kv;
        kv *= kv;
        //kv * relative other color difference between farther and closer other color above and below
        float dv = kv * calc_dist(nraw[nro(y - 3, x)][hc ^ 2] * nraw[nro(y + 3, x)][hc ^ 2],
                                  nraw[nro(y - 1, x)][hc ^ 2] * nraw[nro(y + 1, x)][hc ^ 2]);

        // ratio of color:green left and right
        float hh1 = 2 * nraw[nro(y, x - 1)][hc] / (nraw[nro(y, x - 2)][1] + nraw[nro(y, x)][1]);
        float hh2 = 2 * nraw[nro(y, x + 1)][hc] / (nraw[nro(y, x + 2)][1] + nraw[nro(y, x)][1]);
        float kh = calc_dist(hh1, hh2) * calc_dist(nraw[nro(y, x)][1] * nraw[nro(y, x)][1],
                                                   (nraw[nro(y, x - 2)][1] * nraw[nro(y, x + 2)][1]));
        kh *= kh;
        kh *= kh;
        kh *= kh;
        //kv * relative other color difference between farther and closer other color left and right
        float dh = kh * calc_dist(nraw[nro(y, x - 3)][hc] * nraw[nro(y, x + 3)][hc],
                                  nraw[nro(y, x - 1)][hc] * nraw[nro(y, x + 1)][hc]);
        float e = calc_dist(dh, dv);
        // assign to code depending on strength of difference
        // if horizontal is much larger than vertical: 5
        // if H is larger than V: 4
        // if V in much larger than H: 3
        // if V is larger than H: 2
        char d = dh < dv ? (e > Tg() ? HORSH : HOR) : (e > Tg() ? VERSH : VER);
        return d;
    }

    // pick diagonal at red/blue pixel
    int get_diag_grb(int x, int y, int kc) {
        // similar to hv but without extra width on kernel (only looks to nearest neighbor because green is already there)
        // happens after initial green interpolation step
        // ratio of green:color upper-left (hlu) and lower-right (hrd)
        // SW: I don't think these values are actually used as they will cancel at e = ...
        // maybe there is supposed to be an additional calculation of lower-left and upper-right to weight druld?
        // Edit: added hld, hru below to match pattern in get_hv_grb
        // An alternative would be to get rid of hlu,hrd,hld,hru. Either way, no longer sure that T() is scaled correctly?
        float hlu = nraw[nro(y - 1, x - 1)][1] / nraw[nro(y - 1, x - 1)][kc];
        float hrd = nraw[nro(y + 1, x + 1)][1] / nraw[nro(y + 1, x + 1)][kc];
        // change in greenness up-left to low-right * relative green difference center:UL/LR
        float dlurd = calc_dist(hlu, hrd) *
                      calc_dist(nraw[nro(y - 1, x - 1)][1] * nraw[nro(y + 1, x + 1)][1],
                                nraw[nro(y, x)][1] * nraw[nro(y, x)][1]);

        float hld = nraw[nro(y + 1, x - 1)][1] / nraw[nro(y + 1, x - 1)][kc];
        float hru = nraw[nro(y - 1, x + 1)][1] / nraw[nro(y - 1, x + 1)][kc];
        // change in greenness up-right to low-left * relative green difference center:UR/LL
        float druld = calc_dist(hld, hru) *
                      calc_dist(nraw[nro(y - 1, x + 1)][1] * nraw[nro(y + 1, x - 1)][1],
                                nraw[nro(y, x)][1] * nraw[nro(y, x)][1]);
        float e = calc_dist(dlurd, druld);
        // assign to code depending on strength of difference
        // if + (ruld) is much larger than - (lurd): 24
        // if + is larger than -: 16
        // if - in much larger than +: 40
        // if - is larger than +: 32
        char d = druld < dlurd ? (e > T() ? RULDSH : RULD) : (e > T() ? LURDSH : LURD);
        return d;
    }

    // pick diagonal at green pixel
    int get_diag_rbg(int x, int y, int /* hc */) {
        // relative green difference center:UL/LR
        float dlurd = calc_dist(nraw[nro(y - 1, x - 1)][1] * nraw[nro(y + 1, x + 1)][1],
                                nraw[nro(y, x)][1] * nraw[nro(y, x)][1]);
        // relative green difference center:UR/LL
        float druld = calc_dist(nraw[nro(y - 1, x + 1)][1] * nraw[nro(y + 1, x - 1)][1],
                                nraw[nro(y, x)][1] * nraw[nro(y, x)][1]);
        float e = calc_dist(dlurd, druld);
        // assign to code depending on strength of difference
        // if + (ruld) is much larger than - (lurd): 24
        // if + is larger than -: 16
        // if - in much larger than +: 40
        // if - is larger than +: 32
        char d = druld < dlurd ? (e > T() ? RULDSH : RULD) : (e > T() ? LURDSH : LURD);
        return d;
    }

    static inline float scale_over(float ec, float base) {
        float s = base * .4;
        float o = ec - base;
        return base + sqrt(s * (o + s)) - s;
    }

    static inline float scale_under(float ec, float base) {
        float s = base * .6;
        float o = base - ec;
        return base - sqrt(s * (o + s)) + s;
    }

    ~DHT();
    DHT(pfs::Array2D *Xj, pfs::Array2D *Yj, pfs::Array2D *Zj); //SW base on pfs::Array

    // main executions steps
    void make_hv_dirs();
    void make_greens();
    void make_diag_dirs();
    void make_rb();
    void copy_to_image();

    // internal functions
    void refine_hv_dirs(int i, int js);
    void refine_ihv_dirs(int i);
    void refine_idiag_dirs(int i);
    void make_hv_dline(int i);
    void make_diag_dline(int i);
    void make_gline(int i);
    void make_rbdiag(int i);
    void make_rbhv(int i);

//  void hide_hots(); // thresholds don't make sense with arbitrary hdr scaling
//  void restore_hots();
    int COLOR(int i, int j);
    pfs::Array2D *imgdata[3];
};

typedef float float3[3];

/*
 * информация о цветах копируется во float в общем то с одной целью -- чтобы
 * вместо 0 можно было писать 0.5, что больше подходит для вычисления яркостной
 * разницы. причина: в целых числах разница в 1 стоп составляет ряд 8,4,2,1,0 --
 * последнее число должно быть 0.5, которое непредствамио в целых числах. так же
 * это изменение позволяет не думать о специальных случаях деления на ноль.
 *
 * альтернативное решение: умножить все данные на 2 и в младший бит внести 1.
 * правда, всё равно придётся следить, чтобы при интерпретации зелёного цвета не
 * получился 0 при округлении, иначе проблема при интерпретации синих и красных.
 *
 */

// deepl translation:
/*
 * color information is copied to float for one purpose in general -- so that
 * instead of 0, you can write 0.5, which is more suitable for calculating brightness differences.
 * difference. reason: in integers the difference of 1 stop is a series of 8,4,2,1,0 --
 * the last number should be 0.5, which is unpredictable in integers. also
 * this change allows us not to think about special cases of division by zero.
 *
 * an alternative solution is to multiply all the data by 2 and put 1 in the low bit.
 * However, you still have to make sure that the green color is not interpreted as a 0 when rounding.
 * * you have to make sure that you don't get a 0 when rounding, otherwise you'll have a problem when
 * interpreting blue and red colors.
 *
 * SW Note: with HDR input everything is assumed float, and values are allowed to go to zero
 *
 */

// cant use dcraw codes because image is cropped. this is reliable as long as r0 and g0 (see next functions)
// are properly set
int DHT::COLOR(int i, int j) {
    // is green
    if (j % 2 == ((g0 + i) % 2))
        return 1;
    // is red row
    if (i % 2 == r0)
        return 0;
    // is blue row
    return 2;
}


// assumes RGGB bayer grid, other sensor configs will not work
// use with red or blue
// is this color channel in the first (returns 0) or second row (returns 1)?
// requires non-zero values as some point (and also no data corruption (non-zero channel values off grid)
int first_non_zero_row(pfs::Array2D *X) {
    for (int i = 0; i < X->getRows(); i++) {
        for (int j = 0; j < X->getCols(); j++)
            if ((*X)(j, i) > 1e-6)
                return i % 2;
    }
    return 0;
}

// assumes X is green channel
// is the first pixel (0,0) green (returns 0) or not (returns 1, implying second pixel in first row is green)
// requires non-zero values as some point (and also no data corruption (non-zero channel values off grid)
int first_non_zero(pfs::Array2D *X) {
    for (int i = 0; i < X->getRows(); i++) {
        for (int j = 0; j < X->getCols(); j++)
            if ((*X)(j, i) > 1e-6)
                return (j + (i % 2)) % 2;
    }
    return 0;
}

// class initialization
DHT::DHT(pfs::Array2D *Xj, pfs::Array2D *Yj, pfs::Array2D *Zj) {
    imgdata[0] = Xj;
    imgdata[1] = Yj;
    imgdata[2] = Zj;
    iwidth = imgdata[0]->getCols();
    iheight = imgdata[0]->getRows();
    g0 = first_non_zero(Yj);
    r0 = first_non_zero_row(Xj);

    nr_height = iheight + nr_topmargin * 2;
    nr_width = iwidth + nr_leftmargin * 2;
    nraw = (float3 *) malloc(nr_height * nr_width * sizeof(float3));
    ndir = (char *) calloc(nr_height * nr_width, 1);
    channel_maximum[0] = channel_maximum[1] = channel_maximum[2] = 0;
    channel_minimum[0] = (*imgdata[0])(0);
    channel_minimum[1] = (*imgdata[1])(0);
    channel_minimum[2] = (*imgdata[2])(0);
    //SW set base level to FLOOR, this value is subtracted back at the end
    for (int i = 0; i < nr_height * nr_width; ++i)
        nraw[i][0] = nraw[i][1] = nraw[i][2] = FLOOR;
    // read in raw data and record min/max
    for (int i = 0; i < iheight; ++i) {
        int col_cache[48];
        // use a cache of indices for color for efficiency
        for (int j = 0; j < 48; ++j) {
            int l = COLOR(i, j);
            col_cache[j] = l;
        }
        for (int j = 0; j < iwidth; ++j) {
            int l = col_cache[j % 48];
            float c = (*imgdata[l])(j, i) + FLOOR;
            if (c != 0) {
                if (channel_maximum[l] < c)
                    channel_maximum[l] = c;
                if (channel_minimum[l] > c)
                    channel_minimum[l] = c;
                nraw[nro(i + nr_topmargin, j + nr_leftmargin)][l] = c;
            }
        }
    }
}

// step 1: add interpolation codes for orthogonal interpolation
void DHT::make_hv_dirs() {
    // assigns values 2-5 to DHT::ndir
    for (int i = 0; i < iheight; ++i) {
        make_hv_dline(i);
    }
    for (int i = 0; i < iheight; ++i) {
        refine_hv_dirs(i, i & 1);
    }
    for (int i = 0; i < iheight; ++i) {
        refine_hv_dirs(i, (i & 1) ^ 1);
    }
    for (int i = 0; i < iheight; ++i) {
        refine_ihv_dirs(i);
    }
    // at this point, values are all still 2-5 but isolated values are cleaned up
}

// step 1a: encode directional gradients
void DHT::make_hv_dline(int i) {
    int js = COLOR(i, 0) & 1;
    int kc = COLOR(i, js);
    /*
     * js -- начальная х-координата, которая попадает мимо известного зелёного
     * kc -- известный цвет в точке интерполирования
     *
     * * js -- an initial x-coordinate that falls past the known green
     * kc -- known color at the interpolation point
     */
    for (int j = 0; j < iwidth; j++) {
        int x = j + nr_leftmargin;
        int y = i + nr_topmargin;
        char d = 0;
        if ((j & 1) == js) {
            // red/blue pixel
            d = get_hv_grb(x, y, kc);
        } else {
            // green pixel
            d = get_hv_rbg(x, y, kc);
        }
        // assign to code depending on strength of difference
        // if horizontal is much larger than vertical: 5
        // if H is larger than V: 4
        // if V in much larger than H: 3
        // if V is larger than H: 2
        ndir[nro(y, x)] |= d;
    }
}

// step 1b,c
void DHT::refine_hv_dirs(int i, int js) {
    for (int j = js; j < iwidth; j += 2) {
        int x = j + nr_leftmargin;
        int y = i + nr_topmargin;
        // skip extreme gradients (when horizontal or vertical dominates)
        if (ndir[nro(y, x)] & HVSH)
            continue;
        // check how many of UDLR are vertical (0-16)
        int nv = (ndir[nro(y - 1, x)] & VER) + (ndir[nro(y + 1, x)] & VER) +
                 (ndir[nro(y, x - 1)] & VER) + (ndir[nro(y, x + 1)] & VER);
        // check how many of UDLR are horizontal (0-8)
        int nh = (ndir[nro(y - 1, x)] & HOR) + (ndir[nro(y + 1, x)] & HOR) +
                 (ndir[nro(y, x - 1)] & HOR) + (ndir[nro(y, x + 1)] & HOR);
        // check if current pixel matches an adjacent (in direction of current pixel)
        bool codir = (ndir[nro(y, x)] & VER) ? ((ndir[nro(y - 1, x)] & VER) ||
                                                (ndir[nro(y + 1, x)] & VER)) : (
                             (ndir[nro(y, x - 1)] & HOR) || (ndir[nro(y, x + 1)] & HOR));
        // scale to 0-4
        nv /= VER;
        nh /= HOR;
        // current pixel is vertical, at least 3 are horizontal, and no codirection
        if ((ndir[nro(y, x)] & VER) && (nh > 2 && !codir)) {
            // change 4,5 -> 2,3 (to horizontal to match neighborhood)
            ndir[nro(y, x)] &= ~VER;
            ndir[nro(y, x)] |= HOR;
        }
        // current pixel is horizontal, at least 3 are vertical, and no codirection
        if ((ndir[nro(y, x)] & HOR) && (nv > 2 && !codir)) {
            // change 2,3 -> 4,5 (to vertical to match neighborhood)
            ndir[nro(y, x)] &= ~HOR;
            ndir[nro(y, x)] |= VER;
        }
    }
}

// step 1d: same as above but only refines if all neighbors are different
void DHT::refine_ihv_dirs(int i) {
    for (int j = 0; j < iwidth; j++) {
        int x = j + nr_leftmargin;
        int y = i + nr_topmargin;
        if (ndir[nro(y, x)] & HVSH)
            continue;
        int nv = (ndir[nro(y - 1, x)] & VER) + (ndir[nro(y + 1, x)] & VER) +
                 (ndir[nro(y, x - 1)] & VER) + (ndir[nro(y, x + 1)] & VER);
        int nh = (ndir[nro(y - 1, x)] & HOR) + (ndir[nro(y + 1, x)] & HOR) +
                 (ndir[nro(y, x - 1)] & HOR) + (ndir[nro(y, x + 1)] & HOR);
        nv /= VER;
        nh /= HOR;
        if ((ndir[nro(y, x)] & VER) && nh > 3) {
            ndir[nro(y, x)] &= ~VER;
            ndir[nro(y, x)] |= HOR;
        }
        if ((ndir[nro(y, x)] & HOR) && nv > 3) {
            ndir[nro(y, x)] &= ~HOR;
            ndir[nro(y, x)] |= VER;
        }
    }
}

/*
 * вычисление недостающих зелёных точек.
 * calculating missing green dots.
 */
// step 2: interpolate green
void DHT::make_greens() {
    for (int i = 0; i < iheight; ++i) {
        make_gline(i);
    }
}

// step 2a: interpolate each row
void DHT::make_gline(int i) {
    // based on the data in ndir, perform initial interpolation of green channel
    // choose better direction (up/down left/right) and average look-ahead and look-behind
    // weighted by the ratio between relative difference of the look-ahead/behind to the current pixel
    int js = COLOR(i, 0) & 1;
    int kc = COLOR(i, js);
    /*
     * js -- начальная х-координата, которая попадает мимо известного зелёного
     * kc -- известный цвет в точке интерполирования
     *
     * js -- an initial x-coordinate that falls past the known green
     * kc -- known color at the interpolation point
     */
    for (int j = js; j < iwidth; j += 2) {
        int x = j + nr_leftmargin;
        int y = i + nr_topmargin;
        int dx, dy, dx2, dy2;
        float h1, h2;
        // h1, h2: interpolate vertically
        if (ndir[nro(y, x)] & VER) {
            dx = dx2 = 0;
            dy = -1;
            dy2 = 1;
            // green value scaled by average of current pixel color and two away (above, below)
            h1 = 2 * nraw[nro(y - 1, x)][1] / (nraw[nro(y - 2, x)][kc] + nraw[nro(y, x)][kc]);
            h2 = 2 * nraw[nro(y + 1, x)][1] / (nraw[nro(y + 2, x)][kc] + nraw[nro(y, x)][kc]);
        } else {
            dy = dy2 = 0;
            dx = 1;
            dx2 = -1;
            // green value scaled by average of current pixel color and two away (right, left)
            h1 = 2 * nraw[nro(y, x + 1)][1] / (nraw[nro(y, x + 2)][kc] + nraw[nro(y, x)][kc]);
            h2 = 2 * nraw[nro(y, x - 1)][1] / (nraw[nro(y, x - 2)][kc] + nraw[nro(y, x)][kc]);
        }
        // weights left/right or up/down (which is closer to current pixel
        float b1 = 1 / calc_dist(nraw[nro(y, x)][kc], nraw[nro(y + dy * 2, x + dx * 2)][kc]);
        float b2 = 1 / calc_dist(nraw[nro(y, x)][kc], nraw[nro(y + dy2 * 2, x + dx2 * 2)][kc]);
        // squared
        b1 *= b1;
        b2 *= b2;
        // current pixel value scaled by h1, h2 according to weighting b1,b2. since h = g/c, this restores the green
        // value
        float eg = nraw[nro(y, x)][kc] * (b1 * h1 + b2 * h2) / (b1 + b2);
        float min, max;
        min = MIN(nraw[nro(y + dy, x + dx)][1], nraw[nro(y + dy2, x + dx2)][1]);
        max = MAX(nraw[nro(y + dy, x + dx)][1], nraw[nro(y + dy2, x + dx2)][1]);
        min /= 1.2f;
        max *= 1.2f;
        // if value is well outsides bounds of two adjacent green pixels, scale back towards center
        if (eg < min)
            eg = scale_under(eg, min);
        else if (eg > max)
            eg = scale_over(eg, max);
        // don't let this value exceed the maximum value anywhere in the raw image
        if (eg > channel_maximum[1])
            eg = channel_maximum[1];
        else if (eg < channel_minimum[1])
            eg = channel_minimum[1];
        // at this point greens have been assigned everywhere in the image based on initial interpolation
        nraw[nro(y, x)][1] = eg;
    }
}

// step 3: add interpolation codes for diagonal interpolation
void DHT::make_diag_dirs() {
    for (int i = 0; i < iheight; ++i) {
        make_diag_dline(i);
    }
    for (int i = 0; i < iheight; ++i) {
        refine_idiag_dirs(i);
    }
    // at this point, values are 2-5 + 16,24,32,40
}

// step 3a:
void DHT::make_diag_dline(int i) {
    int js = COLOR(i, 0) & 1;
    int kc = COLOR(i, js);
    /*
     * js -- начальная х-координата, которая попадает мимо известного зелёного
     * kc -- известный цвет в точке интерполирования
     *
     * * js -- an initial x-coordinate that falls past the known green
     * kc -- known color at the interpolation point
     */
    for (int j = 0; j < iwidth; j++) {
        int x = j + nr_leftmargin;
        int y = i + nr_topmargin;
        char d = 0;
        if ((j & 1) == js) {
            // red/blue pixel
            d = get_diag_grb(x, y, kc);
        } else {
            // green pixel
            d = get_diag_rbg(x, y, kc);
        }
        // assign to code depending on strength of difference
        // if + (ruld) is much larger than - (lurd): 24
        // if + is larger than -: 16
        // if - in much larger than +: 40
        // if - is larger than +: 32
        // these values add onto current 2-5 codes from HV interpolation
        ndir[nro(y, x)] |= d;
    }
}

// step 3b:
void DHT::refine_idiag_dirs(int i) {
    for (int j = 0; j < iwidth; j++) {
        int x = j + nr_leftmargin;
        int y = i + nr_topmargin;
        // skip extreme gradients
        if (ndir[nro(y, x)] & DIASH)
            continue;
        // check how many of surrounding 8 pixels are LURD
        int nv = (ndir[nro(y - 1, x)] & LURD) + (ndir[nro(y + 1, x)] & LURD) +
                 (ndir[nro(y, x - 1)] & LURD) + (ndir[nro(y, x + 1)] & LURD) +
                 (ndir[nro(y - 1, x - 1)] & LURD) + (ndir[nro(y - 1, x + 1)] & LURD) +
                 (ndir[nro(y + 1, x - 1)] & LURD) + (ndir[nro(y + 1, x + 1)] & LURD);
        // check how many of surrounding 8 pixels are RULD
        int nh = (ndir[nro(y - 1, x)] & RULD) + (ndir[nro(y + 1, x)] & RULD) +
                 (ndir[nro(y, x - 1)] & RULD) + (ndir[nro(y, x + 1)] & RULD) +
                 (ndir[nro(y - 1, x - 1)] & RULD) + (ndir[nro(y - 1, x + 1)] & RULD) +
                 (ndir[nro(y + 1, x - 1)] & RULD) + (ndir[nro(y + 1, x + 1)] & RULD);
        nv /= LURD;
        nh /= RULD;
        if ((ndir[nro(y, x)] & LURD) && nh > 7) {
            // change LURD to RULD if all neighbors are RULD
            ndir[nro(y, x)] &= ~LURD;
            ndir[nro(y, x)] |= RULD;
        }
        if ((ndir[nro(y, x)] & RULD) && nv > 7) {
            // change RULD to LURD if all neighbors are LURD
            ndir[nro(y, x)] &= ~RULD;
            ndir[nro(y, x)] |= LURD;
        }
    }
}

// step 4:
void DHT::make_rb() {
    for (int i = 0; i < iheight; ++i) {
        make_rbdiag(i);
    }
    for (int i = 0; i < iheight; ++i) {
        make_rbhv(i);
    }
}
/*
 * интерполяция красных и синих.
 *
 * сначала интерполируются недостающие цвета, по диагональным направлениям от
 * которых находятся известные, затем ситуация сводится к тому как
 * интерполировались зелёные.
 */

/*
 * interpolation of reds and blues.
 *
 * first the missing colors are interpolated, along the diagonal directions from
 * from which the known colors are found, then the situation is reduced to the way in which
 * interpolated green.
 */
// step 4a:
void DHT::make_rbdiag(int i) {
    int js = COLOR(i, 0) & 1;
    int uc = COLOR(i, js);
    int cl = uc ^ 2;
    /*
     * js -- начальная х-координата, которая попадает на уже интерполированный
     * зелёный al -- известный цвет (кроме зелёного) в точке интерполирования cl
     * -- неизвестный цвет
     */

    /*
     * js -- initial x-coordinate, which falls on the already interpolated one green
     * uc -- known color (other than green) at the interpolation point
     * cl -- unknown color
     */
    for (int j = js; j < iwidth; j += 2) {
        int x = j + nr_leftmargin;
        int y = i + nr_topmargin;
        int dx, dy, dx2, dy2;
        // interpolate LURD
        if (ndir[nro(y, x)] & LURD) {
            dx = -1;
            dx2 = 1;
            dy = -1;
            dy2 = 1;
        } else {
            dx = -1;
            dx2 = 1;
            dy = 1;
            dy2 = -1;
        }
        // green value scaled by average of current pixel color and  along chosen diagonal
        float g1 = 1 / calc_dist(nraw[nro(y, x)][1], nraw[nro(y + dy, x + dx)][1]);
        float g2 = 1 / calc_dist(nraw[nro(y, x)][1], nraw[nro(y + dy2, x + dx2)][1]);
        // cubed
        g1 *= g1 * g1;
        g2 *= g2 * g2;

        float eg;
        // current pixel value (green) * weighted average of diagonal up and down colorness (relative to green),
        // this cancels the green values restoring
        // the energy level of the red/blue interpolated value.
        // G*(g1*R_a/G_a + g2*R_b/G_b)/(g1+g2)
        eg = nraw[nro(y, x)][1] * (g1 * nraw[nro(y + dy, x + dx)][cl] / nraw[nro(y + dy, x + dx)][1] +
                                   g2 * nraw[nro(y + dy2, x + dx2)][cl] /
                                   nraw[nro(y + dy2, x + dx2)][1]) / (g1 + g2);
        float min, max;
        min = MIN(nraw[nro(y + dy, x + dx)][cl], nraw[nro(y + dy2, x + dx2)][cl]);
        max = MAX(nraw[nro(y + dy, x + dx)][cl], nraw[nro(y + dy2, x + dx2)][cl]);
        min /= 1.2f;
        max *= 1.2f;
        // if value is well outsides bounds of two adjacent pixels, scale back towards center
        if (eg < min)
            eg = scale_under(eg, min);
        else if (eg > max)
            eg = scale_over(eg, max);
        // don't let this value exceed the maximum value anywhere in the raw image
        if (eg > channel_maximum[cl])
            eg = channel_maximum[cl];
        else if (eg < channel_minimum[cl])
            eg = channel_minimum[cl];
        // at this point color is now a complete checkerboard, opposite the original green raw values
        nraw[nro(y, x)][cl] = eg;
    }
}

/*
 * интерполяция красных и синих в точках где был известен только зелёный,
 * направления горизонтальные или вертикальные
 *
 * interpolation of red and blue at points where only green was known, directions horizontal or vertical
 */
// step 4b: this fills in red and blue at the same time
void DHT::make_rbhv(int i) {
    int js = (COLOR(i, 0) & 1) ^ 1; // first raw green pixel
    for (int j = js; j < iwidth; j += 2) {
        int x = j + nr_leftmargin;
        int y = i + nr_topmargin;
        /*
         * поскольку сверху-снизу и справа-слева уже есть все необходимые красные и
         * синие, то можно выбрать наилучшее направление исходя из информации по
         * обоим цветам.
         */
        /*
         * since top-bottom and right-left already have all the necessary red and
         * blue, you can choose the best direction based on the information on
         * both colors.
         */
        int dx, dy, dx2, dy2;
        // choose vertical or horizontal
        if (ndir[nro(y, x)] & VER) {
            dx = dx2 = 0;
            dy = -1;
            dy2 = 1;
        } else {
            dy = dy2 = 0;
            dx = 1;
            dx2 = -1;
        }
        // weights to left-right or up down
        float g1 = 1 / calc_dist(nraw[nro(y, x)][1], nraw[nro(y + dy, x + dx)][1]);
        float g2 = 1 / calc_dist(nraw[nro(y, x)][1], nraw[nro(y + dy2, x + dx2)][1]);
        // squared
        g1 *= g1;
        g2 *= g2;
        float eg_r, eg_b;
        eg_r = nraw[nro(y, x)][1] *
               (g1 * nraw[nro(y + dy, x + dx)][0] / nraw[nro(y + dy, x + dx)][1] +
                g2 * nraw[nro(y + dy2, x + dx2)][0] / nraw[nro(y + dy2, x + dx2)][1]) / (g1 + g2);
        eg_b = nraw[nro(y, x)][1] *
               (g1 * nraw[nro(y + dy, x + dx)][2] / nraw[nro(y + dy, x + dx)][1] +
                g2 * nraw[nro(y + dy2, x + dx2)][2] / nraw[nro(y + dy2, x + dx2)][1]) / (g1 + g2);
        float min_r, max_r;
        min_r = MIN(nraw[nro(y + dy, x + dx)][0], nraw[nro(y + dy2, x + dx2)][0]);
        max_r = MAX(nraw[nro(y + dy, x + dx)][0], nraw[nro(y + dy2, x + dx2)][0]);
        float min_b, max_b;
        min_b = MIN(nraw[nro(y + dy, x + dx)][2], nraw[nro(y + dy2, x + dx2)][2]);
        max_b = MAX(nraw[nro(y + dy, x + dx)][2], nraw[nro(y + dy2, x + dx2)][2]);
        min_r /= 1.2f;
        max_r *= 1.2f;
        min_b /= 1.2f;
        max_b *= 1.2f;

        if (eg_r < min_r)
            eg_r = scale_under(eg_r, min_r);
        else if (eg_r > max_r)
            eg_r = scale_over(eg_r, max_r);
        if (eg_b < min_b)
            eg_b = scale_under(eg_b, min_b);
        else if (eg_b > max_b)
            eg_b = scale_over(eg_b, max_b);

        if (eg_r > channel_maximum[0])
            eg_r = channel_maximum[0];
        else if (eg_r < channel_minimum[0])
            eg_r = channel_minimum[0];
        if (eg_b > channel_maximum[2])
            eg_b = channel_maximum[2];
        else if (eg_b < channel_minimum[2])
            eg_b = channel_minimum[2];
        nraw[nro(y, x)][0] = eg_r;
        nraw[nro(y, x)][2] = eg_b;
    }
}



/*
 * перенос изображения в выходной массив
 *
 * transfer image to the output array
 */
// step 5
void DHT::copy_to_image() {
    int pidx;
    for (int i = 0; i < iheight; ++i)
        for (int j = 0; j < iwidth; ++j) {
            pidx = nro(i + nr_topmargin, j + nr_leftmargin);
            for (int cc = 0; cc < 3; cc++)
                (*imgdata[cc])(j, i) = nraw[pidx][cc] - FLOOR;
        }
}

DHT::~DHT() {
    free(nraw);
    free(ndir);
}

void dht_interpolate(pfs::Array2D *Xj, pfs::Array2D *Yj, pfs::Array2D *Zj) {
    DHT dht(Xj, Yj, Zj);
    dht.make_hv_dirs();
    dht.make_greens();
    dht.make_diag_dirs();
    dht.make_rb();
    dht.copy_to_image();
}
