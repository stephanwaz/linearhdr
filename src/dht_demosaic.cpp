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

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

static inline float calc_dist(float c1, float c2)
{
  return c1 > c2 ? c1 / c2 : c2 / c1;
}

struct DHT
{
  int nr_height, nr_width, iwidth, iheight; //SW
  static const int nr_topmargin = 4, nr_leftmargin = 4;
  float (*nraw)[3];
  float channel_maximum[3]; //SW change to float
  float channel_minimum[3];
//  LibRaw &libraw; //SW don't use
  enum
  {
    HVSH = 1,
    HOR = 2,
    VER = 4,
    HORSH = HOR | HVSH,
    VERSH = VER | HVSH,
    DIASH = 8,
    LURD = 16,
    RULD = 32,
    LURDSH = LURD | DIASH,
    RULDSH = RULD | DIASH,
    HOT = 64
  };
  static inline float Thot(void) throw() { return 64.0f; }
  static inline float Tg(void) throw() { return 256.0f; }
  static inline float T(void) throw() { return 1.4f; }
  char *ndir;
  inline int nr_offset(int row, int col) throw()
  {
    return (row * nr_width + col);
  }
  int get_hv_grb(int x, int y, int kc)
  {
    float hv1 = 2 * nraw[nr_offset(y - 1, x)][1] /
                (nraw[nr_offset(y - 2, x)][kc] + nraw[nr_offset(y, x)][kc]);
    float hv2 = 2 * nraw[nr_offset(y + 1, x)][1] /
                (nraw[nr_offset(y + 2, x)][kc] + nraw[nr_offset(y, x)][kc]);
    float kv = calc_dist(hv1, hv2) *
               calc_dist(nraw[nr_offset(y, x)][kc] * nraw[nr_offset(y, x)][kc],
                         (nraw[nr_offset(y - 2, x)][kc] *
                          nraw[nr_offset(y + 2, x)][kc]));
    kv *= kv;
    kv *= kv;
    kv *= kv;
    float dv =
        kv *
        calc_dist(nraw[nr_offset(y - 3, x)][1] * nraw[nr_offset(y + 3, x)][1],
                  nraw[nr_offset(y - 1, x)][1] * nraw[nr_offset(y + 1, x)][1]);
    float hh1 = 2 * nraw[nr_offset(y, x - 1)][1] /
                (nraw[nr_offset(y, x - 2)][kc] + nraw[nr_offset(y, x)][kc]);
    float hh2 = 2 * nraw[nr_offset(y, x + 1)][1] /
                (nraw[nr_offset(y, x + 2)][kc] + nraw[nr_offset(y, x)][kc]);
    float kh = calc_dist(hh1, hh2) *
               calc_dist(nraw[nr_offset(y, x)][kc] * nraw[nr_offset(y, x)][kc],
                         (nraw[nr_offset(y, x - 2)][kc] *
                          nraw[nr_offset(y, x + 2)][kc]));
    kh *= kh;
    kh *= kh;
    kh *= kh;
    float dh =
        kh *
        calc_dist(nraw[nr_offset(y, x - 3)][1] * nraw[nr_offset(y, x + 3)][1],
                  nraw[nr_offset(y, x - 1)][1] * nraw[nr_offset(y, x + 1)][1]);
    float e = calc_dist(dh, dv);
    char d = dh < dv ? (e > Tg() ? HORSH : HOR) : (e > Tg() ? VERSH : VER);
    return d;
  }
  int get_hv_rbg(int x, int y, int hc)
  {
    float hv1 = 2 * nraw[nr_offset(y - 1, x)][hc ^ 2] /
                (nraw[nr_offset(y - 2, x)][1] + nraw[nr_offset(y, x)][1]);
    float hv2 = 2 * nraw[nr_offset(y + 1, x)][hc ^ 2] /
                (nraw[nr_offset(y + 2, x)][1] + nraw[nr_offset(y, x)][1]);
    float kv = calc_dist(hv1, hv2) *
               calc_dist(nraw[nr_offset(y, x)][1] * nraw[nr_offset(y, x)][1],
                         (nraw[nr_offset(y - 2, x)][1] *
                          nraw[nr_offset(y + 2, x)][1]));
    kv *= kv;
    kv *= kv;
    kv *= kv;
    float dv = kv * calc_dist(nraw[nr_offset(y - 3, x)][hc ^ 2] *
                                  nraw[nr_offset(y + 3, x)][hc ^ 2],
                              nraw[nr_offset(y - 1, x)][hc ^ 2] *
                                  nraw[nr_offset(y + 1, x)][hc ^ 2]);
    float hh1 = 2 * nraw[nr_offset(y, x - 1)][hc] /
                (nraw[nr_offset(y, x - 2)][1] + nraw[nr_offset(y, x)][1]);
    float hh2 = 2 * nraw[nr_offset(y, x + 1)][hc] /
                (nraw[nr_offset(y, x + 2)][1] + nraw[nr_offset(y, x)][1]);
    float kh = calc_dist(hh1, hh2) *
               calc_dist(nraw[nr_offset(y, x)][1] * nraw[nr_offset(y, x)][1],
                         (nraw[nr_offset(y, x - 2)][1] *
                          nraw[nr_offset(y, x + 2)][1]));
    kh *= kh;
    kh *= kh;
    kh *= kh;
    float dh =
        kh * calc_dist(
                 nraw[nr_offset(y, x - 3)][hc] * nraw[nr_offset(y, x + 3)][hc],
                 nraw[nr_offset(y, x - 1)][hc] * nraw[nr_offset(y, x + 1)][hc]);
    float e = calc_dist(dh, dv);
    char d = dh < dv ? (e > Tg() ? HORSH : HOR) : (e > Tg() ? VERSH : VER);
    return d;
  }
  int get_diag_grb(int x, int y, int kc)
  {
    float hlu =
        nraw[nr_offset(y - 1, x - 1)][1] / nraw[nr_offset(y - 1, x - 1)][kc];
    float hrd =
        nraw[nr_offset(y + 1, x + 1)][1] / nraw[nr_offset(y + 1, x + 1)][kc];
    float dlurd =
        calc_dist(hlu, hrd) *
        calc_dist(nraw[nr_offset(y - 1, x - 1)][1] *
                      nraw[nr_offset(y + 1, x + 1)][1],
                  nraw[nr_offset(y, x)][1] * nraw[nr_offset(y, x)][1]);
    float druld =
        calc_dist(hlu, hrd) *
        calc_dist(nraw[nr_offset(y - 1, x + 1)][1] *
                      nraw[nr_offset(y + 1, x - 1)][1],
                  nraw[nr_offset(y, x)][1] * nraw[nr_offset(y, x)][1]);
    float e = calc_dist(dlurd, druld);
    char d =
        druld < dlurd ? (e > T() ? RULDSH : RULD) : (e > T() ? LURDSH : LURD);
    return d;
  }
  int get_diag_rbg(int x, int y, int /* hc */)
  {
    float dlurd = calc_dist(
        nraw[nr_offset(y - 1, x - 1)][1] * nraw[nr_offset(y + 1, x + 1)][1],
        nraw[nr_offset(y, x)][1] * nraw[nr_offset(y, x)][1]);
    float druld = calc_dist(
        nraw[nr_offset(y - 1, x + 1)][1] * nraw[nr_offset(y + 1, x - 1)][1],
        nraw[nr_offset(y, x)][1] * nraw[nr_offset(y, x)][1]);
    float e = calc_dist(dlurd, druld);
    char d =
        druld < dlurd ? (e > T() ? RULDSH : RULD) : (e > T() ? LURDSH : LURD);
    return d;
  }
  static inline float scale_over(float ec, float base)
  {
    float s = base * .4;
    float o = ec - base;
    return base + sqrt(s * (o + s)) - s;
  }
  static inline float scale_under(float ec, float base)
  {
    float s = base * .6;
    float o = base - ec;
    return base - sqrt(s * (o + s)) + s;
  }
  ~DHT();
  DHT(pfs::Array2D *Xj, pfs::Array2D *Yj, pfs::Array2D *Zj); //SW base on pfs::Array
  void copy_to_image();
  void make_greens();
  void make_diag_dirs();
  void make_hv_dirs();
  void refine_hv_dirs(int i, int js);
  void refine_ihv_dirs(int i);
  void refine_idiag_dirs(int i);
  void illustrate_dline(int i);
  void make_hv_dline(int i);
  void make_diag_dline(int i);
  void make_gline(int i);
  void make_rbdiag(int i);
  void make_rbhv(int i);
  void make_rb();
  void hide_hots();
  void restore_hots();
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
 * * you have to make sure that you don't get a 0 when rounding, otherwise you'll have a problem when interpreting blue and red colors.
 *
 */
int DHT::COLOR(int i, int j){
    int ci = 0;
    float mv = (*imgdata[0])(j, i);
    if ((*imgdata[1])(j, i) > mv){
        ci = 1;
        mv = (*imgdata[1])(j, i);
    };
    if ((*imgdata[2])(j, i) > mv)
        ci = 2;
    return ci;
}


DHT::DHT(pfs::Array2D *Xj, pfs::Array2D *Yj, pfs::Array2D *Zj)
{
  imgdata[0] = Xj;
  imgdata[1] = Yj;
  imgdata[2] = Zj;
  iwidth = imgdata[0]->getCols();
  iheight = imgdata[0]->getRows();
  nr_height = iheight + nr_topmargin * 2;
  nr_width = iwidth + nr_leftmargin * 2;
  nraw = (float3 *)malloc(nr_height * nr_width * sizeof(float3));

  ndir = (char *)calloc(nr_height * nr_width, 1);
  channel_maximum[0] = channel_maximum[1] = channel_maximum[2] = 0;
  channel_minimum[0] = (*imgdata[0])(0);
  channel_minimum[1] = (*imgdata[1])(0);
  channel_minimum[2] = (*imgdata[2])(0);
  for (int i = 0; i < nr_height * nr_width; ++i)
    nraw[i][0] = nraw[i][1] = nraw[i][2] = 0.5;
  for (int i = 0; i < iheight; ++i)
  {
    int col_cache[48];
    for (int j = 0; j < 48; ++j)
    {
      int l = COLOR(i, j);
//      if (l == 3)
//        l = 1;
      col_cache[j] = l;
    }
    for (int j = 0; j < iwidth; ++j)
    {
      int l = col_cache[j % 48];
      float c = (*imgdata[l])(j, i);
      if (c != 0)
      {
        if (channel_maximum[l] < c)
          channel_maximum[l] = c;
        if (channel_minimum[l] > c)
          channel_minimum[l] = c;
        nraw[nr_offset(i + nr_topmargin, j + nr_leftmargin)][l] = (float)c;
      }
    }
  }
//  channel_minimum[0] += 1e-5;
//  channel_minimum[1] += 1e-5;
//  channel_minimum[2] += 1e-5;
}

void DHT::hide_hots()
{
  //SW int iwidth = libraw.imgdata.sizes.iwidth;
#if defined(LIBRAW_USE_OPENMP)
#pragma omp parallel for schedule(guided) firstprivate(iwidth)
#endif
  for (int i = 0; i < iheight; ++i)
  {
    int js = COLOR(i, 0) & 1;
    int kc = COLOR(i, js);
    /*
     * js -- начальная х-координата, которая попадает мимо известного зелёного
     * kc -- известный цвет в точке интерполирования
     */
      /*
       * js -- an initial x-coordinate that gets past the known green
       * kc -- known color at the interpolation point
       */
    for (int j = js; j < iwidth; j += 2)
    {
      int x = j + nr_leftmargin;
      int y = i + nr_topmargin;
      float c = nraw[nr_offset(y, x)][kc];
      if ((c > nraw[nr_offset(y, x + 2)][kc] &&
           c > nraw[nr_offset(y, x - 2)][kc] &&
           c > nraw[nr_offset(y - 2, x)][kc] &&
           c > nraw[nr_offset(y + 2, x)][kc] &&
           c > nraw[nr_offset(y, x + 1)][1] &&
           c > nraw[nr_offset(y, x - 1)][1] &&
           c > nraw[nr_offset(y - 1, x)][1] &&
           c > nraw[nr_offset(y + 1, x)][1]) ||
          (c < nraw[nr_offset(y, x + 2)][kc] &&
           c < nraw[nr_offset(y, x - 2)][kc] &&
           c < nraw[nr_offset(y - 2, x)][kc] &&
           c < nraw[nr_offset(y + 2, x)][kc] &&
           c < nraw[nr_offset(y, x + 1)][1] &&
           c < nraw[nr_offset(y, x - 1)][1] &&
           c < nraw[nr_offset(y - 1, x)][1] &&
           c < nraw[nr_offset(y + 1, x)][1]))
      {
        float avg = 0;
        for (int k = -2; k < 3; k += 2)
          for (int m = -2; m < 3; m += 2)
            if (m == 0 && k == 0)
              continue;
            else
              avg += nraw[nr_offset(y + k, x + m)][kc];
        avg /= 8;
        //				float dev = 0;
        //				for (int k = -2; k < 3; k += 2)
        //					for (int l = -2; l < 3; l += 2)
        //						if (k == 0 && l == 0)
        //							continue;
        //						else {
        //							float t = nraw[nr_offset(y + k, x + l)][kc] -
        //avg; 							dev += t * t;
        //						}
        //				dev /= 8;
        //				dev = sqrt(dev);
        if (calc_dist(c, avg) > Thot())
        {
          ndir[nr_offset(y, x)] |= HOT;
          float dv = calc_dist(
              nraw[nr_offset(y - 2, x)][kc] * nraw[nr_offset(y - 1, x)][1],
              nraw[nr_offset(y + 2, x)][kc] * nraw[nr_offset(y + 1, x)][1]);
          float dh = calc_dist(
              nraw[nr_offset(y, x - 2)][kc] * nraw[nr_offset(y, x - 1)][1],
              nraw[nr_offset(y, x + 2)][kc] * nraw[nr_offset(y, x + 1)][1]);
          if (dv > dh)
            nraw[nr_offset(y, x)][kc] = (nraw[nr_offset(y, x + 2)][kc] +
                                         nraw[nr_offset(y, x - 2)][kc]) /
                                        2;
          else
            nraw[nr_offset(y, x)][kc] = (nraw[nr_offset(y - 2, x)][kc] +
                                         nraw[nr_offset(y + 2, x)][kc]) /
                                        2;
        }
      }
    }
    for (int j = js ^ 1; j < iwidth; j += 2)
    {
      int x = j + nr_leftmargin;
      int y = i + nr_topmargin;
      float c = nraw[nr_offset(y, x)][1];
      if ((c > nraw[nr_offset(y, x + 2)][1] &&
           c > nraw[nr_offset(y, x - 2)][1] &&
           c > nraw[nr_offset(y - 2, x)][1] &&
           c > nraw[nr_offset(y + 2, x)][1] &&
           c > nraw[nr_offset(y, x + 1)][kc] &&
           c > nraw[nr_offset(y, x - 1)][kc] &&
           c > nraw[nr_offset(y - 1, x)][kc ^ 2] &&
           c > nraw[nr_offset(y + 1, x)][kc ^ 2]) ||
          (c < nraw[nr_offset(y, x + 2)][1] &&
           c < nraw[nr_offset(y, x - 2)][1] &&
           c < nraw[nr_offset(y - 2, x)][1] &&
           c < nraw[nr_offset(y + 2, x)][1] &&
           c < nraw[nr_offset(y, x + 1)][kc] &&
           c < nraw[nr_offset(y, x - 1)][kc] &&
           c < nraw[nr_offset(y - 1, x)][kc ^ 2] &&
           c < nraw[nr_offset(y + 1, x)][kc ^ 2]))
      {
        float avg = 0;
        for (int k = -2; k < 3; k += 2)
          for (int m = -2; m < 3; m += 2)
            if (k == 0 && m == 0)
              continue;
            else
              avg += nraw[nr_offset(y + k, x + m)][1];
        avg /= 8;
        //				float dev = 0;
        //				for (int k = -2; k < 3; k += 2)
        //					for (int l = -2; l < 3; l += 2)
        //						if (k == 0 && l == 0)
        //							continue;
        //						else {
        //							float t = nraw[nr_offset(y + k, x + l)][1] -
        //avg; 							dev += t * t;
        //						}
        //				dev /= 8;
        //				dev = sqrt(dev);
        if (calc_dist(c, avg) > Thot())
        {
          ndir[nr_offset(y, x)] |= HOT;
          float dv = calc_dist(
              nraw[nr_offset(y - 2, x)][1] * nraw[nr_offset(y - 1, x)][kc ^ 2],
              nraw[nr_offset(y + 2, x)][1] * nraw[nr_offset(y + 1, x)][kc ^ 2]);
          float dh = calc_dist(
              nraw[nr_offset(y, x - 2)][1] * nraw[nr_offset(y, x - 1)][kc],
              nraw[nr_offset(y, x + 2)][1] * nraw[nr_offset(y, x + 1)][kc]);
          if (dv > dh)
            nraw[nr_offset(y, x)][1] =
                (nraw[nr_offset(y, x + 2)][1] + nraw[nr_offset(y, x - 2)][1]) /
                2;
          else
            nraw[nr_offset(y, x)][1] =
                (nraw[nr_offset(y - 2, x)][1] + nraw[nr_offset(y + 2, x)][1]) /
                2;
        }
      }
    }
  }
}

void DHT::restore_hots()
{
  //SW int iwidth = libraw.imgdata.sizes.iwidth;
#if defined(LIBRAW_USE_OPENMP)
#ifdef _MSC_VER
#pragma omp parallel for firstprivate(iwidth)
#else
#pragma omp parallel for schedule(guided) firstprivate(iwidth) collapse(2)
#endif
#endif
  for (int i = 0; i < iheight; ++i)
  {
    for (int j = 0; j < iwidth; ++j)
    {
      int x = j + nr_leftmargin;
      int y = i + nr_topmargin;
      if (ndir[nr_offset(y, x)] & HOT)
      {
        int l = COLOR(i, j);
        nraw[nr_offset(i + nr_topmargin, j + nr_leftmargin)][l] = (*imgdata[l])(j, i);
      }
    }
  }
}

void DHT::make_diag_dirs()
{
#if defined(LIBRAW_USE_OPENMP)
#pragma omp parallel for schedule(guided)
#endif
  for (int i = 0; i < iheight; ++i)
  {
    make_diag_dline(i);
  }
//#if defined(LIBRAW_USE_OPENMP)
//#pragma omp parallel for schedule(guided)
//#endif
//	for (int i = 0; i < libraw.imgdata.sizes.iheight; ++i) {
//		refine_diag_dirs(i, i & 1);
//	}
//#if defined(LIBRAW_USE_OPENMP)
//#pragma omp parallel for schedule(guided)
//#endif
//	for (int i = 0; i < libraw.imgdata.sizes.iheight; ++i) {
//		refine_diag_dirs(i, (i & 1) ^ 1);
//	}
#if defined(LIBRAW_USE_OPENMP)
#pragma omp parallel for schedule(guided)
#endif
  for (int i = 0; i < iheight; ++i)
  {
    refine_idiag_dirs(i);
  }
}

void DHT::make_hv_dirs()
{
#if defined(LIBRAW_USE_OPENMP)
#pragma omp parallel for schedule(guided)
#endif
  for (int i = 0; i < iheight; ++i)
  {
    make_hv_dline(i);
  }
#if defined(LIBRAW_USE_OPENMP)
#pragma omp parallel for schedule(guided)
#endif
  for (int i = 0; i < iheight; ++i)
  {
    refine_hv_dirs(i, i & 1);
  }
#if defined(LIBRAW_USE_OPENMP)
#pragma omp parallel for schedule(guided)
#endif
  for (int i = 0; i < iheight; ++i)
  {
    refine_hv_dirs(i, (i & 1) ^ 1);
  }
#if defined(LIBRAW_USE_OPENMP)
#pragma omp parallel for schedule(guided)
#endif
  for (int i = 0; i < iheight; ++i)
  {
    refine_ihv_dirs(i);
  }
}

void DHT::refine_hv_dirs(int i, int js)
{
  //SW int iwidth = libraw.imgdata.sizes.iwidth;
  for (int j = js; j < iwidth; j += 2)
  {
    int x = j + nr_leftmargin;
    int y = i + nr_topmargin;
    if (ndir[nr_offset(y, x)] & HVSH)
      continue;
    int nv =
        (ndir[nr_offset(y - 1, x)] & VER) + (ndir[nr_offset(y + 1, x)] & VER) +
        (ndir[nr_offset(y, x - 1)] & VER) + (ndir[nr_offset(y, x + 1)] & VER);
    int nh =
        (ndir[nr_offset(y - 1, x)] & HOR) + (ndir[nr_offset(y + 1, x)] & HOR) +
        (ndir[nr_offset(y, x - 1)] & HOR) + (ndir[nr_offset(y, x + 1)] & HOR);
    bool codir = (ndir[nr_offset(y, x)] & VER)
                     ? ((ndir[nr_offset(y - 1, x)] & VER) ||
                        (ndir[nr_offset(y + 1, x)] & VER))
                     : ((ndir[nr_offset(y, x - 1)] & HOR) ||
                        (ndir[nr_offset(y, x + 1)] & HOR));
    nv /= VER;
    nh /= HOR;
    if ((ndir[nr_offset(y, x)] & VER) && (nh > 2 && !codir))
    {
      ndir[nr_offset(y, x)] &= ~VER;
      ndir[nr_offset(y, x)] |= HOR;
    }
    if ((ndir[nr_offset(y, x)] & HOR) && (nv > 2 && !codir))
    {
      ndir[nr_offset(y, x)] &= ~HOR;
      ndir[nr_offset(y, x)] |= VER;
    }
  }
}

void DHT::refine_ihv_dirs(int i)
{
  //SW int iwidth = libraw.imgdata.sizes.iwidth;
  for (int j = 0; j < iwidth; j++)
  {
    int x = j + nr_leftmargin;
    int y = i + nr_topmargin;
    if (ndir[nr_offset(y, x)] & HVSH)
      continue;
    int nv =
        (ndir[nr_offset(y - 1, x)] & VER) + (ndir[nr_offset(y + 1, x)] & VER) +
        (ndir[nr_offset(y, x - 1)] & VER) + (ndir[nr_offset(y, x + 1)] & VER);
    int nh =
        (ndir[nr_offset(y - 1, x)] & HOR) + (ndir[nr_offset(y + 1, x)] & HOR) +
        (ndir[nr_offset(y, x - 1)] & HOR) + (ndir[nr_offset(y, x + 1)] & HOR);
    nv /= VER;
    nh /= HOR;
    if ((ndir[nr_offset(y, x)] & VER) && nh > 3)
    {
      ndir[nr_offset(y, x)] &= ~VER;
      ndir[nr_offset(y, x)] |= HOR;
    }
    if ((ndir[nr_offset(y, x)] & HOR) && nv > 3)
    {
      ndir[nr_offset(y, x)] &= ~HOR;
      ndir[nr_offset(y, x)] |= VER;
    }
  }
}
void DHT::make_hv_dline(int i)
{
  //SW int iwidth = libraw.imgdata.sizes.iwidth;
  int js = COLOR(i, 0) & 1;
  int kc = COLOR(i, js);
  /*
   * js -- начальная х-координата, которая попадает мимо известного зелёного
   * kc -- известный цвет в точке интерполирования
   */
  for (int j = 0; j < iwidth; j++)
  {
    int x = j + nr_leftmargin;
    int y = i + nr_topmargin;
    char d = 0;
    if ((j & 1) == js)
    {
      d = get_hv_grb(x, y, kc);
    }
    else
    {
      d = get_hv_rbg(x, y, kc);
    }
    ndir[nr_offset(y, x)] |= d;
  }
}

void DHT::make_diag_dline(int i)
{
  //SW int iwidth = libraw.imgdata.sizes.iwidth;
  int js = COLOR(i, 0) & 1;
  int kc = COLOR(i, js);
  /*
   * js -- начальная х-координата, которая попадает мимо известного зелёного
   * kc -- известный цвет в точке интерполирования
   */
  for (int j = 0; j < iwidth; j++)
  {
    int x = j + nr_leftmargin;
    int y = i + nr_topmargin;
    char d = 0;
    if ((j & 1) == js)
    {
      d = get_diag_grb(x, y, kc);
    }
    else
    {
      d = get_diag_rbg(x, y, kc);
    }
    ndir[nr_offset(y, x)] |= d;
  }
}

void DHT::refine_idiag_dirs(int i)
{
  //SW int iwidth = libraw.imgdata.sizes.iwidth;
  for (int j = 0; j < iwidth; j++)
  {
    int x = j + nr_leftmargin;
    int y = i + nr_topmargin;
    if (ndir[nr_offset(y, x)] & DIASH)
      continue;
    int nv = (ndir[nr_offset(y - 1, x)] & LURD) +
             (ndir[nr_offset(y + 1, x)] & LURD) +
             (ndir[nr_offset(y, x - 1)] & LURD) +
             (ndir[nr_offset(y, x + 1)] & LURD) +
             (ndir[nr_offset(y - 1, x - 1)] & LURD) +
             (ndir[nr_offset(y - 1, x + 1)] & LURD) +
             (ndir[nr_offset(y + 1, x - 1)] & LURD) +
             (ndir[nr_offset(y + 1, x + 1)] & LURD);
    int nh = (ndir[nr_offset(y - 1, x)] & RULD) +
             (ndir[nr_offset(y + 1, x)] & RULD) +
             (ndir[nr_offset(y, x - 1)] & RULD) +
             (ndir[nr_offset(y, x + 1)] & RULD) +
             (ndir[nr_offset(y - 1, x - 1)] & RULD) +
             (ndir[nr_offset(y - 1, x + 1)] & RULD) +
             (ndir[nr_offset(y + 1, x - 1)] & RULD) +
             (ndir[nr_offset(y + 1, x + 1)] & RULD);
    nv /= LURD;
    nh /= RULD;
    if ((ndir[nr_offset(y, x)] & LURD) && nh > 7)
    {
      ndir[nr_offset(y, x)] &= ~LURD;
      ndir[nr_offset(y, x)] |= RULD;
    }
    if ((ndir[nr_offset(y, x)] & RULD) && nv > 7)
    {
      ndir[nr_offset(y, x)] &= ~RULD;
      ndir[nr_offset(y, x)] |= LURD;
    }
  }
}

/*
 * вычисление недостающих зелёных точек.
 */
void DHT::make_greens()
{
#if defined(LIBRAW_USE_OPENMP)
#pragma omp parallel for schedule(guided)
#endif
  for (int i = 0; i < iheight; ++i)
  {
    make_gline(i);
  }
}

void DHT::make_gline(int i)
{
  //SW int iwidth = libraw.imgdata.sizes.iwidth;
  int js = COLOR(i, 0) & 1;
  int kc = COLOR(i, js);
  /*
   * js -- начальная х-координата, которая попадает мимо известного зелёного
   * kc -- известный цвет в точке интерполирования
   */
  for (int j = js; j < iwidth; j += 2)
  {
    int x = j + nr_leftmargin;
    int y = i + nr_topmargin;
    int dx, dy, dx2, dy2;
    float h1, h2;
    if (ndir[nr_offset(y, x)] & VER)
    {
      dx = dx2 = 0;
      dy = -1;
      dy2 = 1;
      h1 = 2 * nraw[nr_offset(y - 1, x)][1] /
           (nraw[nr_offset(y - 2, x)][kc] + nraw[nr_offset(y, x)][kc]);
      h2 = 2 * nraw[nr_offset(y + 1, x)][1] /
           (nraw[nr_offset(y + 2, x)][kc] + nraw[nr_offset(y, x)][kc]);
    }
    else
    {
      dy = dy2 = 0;
      dx = 1;
      dx2 = -1;
      h1 = 2 * nraw[nr_offset(y, x + 1)][1] /
           (nraw[nr_offset(y, x + 2)][kc] + nraw[nr_offset(y, x)][kc]);
      h2 = 2 * nraw[nr_offset(y, x - 1)][1] /
           (nraw[nr_offset(y, x - 2)][kc] + nraw[nr_offset(y, x)][kc]);
    }
    float b1 = 1 / calc_dist(nraw[nr_offset(y, x)][kc],
                             nraw[nr_offset(y + dy * 2, x + dx * 2)][kc]);
    float b2 = 1 / calc_dist(nraw[nr_offset(y, x)][kc],
                             nraw[nr_offset(y + dy2 * 2, x + dx2 * 2)][kc]);
    b1 *= b1;
    b2 *= b2;
    float eg = nraw[nr_offset(y, x)][kc] * (b1 * h1 + b2 * h2) / (b1 + b2);
    float min, max;
    min = MIN(nraw[nr_offset(y + dy, x + dx)][1],
              nraw[nr_offset(y + dy2, x + dx2)][1]);
    max = MAX(nraw[nr_offset(y + dy, x + dx)][1],
              nraw[nr_offset(y + dy2, x + dx2)][1]);
    min /= 1.2f;
    max *= 1.2f;
    if (eg < min)
      eg = scale_under(eg, min);
    else if (eg > max)
      eg = scale_over(eg, max);
    if (eg > channel_maximum[1])
      eg = channel_maximum[1];
    else if (eg < channel_minimum[1])
      eg = channel_minimum[1];
    nraw[nr_offset(y, x)][1] = eg;
  }
}

void DHT::illustrate_dline(int i)
{
  //SW int iwidth = libraw.imgdata.sizes.iwidth;
  for (int j = 0; j < iwidth; j++)
  {
    int x = j + nr_leftmargin;
    int y = i + nr_topmargin;
    nraw[nr_offset(y, x)][0] = nraw[nr_offset(y, x)][1] =
        nraw[nr_offset(y, x)][2] = 0.5;
    int l = ndir[nr_offset(y, x)] & 8;
    // l >>= 3; // WTF?
    l = 1;
    if (ndir[nr_offset(y, x)] & HOT)
      nraw[nr_offset(y, x)][0] =
          l * channel_maximum[0] / 4 + channel_maximum[0] / 4;
    else
      nraw[nr_offset(y, x)][2] =
          l * channel_maximum[2] / 4 + channel_maximum[2] / 4;
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

void DHT::make_rbdiag(int i)
{
  //SW int iwidth = libraw.imgdata.sizes.iwidth;
  int js = COLOR(i, 0) & 1;
  int uc = COLOR(i, js);
  int cl = uc ^ 2;
  /*
   * js -- начальная х-координата, которая попадает на уже интерполированный
   * зелёный al -- известный цвет (кроме зелёного) в точке интерполирования cl
   * -- неизвестный цвет
   */
  for (int j = js; j < iwidth; j += 2)
  {
    int x = j + nr_leftmargin;
    int y = i + nr_topmargin;
    int dx, dy, dx2, dy2;
    if (ndir[nr_offset(y, x)] & LURD)
    {
      dx = -1;
      dx2 = 1;
      dy = -1;
      dy2 = 1;
    }
    else
    {
      dx = -1;
      dx2 = 1;
      dy = 1;
      dy2 = -1;
    }
    float g1 = 1 / calc_dist(nraw[nr_offset(y, x)][1],
                             nraw[nr_offset(y + dy, x + dx)][1]);
    float g2 = 1 / calc_dist(nraw[nr_offset(y, x)][1],
                             nraw[nr_offset(y + dy2, x + dx2)][1]);
    g1 *= g1 * g1;
    g2 *= g2 * g2;

    float eg;
    eg = nraw[nr_offset(y, x)][1] *
         (g1 * nraw[nr_offset(y + dy, x + dx)][cl] /
              nraw[nr_offset(y + dy, x + dx)][1] +
          g2 * nraw[nr_offset(y + dy2, x + dx2)][cl] /
              nraw[nr_offset(y + dy2, x + dx2)][1]) /
         (g1 + g2);
    float min, max;
    min = MIN(nraw[nr_offset(y + dy, x + dx)][cl],
              nraw[nr_offset(y + dy2, x + dx2)][cl]);
    max = MAX(nraw[nr_offset(y + dy, x + dx)][cl],
              nraw[nr_offset(y + dy2, x + dx2)][cl]);
    min /= 1.2f;
    max *= 1.2f;
    if (eg < min)
      eg = scale_under(eg, min);
    else if (eg > max)
      eg = scale_over(eg, max);
    if (eg > channel_maximum[cl])
      eg = channel_maximum[cl];
    else if (eg < channel_minimum[cl])
      eg = channel_minimum[cl];
    nraw[nr_offset(y, x)][cl] = eg;
  }
}

/*
 * интерполяция красных и синих в точках где был известен только зелёный,
 * направления горизонтальные или вертикальные
 *
 * interpolation of red and blue at points where only green was known, directions horizontal or vertical
 */

void DHT::make_rbhv(int i)
{
  //SW int iwidth = libraw.imgdata.sizes.iwidth;
  int js = (COLOR(i, 0) & 1) ^ 1;
  for (int j = js; j < iwidth; j += 2)
  {
    int x = j + nr_leftmargin;
    int y = i + nr_topmargin;
    /*
     * поскольку сверху-снизу и справа-слева уже есть все необходимые красные и
     * синие, то можно выбрать наилучшее направление исходя из информации по
     * обоим цветам.
     */
    int dx, dy, dx2, dy2;
    if (ndir[nr_offset(y, x)] & VER)
    {
      dx = dx2 = 0;
      dy = -1;
      dy2 = 1;
    }
    else
    {
      dy = dy2 = 0;
      dx = 1;
      dx2 = -1;
    }
    float g1 = 1 / calc_dist(nraw[nr_offset(y, x)][1],
                             nraw[nr_offset(y + dy, x + dx)][1]);
    float g2 = 1 / calc_dist(nraw[nr_offset(y, x)][1],
                             nraw[nr_offset(y + dy2, x + dx2)][1]);
    g1 *= g1;
    g2 *= g2;
    float eg_r, eg_b;
    eg_r = nraw[nr_offset(y, x)][1] *
           (g1 * nraw[nr_offset(y + dy, x + dx)][0] /
                nraw[nr_offset(y + dy, x + dx)][1] +
            g2 * nraw[nr_offset(y + dy2, x + dx2)][0] /
                nraw[nr_offset(y + dy2, x + dx2)][1]) /
           (g1 + g2);
    eg_b = nraw[nr_offset(y, x)][1] *
           (g1 * nraw[nr_offset(y + dy, x + dx)][2] /
                nraw[nr_offset(y + dy, x + dx)][1] +
            g2 * nraw[nr_offset(y + dy2, x + dx2)][2] /
                nraw[nr_offset(y + dy2, x + dx2)][1]) /
           (g1 + g2);
    float min_r, max_r;
    min_r = MIN(nraw[nr_offset(y + dy, x + dx)][0],
                nraw[nr_offset(y + dy2, x + dx2)][0]);
    max_r = MAX(nraw[nr_offset(y + dy, x + dx)][0],
                nraw[nr_offset(y + dy2, x + dx2)][0]);
    float min_b, max_b;
    min_b = MIN(nraw[nr_offset(y + dy, x + dx)][2],
                nraw[nr_offset(y + dy2, x + dx2)][2]);
    max_b = MAX(nraw[nr_offset(y + dy, x + dx)][2],
                nraw[nr_offset(y + dy2, x + dx2)][2]);
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
    nraw[nr_offset(y, x)][0] = eg_r;
    nraw[nr_offset(y, x)][2] = eg_b;
  }
}

void DHT::make_rb()
{
#if defined(LIBRAW_USE_OPENMP)
#pragma omp barrier
#pragma omp parallel for schedule(guided)
#endif
  for (int i = 0; i < iheight; ++i)
  {
    make_rbdiag(i);
  }
#if defined(LIBRAW_USE_OPENMP)
#pragma omp barrier
#pragma omp parallel for schedule(guided)
#endif
  for (int i = 0; i < iheight; ++i)
  {
    make_rbhv(i);
  }
}

/*
 * перенос изображения в выходной массив
 *
 * transfer image to the output array
 */
void DHT::copy_to_image()
{
  //SW int iwidth = libraw.imgdata.sizes.iwidth;
#if defined(LIBRAW_USE_OPENMP)
#ifdef _MSC_VER
#pragma omp parallel for
#else
#pragma omp parallel for schedule(guided) collapse(2)
#endif
#endif
  for (int i = 0; i < iheight; ++i)
  {
    for (int j = 0; j < iwidth; ++j)
    {
        (*imgdata[0])(j, i) = nraw[nr_offset(i + nr_topmargin, j + nr_leftmargin)][0];
        (*imgdata[1])(j, i) = nraw[nr_offset(i + nr_topmargin, j + nr_leftmargin)][1];
        (*imgdata[2])(j, i) = nraw[nr_offset(i + nr_topmargin, j + nr_leftmargin)][2];
    }
  }
}

DHT::~DHT()
{
  free(nraw);
  free(ndir);
}

void dht_interpolate(pfs::Array2D *Xj, pfs::Array2D *Yj, pfs::Array2D *Zj)
{
    DHT dht(Xj, Yj, Zj);
    dht.hide_hots();
    dht.make_hv_dirs();
    //	dht.illustrate_dirs();
    dht.make_greens();
    dht.make_diag_dirs();
    //	dht.illustrate_dirs();
    dht.make_rb();
    dht.restore_hots();
    dht.copy_to_image();
}
