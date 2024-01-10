/* -*- C++ -*-
 * Copyright 2023 Stephen Wasilewski
 * Modified from:
 *
 * File: dcraw_emu.cpp
 * Copyright 2008-2021 LibRaw LLC (info@libraw.org)
 * Created: Sun Mar 23,   2008
 *
 * LibRaw simple C++ API sample: almost complete dcraw emulator
 *

LibRaw is free software; you can redistribute it and/or modify
it under the terms of the one of two licenses as you choose:

1. GNU LESSER GENERAL PUBLIC LICENSE version 2.1
   (See file LICENSE.LGPL provided in LibRaw distribution archive for details).

2. COMMON DEVELOPMENT AND DISTRIBUTION LICENSE (CDDL) Version 1.0
   (See file LICENSE.CDDL provided in LibRaw distribution archive for details).
*/

/*
 * Copyright (c) 2023 Stephen Wasilewski, EPFL
 * primary modification is the elimation of some options and a change in
 * the default arguments targetting hdr generation from raw pixel values
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


#ifdef _MSC_VER
// suppress sprintf-related warning. sprintf() is permitted in sample code
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "libraw/libraw.h"

#ifndef LIBRAW_WIN32_CALLS
#include <sys/mman.h>
#include <sys/time.h>
#include <unistd.h>
#else
#include <io.h>
#endif
#include <fcntl.h>
#include <sys/stat.h>
#include <iostream>

#ifdef LIBRAW_WIN32_CALLS
#define snprintf _snprintf
#include <windows.h>
#else
#define O_BINARY 0
#endif

#ifdef USE_DNGSDK
#include "dng_host.h"
#include "dng_negative.h"
#include "dng_simple_image.h"
#include "dng_info.h"
#endif

void usage(const char *prog)
{
  printf("rawconvert: dcraw_emu fork that ensures consistent scaling of raw output as raw RGB\n");
  printf("IMPORTANT! you must set -S and -k for reliable results! use:\n");
  printf("\texiftool -AverageBlackLevel -LinearityUpperMargin (or equivalent metadata for your camera\n\n");
  printf("Usage:  %s [OPTION]... [FILE]...\n", prog);
  printf("-v        Verbose: print progress messages (repeated -v will add "
         "verbosity)\n"
         "-P <file> Fix the dead pixels listed in this file\n"
         "-K <file> Subtract dark frame (16-bit raw PGM)\n"
         "-k <num>  Set the darkness level\n"
         "-S <num>  Set the saturation level\n"
         "-T        Write TIFF instead of PPM (default)\n"
         "-p        Write PPM instead of TIFF\n"
         "-G        Use green_matching() filter\n"
         "-r <r g b g> Set custom white balance (used to normalize cam_xyz matrix)\n"
         "-B <x y w h> use cropbox\n"
         "-Z <suf>  Output filename generation rules\n"
         "          .suf => append .suf to input name, keeping existing suffix "
         "too\n"
         "           suf => replace input filename last extension\n"
         "          - => output to stdout\n"
         "          filename.suf => output to filename.suf\n"
         "-disinterp Do not run interpolation step\n"
         "-h         Half-size color image\n"
         "-q N      Set the interpolation quality:\n"
         "          0 - linear, 1 - VNG, 2 - PPG, 3 - AHD, 4 - DCB\n"
         "          11 - DHT, 12 - AAHD\n"
         "-identify skip all processing and simply print XYZ->CamRGB matrix\n"
  );
  exit(1);
}

static int verbosity = 0;
int cnt = 0;
int my_progress_callback(void *d, enum LibRaw_progress p, int iteration,
                         int expected)
{
  char *passed = (char *)(d ? d : "default string"); // data passed to callback
                                                     // at set_callback stage

  if (verbosity > 2) // verbosity set by repeat -v switches
  {
    printf("CB: %s  pass %d of %d (data passed=%s)\n", libraw_strprogress(p),
           iteration, expected, passed);
  }
  else if (iteration == 0) // 1st iteration of each step
    printf("Starting %s (expecting %d iterations)\n", libraw_strprogress(p),
           expected);
  else if (iteration == expected - 1)
    printf("%s finished\n", libraw_strprogress(p));

  ///    if(++cnt>10) return 1; // emulate user termination on 10-th callback
  ///    call

  return 0; // always return 0 to continue processing
}

// timer
#ifndef LIBRAW_WIN32_CALLS
static struct timeval start, end;
void timerstart(void) { gettimeofday(&start, NULL); }
void timerprint(const char *msg, const char *filename)
{
  gettimeofday(&end, NULL);
  float msec = (end.tv_sec - start.tv_sec) * 1000.0f +
               (end.tv_usec - start.tv_usec) / 1000.0f;
  printf("Timing: %s/%s: %6.3f msec\n", filename, msg, msec);
}
#else
LARGE_INTEGER start;
void timerstart(void) { QueryPerformanceCounter(&start); }
void timerprint(const char *msg, const char *filename)
{
  LARGE_INTEGER unit, end;
  QueryPerformanceCounter(&end);
  QueryPerformanceFrequency(&unit);
  float msec = (float)(end.QuadPart - start.QuadPart);
  msec /= (float)unit.QuadPart / 1000.0f;
  printf("Timing: %s/%s: %6.3f msec\n", filename, msg, msec);
}

#endif

struct file_mapping
{
	void *map;
	INT64 fsize;
#ifdef LIBRAW_WIN32_CALLS
	HANDLE fd, fd_map;
	file_mapping() : map(0), fsize(0), fd(INVALID_HANDLE_VALUE), fd_map(INVALID_HANDLE_VALUE){}
#else
	int  fd;
	file_mapping() : map(0), fsize(0), fd(-1){}
#endif
};

void create_mapping(struct file_mapping& data, const std::string& fn)
{
#ifdef LIBRAW_WIN32_CALLS
	std::wstring fpath(fn.begin(), fn.end());
	if ((data.fd = CreateFileW(fpath.c_str(), GENERIC_READ, FILE_SHARE_READ, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0)) == INVALID_HANDLE_VALUE) return;
	LARGE_INTEGER fs;
	if (!GetFileSizeEx(data.fd, &fs)) return;
	data.fsize = fs.QuadPart;
	if ((data.fd_map = ::CreateFileMapping(data.fd, 0, PAGE_READONLY, fs.HighPart, fs.LowPart, 0)) == INVALID_HANDLE_VALUE) return;
	data.map = MapViewOfFile(data.fd_map, FILE_MAP_READ, 0, 0, data.fsize);
#else
	struct stat stt;
	if ((data.fd = open(fn.c_str(), O_RDONLY)) < 0) return;
	if (fstat(data.fd, &stt) != 0) return;
	data.fsize = stt.st_size;
	data.map = mmap(0, data.fsize, PROT_READ | PROT_WRITE, MAP_PRIVATE, data.fd, 0);
	return;
#endif
}

void close_mapping(struct file_mapping& data)
{
#ifdef LIBRAW_WIN32_CALLS
	if (data.map) UnmapViewOfFile(data.map);
	if (data.fd_map != INVALID_HANDLE_VALUE) CloseHandle(data.fd_map);
	if (data.fd != INVALID_HANDLE_VALUE) CloseHandle(data.fd);
	data.map = 0;
	data.fsize = 0;
	data.fd = data.fd_map = INVALID_HANDLE_VALUE;
#else
	if (data.map)
		munmap(data.map, data.fsize);
	if (data.fd >= 0)
		close(data.fd);
	data.map = 0;
	data.fsize = 0;
	data.fd = -1;
#endif
}


int main(int argc, char *argv[])
{
  if (argc == 1)
    usage(argv[0]);

  LibRaw RawProcessor;
  int i, arg, c, ret;
  int identify = 0;
  char opm, opt, *cp, *sp;
  char *outext = NULL;
#ifdef OUT
#undef OUT
#endif
#define OUT RawProcessor.imgdata.params
  OUT.user_mul[1] = OUT.user_mul[3] = 1;
  OUT.user_mul[0] = OUT.user_mul[2] = 1;
  OUT.output_tiff = 1;
  OUT.user_qual = 11;
  argv[argc] = (char *)"";
  for (arg = 1; (((opm = argv[arg][0]) - 2) | 2) == '+';)
  {
    char *optstr = argv[arg];
    opt = argv[arg++][1];
    if ((cp = strchr(sp = (char *)"cnbrkStqmHABCgU", opt)) != 0)
      for (i = 0; i < "111411111144221"[cp - sp] - '0'; i++)
        if (!isdigit(argv[arg + i][0]) && !optstr[2])
        {
          fprintf(stderr, "Non-numeric argument to \"-%c\"\n", opt);
          return 1;
        }
    switch (opt)
    {
    case 'v':
      verbosity++;
      break;
    case 'i':
      identify = 1;
      break;
    case 'P':
      OUT.bad_pixels = argv[arg++];
      break;
    case 'K':
      OUT.dark_frame = argv[arg++];
      break;
    case 'r':
      for (c = 0; c < 4; c++)
        OUT.user_mul[c] = (float)atof(argv[arg++]);
      break;
    case 'q':
        OUT.user_qual = atoi(argv[arg++]);
        break;
    case 'k':
      OUT.user_black = atoi(argv[arg++]);
      break;
    case 'S':
      OUT.user_sat = atoi(argv[arg++]);
      break;
    case 'B':
      for (c = 0; c < 4; c++)
        OUT.cropbox[c] = atoi(argv[arg++]);
      break;
    case 'T':
      OUT.output_tiff = 1;
      break;
    case 'p':
        OUT.output_tiff = 0;
        break;
    case 'G':
        OUT.green_matching = 1;
        break;
    case 'h':
        OUT.half_size = 1;
        break;
    case 'Z':
      outext = strdup(argv[arg++]);
      break;
    case 'd':
        if (!strcmp(optstr, "-disinterp"))
            OUT.no_interpolation = 1;
        else
            fprintf(stderr, "Unknown option \"%s\".\n", argv[arg - 1]);
        break;
    default:
      fprintf(stderr, "Unknown option \"-%c\".\n", opt);
      break;
    }
  }
  OUT.adjust_maximum_thr = 0;
  OUT.output_color = 0;
  OUT.gamm[0] = OUT.gamm[1] = OUT.no_auto_bright = 1;
  OUT.output_bps = 16;

#ifndef LIBRAW_WIN32_CALLS
  putenv((char *)"TZ=UTC"); // dcraw compatibility, affects TIFF datestamp field
#else
  _putenv(
      (char *)"TZ=UTC"); // dcraw compatibility, affects TIFF datestamp field
#endif
#define P1 RawProcessor.imgdata.idata
#define S RawProcessor.imgdata.sizes
#define C RawProcessor.imgdata.color
#define T RawProcessor.imgdata.thumbnail

if (verbosity > 1)
RawProcessor.set_progress_handler(my_progress_callback,
                                  (void *)"Sample data passed");
#ifdef LIBRAW_USE_OPENMP
  if (verbosity)
    printf("Using %d threads\n", omp_get_max_threads());
#endif

  int done = 0;
  int total = argc - arg;
  for (; arg < argc; arg++)
  {
    char outfn[1024];

    if (verbosity)
      printf("Processing file %s\n", argv[arg]);


    ret = RawProcessor.open_file(argv[arg]);

    if (identify){
      fprintf(stdout, "RGBG multipliers:\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n", OUT.user_mul[0], OUT.user_mul[1], OUT.user_mul[2], OUT.user_mul[3]);
      fprintf(stdout, "XYZ->CamRGB:");
      for (int r = 0; r < P1.colors; r++)
          fprintf(stdout, "\t%6.4f\t%6.4f\t%6.4f", C.cam_xyz[r][0] * OUT.user_mul[r], C.cam_xyz[r][1] * OUT.user_mul[r], C.cam_xyz[r][2] * OUT.user_mul[r]);
      fprintf(stdout, "\nD65_multips:");
      for (int c = 0; c < P1.colors; c++)
          fprintf(stdout, "\t%f", C.pre_mul[c]);
      fprintf(stdout, "\n");
      return 0;
    }

    if (ret != LIBRAW_SUCCESS)
    {
    fprintf(stderr, "Cannot open %s: %s\n", argv[arg],
            libraw_strerror(ret));
    continue; // no recycle b/c open_file will recycle itself
    }


    if ((ret = RawProcessor.unpack()) != LIBRAW_SUCCESS)
    {
      fprintf(stderr, "Cannot unpack %s: %s\n", argv[arg],
              libraw_strerror(ret));
      continue;
    }

    if (OUT.user_black < 0) {
        std::cerr << "WARNING! -k not set, computed black point is: " << C.black << std::endl;
    }

    if (LIBRAW_SUCCESS != (ret = RawProcessor.dcraw_process()))
    {
      fprintf(stderr, "Cannot do postprocessing on %s: %s\n", argv[arg],
              libraw_strerror(ret));
      if (LIBRAW_FATAL_ERROR(ret))
        continue;
    }

      if (OUT.user_sat < 0) {
          std::cerr << "WARNING! -S not set, computed white point is: " << C.maximum << std::endl;
      }


    if (!outext)
      snprintf(outfn, sizeof(outfn), "%s.%s", argv[arg],
               OUT.output_tiff ? "tiff" : (P1.colors > 1 ? "ppm" : "pgm"));
    else if (!strcmp(outext, "-"))
      snprintf(outfn, sizeof(outfn), "-");
    else
    {
      if (*outext == '.') // append
        snprintf(outfn, sizeof(outfn), "%s%s", argv[arg], outext);
      else if (strchr(outext, '.') && *outext != '.') // dot is not 1st char
        strncpy(outfn, outext, sizeof(outfn));
      else
      {
        strncpy(outfn, argv[arg], sizeof(outfn));
        if (strlen(outfn) > 0)
        {
          char *lastchar = outfn + strlen(outfn); // points to term 0
          while (--lastchar > outfn)
          {
            if (*lastchar == '/' || *lastchar == '\\')
              break;
            if (*lastchar == '.')
            {
              *lastchar = 0;
              break;
            };
          }
        }
        strncat(outfn, ".", sizeof(outfn) - strlen(outfn) - 1);
        strncat(outfn, outext, sizeof(outfn) - strlen(outfn) - 1);
      }
    }

    if (verbosity)
    {
      printf("Writing file %s\n", outfn);
    }

    if (LIBRAW_SUCCESS != (ret = RawProcessor.dcraw_ppm_tiff_writer(outfn)))
      fprintf(stderr, "Cannot write %s: %s\n", outfn, libraw_strerror(ret));
    else
      done++;

	RawProcessor.recycle(); // just for show this call

  }
#ifdef USE_DNGSDK
  if (dnghost)
    delete dnghost;
#endif
  if (total == 0)
    return 1;
  if (done < total)
    return 2;
  return 0;
}
