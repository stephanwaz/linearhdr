/* This file includes code modified from the PFSTOOLS package.
 * Copyright (C) 2003,2004 Rafal Mantiuk and Grzegorz Krawczyk
 *
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

#include <iostream>
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <sstream>
#include <rgbeio.h>
#include <linearhdr.h>
#include <pfs.h>

#define PROG_NAME "convertcolor"

class QuietException 
{
};

void printHelp()
{
  fprintf( stderr, PROG_NAME " convert pfs stream to destination, options: \n"
                             "\t[--colorspace, -c]: output colorspace. choices: XYZ, RGB, SRGB, YUV, Yxy, PQYCbCr2020, YCbCr709, HLGYCbCr2020, RGB2020\n"
                             "\t[--hdr, -H]: output rgbe image format instead of pfs stream\n"
                             "\t[--tsv, -t]: output tsv data instead of pfs stream\n"
                             "\t[--raw, -r]: if --hdr do not apply 1/179 multiplier\n"
                             "\t[--demosaic, -d]: assumes bayer grid input, apply DHT demosaicing and then color transform. use XYZ to pass raw color out.\n"
                             "\t[--mtx, -m '<val> <val> <val> <val> <val> <val> <val> <val> <val>']: input to output color transform matrix, overrides --colorspace give in row major order\n"
                             "\t[--verbose, -v]: verbose messages\n\n");
}

bool verbose = false;

void ConvertColor( int argc, char* argv[] ) {
    pfs::DOMIO pfsio;

    int c;
    std::string csout = "RGB";
    pfs::ColorSpace cs_out;;
    bool data = false;
    bool hdr = false;
    bool radiance = true;
    bool usemtx = false;
    bool demosaic = false;
    float rgb2rgb[3][3] = {{0.0, 0.0, 0.0},
                           {0.0, 0.0, 0.0},
                           {0.0, 0.0, 0.0}};

    int width, height, size;

    static struct option cmdLineOptions[] = {
    { "help", no_argument, NULL, 'h' },
    { "colorspace", required_argument, NULL, 'c' },
    { "tsv", no_argument, NULL, 't' },
    { "hdr", no_argument, NULL, 'H' },
    { "raw", no_argument, NULL, 'r' },
    { "mtx", required_argument, NULL, 'm' },
    { "demosaic", no_argument, NULL, 'd' },
    { "verbose", no_argument, NULL, 'v' },
    { NULL, 0, NULL, 0 }
    };

    int optionIndex = 0;
    std::stringstream k; //to read in multivalue arguments
    int ci = 0;
    while ((c = getopt_long(argc, argv, "Hhrvdc:m:", cmdLineOptions, &optionIndex)) != -1) {
        ci++;
        if (strlen(argv[ci]) > 2 && argv[ci][1] != '-'){
            char message[100];
            snprintf(message, 100, "bad option : %s, all long options need --", argv[ci]);
            throw pfs::Exception(message);
        }
        switch( c ) {
            case 'h':
              printHelp();
              throw QuietException();
            case 'H':
                data = true;
                hdr = true;
                break;
            case 'c':
                ci++;
                csout.assign(optarg);
                break;
            case 't':
                data = true;
                break;
            case 'v':
                verbose = true;
                break;
            case 'r':
                radiance = false;
                break;
            case 'd':
                demosaic = true;
                break;
            case 'm':
                ci++;
                k.str("");
                k.clear();
                k << optarg;
                k >> rgb2rgb[0][0] >> rgb2rgb[0][1]  >> rgb2rgb[0][2];
                k >> rgb2rgb[1][0] >> rgb2rgb[1][1]  >> rgb2rgb[1][2];
                k >> rgb2rgb[2][0] >> rgb2rgb[2][1]  >> rgb2rgb[2][2];
                usemtx = true;
                break;
            default:
              throw QuietException();
        }
    }

    char a = csout.at(1);
    uint la = csout.length();

    switch( a ){
        case 'G':
            if (la == 3)
                cs_out = pfs::CS_RGB;
            else
                cs_out = pfs::CS_RGB2020;
            break;
        case 'Y':
            cs_out = pfs::CS_XYZ;
            break;
        case 'R':
            cs_out = pfs::CS_SRGB;
            break;
        case 'U':
            cs_out = pfs::CS_YUV;
            break;
        case 'x':
            cs_out = pfs::CS_Yxy;
            break;
        case 'Q':
            cs_out = pfs::CS_PQYCbCr2020;
            break;
        case 'C':
            cs_out = pfs::CS_YCbCr709;
            break;
        case 'L':
            cs_out = pfs::CS_HLGYCbCr2020;
            break;
        default:
            throw QuietException();
    }


  
  while( true ) {
      pfs::Frame *frame = pfsio.readFrame( stdin );
      if( frame == NULL ) break; // No more frames


      pfs::Channel *X, *Y, *Z;
      frame->getXYZChannels( X, Y, Z );

      width = Y->getCols();
      height = Y->getRows();
      size = width * height;

      pfs::Array2D *xyz[3] = {X, Y, Z};

    if (demosaic) {
        int g0 = first_non_zero(Y);
        int r0 = first_non_zero_row(X);
        VERBOSE_STR << "go: " << g0 << " r0: " << r0 << std::endl;
        dht_interpolate(X, Y, Z);
    }

    if (usemtx) {
        VERBOSE_STR << rgb2rgb[0][0] << " " << rgb2rgb[0][1] << " " << rgb2rgb[0][2] << std::endl;
        VERBOSE_STR << rgb2rgb[1][0] << " " << rgb2rgb[1][1] << " " << rgb2rgb[1][2] << std::endl;
        VERBOSE_STR << rgb2rgb[2][0] << " " << rgb2rgb[2][1] << " " << rgb2rgb[2][2] << std::endl;
        for (int j = 0; j < size; j++)
            apply_color_transform(j, xyz, rgb2rgb);

    } else {
        pfs::transformColorSpace( pfs::CS_XYZ, X, Y, Z, cs_out, X, Y, Z );
    }

    if (data) {
        if (hdr) {
            std::stringstream header;
            RGBEWriter writer( stdout, radiance );
            header << "RGB_CHANNELS= " << csout << std::endl;
            std::string hstring = header.str();
            writer.writeImage( X, Y, Z, hstring );
        } else {
            for (int i = 0; i < size; i++)
                std::cout << (*X)(i) << "\t" << (*Y)(i) << "\t" << (*Z)(i) << std::endl;
        }

    } else {
        pfsio.writeFrame( frame, stdout );
    }
    pfsio.freeFrame( frame );

  }
}

int main( int argc, char* argv[] )
{
  try {
      ConvertColor( argc, argv );
  }
  catch( pfs::Exception ex ) {
    fprintf( stderr, PROG_NAME " error: %s\n", ex.getMessage() );
    return EXIT_FAILURE;
  }        
  catch( QuietException  ex ) {
    return EXIT_FAILURE;
  }        
  return EXIT_SUCCESS;
}
