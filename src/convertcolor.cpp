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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <sstream>
#include <rgbeio.h>

#include <pfs.h>

#define PROG_NAME "convertcolor"

class QuietException 
{
};

void printHelp()
{
  fprintf( stderr, PROG_NAME " convert pfs stream from XYZ to dest, options are: \n"
                             "XYZ, RGB, SRGB, YUV, Yxy, PQYCbCr2020, YCbCr709, HLGYCbCr2020, RGB2020\n" );
}

void ConvertColor( int argc, char* argv[] ) {
    pfs::DOMIO pfsio;


    bool verbose = false;
    int c;
    std::string csout = "RGB";
    pfs::ColorSpace cs_out;;
    bool data = false;
    bool hdr = false;
    bool radiance = true;

    static struct option cmdLineOptions[] = {
    { "help", no_argument, NULL, 'h' },
    { "verbose", no_argument, NULL, 'v' },
    { "colorspace", required_argument, NULL, 'c' },
    { "data", no_argument, NULL, 'd' },
    { "hdr", no_argument, NULL, 'H' },
    { "raw", no_argument, NULL, 'r' },
    { NULL, 0, NULL, 0 }
    };

    int optionIndex = 0;
    while ((c = getopt_long(argc, argv, "Hhrvdc:", cmdLineOptions, &optionIndex)) != -1) {

        switch( c ) {
            case 'h':
              printHelp();
              throw QuietException();
            case 'v':
              verbose = true;
              break;
            case 'H':
                data = true;
                hdr = true;
                break;
            case 'c':
              csout.assign(optarg);
              break;
            case 'd':
                data = true;
                break;
            case 'r':
                radiance = false;
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

//    VERBOSE_STR << cs_out << std::endl;


  
  while( true ) {
    pfs::Frame *frame = pfsio.readFrame( stdin );
    if( frame == NULL ) break; // No more frames


    pfs::Channel *X, *Y, *Z;
    frame->getXYZChannels( X, Y, Z );

    pfs::transformColorSpace( pfs::CS_XYZ, X, Y, Z, cs_out, X, Y, Z );
    if (data) {
        if (hdr) {
            std::stringstream header;
            RGBEWriter writer( stdout, radiance );
            header << "RGB_CHANNELS= " << csout << std::endl;
            std::string hstring = header.str();
            writer.writeImage( X, Y, Z, hstring );
        } else {
            int width = Y->getCols();
            int height = Y->getRows();
            int size = width * height;
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
