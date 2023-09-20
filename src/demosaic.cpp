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

#define PROG_NAME "demosaic"


class QuietException 
{
};

void printHelp()
{
  fprintf( stderr, PROG_NAME " for testing demosaic" );
}


//void AHD(float *sensor2xyz) {
//    /* This function is based on the AHD code from dcraw, which in turn
//       builds on work by Keigo Hirakawa, Thomas Parks, and Paul Lee. */
//    const int G = 1, tsize = 256;
//
//    struct DemosaicBuffer {
//        /* Horizontally and vertically interpolated sensor colors */
//        float3 rgb[2][tsize][tsize];
//
//        /* CIElab color values */
//        float3 cielab[2][tsize][tsize];
//
//        /* Homogeneity map */
//        uint8_t homo[2][tsize][tsize];
//    };
//
//    /* Temporary tile storage */
//    DemosaicBuffer *buffers = new DemosaicBuffer[omp_get_max_threads()];
//
//    cout << "AHD demosaicing .." << endl;
//
//    /* Allocate a big buffer for the interpolated colors */
//    image_demosaiced = new float3[width*height];
//
//    size_t offset = 0;
//    float maxvalue = 0;
//    for (size_t y=0; y<height; ++y) {
//        for (size_t x=0; x<width; ++x) {
//            float value = image_merged[offset];
//            image_demosaiced[offset][fc(x, y)] = value;
//            if (value > maxvalue)
//                maxvalue = value;
//            offset++;
//        }
//    }
//
//    /* The AHD implementation below doesn't interpolate colors on a 5-pixel wide
//       boundary region -> use a naive averaging method on this region instead. */
//    const size_t border = 5;
//    for (size_t y=0; y<height; ++y) {
//        for (size_t x=0; x<width; ++x) {
//            if (x == border && y >= border && y < height-border)
//                x = width-border; /* Jump over the center part of the image */
//
//            float binval[3] = {0, 0, 0};
//            int bincount[3] = {0, 0, 0};
//
//            for (size_t ys=y-1; ys != y+2; ++ys) {
//                for (size_t xs=x-1; xs != x+2; ++xs) {
//                    if (ys < height && xs < width) {
//                        int col = fc(xs, ys);
//                        binval[col] += image_demosaiced[ys*width+xs][col];
//                        ++bincount[col];
//                    }
//                }
//            }
//
//            int col = fc(x, y);
//            for (int c=0; c<3; ++c) {
//                if (col != c)
//                    image_demosaiced[y*width+x][c] = bincount[c] ? (binval[c]/bincount[c]) : 1.0f;
//            }
//        }
//    }
//
//    /* Matrix that goes from sensor to normalized XYZ tristimulus values */
//    float sensor2xyz_n[3][3], sensor2xyz_n_maxvalue = 0;
//    const float d65_white[3] = { 0.950456, 1, 1.088754 };
//    for (int i=0; i<3; ++i) {
//        for (int j=0; j<3; ++j) {
//            sensor2xyz_n[i][j] = sensor2xyz[i*3+j] / d65_white[i];
//            sensor2xyz_n_maxvalue = std::max(sensor2xyz_n_maxvalue, sensor2xyz_n[i][j]);
//        }
//    }
//
//    /* Scale factor that is guaranteed to push XYZ values into the range [0, 1] */
//    float scale = 1.0 / (maxvalue * sensor2xyz_n_maxvalue);
//
//    /* Precompute a table for the nonlinear part of the CIELab conversion */
//    const int cielab_table_size = 0xFFFF;
//    float cielab_table[cielab_table_size];
//    for (int i=0; i<cielab_table_size; ++i) {
//        float r = i * 1.0f / (cielab_table_size-1);
//        cielab_table[i] = r > 0.008856 ? std::pow(r, 1.0f / 3.0f) : 7.787f*r + 4.0f/29.0f;
//    }
//
//    /* Process the image in tiles */
//    std::vector<std::pair<size_t, size_t>> tiles;
//    for (size_t top = 2; top < height - 5; top += tsize - 6)
//        for (size_t left = 2; left < width - 5; left += tsize - 6)
//            tiles.push_back(std::make_pair(left, top));
//
//#pragma omp parallel for /* Parallelize over tiles */
//    for (int tile=0; tile<tiles.size(); ++tile) {
//        DemosaicBuffer &buf = buffers[omp_get_thread_num()];
//        size_t left = tiles[tile].first, top = tiles[tile].second;
//
//        for (size_t y=top; y<top+tsize && y<height-2; ++y) {
//            /* Interpolate green horizontally and vertically, starting
//               at the first position where it is missing */
//            size_t x = left + (fc(left, y) & 1), color = fc(x, y);
//
//            for (; x<left+tsize && x<width-2; x += 2) {
//                float3 *pix = image_demosaiced + y*width + x;
//
//                float interp_h = 0.25f * ((pix[-1][G] + pix[0][color] + pix[1][G]) * 2
//                                          - pix[-2][color] - pix[2][color]);
//                float interp_v = 0.25*((pix[-width][G] + pix[0][color] + pix[width][G]) * 2
//                                       - pix[-2*width][color] - pix[2*width][color]);
//
//                /* Don't allow the interpolation to create new local maxima / minima */
//                buf.rgb[0][y-top][x-left][G] = clamp(interp_h, pix[-1][G], pix[1][G]);
//                buf.rgb[1][y-top][x-left][G] = clamp(interp_v, pix[-width][G], pix[width][G]);
//            }
//        }
//
//        /* Interpolate red and blue, and convert to CIELab */
//        for (int dir=0; dir<2; ++dir) {
//            for (size_t y=top+1; y<top+tsize-1 && y<height-3; ++y) {
//                for (size_t x = left+1; x<left+tsize-1 && x<width-3; ++x) {
//                    float3 *pix = image_demosaiced + y*width + x;
//                    float3 *interp = &buf.rgb[dir][y-top][x-left];
//                    float3 *lab = &buf.cielab[dir][y-top][x-left];
//
//                    /* Determine the color at the current pixel */
//                    int color = fc(x, y);
//
//                    if (color == G) {
//                        color = fc(x, y+1);
//                        /* Interpolate both red and green */
//                        interp[0][2-color] = std::max(0.0f, pix[0][G] + (0.5f*(
//                                pix[-1][2-color] + pix[1][2-color] - interp[-1][G] - interp[1][G])));
//
//                        interp[0][color] = std::max(0.0f,  pix[0][G] + (0.5f*(
//                                pix[-width][color] + pix[width][color] - interp[-tsize][1] - interp[tsize][1])));
//                    } else {
//                        /* Interpolate the other color */
//                        color = 2 - color;
//                        interp[0][color] = std::max(0.0f, interp[0][G] + (0.25f * (
//                                pix[-width-1][color] + pix[-width+1][color]
//                                + pix[+width-1][color] + pix[+width+1][color]
//                                - interp[-tsize-1][G] - interp[-tsize+1][G]
//                                - interp[+tsize-1][G] - interp[+tsize+1][G])));
//                    }
//
//                    /* Forward the color at the current pixel with out modification */
//                    color = fc(x, y);
//                    interp[0][color] = pix[0][color];
//
//                    /* Convert to CIElab */
//                    float xyz[3] = { 0, 0, 0 };
//                    for (int i=0; i<3; ++i)
//                        for (int j=0; j<3; ++j)
//                            xyz[i] += sensor2xyz_n[i][j] * interp[0][j];
//
//                    for (int i=0; i<3; ++i)
//                        xyz[i] = cielab_table[std::max(0, std::min(cielab_table_size-1,
//                                                                   (int) (xyz[i] * scale * cielab_table_size)))];
//
//                    lab[0][0] = (116.0f * xyz[1] - 16);
//                    lab[0][1] = 500.0f * (xyz[0] - xyz[1]);
//                    lab[0][2] = 200.0f * (xyz[1] - xyz[2]);
//                }
//            }
//        }
//
//        /*  Build homogeneity maps from the CIELab images: */
//        const int offset_table[4] = { -1, 1, -tsize, tsize };
//        memset(buf.homo, 0, 2*tsize*tsize);
//        for (size_t y=top+2; y < top+tsize-2 && y < height-4; ++y) {
//            for (size_t x=left+2; x< left+tsize-2 && x < width-4; ++x) {
//                float ldiff[2][4], abdiff[2][4];
//
//                for (int dir=0; dir < 2; dir++) {
//                    float3 *lab = &buf.cielab[dir][y-top][x-left];
//
//                    for (int i=0; i < 4; i++) {
//                        int offset = offset_table[i];
//
//                        /* Luminance and chromaticity differences in 4 directions,
//                           for each of the two interpolated images */
//                        ldiff[dir][i] = std::abs(lab[0][0] - lab[offset][0]);
//                        abdiff[dir][i] = square(lab[0][1] - lab[offset][1])
//                                         + square(lab[0][2] - lab[offset][2]);
//                    }
//                }
//
//                float leps  = std::min(std::max(ldiff[0][0], ldiff[0][1]),
//                                       std::max(ldiff[1][2], ldiff[1][3]));
//                float abeps = std::min(std::max(abdiff[0][0], abdiff[0][1]),
//                                       std::max(abdiff[1][2], abdiff[1][3]));
//
//                /* Count the number directions in which the above thresholds can
//                   be maintained, for each of the two interpolated images */
//                for (int dir=0; dir < 2; dir++)
//                    for (int i=0; i < 4; i++)
//                        if (ldiff[dir][i] <= leps && abdiff[dir][i] <= abeps)
//                            buf.homo[dir][y-top][x-left]++;
//            }
//        }
//
//        /*  Combine the most homogenous pixels for the final result */
//        for (size_t y=top+3; y < top+tsize-3 && y < height-5; ++y) {
//            for (size_t x=left+3; x < left+tsize-3 && x < width-5; ++x) {
//                /* Look, which of the to images is more homogeneous in a 3x3 neighborhood */
//                int hm[2] = {0, 0};
//                for (int dir=0; dir < 2; dir++)
//                    for (size_t i=y-top-1; i <= y-top+1; i++)
//                        for (size_t j=x-left-1; j <= x-left+1; j++)
//                            hm[dir] += buf.homo[dir][i][j];
//
//                if (hm[0] != hm[1]) {
//                    /* One of the images was more homogeneous */
//                    for (int col=0; col<3; ++col)
//                        image_demosaiced[y*width+x][col] = buf.rgb[hm[1] > hm[0] ? 1 : 0][y-top][x-left][col];
//                } else {
//                    /* No clear winner, blend */
//                    for (int col=0; col<3; ++col)
//                        image_demosaiced[y*width+x][col] = 0.5f*(buf.rgb[0][y-top][x-left][col]
//                                                                 + buf.rgb[1][y-top][x-left][col]);
//                }
//            }
//        }
//    }
//
//    delete[] buffers;
//    delete[] image_merged;
//    image_merged = NULL;
//}

void demosaic( int argc, char* argv[] ) {
    pfs::DOMIO pfsio;


    bool verbose = false;
    int c;


    static struct option cmdLineOptions[] = {
    { "help", no_argument, NULL, 'h' },
    { "verbose", no_argument, NULL, 'v' },
    { NULL, 0, NULL, 0 }
    };

    int optionIndex = 0;
    while ((c = getopt_long(argc, argv, "hv", cmdLineOptions, &optionIndex)) != -1) {

        switch( c ) {
            case 'h':
              printHelp();
              throw QuietException();
            case 'v':
              verbose = true;
              break;
            default:
              throw QuietException();
        }
    }

  
  while( true ) {
    pfs::Frame *frame = pfsio.readFrame( stdin );
    if( frame == NULL ) break; // No more frames


    pfs::Channel *X, *Y, *Z;
    frame->getXYZChannels( X, Y, Z );
    pfs::transformColorSpace( pfs::CS_XYZ, X, Y, Z, pfs::CS_RGB, X, Y, Z );



    RGBEWriter writer( stdout, true );
    writer.writeImage( X, Y, Z, "" );

    pfsio.freeFrame( frame );

  }
}

int main( int argc, char* argv[] )
{
  try {
      demosaic( argc, argv );
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
