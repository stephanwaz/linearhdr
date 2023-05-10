/**
 * simplified hdr generation assuming raw linear response, fixes possible bugs in old code that lost absolute 
 * calibration
 * @author Stephen Wasilewski stephen.wasilewski@epfl.ch 
 *
 * Forked from: version pfstools 2.2.0:
 *
 * @brief Create an HDR image or calibrate a response curve from a set
 * of differently exposed images supplied in PFS stream
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
 * @author Grzegorz Krawczyk, <krawczyk@mpi-sb.mpg.de>
 * @author Ivo Ihrke, <ihrke@mmci.uni-saarland.de>
 *
 * $Id: pfshdrcalibrate.cpp,v 1.16 2011/02/24 17:35:59 ihrke Exp $
 */

#include <config.h>

#include <iostream>
#include <vector>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <getopt.h>

#include <pfs.h>

#include <responses.h>
#include <linearhdr.h>
#include <fstream>
#include <sstream>
#include <rgbeio.h>

#ifdef NETPBM_FOUND
#include <ppmio.h>
#endif

using namespace std;

#define PROG_NAME "linearhdr"

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
    return min(max(a, b), c);
}

inline float max(float a, float b, float c, float d) {
    return max(max(max(a, b), c), d);
}

struct FrameInfo {float etime; float iso; float aperture; float factor;};

inline FrameInfo frame_info(pfs::Frame *frame, float scale){
    FrameInfo result = {1.f, 100.f, 1.f, 1.f};

    const char *etime_str =
            frame->getTags()->getString("exposure_time");

    result.etime = atof(etime_str);

    const char *iso_str =
            frame->getTags()->getString("ISO");
    if (iso_str != nullptr)
        result.iso = atof(iso_str);

    const char *aperture_str = frame->getTags()->getString("aperture");
    if (aperture_str != nullptr)
        result.aperture = atof(aperture_str);

    result.factor = scale * 100.0f * result.aperture * result.aperture / ( result.iso * result.etime );
    return result;

}

FrameInfo correct_exposure(FrameInfo info) {
    FrameInfo result = {1.f, 100.f, 1.f, 1.f};
    result.etime = std::pow(2, -std::round(std::log2(1/info.etime)*3)/3);
    result.aperture = std::pow(2, std::round(std::log2(info.aperture* info.aperture)*3)/6);
    result.iso = info.iso;
    return result;
}

//---------------------------------------------------
//--- standard PFS stuff
bool verbose = false;

class QuietException {
};


void printHelp() {
    fprintf(stderr, PROG_NAME " [Options] [exposure_list]\n"
                    "Options:\n"
                    "\t[--saturation-offset, -o <val>]: exclude images within <val> of 1 default=0.2\n"
                    "\t[--range, -r <val>]: dynamic range of single raw exposure, used to set lower cutoff,"
                    "\n\t\tgive as power of 2 default=6.64386\n"
                    "\t[--deghosting, -d <val>]: relative difference for outlier detection when less than 1,"
                    "\n\t\totherwise absolute difference (good for clouds) default=OFF\n"
                    "\t[--tsv, -t]: output raw data as tsv, exposures seperated by extra linebreak,"
                    "\n\t\tdo not use with large files!\n"
                    "\t[--scale, -s <val>]: absolute scaling for hdr (eg ND filter, known response, etc.) default=1.0\n"
                    "\t\tuse linearhdr_calibrate to calculate\n"
                    "\t[--use-yxy, -X]: merge hdr in Yxy space instead of RGB\n"
                    "\t[--cull, -c]: throw away extra exposures that are not needed to keep output in range\n"
                    "\t[--rgbe, -R]: output radiance rgbe (default)\n"
                    "\t[--pfs, -P]: output pfs stream\n"
                    "\t[--exact, -e]: input camera values interpreted as exact (default=False)\n"
                    "\t[--nominal, -n]: input camera values interpreted as nominal (default=True)\n"
                    "\t[--verbose, -v]\n\t[--help]\n\n"
                    "If exposure_list is given, images are read from file formatted as:\n"
                    "\t<image1.ppm> <iso> <aperture> <exposure_time>\n"
                    "\t<image2.ppm> <iso> <aperture> <exposure_time>\n\t...\n\n"
                    "list should be sorted by longest exposure time to shortest (only critical if --cull)\n"
                    "else, program expects a sequence of images (formatted as from pfsin on the stdin),\n"
                    "use/see 'linearhdr_make_list' for an example.\n"
                    "By default, linearhdr expects nominal aperture and shutter speed.\n"
                    "If using pfsinme, note that nominal camera values are manipulated by dcraw\n"
                    " (but with less accuracy) so make sure to use the --exact flag so shutter \n"
                    "and aperture are not double corrected.\n\n");
}

void pfshdrraw(int argc, char *argv[]) {


    pfs::DOMIO pfsio;

    /* defaults */

    float opt_saturation_offset_perc = 0.01;
    float opt_black_offset_perc = 0.01;
    float range = 6.64386;
    float opt_deghosting = -1;
    float opt_scale = 1.0f;
    bool yxy = false;
    bool keep = true;
    int lead_channel = -1;
    bool dataonly = false;
    bool rgbe = true;
    bool nominal = true;

    /* helper */
    int c;

    static struct option cmdLineOptions[] = {
            {"help",       no_argument,       nullptr, 'h'},
            {"verbose",    no_argument,       nullptr, 'v'},
            {"use-yxy",    no_argument,       nullptr, 'X'},
            {"cull",    no_argument,       nullptr, 'c'},
            {"rgbe",    no_argument,       nullptr, 'R'},
            {"pfs",    no_argument,       nullptr, 'P'},
            {"exact",    no_argument,       nullptr, 'e'},
            {"nominal",    no_argument,       nullptr, 'n'},
            {"deghosting", required_argument, nullptr, 'd'},
            {"tsv", no_argument, nullptr, 't'},
            { "saturation-offset", required_argument, nullptr, 'o' },
            { "range", required_argument, nullptr, 'r' },
            { "scale", required_argument, nullptr, 's' },
            {nullptr, 0,                         nullptr, 0}
    };

    int optionIndex = 0;
    while ((c = getopt_long(argc, argv, "hnevud:s:r:X:o:", cmdLineOptions, &optionIndex)) != -1) {

        switch (c) {
            /* help */
            case 'h':
                printHelp();
                throw QuietException();
            /* verbose */
            case 'v':
                verbose = true;
                break;
            case 'n':
                nominal = true;
                break;
            case 'e':
                nominal = false;
                break;
            case 'X':
                yxy = true;
                lead_channel = 0;
                break;
            case 'o':
                opt_saturation_offset_perc = atof(optarg);
                if( opt_saturation_offset_perc < 0 || opt_saturation_offset_perc > 0.25 )
                    throw pfs::Exception("saturation offset should be between 0 and 0.25");
                break;
            case 'r':
                range = atof(optarg);
                if( range < 2)
                    throw pfs::Exception("range should be greater than 2");
                opt_black_offset_perc = pow(2, -range);
                break;
            case 's':
                opt_scale = atof(optarg);
                if( opt_scale <= 0)
                    throw pfs::Exception("scale must be positive");
                break;
            case 'd':
                opt_deghosting = atof(optarg);
                if (opt_deghosting < 0)
                    throw pfs::Exception("deghosting threshold should be >0");
                break;
            case 't':
                dataonly = true;
                break;
            case 'c':
                keep = false;
                break;
            case 'R':
                rgbe = true;
                break;
            case 'P':
                rgbe = false;
                break;
            default:
                throw QuietException();
        }
    }
    bool fromfile = argv[optind] != nullptr;

    std::ifstream infile(argv[optind]);


    int frame_no = 0;
    int width = 0;
    int height = 0;
    int size = 0;
    float fmin, fmax;
    float gmin = 1e30;
    float gmax = 1e-30;
    bool oorange = true;
    float pmax;

    // collected exposures
    ExposureList imgsR;
    ExposureList imgsG;
    ExposureList imgsB;

    pfs::ColorSpace hdr_target = yxy ? pfs::CS_Yxy : pfs::CS_RGB;

    while (true) {
        pmax = 0;
        pfs::Frame *iframe = nullptr;
        FrameInfo info = {1.f, 100.f, 1.f, 1.f};

        if (!fromfile){
            //--- read frames from pfs stream
            iframe = pfsio.readFrame(stdin);
            if (iframe == nullptr)
                break; // No more frames
        } else {
            //--- read frames from list by parsing file
            std::string line, framefile;
            if (std::getline(infile, line)) {
                std::istringstream iss(line);
                if (!(iss >> framefile >> info.iso >> info.aperture >> info.etime)) { continue;}
                if (nominal)
                    info = correct_exposure(info);
                info.factor = opt_scale * 100.0f * info.aperture * info.aperture / ( info.iso * info.etime );

                FILE *fh = fopen( framefile.c_str(), "rb");
                pfs::FrameFile ff = pfs::FrameFile( fh, framefile.c_str());
                #ifdef NETPBM_FOUND
                PPMReader reader( PROG_NAME, ff.fh );
                iframe = pfsio.createFrame( reader.getWidth(), reader.getHeight() );
                pfs::Channel *X, *Y, *Z;
                iframe->createXYZChannels( X, Y, Z );

                //Store sRGB data temporarily in XYZ channels
                reader.readImage( X, Y, Z );
                pfs::transformColorSpace(pfs::CS_RGB, X, Y, Z, hdr_target, X, Y, Z);
                #else
                throw pfs::Exception("linearhdr compiled without ppm support use pfsin to load files");
                #endif

            } else
                break; //no more lines
        }

        if (keep || oorange) {

            pfs::Channel *X = nullptr;
            pfs::Channel *Y = nullptr;
            pfs::Channel *Z = nullptr;

            iframe->getXYZChannels(X, Y, Z);

            if (!fromfile){
                info = frame_info(iframe, opt_scale);
                pfs::transformColorSpace(pfs::CS_XYZ, X, Y, Z, hdr_target, X, Y, Z);
            }


            if (X == nullptr || Y == nullptr || Z == nullptr)
                throw pfs::Exception("missing XYZ channels in the PFS stream (try to preview your files using pfsview)");

            // frame size
            width = Y->getCols();
            height = Y->getRows();
            size = width * height;

            if (dataonly){
                fmax = info.factor;
                for (int i = 0; i < size; i++) {
                    std::cout << (*X)(i) << "\t" << (*Y)(i) << "\t" << (*Z)(i) << "\t";
                    if (yxy) {
                        pmax = max((*X)(i), pmax);
                        float below = (*X)(i) < opt_black_offset_perc;
                        float above = (*X)(i) > 1 - opt_saturation_offset_perc;
                        std::cout << (*X)(i) * fmax << "\t" << (*Y)(i) << "\t" << (*Z)(i) << "\t" << (*X)(i) * fmax << "\t" << below << "\t" << above << std::endl;
                    } else {
                        pmax = max((*X)(i), (*Y)(i), (*Z)(i), pmax);
                        float below = min((*X)(i), (*Y)(i), (*Z)(i)) < opt_black_offset_perc;
                        float above = max((*X)(i), (*Y)(i), (*Z)(i)) > 1 - opt_saturation_offset_perc;
                        float lum = (0.265074126*(*X)(i) + 0.670114631*(*Y)(i) + 0.064811243*(*Z)(i)) * fmax;
                        std::cout << (*X)(i) * fmax << "\t" << (*Y)(i) * fmax << "\t" << (*Z)(i) * fmax << "\t" << lum << "\t" << below << "\t" << above << std::endl;
                    }
                }
                std::cout  << std::endl;
            } else {
                Exposure eR, eG, eB;
                eR.iso = eG.iso = eB.iso = info.iso;
                eR.aperture = eG.aperture = eB.aperture = info.aperture;
                eR.exposure_time = eG.exposure_time = eB.exposure_time = info.etime;
                eR.yi = new pfs::Array2DImpl(width, height);
                eG.yi = new pfs::Array2DImpl(width, height);
                eB.yi = new pfs::Array2DImpl(width, height);

                if (eR.yi == nullptr || eG.yi == nullptr || eB.yi == nullptr)
                    throw pfs::Exception("could not allocate memory for source exposure");

                for (int i = 0; i < size; i++) {

                    (*eR.yi)(i) = (*X)(i);
                    (*eG.yi)(i) = (*Y)(i);
                    (*eB.yi)(i) = (*Z)(i);
                    if (yxy)
                        pmax = max((*X)(i), pmax);
                    else
                        pmax = max((*X)(i), (*Y)(i), (*Z)(i), pmax);
                }

                // add to exposures list
                imgsR.push_back(eR);
                imgsG.push_back(eG);
                imgsB.push_back(eB);
            }
            fmax = info.factor * (1 - opt_saturation_offset_perc);
            fmin = fmax * opt_black_offset_perc;
            gmax = max(gmax, fmax);
            gmin = min(gmin, fmin);

            frame_no++;
            VERBOSE_STR << "frame #" << frame_no << ", min:" << fmin << ", max:" << fmax <<  endl;
        }

        oorange = oorange && (pmax > 1 - opt_saturation_offset_perc);
        pfsio.freeFrame(iframe);
    }
    if (dataonly)
        return;
    if (frame_no < 1)
        throw pfs::Exception("at least one image required for calibration (check paths in hdrgen script?)");

    VERBOSE_STR << "using " << frame_no  << " frames, range min:" << gmin << ", max:" << gmax <<  endl;

    if (oorange) {
        VERBOSE_STR << "Warning: some pixels out of range..."  <<  endl;
    }


    // create channels for output
    pfs::Frame *frame = pfsio.createFrame(width, height);

    pfs::Channel *Xj = nullptr;
    pfs::Channel *Yj = nullptr;
    pfs::Channel *Zj = nullptr;

    frame->createXYZChannels(Xj, Yj, Zj);

    /* counter for saturated pixels */
    long sp;

    pfs::Array2D *RGB_out[3] = {Xj, Yj, Zj};
    const ExposureList *exposures[] = {&imgsR, &imgsG, &imgsB};

    VERBOSE_STR << "applying response..." << endl;
    sp = linear_Response(RGB_out, exposures, opt_saturation_offset_perc, opt_black_offset_perc,
                         opt_deghosting, opt_scale, lead_channel);

    if (sp > 0) {
        float perc = ceilf(100.0f * sp / size);
        VERBOSE_STR << "saturated pixels found in " << perc << "% of the image!" << endl;
    }
    if (rgbe){
        RGBEWriter writer( stdout, true );
        pfs::transformColorSpace( hdr_target, Xj, Yj, Zj, pfs::CS_RGB, Xj, Yj, Zj );
        writer.writeImage( Xj, Yj, Zj );
    } else {
        pfs::transformColorSpace(hdr_target, Xj, Yj, Zj, pfs::CS_XYZ, Xj, Yj, Zj);
        pfsio.writeFrame(frame, stdout);
    }

    // clean up memory
    pfsio.freeFrame(frame);

    for (int i = 0; i < imgsR.size(); i++) {
        delete imgsR[i].yi;
        delete imgsG[i].yi;
        delete imgsB[i].yi;
    }
}


int main(int argc, char *argv[]) {
    try {
        pfshdrraw(argc, argv);
    }

    catch (pfs::Exception ex) {
        fprintf(stderr, PROG_NAME " error: %s\n", ex.getMessage());
        return EXIT_FAILURE;
    }

    catch (QuietException ex) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
