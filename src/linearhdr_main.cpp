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
 * This file is derived from a part of PFS CALIBRATION package.
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
#include <iomanip>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <getopt.h>

#include <pfs.h>

#include <linearhdr.h>
#include <fstream>
#include <sstream>
#include <rgbeio.h>
#include <hdrtiffio.h>

using namespace std;

#define PROG_NAME "linearhdr"
#define PROG_VERSION "0.2.1 (compiled on: " __DATE__ " @ " __TIME__ ")"

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
                    "Version: " PROG_VERSION "\n"
                    "Options:\n"
                    "\t[--verbose, -v]\n\t[--help]\n\n"
                    "\tInput Parameters:\n"
                    "\t[--exact, -E]: input camera values interpreted as exact (default=True)\n"
                    "\t[--nominal, -N]: input camera values interpreted as nominal (default=False)\n"
                    "\t[--crop, -c <x y w h> use cropbox x and y given from upper left\n"
                    "\t[--bayer, -B]: toggle  mosaic input (rawconvert -disinterp), default is true\n\n"
                    "\tCorrections:\n"
                    "\t[--oor-low, -m <val>]: value to use for out of range low, default from data\n"
                    "\t[--oor-high, -x <val>]: value to use for out of range high, default from data\n"
                    "\t[--rgbs, -k '<val> <val> <val>']: rgb channel calibration default=1.0 1.0 1.0\n"
                    "\t\toverriden by RGBcalibration in header line. applies to output colorspace (after Camera2RGB in header line)\n"
                    "\t[--scale, -s <val>]: absolute scaling for hdr (eg ND filter, known response, etc.) default=1.0\n"
                    "\t[--efc, -C '<val> <val> <val> <val>']: electronic front curtain shutter correction, give as m a b c\n"
                    "\t\tto correct on image height according to the function: y = 1/((x/c+a*t)/(1+a*t))^b where t is\n"
                    "\t\tthe (mechanical shutter) corrected exposure time. 'm' is the maximum exposure time to which these coefficients apply.\n"
                    "\t\tgive multiple times, starting with the shortest maximum time, to apply up to three ranges of efcs.\n"
                    "\t\tNote that if three sets are given, the last coeffs will apply to all exposure times regardless of 'm'\n"
                    "\t\tt*a*x**2 + t*b*x + 1+t*c\n\n"
                    "\tMerge Settings:\n"
                    "\t[--me, -M]: merge each, weight each color channel independently pixel wise. fastest performance, but most artifacts\n"
                    "\t[--ub, -U]: use best, pick exposure with best weight instead of weighted average\n"
                    "\t[--usemax, -Z]: use highest recorded value (only with -B and not -M) for checking consistency\n"
                    "\t[--usemin, -Y]: use lowest recorded value (only with -B and not -M) (ignored if -Z) for checking consistency\n"
                    "\t[--onepix, -O]: only use a one row buffer for saturation test, instead of the default 2 (only with -B and not -M)\n"
                    "\t\tfaster performance but more artifacts\n"
                    "\t[--saturation-offset, -o <val>]: exclude images within <val> of 1 default=0.01\n"
                    "\t[--range, -r <val>]: lower range of single raw exposure, used to set lower cutoff default=0.0,"
                    "\n\t\tgive as value between 0 and 0.25, default=0.01\n"
                    "\t[--ignore, -G]: ignore all camera data, only use with single frame or all with same exposure\n\n"
                    "\tPostprocessing:\n"
                    "\t[--debayer, -D]: interpolate hdr output, overrides --bayer, but expects same input (rawconvert -disinterp)\n"
                    "\t[--filter, -F]: toggle rank filter (only applies with --debayer), filters to 25th and 75th percentile of surrounding pixels. Default=true\n\n"
                    "\tOutput:\n"
                    "\t[--tsv, -T]: output data as tsv (with header for pvalue -r)\n"
                    "\t[--rgbe, -R]: output radiance rgbe (default)\n"
                    "\t[--pfs, -P]: output pfs stream\n\n"
                    "images are read from file formatted as:\n"
                    "\t<image1.tiff> <iso> <aperture> <exposure_time>\n"
                    "\t<image2.tiff> <iso> <aperture> <exposure_time>\n\t...\n\n"
                    "use/see 'pylinearhdr makelist' for an example.\n"
                    "By default, linearhdr expects exact aperture and shutter speed. so if you are making this list manually"
                    " using nominal values, be sure to use --nominal to better estimate exposure. Note the is generally"
                    " not very reliable as fast exposure times / small apertures can be dramatically different from the "
                    "default correction.\n\n");
}

void linearhdr_main(int argc, char *argv[]) {

    std::stringstream header;


    pfs::DOMIO pfsio;

    /* defaults */

    float sat_off = 0.01;
    float blk_off = 0.00;
    float scale = 1.0f;
    bool tsv = false;
    bool rgbe = true;
    bool nominal = false;
    bool isbayer = true;
    bool demosaic = false;
    bool median = true;
    bool ignore = false;
    bool mergeeach = false;
    bool usebest = false;
    bool usemax = false;
    bool usemin = false;
    bool onep =false;
    float oor_high = -1;
    float oor_low = -1;
    int crop[4] = {0, 0, 0, 0};
    float rgb_corr[3][3] = {{1.0, 0.0, 0.0},
                            {0.0, 1.0, 0.0},
                            {0.0, 0.0, 1.0}};
    float vlambda[3] = {0.333333, 0.333334, 0.333333};
    float wsp[3] = {1, 1, 1};
    float rgbcal[3] = {1.0, 1.0, 1.0};
    float efc[3][4] = {{0.0, 0.0, 0.0, 0.0},
                       {0.0, 0.0, 0.0, 0.0},
                       {0.0, 0.0, 0.0, 0.0}};
    int efci = 0;
    /* helper */
    int c;

    static struct option cmdLineOptions[] = {
            {"help",       no_argument,       nullptr, 'h'},
            {"verbose",    no_argument,       nullptr, 'v'},
            {"rgbe",    no_argument,       nullptr, 'R'},
            {"pfs",    no_argument,       nullptr, 'P'},
            {"exact",    no_argument,       nullptr, 'E'},
            {"bayer",    no_argument,       nullptr, 'B'},
            {"crop",    required_argument,       nullptr, 'c'},
            {"debayer",    no_argument,       nullptr, 'D'},
            {"nominal",    no_argument,       nullptr, 'N'},
            {"tsv", no_argument, nullptr, 'T'},
            { "saturation-offset", required_argument, nullptr, 'o' },
            { "range", required_argument, nullptr, 'r' },
            { "scale", required_argument, nullptr, 's' },
            { "rgbs", required_argument, nullptr, 'k' },
            { "efc", required_argument, nullptr, 'e' },
            { "oor-low", required_argument, nullptr, 'm' },
            { "oor-high", required_argument, nullptr, 'x' },
            { "oob-low", required_argument, nullptr, 'm' },
            { "oob-high", required_argument, nullptr, 'x' },
            { "ignore", no_argument, nullptr, 'G' },
            { "filter", no_argument, nullptr, 'F' },
            { "me", no_argument, nullptr, 'M' },
            { "ub", no_argument, nullptr, 'U' },
            { "usemax", no_argument, nullptr, 'Z' },
            { "usemin", no_argument, nullptr, 'Y' },
            { "onepix", no_argument, nullptr, 'O' },
            {nullptr, 0,                         nullptr, 0}
    };

    std::stringstream k; //to read in multivalue arguments
    int optionIndex = 0;
    int ci = 0;
    while ((c = getopt_long(argc, argv, "hvBDEFGMNOPRTUYZc:e:k:m:o:r:s:x:", cmdLineOptions, &optionIndex)) != -1) {
        ci++;
        if (strlen(argv[ci]) > 2 && argv[ci][1] != '-'){
            char message[100];
            snprintf(message, 100, "bad option : %s, all long options need --", argv[ci]);
            throw pfs::Exception(message);
        }
        switch (c) {
            /* help */
            case 'h':
                printHelp();
                throw QuietException();
            /* verbose */
            case 'v':
                verbose = true;
                break;
            case 'B':
                isbayer = not isbayer;
                break;
            case 'D':
                demosaic = true;
                break;
            case 'E':
                nominal = false;
                break;
            case 'F':
                median = not median;
                break;
            case 'G':
                ignore = true;
                break;
            case 'M':
                mergeeach = true;
                break;
            case 'N':
                nominal = true;
                break;
            case 'O':
                onep = true;
                break;
            case 'P':
                rgbe = false;
                break;
            case 'R':
                rgbe = true;
                break;
            case 'T':
                tsv = true;
                break;
            case 'U':
                usebest = true;
                break;
            case 'Y':
                usemin = true;
                break;
            case 'Z':
                usemax = true;
                break;
            case 'c':
                ci++;
                k.str("");
                k.clear();
                k << optarg;
                k >> crop[0] >> crop[1]  >> crop[2] >> crop[3];
                break;
            case 'e':
                ci++;
                k.str("");
                k.clear();
                k << optarg;
                k >> efc[efci][0] >> efc[efci][1] >> efc[efci][2] >> efc[efci][3];
                efci++;
                break;
            case 'k':
                ci++;
                k.str("");
                k.clear();
                k << optarg;
                k >> rgbcal[0] >> rgbcal[1] >> rgbcal[2];
                break;
            case 'm':
                ci++;
                oor_low = atof(optarg);
                break;
            case 'o':
                ci++;
                sat_off = atof(optarg);
//                if( opt_saturation_offset_perc < 0 || opt_saturation_offset_perc > 0.25 )
//                    throw pfs::Exception("saturation offset should be between 0 and 0.25");
                break;
            case 'r':
                ci++;
                blk_off = atof(optarg);
//                if( opt_black_offset_perc < 0 || opt_black_offset_perc > 0.25 )
//                    throw pfs::Exception("saturation offset should be between 0 and 0.25");
                break;
            case 's':
                ci++;
                scale = atof(optarg);
                if( scale <= 0)
                    throw pfs::Exception("scale must be positive");
                break;
            case 'x':
                ci++;
                oor_high = atof(optarg);
                break;
            default:
                throw QuietException();
        }
    }

    if (argv[optind] == nullptr){
        printHelp();
        throw QuietException();
    }

    string cprefix = tsv ? "" : "";
    header << cprefix << PROG_NAME << "_VERSION= " << PROG_VERSION << endl;
    header << cprefix << argv[0];
    for (int i = 1; i < argc-1; i++){
        header << " " << argv[i];
    }
    header << endl;

    std::ifstream infile(argv[optind]);
    isbayer = isbayer && not demosaic;
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

    // read through specification file loading header information and image frames
    while (true) {
        pmax = 0;
        pfs::Frame *iframe = nullptr;
        pfs::Frame *outframe = nullptr;
        FrameInfo info = {1.f, 100.f, 1.f, 1.f};

        //--- read frames from list by parsing file
        std::string line, framefile;
        if (std::getline(infile, line)) {
            std::istringstream iss(line);
            if (!(iss >> framefile)){
                continue;
            }
            if (!(iss >> info.iso >> info.aperture >> info.etime)) {
                if (framefile.at(0) == '#') {
                    std::string comment = iss.str();
                    const uint begin =  comment.find_first_not_of("# \t");
                    const uint equal = comment.find_first_of('=');
                    if (!isbayer && comment.substr(begin, equal - begin) == "Camera2RGB"){
                        istringstream ss(comment.substr(equal+1, comment.size()));
                        ss >> rgb_corr[0][0] >> rgb_corr[0][1] >> rgb_corr[0][2];
                        ss >> rgb_corr[1][0] >> rgb_corr[1][1] >> rgb_corr[1][2];
                        ss >> rgb_corr[2][0] >> rgb_corr[2][1] >> rgb_corr[2][2];
                    }
                    if (demosaic && comment.substr(begin, equal - begin) == "RGBcalibration"){
                        istringstream ss(comment.substr(equal+1, comment.size()));
                        ss >> rgbcal[0] >> rgbcal[1] >> rgbcal[2];
                    }
                    if ( comment.substr(begin, equal - begin) == "LuminanceRGB"){
                        istringstream ss(comment.substr(equal+1, comment.size()));
                        ss >> vlambda[0] >> vlambda[1] >> vlambda[2];
                    }
                    if ( comment.substr(begin, equal - begin) == "WhiteSaturation"){
                        istringstream ss(comment.substr(equal+1, comment.size()));
                        ss >> wsp[0] >> wsp[1] >> wsp[2];
                    }
                    header  << cprefix << comment.substr(begin, comment.size()) << endl;
                }
                continue;
            }
            if (nominal)
                info = correct_exposure(info);
            if (ignore) {
                info.etime = 1;
                info.iso = 100;
                info.aperture = 1;
            }
            info.factor = scale * 100.0f * info.aperture * info.aperture / ( info.iso * info.etime );

            FILE *fh = fopen( framefile.c_str(), "rb");
            pfs::FrameFile ff = pfs::FrameFile( fh, framefile.c_str());
            HDRTiffReader reader( ff.fileName);
            pfs::Channel *X, *Y, *Z;
            if (crop[2] == 0 || crop[3] == 0){
                iframe = pfsio.createFrame( reader.getWidth(), reader.getHeight() );
                iframe->createXYZChannels( X, Y, Z );
                reader.readImage( X, Y, Z );
            } else {
                outframe = pfsio.createFrame( reader.getWidth(), reader.getHeight() );
                pfs::Channel *Xi, *Yi, *Zi;
                outframe->createXYZChannels( Xi, Yi, Zi );
                reader.readImage( Xi, Yi, Zi );
                iframe = pfsio.createFrame(crop[2], crop[3]);
                iframe->createXYZChannels( X, Y, Z );
                for (int i=0; i < crop[3]; i++)
                    for (int j=0; j< crop[2]; j++) {
                        (*X)(j, i) = (*Xi)(j + crop[0],i+crop[1]);
                        (*Y)(j, i) = (*Yi)(j + crop[0],i+crop[1]);
                        (*Z)(j, i) = (*Zi)(j + crop[0],i+crop[1]);
                    }
            }



        } else
            break; //no more lines

        pfs::Channel *X = nullptr;
        pfs::Channel *Y = nullptr;
        pfs::Channel *Z = nullptr;

        iframe->getXYZChannels(X, Y, Z);

        if (X == nullptr || Y == nullptr || Z == nullptr)
            throw pfs::Exception("missing XYZ channels in the PFS stream (try to preview your files using pfsview)");

        // frame size
        width = Y->getCols();
        height = Y->getRows();
        size = width * height;

        Exposure eR, eG, eB;
        eR.iso = info.iso;
        eG.iso = info.iso;
        eB.iso = info.iso;
        eR.aperture = eG.aperture = eB.aperture = info.aperture;
        eR.exposure_time = eG.exposure_time = eB.exposure_time = info.etime;
        eR.yi = new pfs::Array2DImpl(width, height);
        eG.yi = new pfs::Array2DImpl(width, height);
        eB.yi = new pfs::Array2DImpl(width, height);

        if (eR.yi == nullptr || eG.yi == nullptr || eB.yi == nullptr)
            throw pfs::Exception("could not allocate memory for source exposure");

        // choose right efc coefficients
        for (int i = 0; i < 3; i++) {
            efci = i;
            if (info.etime < efc[0][i])
                break;
        }

        int s;
        for (int i = 0; i < height; i++)
            for (int j = 0; j < width; j++) {
                s = j + i * width;
                (*eR.yi)(s) = (*X)(s);
                (*eG.yi)(s) = (*Y)(s);
                (*eB.yi)(s) = (*Z)(s);
                pmax = max((*eR.yi)(s), (*eG.yi)(s), (*eB.yi)(s), pmax);
        }

        // add to exposures list
        imgsR.push_back(eR);
        imgsG.push_back(eG);
        imgsB.push_back(eB);

        float calcfac = 0.0;
        for (int m = 0; m < 3; m++) {
            calcfac += vlambda[m] * rgbcal[m] * (rgb_corr[m][0] + rgb_corr[m][1] + rgb_corr[m][2]);
        }

        fmax = info.factor * (1 - sat_off) * calcfac;
        fmin = fmax * max(blk_off, 0.005);
        gmax = max(gmax, fmax);
        gmin = min(gmin, fmin);

        frame_no++;
        VERBOSE_STR << "frame #" << frame_no << ", min:" << fmin << ", max:" << fmax <<  endl;

        oorange = oorange && (pmax > 1 - sat_off);
        pfsio.freeFrame(iframe);
    }

    if (frame_no < 1)
        throw pfs::Exception("at least one image required for calibration");


    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++) {
            rgb_corr[i][j] = rgb_corr[i][j] * rgbcal[i];
        }
    }

    header  << cprefix << "HDR_SEQUENCE_COUNT= " << frame_no << endl;
    header  << cprefix << "HDR_VALID_RANGE= " << std::setprecision(5) << gmin << "-" << gmax << " cd/m^2" << endl;
    VERBOSE_STR << "using " << frame_no  << std::setprecision(5) << " frames, range min:" << gmin << ", max:" << gmax <<  endl;

    if (oorange) {
        VERBOSE_STR << "Warning: some pixels out of range..."  <<  endl;
    }

    // create channels for output
    pfs::Frame *frame = pfsio.createFrame(width, height);

    pfs::Channel *Xj = nullptr;
    pfs::Channel *Yj = nullptr;
    pfs::Channel *Zj = nullptr;

    frame->createXYZChannels(Xj, Yj, Zj);


    pfs::Array2D *RGB_out[3] = {Xj, Yj, Zj};
    const ExposureList *exposures[] = {&imgsR, &imgsG, &imgsB};

    VERBOSE_STR << "merging hdr..." << endl;
    long saturated_pixels, under_pixels;

    tie(saturated_pixels, under_pixels) = linear_response(RGB_out, exposures, sat_off,
                         blk_off, scale, rgb_corr, wsp, efc,
                         oor_high, oor_low, isbayer, demosaic, mergeeach, usebest, median, usemax, usemin, onep);

    if (under_pixels > 0) {
        header  << cprefix << "UNDEREXPOSED_PIXEL_CNT= " << under_pixels << endl;
        float perc = ceilf(100.0f * under_pixels / size);
        VERBOSE_STR << PROG_NAME << ": " << "under-exposed pixels found in " << perc << "% (" << under_pixels
                    << " pixels) of the image!" << endl;
    }

    if (saturated_pixels > 0) {
        header  << cprefix << "OVEREXPOSED_PIXEL_CNT= " << saturated_pixels << endl;
        float perc = ceilf(100.0f * saturated_pixels / size);
        // this one might be important, so always report saturated pixels regardless of --verbose
        std::cerr << PROG_NAME << ": " << "saturated pixels found in " << perc << "% (" << saturated_pixels
                  << " pixels) of the image!" << endl;
    }

    if (tsv) {
        std::string hstring = header.str().substr(0,-1);
        std::cout << "#?RADIANCE"  << std::endl;
        std::cout << hstring << "NCOMP=3" << std::endl << "FORMAT=ascii" << std::endl << std::endl;
        std::cout << "-Y " << height << "     +X " << width << std::endl;
        int s;
        std::cout.precision(10);
        for (int i = 0; i < height; i++)
            for (int j = 0; j < width; j++) {
                s = j + i * width;
                std::cout << std::scientific << (*Xj)(s) << "\t" << (*Yj)(s) << "\t" << (*Zj)(s) << std::endl;
            }

    } else if (rgbe){
        RGBEWriter writer( stdout, not ignore );
        std::string hstring = header.str().substr(0,-1);
        writer.writeImage( Xj, Yj, Zj, hstring );
    } else {
        if (demosaic)
            pfs::transformColorSpace(pfs::CS_RGB, Xj, Yj, Zj, pfs::CS_XYZ, Xj, Yj, Zj);
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
        linearhdr_main(argc, argv);
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
