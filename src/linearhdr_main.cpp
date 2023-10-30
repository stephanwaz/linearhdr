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
#include <vector>
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
#define PROG_VERSION "0.1.2"

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
                    "\t[--range, -r <val>]: lower range of single raw exposure, used to set lower cutoff,"
                    "\n\t\tgive as value between 0 and 0.25, default=0.01\n"
                    "\t[--deghosting, -d <val>]: relative difference for outlier detection when less than 1,"
                    "\n\t\totherwise absolute difference (good for clouds) default=OFF\n"
                    "\t[--tsv, -t]: output raw data as tsv, exposures seperated by extra linebreak,"
                    "\n\t\tdo not use with large files!\n"
                    "\t[--oor-low, -m <val>]: value to use for out of range low, default from data\n"
                    "\t[--oor-high, -x <val>]: value to use for out of range high, default from data\n"
                    "\t[--rgbs, -k '<val> <val> <val>']: rgb channel calibration default=1.0 1.0 1.0\n"
                    "\t\toverriden by RGBcalibration in header line applies to output colorspace (after Camera2RGB in header line)\n"
                    "\t[--scale, -s <val>]: absolute scaling for hdr (eg ND filter, known response, etc.) default=1.0\n"
                    "\t\tuse linearhdr_calibrate to calculate\n"
                    "\t[--efc, -C '<val> <val> <val> <val>']: electronic front curtain shutter correction, give as m a b c\n"
                    "\t\t to correct on image height according to the function: y = 1/((x/c+a*t)/(1+a*t))^b where t is\n"
                    "\t\t the (mechanical shutter) corrected exposure time. 'm' is the maximum exposure time to which these coefficients apply.\n"
                    "\tgive multiple times, starting with the shortest maximum time, to apply up to three ranges of efsc.\n"
                    "\t Note that if three sets are given, the last coeffs will apply to all exposure times regardless of 'm'\n"
                    "\t[--rgbe, -R]: output radiance rgbe (default)\n"
                    "\t[--bayer, -B]: expect mosaic input (rawconvert --disinterp) ignores color correction in exposure_list header\n"
                    "\t[--debayer, -D]: interpolate hdr output, overrides --bayer, but expects same input (rawconvert --disinterp)\n"
                    "\t[--pfs, -P]: output pfs stream\n"
                    "\t[--exact, -e]: input camera values interpreted as exact (default=True)\n"
                    "\t[--nominal, -n]: input camera values interpreted as nominal (default=False)\n"
                    "\t[--verbose, -v]\n\t[--help]\n\n"
                    "images are read from file formatted as:\n"
                    "\t<image1.tiff> <iso> <aperture> <exposure_time>\n"
                    "\t<image2.tiff> <iso> <aperture> <exposure_time>\n\t...\n\n"
                    "list should be sorted by longest exposure time to shortest (only critical if --cull)\n"
                    "use/see 'pylinearhdr makelist' for an example.\n"
                    "By default, linearhdr expects exact aperture and shutter speed. so if you are making this list manually"
                    " using nominal values, be sure to use --nominal to better estimate exposure. Note the is generally"
                    " not very reliable as fast exposure times / small apertures can be dramatically different from the "
                    "default correction.\n\n");
}

void linearhdr_main(int argc, char *argv[]) {

    std::stringstream header;
    header << PROG_NAME << "_VERSION= " << PROG_VERSION << endl;
    header << argv[0];
    for (int i = 1; i < argc-1; i++){
        header << " " << argv[i];
    }
    header << endl;


    pfs::DOMIO pfsio;

    /* defaults */

    float opt_saturation_offset_perc = 0.2;
    float opt_black_offset_perc = 0.01;
    float opt_deghosting = -1;
    float opt_scale = 1.0f;
    bool dataonly = false;
    bool rgbe = true;
    bool nominal = false;
    bool isbayer = false;
    bool demosaic = false;
    float oor_high = -1;
    float oor_low = -1;
    float rgb_corr[3][3] = {{1.0, 0.0, 0.0},
                            {0.0, 1.0, 0.0},
                            {0.0, 0.0, 1.0}};
    float vlambda[3] = {0.333333, 0.333334, 0.333333};
    float rgbcal[3] = {1.0, 1.0, 1.0};
    float efc[3][3] = {{0.0, 0.0, 0.0},
                       {0.0, 0.0, 0.0},
                       {0.0, 0.0, 0.0}};
    float efcr[3] = {100.0, 100.0, 100.0};
    int efci = 0;
    /* helper */
    int c;

    static struct option cmdLineOptions[] = {
            {"help",       no_argument,       nullptr, 'h'},
            {"verbose",    no_argument,       nullptr, 'v'},
            {"rgbe",    no_argument,       nullptr, 'R'},
            {"pfs",    no_argument,       nullptr, 'P'},
            {"exact",    no_argument,       nullptr, 'e'},
            {"bayer",    no_argument,       nullptr, 'B'},
            {"debayer",    no_argument,       nullptr, 'D'},
            {"nominal",    no_argument,       nullptr, 'n'},
            {"deghosting", required_argument, nullptr, 'd'},
            {"tsv", no_argument, nullptr, 't'},
            { "saturation-offset", required_argument, nullptr, 'o' },
            { "range", required_argument, nullptr, 'r' },
            { "scale", required_argument, nullptr, 's' },
            { "rgbs", required_argument, nullptr, 'k' },
            { "efc", required_argument, nullptr, 'C' },
            { "oor-low", required_argument, nullptr, 'm' },
            { "oor-high", required_argument, nullptr, 'x' },
            {nullptr, 0,                         nullptr, 0}
    };

    std::stringstream k; //to read in multivalue arguments
    int optionIndex = 0;
    while ((c = getopt_long(argc, argv, "hnevuBDRd:s:r:o:m:x:k:C:", cmdLineOptions, &optionIndex)) != -1) {
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
            case 'B':
                isbayer = true;
                break;
            case 'D':
                isbayer = false;
                demosaic = true;
                break;
            case 'x':
                oor_high = atof(optarg);
                break;
            case 'm':
                oor_low = atof(optarg);
                break;
            case 'o':
                opt_saturation_offset_perc = atof(optarg);
                if( opt_saturation_offset_perc < 0 || opt_saturation_offset_perc > 0.25 )
                    throw pfs::Exception("saturation offset should be between 0 and 0.25");
                break;
            case 'r':
                opt_black_offset_perc = atof(optarg);
                if( opt_black_offset_perc < 0 || opt_black_offset_perc > 0.25 )
                    throw pfs::Exception("saturation offset should be between 0 and 0.25");
                break;
            case 's':
                opt_scale = atof(optarg);
                if( opt_scale <= 0)
                    throw pfs::Exception("scale must be positive");
                break;
            case 'k':
                k.str("");
                k.clear();
                k << optarg;
                k >> rgbcal[0] >> rgbcal[1] >> rgbcal[2];
                break;
            case 'C':
                k.str("");
                k.clear();
                k << optarg;
                k >> efcr[efci] >> efc[efci][0] >> efc[efci][1] >> efc[efci][2];
                efci++;
                break;
            case 'd':
                opt_deghosting = atof(optarg);
                if (opt_deghosting < 0)
                    throw pfs::Exception("deghosting threshold should be >0");
                break;
            case 't':
                dataonly = true;
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

    if (argv[optind] == nullptr){
        printHelp();
        throw QuietException();
    }

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


    // read through specification file loading header information and image frames
    while (true) {
        pmax = 0;
        pfs::Frame *iframe = nullptr;
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
                    if (!isbayer && comment.substr(begin, equal - begin) == "RGBcalibration"){
                        istringstream ss(comment.substr(equal+1, comment.size()));
                        ss >> rgbcal[0] >> rgbcal[1] >> rgbcal[2];
                    }
                    if ( comment.substr(begin, equal - begin) == "LuminanceRGB"){
                        istringstream ss(comment.substr(equal+1, comment.size()));
                        ss >> vlambda[0] >> vlambda[1] >> vlambda[2];
                    }
                    header << comment.substr(begin, comment.size()) << endl;
                }
                continue;
            }
            if (nominal)
                info = correct_exposure(info);
            info.factor = opt_scale * 100.0f * info.aperture * info.aperture / ( info.iso * info.etime );

            FILE *fh = fopen( framefile.c_str(), "rb");
            pfs::FrameFile ff = pfs::FrameFile( fh, framefile.c_str());
            HDRTiffReader reader( ff.fileName);
            iframe = pfsio.createFrame( reader.getWidth(), reader.getHeight() );
            pfs::Channel *X, *Y, *Z;
            iframe->createXYZChannels( X, Y, Z );
            reader.readImage( X, Y, Z );

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

        // for direct standard output (no merging)
        if (dataonly){
            fmax = info.factor;
            float r, g, b;
            for (int i = 0; i < size; i++) {
                r = ((*X)(i) * rgb_corr[0][0] + (*Y)(i) * rgb_corr[0][1] + (*Z)(i) * rgb_corr[0][2]) * rgbcal[0];
                g = ((*X)(i) * rgb_corr[1][0] + (*Y)(i) * rgb_corr[1][1] + (*Z)(i) * rgb_corr[1][2]) * rgbcal[1];
                b = ((*X)(i) * rgb_corr[2][0] + (*Y)(i) * rgb_corr[2][1] + (*Z)(i) * rgb_corr[2][2]) * rgbcal[2];
                std::cout << (*X)(i) << "\t" << (*Y)(i) << "\t" << (*Z)(i) << "\t";

                pmax = max((*X)(i), (*Y)(i), (*Z)(i), pmax);
                float below = min((*X)(i), (*Y)(i), (*Z)(i)) < opt_black_offset_perc;
                float above = max((*X)(i), (*Y)(i), (*Z)(i)) > 1 - opt_saturation_offset_perc;
                float lum = (vlambda[0]*r + vlambda[1]*g + vlambda[2]*b) * fmax;
                std::cout << r * fmax << "\t" << g * fmax << "\t" << b * fmax << "\t" << lum << "\t" << below << "\t" << above << std::endl;
            }
            std::cout  << std::endl;
        } else {
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
                if (info.etime < efcr[i])
                    break;
            }
            int s;
            float efcc = 1.0;
            for (int i = 0; i < height; i++)
                for (int j = 0; j < width; j++) {
                    s = j + i * width;
                    // apply electronic front curtain shutter correction
                    if (efc[efci][0] > 0){
                        efcc = 1/pow(((height-i) / efc[efci][2] + efc[efci][0] * info.etime)/(1 + efc[efci][0] * info.etime), efc[efci][1]);
                    }
                    (*eR.yi)(s) = (*X)(s) * efcc;
                    (*eG.yi)(s) = (*Y)(s) * efcc;
                    (*eB.yi)(s) = (*Z)(s) * efcc;
                    pmax = max((*X)(s), (*Y)(s), (*Z)(s), pmax);
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

        oorange = oorange && (pmax > 1 - opt_saturation_offset_perc);
        pfsio.freeFrame(iframe);
    }

    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++) {
            rgb_corr[i][j] = rgb_corr[i][j] * rgbcal[i];
        }
    }

    if (dataonly)
        return;
    if (frame_no < 1)
        throw pfs::Exception("at least one image required for calibration (check paths in hdrgen script?)");

    header << "HDR_SEQUENCE_COUNT= " << frame_no << endl;
    header << "HDR_VALID_RANGE= " << std::setprecision(3) << gmin << "-" << gmax << " cd/m^2" << endl;
    VERBOSE_STR << "using " << frame_no  << std::setprecision(3) << " frames, range min:" << gmin << ", max:" << gmax <<  endl;

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
    sp = linear_response(RGB_out, exposures, opt_saturation_offset_perc,
                         opt_black_offset_perc, opt_deghosting,
                         opt_scale, vlambda, rgb_corr, oor_high, oor_low, isbayer, demosaic);

    if (sp > 0) {
        float perc = ceilf(100.0f * sp / size);
        // this one might be important, so always report saturated pixels regardless of --verbose
        std::cerr << PROG_NAME << ": " << "saturated pixels found in " << perc << "% (" << sp
                  << " pixels) of the image!" << endl;
    }

    if (rgbe){
        RGBEWriter writer( stdout, true );
        std::string hstring = header.str().substr(0,-1);
        writer.writeImage( Xj, Yj, Zj, hstring );
    } else {
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
