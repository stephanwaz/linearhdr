=========
linearhdr
=========

provides a simple command line tool for merging camera raw files
(after conversion with dcraw_emu) into HDR images.

Dependencies
------------

For usage:

    - dcraw_emu : https://www.libraw.org/download
    - exiftool : https://exiftool.org/
    - radiance (optional) : https://github.com/LBNL-ETA/Radiance/releases
    - pfstools (optional) : https://pfstools.sourceforge.net/pfstmo.html

For compiling:

    - cmake : https://cmake.org/download/ (also via pip, macports, etc...)
    - netpbm (optional if you have pfstools) : https://netpbm.sourceforge.net/getting_netpbm.php


Installation
------------

From a command line in this directory::

    mkdir build
    cd build
    cmake ..
    ccmake .. # to check/set install location (default is build/bin), <enter> to edit, then 'c' then 'g'
    make install

make sure that your install location is in you $PATH

This will install three programs::

    linearhdr: the main c++ executable
    linearhdr_make_list: a python helper script to generate input sequence files
    linearhdr_calibrate: a python script for determining a camera/lens/whitebal scale factor (linearhdr -s)

Calibration
-----------
Camera raw files have a broadly linear response, and because of reciprocity calculating absolute luminance across
ISO, Aperture, and shutter speed are possible across this linear range so long as a single scaling factor
(accounting for differences in ISO digital film standards, lens transmission, etc..) is known.

The linearhdr_calibrate script will determine this scaling factor and output data to check the valid upper limit of
the linear range (the -s parameter of linearhdr). Test images can be used to adjust -r, if this value is set too high
 (meaning low values are included) the result is noise in the output). For standard full stop hdr sequences, there is
little reason to adjust the default, as it provides more than enough dynamic range for sufficient overlap between
exposures.

To capture a calibration sequence, you either need a luminance meter, or a source with a known brightness (in cd/m^2).
You can find this information for many smart phones. Using published values will not be as reliable as a measurment,
but will get HDR images within a reasonable range. For example my iphone 12 mini has a measured brightness,
achieved with Konica-Minolta LS 110 by:

    1. disabling settings > display > truetone
    2. disabling settings > accessibility > display & text size > auto brightness
    3. setting brightness to max
    4. viewing an all white picture
    5. in a dark room
    5. measuring the center of the phone from a perpendicular vantage at least 1 m away.

of 587 cd/m^2. This is compared to review values of 627. This discrepancy is likely due to the age of my phone, and that
it has a scree protector, but without a luminance meter, using this value I would be within 10% of truth.

For the image sequence:

    1. place the phone (or other similar light source, monitor, etc..) in a dark room,
    2. as close (and centered) to the camera as possible (so the camera is focused)
    3. with ISO 100 and a middle aperture for your lens take a complete sequence of shutter speeds
    4. start at nearly black
    5. go to fully white
    6. increment in 1/3 stops (usually one rotation of the wheel).
    7. make sure you are capturing raw images

To use linearhdr_calibrate:

    1. identify the pixel coordinate of a 100-300 by 100-300 pixel area covering the perpendicular view to the phone
    2. mark down the upper left corner and check that the pixel dimensions will only see screen.
    3. you can do this with pfsin img.raw | pfsview, photoshop, or other image software.

the command for linearhdr_calibrate assuming CR3 files and pixel coordinates of left=2000, top==1000::

    linearhdr_calibrate sequence/*.raw 2000 1000 100 100

This will output for every frame the raw r g b channels, there exposure compensated luminance values (r g b l),
where l is photopic luminance, and two columns indicating with 0 or 1 whether the raw values are out of range. For
values between zero and one, this means that some pixels are out of range. The last lines will print the average and
the min and max taken from the rows where both the lower and upper ranges are 0. if the variance between max and min
is too much, examinine the data and determine if changing the cutoffs from the top (by adjusting -o) or the bottom
(by adjusting -r) would improve the result. give these arguments in quotes::

    linearhdr_calibrate '-o .2' sequence/*.raw 2000 1000 100 100

Once satisfied with the average take your measured value (ref) to calculate your camera's calibration::

    s = ref/camera_avg

always give this as an argument to linearhdr, or correct the output by multiplying it by this scale factor.

Usage
-----

For best results capture tripod mounted sequences with shutter speed varying by
one full stop (3 clicks) between frames, beginning with no white pixels
(or upper limit found in calibration) and ending with no black pixels. Most dSLR cameras have
a histogram display with the image preview to aid with this. ISO and aperture should be kept
constant, although in theory these will be properly compensated for. White balance should also
be held constant with any pre-calibration values.

linearhdr --help::

    linearhdr [Options] [exposure_list]
    Options:
        [--saturation-offset, -o <val>]: exclude images within <val> of 1 default=0.2
        [--range, -r <val>]: dynamic range of single raw exposure, used to set lower cutoff,
            give as power of 2 default=6.64386
        [--deghosting, -d <val>]: relative difference for outlier detection when less than 1,
            otherwise absolute difference (good for clouds) default=OFF
        [--tsv, -t]: output raw data as tsv, exposures seperated by extra linebreak,
            do not use with large files!
        [--scale, -s <val>]: absolute scaling for hdr (eg ND filter, known response, etc.) default=1.0
        [--use-yxy, -X]: merge hdr in Yxy space instead of RGB
        [--cull, -c]: reduce number of input exposures used
        [--rgbe, -R]: output radiance rgbe (default)
        [--pfs, -P]: output pfs stream
        [--exact, -e]: input camera values interpreted as exact (default=False)
        [--nominal, -n]: input camera values interpreted as nominal (default=True)
        [--verbose, -v]
        [--help]

    If exposure_list is given, images are read from file formatted as:
        <image1.ppm> <iso> <aperture> <exposure_time>
        <image2.ppm> <iso> <aperture> <exposure_time>
        ...

        list should be sorted by longest exposure time to shortest (only critical if --cull)
        else, program expects a sequence of images (formatted as from pfsin on the stdin),
    use/see 'make_hdr_list' for an example.
    By default, linearhdr expects nominal aperture and shutter speed.If using pfsinme, note that nominal camera values are manipulated by dcraw (but with less accuracy)
     so make sure to use the --exact flag so shutter and aperture are not double corrected.