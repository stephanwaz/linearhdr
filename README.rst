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
    convertcolor: a utility compatible with pfstools
    linearhdr_extract: a shell script for processing image regions


Additional optional python scripts can be installed with::

    pip3 install .

this will install::

    pylinearhdr

see pylinearhdr --help for all options

Calibration
-----------
Camera raw files have a broadly linear response, and because of reciprocity calculating absolute luminance across
ISO, Aperture, and shutter speed is possible across this linear range so long as a single scaling factor
(accounting for differences in ISO digital film standards, lens transmission, etc..) is known.

The linearhdr_calibrate script will determine this scaling factor and output data to check the valid upper limit of
the linear range (the -s parameter of linearhdr). Test images can be used to adjust -r, if this value is set too high
(meaning low values are included) the result is noise in the output. For standard full stop hdr sequences, there is
little reason to adjust the default, as it provides more than enough dynamic range for sufficient overlap between
exposures.

To capture a calibration sequence, you either need a luminance meter, or a source with a known brightness (in cd/m^2).
You can find this information for many smart phones. Using published values will not be as reliable as a measurment,
but will get HDR images within a reasonable range. For example my iphone 12 mini has a measured brightness of 587 cd/m^2,
achieved with Konica-Minolta LS 110 by:

    1. disabling: settings > display & brightness > truetone
    2. disabling: settings > accessibility > display & text size > auto brightness
    3. set to Never: settings > display & brightness > auto-lock
    4. setting brightness to max
    5. viewing an all white picture
    6. in a dark room
    7. measuring the center of the phone from a perpendicular vantage at least 1 m away
       (or minimum focal distance of luminance meter).

This is compared to review values of 627. This discrepancy is likely due to the age of my phone, and that
it has a scree protector, but without a luminance meter, using this value I would be within 10%.

For the image sequence:

    1. place the phone (or other similar light source, monitor, etc..) in a dark room,
    2. as close (and centered) to the camera as possible (so the camera is focused)
    3. with ISO 100 and a middle aperture for your lens take a complete sequence of shutter speeds
    4. start at nearly black
    5. go to fully white
    6. increment in 1/3 stops (usually one click of the wheel).
    7. make sure you are capturing raw images

To use linearhdr_calibrate:

    1. identify the pixel coordinate of a 100 to 300 by 100 to 300 pixel area covering the perpendicular view to the phone
    2. mark down the upper left corner and check that the pixel dimensions will only see screen.
    3. you can do this with pfsin img.raw | pfsview, photoshop, or other image software.

the command for linearhdr_calibrate assuming CR3 files and pixel coordinates of left=2000, top==1000::

    linearhdr_calibrate sequence/*.raw 2000 1000 100 100
    # or
    # pylinearhdr calibrate sequence/*.raw "2000 1000 100 100"

This will output for every frame the raw r g b channels, their exposure compensated luminance values (r g b l),
where l is photopic luminance, and two columns indicating with 0 or 1 whether the raw values are out of range. When
there are values between zero and one, only some pixels are out of range indicating that the view to the source area is
not perfectly consistent. Values close to zero and one are fine, but if most in range exposures are not close to zero,
check that the cropped region is indeed uniform. The last lines will print the average and the min and max taken from
the rows where both the lower and upper ranges are 0. If the variance between max and min is too much, examine the data
and determine if changing the cutoffs from the top (by adjusting -o) or the bottom (by adjusting -r) would improve the
result. give these arguments in quotes::

    linearhdr_calibrate '-o .2' sequence/*.raw 2000 1000 100 100
    # or
    # pylinearhdr calibrate -hdropts '-o .2' sequence/*.raw "2000 1000 100 100"

Once satisfied with the average take your measured value (ref) to calculate your camera's calibration::

    s = ref/camera_avg

always give this as an argument to linearhdr, or correct the output by multiplying it by this scale factor.

Usage
-----

Assuming a folder "HDR" containing a sequence of raw images (change extension to match) and a calibration scale of 1.4 to generate an HDR::

    linearhdr_make_list HDR/*.raw > HDR.txt
    # or
    # pylinearhdr HDR/*.raw > HDR.txt
    linearhdr -s 1.2 HDR.txt > HDR.hdr

For best results capture tripod mounted sequences with shutter speed varying by
one full stop (3 clicks) between frames, beginning with no white pixels
(or upper limit found in calibration) and ending with no black pixels. Most dSLR cameras have
a histogram display with the image preview to aid with this. ISO and aperture should be kept
constant, although in theory these will be properly compensated for. White balance should also
be held constant with any pre-calibration values. Always use the --scale value associated with the
particular camera and lens, as well as the --saturation-offset and --range identified during calibration.

linearhdr --help::

    linearhdr [Options] [exposure_list]
    Options:
        [--saturation-offset, -o <val>]: exclude images within <val> of 1 default=0.2
        [--range, -r <val>]: lower range of single raw exposure, used to set lower cutoff,
            give as value between 0 and 0.25, default=0.01
        [--deghosting, -d <val>]: relative difference for outlier detection when less than 1,
            otherwise absolute difference (good for clouds) default=OFF
        [--tsv, -t]: output raw data as tsv, exposures seperated by extra linebreak,
            do not use with large files!
        [--scale, -s <val>]: absolute scaling for hdr (eg ND filter, known response, etc.) default=1.0
        [--oor-low, -m <val>]: value to use for out of range low, default from data
        [--oor-high, -x <val>]: value to use for out of range high, default from data
        [--scale, -s <val>]: absolute scaling for hdr (eg ND filter, known response, etc.) default=1.0
            use linearhdr_calibrate to calculate
        [--use-yuv, -L]: merge hdr in YUV space instead of RGB (default)
        [--use-rgb, -C]: merge hdr in RGB space instead of YUV
        [--cull, -c]: throw away extra exposures that are not needed to keep output in range
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
    use/see 'linearhdr_make_list' for an example.
    By default, linearhdr expects nominal aperture and shutter speed.
    If using pfsinme, note that nominal camera values are manipulated by dcraw
     (but with less accuracy) so make sure to use the --exact flag so shutter
    and aperture are not double corrected.

The "range" option can be used to set the low end acceptable value, by default and raw values below .01
are counted out of ranges, but for some raw images with higher bit depth there may be useful information in
this low end that could reduce noise. Alternative, low bit depth or less reliable cameras may be too noisy in this
range to provide useful signal:noise ratios. by extending the range parameter, it is possible to build HDR images from
more widely spaced exposures.

The "deghosting" option can attempt to remove moving elements from the sequence. It will use the first image in
the sequence as the reference, assuming the exposure list is order by longest exposure first, this will be the
pixel with the least photon noise. To prioritize a different frame, list that image first in the sequence (but note
that if this is the shortest exposure is not out of range this is incompatible with cull,
as all subsequent exposures will be skipped. The
deghosting works by excluding exposures that deviate from this reference by a given relative factor (when less than 1)
or an absolute factor (when greater than 1). use a relative value to remove object motion (people cars) use an
absolute value to isolate deghosting to the sky (moving clouds/sun).

the "tsv" option is for debugging, raw data analysis and simply dumps the exposure values (raw and compensated) to
the standard output. output columns depend on RGB or Yxy output.
for RGB: R_exp G_exp B_exp R_lum G_lum B_lum lum   below above
for Yxy: Y_exp x     y     Y_lum x     y     Y_lum below above

The "use-yuv" merges hdr in Yuv space, this should not be used for calibration unless the source is perfectly
matched to the white balance of the camera, but does do a better job holding luminance calibration across saturated
colors.

The "use-rgb" merges hdr in rgb space, this should be used for calibration and is appended by default to
linearhdr_calibrate call.

"cull" can be useful when deghosting fails as a way to reduce redundant date in the brightest part of the image
introduced by moving clouds and sun positions, by eliminating exposures that do not add to the usable dynamic range.

"rgbe/pfs" select output format, the last flag takes priority, but "tsv" overrides both.

"exact" directly uses aperture and shutter values for exposure compensation.

"nominal" will correct aperture by F = 2^(round(log2(fn^2) * 3) / 6) and
exposure time by T = 2^(round(log2(1/T) * 3) / 3)

