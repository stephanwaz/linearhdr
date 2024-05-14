=========
linearhdr
=========

provides a command line tool for merging camera raw files into HDR images.

Dependencies
------------

For usage:

    - exiftool : https://exiftool.org/
    - radiance : https://github.com/LBNL-ETA/Radiance/releases
    - python (>= 3.7) :  https://www.python.org/

For compiling:

    - cmake : https://cmake.org/download/ (also via pip, macports, etc...)

Installation
------------
First, clone this pacakge and LibRaw as a submodule::

        git clone https://github.com/stephanwaz/linearhdr.git
        cd linearhdr
        git clone https://github.com/LibRaw/LibRaw.git
        cd LibRaw
        git fetch --all --tags
        git checkout 0.21.1
        cd ..

Then build c++ code command line tools called by pylinearhdr (linearhdr, rawconvert).
From a command line in this directory::

    mkdir build
    cd build
    cmake ..
    ccmake .. # to check/set install location (default is build/bin), <enter> to edit, then 'c' then 'g'
    make install

make sure that your install location is in you $PATH.

Then install pylinearhdr (recommend creating a virtual environment and using pip, or other package manager such as conda)::

    pip3 install .

Upgrading
---------
to upgrade::

    git pull
    cd build
    make clean
    cmake ..
    make install
    cd ..
    pip install . --upgrade

Getting Started
---------------
Primarily all operations can be completed using the pylinearhdr command. See pylinearhdr --help for a list of subcommands. Merging
is done with pylinearhdr run "raw files" > merged.hdr. See pylinearhdr run --help for detailed settings. As an alternative to the run
command, the input list for linearhdr can be generated with pylinearhdr makelist "raw files" > merge_list.txt, which has most of the same
parameters. This allows for a two step process in case you would like to do any manual manipulation to the merge_list file.

Calibration
-----------
One of the main values in using pylinearhdr is that it provides an interface to include shutter speed, aperture, and
color/luminance calibration. These parameters can be saved in a .cfg file and loaded as::

    pylinearhdr -c calibration.cfg run "raw files" > merged.hdr

See::

    pylinearhdr shutter --help
    pylinearhdr aperture --help
    pylinearhdr calibrate --help

for procedures.

License, Copyrights and Acknowledgements
----------------------------------------

Copyright (C) 2023 Stephen Wasilewski, EPFL

This program is free software: you can redistribute it and/or
modify it under the terms of theGNU Lesser General Public License
as published by the Free Software Foundation, either version 3 of
the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.

Linearhdr uses LibRaw:

Copyright (C) 2008-2021 LibRaw LLC
LibRaw uses code from Dave Coffin’s dcraw.c utility (without RESTRICTED/GPL2 code):
Copyright 1997-2018 by Dave Coffin, dcoffin a cybercom o net
LibRaw uses DCB demosaic code by Jaceck Gozdz distributed under BSD license:
Copyright (C) 2010, Jacek Gozdz (mailto:cuniek@kft.umcs.lublin.pl)
LibRaw uses Roland Karlsson’s X3F tools source code, licensed under BSD license:
Copyright (c) 2010, Roland Karlsson (roland@proxel.se)

pfstools:

Copyright (C) 2003,2004 Rafal Mantiuk and Grzegorz Krawczyk

radiance:

Radiance v5.4 Copyright (c) 1990 to 2022, The Regents of the University of
California, through Lawrence Berkeley National Laboratory (subject to receipt
of any required approvals from the U.S. Dept. of Energy).  All rights reserved.
