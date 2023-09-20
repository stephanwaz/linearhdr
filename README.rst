=========
linearhdr
=========

provides a simple command line tool for merging camera raw files into HDR images.

Dependencies
------------

For usage:

    - exiftool : https://exiftool.org/
    - radiance : https://github.com/LBNL-ETA/Radiance/releases

For compiling:

    - cmake : https://cmake.org/download/ (also via pip, macports, etc...)
    - LibRaw::

        git clone https://github.com/LibRaw/LibRaw.git
        cd LibRaw
        git fetch --all --tags
        git checkout 0.21.1


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
