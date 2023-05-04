#!/usr/bin/env python
import os
import sys
import re
import shlex
from math import log2
from subprocess import Popen, PIPE


def info_from_exif(img, correct):
    rawinfo = Popen(shlex.split(f"exiftool -ISO -ShutterSpeed -Aperture -s {img}"), 
                    stdout=PIPE).communicate()[0].decode(sys.stdin.encoding)
    rawinfo = re.split(r'[\n:]', rawinfo)
    iso = int(rawinfo[1].strip())
    if "/" in rawinfo[3]:
        shutter = rawinfo[3].strip().split("/")
        shutter = float(shutter[1])/float(shutter[0])
    else:
        shutter = 1/float(rawinfo[3].strip())
    aperture = float(rawinfo[5].strip())
    if correct:
        shutter = 2**(-round(log2(1/shutter)*3)/3)
        aperture = 2**(round(log2(aperture**2)*3,0)/6)
    return shutter, aperture, iso


def get_raw_frame(img, correct=False, overwrite=False):
    ppm = img + ".ppm"
    if overwrite or not os.path.isfile(ppm):
        Popen(shlex.split(f"dcraw_emu -4 -o 1 -w {img}")).communicate()
    rawinfo = info_from_exif(img, correct)
    return ppm, *rawinfo


def main(*imgs, s=False, c=False, o=False, **kwargs):
    ppms = [get_raw_frame(img, correct=c, overwrite=o) for img in imgs]
    for ppm, sh, ap, iso in sorted(ppms, key=lambda x: x[1]):
        if s:
            print(f"pfsin {ppm} | pfstag --set 'ISO={iso}' --set 'aperture={ap:.03f}' --set 'exposure_time={1/sh:.08f}'")
        else:
            print(f"{ppm} {iso} {ap:.03f} {1/sh:.08f}")


if __name__ == '__main__':
    
    usage = f"""
usage: {sys.argv[0]} [--shell/-s] [--overwrite/-o] [--correct/-c] [--help/-h]  img1 img2 ...

    shell:     output shell file for use with stdin of linearhdr: bash output.sh | linearhdr.
    overwrite: run dcraw_emu even if output file exists.
    correct:   apply correction to nominal aperture and shutter 
               speed values, use with linearhdr --exact.
    help:      print this message.
    """
    
    args = []
    kwargs = dict(s=False, c=False, o=False)
    for a in sys.argv[1:]:
        if a[0] == '-':
            f = a.strip("-")[0]
            if f in kwargs:
                kwargs[f] = True
            elif f == 'h':
                print(usage, file=sys.stderr)
                sys.exit(0)
            else:
                raise ValueError(f"Illegal option {a}")
        elif os.path.isfile(a):
            args.append(a)
        else:
            raise ValueError(f"{a} is not a file")
    main(*args, **kwargs)
