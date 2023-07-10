#!/usr/bin/env python
import os
import sys
import re
import shlex
from math import log2
from subprocess import Popen, PIPE


def info_from_exif(img, correct, times=False):
    exifl = f"-ISO -ShutterSpeed -Aperture"
    if times:
        exifl += " -CreateDate"
    rawinfo = Popen(shlex.split(f"exiftool {exifl} -s {img}"),
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
    if times:
        time = "{}:{}:{}:{}:{}".format(*rawinfo[7:12])
    else:
        time = None
    return shutter, aperture, iso, time


def get_raw_frame(img, correct=False, overwrite=False, listonly=False, crop=None, bad_pixels=None):
    ppm = img + ".ppm"
    if listonly:
        ppm = img
    elif overwrite or not os.path.isfile(ppm):
        cs = ""
        if crop is not None:
            cs = "-B {} {} {} {}".format(*crop)
        if bad_pixels is not None:
            cs += f" -P {bad_pixels}"
        Popen(shlex.split(f"dcraw_emu -4 -o 1 {cs} -Z {ppm} -w {img}")).communicate()
    rawinfo = info_from_exif(img, correct, times=listonly)
    return ppm, *rawinfo


def report(ppms, s=False, l=False, scale=1, sat_w=0.8, sat_b=.01):
    if l:
        print(f"Name Date ISO aperture etime shutter luminance range", file=sys.stderr)
    else:
        ppms = sorted(ppms, key=lambda x: x[1])
    for ppm, sh, ap, iso, time in ppms:
        if l:
            fmax = scale * 100 * ap * ap * sh / iso
            print(f"{ppm} {time} {iso} {ap:.03f} {1/sh:.08f} {sh:.02f} = {sat_b*fmax:.02f} to {sat_w*fmax:.02f}", file=sys.stderr)
        elif s:
            print(f"pfsin {ppm} | pfstag --set 'ISO={iso/scale}' --set 'aperture={ap:.03f}' --set 'exposure_time={1/sh:.08f}'")
        else:
            print(f"{ppm} {iso/scale} {ap:.03f} {1/sh:.08f}")


def main(*imgs, s=False, c=False, o=False, l=False):
    if l:
        s = False
        c = False
        o = False
    ppms = [get_raw_frame(img, correct=c, overwrite=o, listonly=l) for img in imgs]
    report(ppms, s, l)


if __name__ == '__main__':
    
    usage = f"""
usage: {sys.argv[0]} [--shell/-s] [--overwrite/-o] [--correct/-c] [--list/-l] [--help/-h] img1 img2 ...

    shell:     output shell file for use with stdin of linearhdr: bash output.sh | linearhdr.
    overwrite: run dcraw_emu even if output file exists.
    correct:   apply correction to nominal aperture and shutter 
               speed values, use with linearhdr --exact.
    list:      skip execution and just print metadata
    help:      print this message.
    """
    
    args = []
    kwargs = dict(s=False, c=False, o=False, l=False)
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
