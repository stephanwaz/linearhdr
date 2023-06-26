#!/usr/bin/env python
import os
import sys
import re
import shlex
from math import log2
from subprocess import Popen, PIPE


def info_from_exif(img):
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
    shutter = 2**(-round(log2(1/shutter)*3)/3)
    aperture = 2**(round(log2(aperture**2)*3,0)/6)
    return shutter, aperture, iso


def get_raw_frame(img, u, l, w, h, opts, bad_pixels):
    ppm = img + "_calibrate.ppm"
    cs = ""
    if bad_pixels is not None:
        cs += f"-P {bad_pixels}"
    Popen(shlex.split(f"dcraw_emu -4 -o 1 -B {u} {l} {w} {h} {cs} -w -Z {ppm} {img}")).communicate()
    sh, ap, iso = info_from_exif(img)
    txt = img + "_calibrate.txt"
    f = open(txt, 'w')
    print(f"{ppm} {iso} {ap:.03f} {1/sh:.08f}", file=f)
    f.close()
    f = Popen(shlex.split(f"linearhdr {opts} --use-rgb --tsv {txt}"), stdout=PIPE)
    vals = Popen(shlex.split("total -m"), stdin=f.stdout, stdout=PIPE).communicate()[0].split()
    vals = [float(i) for i in vals] + [float(vals[-2]) + float(vals[-1])]
    os.remove(ppm)
    os.remove(txt)
    return ppm, sh, ap, iso, vals


def report(ppms):
    avg = 0
    div = 0
    minv = 1e9
    maxv = 0
    for ppm, sh, ap, iso, rgb in sorted(ppms, key=lambda x: x[1]):
        print(f"{ppm}\t{iso}\t{ap:.02f}\t{1/sh:.10f}\t" + "\t".join([f"{i:.04f}" for i in rgb]))
        if rgb[-2] == 0 and rgb[-1] == 0:
            minv = min(minv, rgb[-4])
            maxv = max(maxv, rgb[-4])
            avg += rgb[-4]
            div += 1
    print("Average value:", avg/div, file=sys.stderr)
    print("min-max for in range exposures:", minv, maxv, file=sys.stderr)
    print("adjust -o and -r settings of linearhdr to change viable range", file=sys.stderr)


def main(*imgs, crop=(0,0,50,50), opts=""):
    ppms = [get_raw_frame(img, *crop, opts=opts) for img in imgs]
    report(ppms)



if __name__ == '__main__':
    
    usage = f"""
usage: {sys.argv[0]} '[linearhdr options]' [--help/-h] img1 img2 ... <left_edge> <top_edge> <width> <height>

    last four arguments identify the crop region for calibration. use linearhdr --help to open calibration instructions
    """
    
    args = []
    opts = ""
    if "-h" in sys.argv or "--help" in sys.argv:
        print(usage, file=sys.stderr)
        file = __file__.rsplit("/", 1)[0] + "/linearhdr_README.rst"
        os.system(f"more {file}")
        sys.exit(0)
    for a in sys.argv[1:-4]:
        if os.path.isfile(a):
            args.append(a)
        elif a.strip()[0] == '-':
            opts = a
        else:
            print(f"{a} is not a file", file=sys.stderr)
            print(usage, file=sys.stderr)
            sys.exit(1)
    try:
        crop = [int(i) for i in sys.argv[-4:]]
    except (ValueError, IndexError, TypeError):
        print(usage, file=sys.stderr)
        sys.exit(1)
    main(*args, crop=crop, opts=opts)
