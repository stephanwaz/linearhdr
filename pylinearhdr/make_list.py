#!/usr/bin/env python
import os
import sys
import re
import shlex
from math import log2
from subprocess import Popen, PIPE


def info_from_exif(img, correct, times=True):
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
        aperture = 2**(round(log2(aperture**2)*3, 0)/6)
    if times:
        time = "{}:{}:{}:{}:{}".format(*rawinfo[7:12])
    else:
        time = None
    return shutter, aperture, iso, time

def header_info(img, fields=("ColorTempAsShot", "Colorspace")):
    fs = " ".join([f"-{i}" for i in fields])
    hinfo = Popen(shlex.split(f"exiftool {fs} -s {img}"),
                    stdout=PIPE).communicate()[0].decode(sys.stdin.encoding)
    hinfo = re.sub(r"\s+:\s+", "= ", hinfo)
    hinfo = "\n".join([f"# HDR_{i.strip()}" for i in hinfo.strip().splitlines(keepends=False)])
    return hinfo


def name_by_exif(img, prefix=None):
    suff = img.rsplit(".", 1)[1]
    fs = "-ISO -Aperture -ShutterSpeed  -ColorTempAsShot"
    hinfo = Popen(shlex.split(f"exiftool {fs} -s {img}"),
                  stdout=PIPE).communicate()[0].decode(sys.stdin.encoding)
    hinfo = re.sub(r"\s+:\s+", " ", hinfo)
    if prefix is None:
        outn = img.rsplit(".", 1)[0]
    elif "/" in img:
        outn = img.rsplit("/", 1)[0] + "/" + prefix
    else:
        outn = prefix
    for i in hinfo.strip().splitlines(keepends=False):
        k, v = i.split()
        if k == "Aperture":
            k = "F-"
            v = str(int(round(float(v))))
        elif k == "ShutterSpeed":
            if "/" in v:
                v = v.split("/")[-1]
            k = "S-"
        elif k == "ColorTempAsShot":
            k = ""
            v += "k"
        outn += "_" + k + v
    return outn + "." + suff


def process_dcraw_opt(val, img, callexif=True):
    try:
        return int(val)
    except ValueError:
        pass
    try:
        vals = " ".join([str(int(i)) for i in val.strip().split()])
    except ValueError:
        pass
    else:
        return f"'{vals}'"
    if callexif:
        rawinfo = Popen(shlex.split(f"exiftool -{val} -s {img}"),
                        stdout=PIPE).communicate()[0].decode(sys.stdin.encoding).split(":", 1)[1].strip()
        return process_dcraw_opt(rawinfo, img, False)
    raise ValueError(f"Bad option given '{val}' could not be processed as a exiftool parameter")


def get_raw_frame(img, correct=False, overwrite=False, listonly=False, crop=None, bad_pixels=None,
                  black="AverageBlackLevel", white="LinearityUpperMargin"):
    black = process_dcraw_opt(black, img)
    white = process_dcraw_opt(white, img)
    ppm = img + ".ppm"
    if listonly:
        ppm = img
    elif overwrite or not os.path.isfile(ppm):
        cs = ""
        if crop is not None:
            cs = "-B {} {} {} {}".format(*crop)
        if bad_pixels is not None:
            cs += f" -P {bad_pixels}"
        Popen(shlex.split(f"dcraw_emu -c 0 -4 -o 0 {cs} -Z {ppm} -k {black} -r 2 1 2 1 -S {white} {img}")).communicate()
    rawinfo = info_from_exif(img, correct)
    return ppm, *rawinfo


def report(ppms, s=False, l=False, scale=1, sat_w=0.8, sat_b=.01, outf=None):
    if outf is None:
        outf = sys.stdout
    if l:
        print(f"Name Date ISO aperture etime shutter luminance range", file=sys.stderr)
    else:
        ppms = sorted(ppms, key=lambda x: x[1])
        if not s:
            capdate = sorted([ppms[0][-1].strip(), ppms[-1][-1].strip()])
            stop = capdate[-1].rsplit(" ", 1)[1]
            print(f"# CAPDATE= {capdate[0]}-{stop}", file=outf)
    for ppm, sh, ap, iso, time in ppms:
        if l:
            fmax = scale * 100 * ap * ap * sh / iso
            print(f"{ppm} {time} {iso} {ap:.03f} {1/sh:.08f} {sh:.02f} = {sat_b*fmax:.02f} to {sat_w*fmax:.02f}", file=sys.stderr)
        elif s:
            print(f"pfsin {ppm} | pfstag --set 'ISO={iso/scale}' --set 'aperture={ap:.03f}' --set 'exposure_time={1/sh:.08f}'", file=outf)
        else:
            print(f"{ppm} {iso/scale} {ap:.03f} {1/sh:.08f}", file=outf)


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
    
    NOTE: use pylinearhdr makelist for more extensive options
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
