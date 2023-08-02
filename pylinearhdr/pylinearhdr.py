from subprocess import Popen, PIPE
import shlex
import sys
import os
import re
from math import log2
from raytools import io
import numpy as np


def camera_raw_values(img):
    ppm = img + "_raw.ppm"
    hdr = img + "_raw.hdr"
    Popen(shlex.split(f"dcraw_emu -disinterp -4 -o 0 -w -Z {ppm} {img}")).communicate()
    f = Popen(shlex.split(f"pfsin {ppm}"), stdout=PIPE)
    Popen(shlex.split(f"pfsout {hdr}"), stdin=f.stdout, stderr=PIPE).communicate()
    imrgb = io.hdr2carray(hdr)
    imar = np.max(imrgb, axis=0)
    os.remove(ppm)
    os.remove(hdr)
    return imar


def get_xyz_cam(img):
    rawinfo = Popen(shlex.split(f"raw-identify -v {img}"),
                    stdout=PIPE).communicate()[0].decode(sys.stdin.encoding)
    a = rawinfo.rsplit("XYZ->CamRGB matrix:", 1)[1].strip().split()[0:9]
    xyz_rgb = np.array([float(i) for i in a]).reshape(3,3)
    return xyz_rgb

def get_cam_rgb(img):
    rawinfo = Popen(shlex.split(f"raw-identify -v {img}"),
                    stdout=PIPE).communicate()[0].decode(sys.stdin.encoding)
    a = re.search(r"Derived D65 multipliers: .*", rawinfo)
    rgb = " ".join(a.group(0).strip().rsplit(None, 3)[-3:])
    return rgb



def info_from_exif(img, correct, times=True, fo=None):
    if fo is None:
        anom = []
        aexc = []
    else:
        anom = [i[0] for i in fo]
        aexc = [i[1] for i in fo]
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
        apmat = np.isclose(aperture, anom)
        if np.any(apmat):
            aperture = np.array(aexc)[apmat][0]
        else:
            aperture = 2**(round(log2(aperture**2)*3,0)/6)
    if times:
        time = "{}:{}:{}:{}:{}".format(*rawinfo[7:12])
    else:
        time = None
    return shutter, aperture, iso, time


def get_info(img, fields):
    fs = " ".join([f"-{i}" for i in fields])
    hinfo = Popen(shlex.split(f"exiftool {fs} -s {img}"),
                  stdout=PIPE).communicate()[0].decode(sys.stdin.encoding)
    return hinfo


def header_info(img, fields=("ColorTempAsShot", "Colorspace")):
    hinfo = get_info(img, fields)
    hinfo = re.sub(r"\s+:\s+", "= ", hinfo)
    hinfo = "\n".join([f"# HDR_{i.strip()}" for i in hinfo.strip().splitlines(keepends=False)])
    return hinfo


def name_by_exif(img, prefix=None):
    suff = img.rsplit(".", 1)[1]
    hinfo = get_info(img, ("ISO",  "Aperture",  "ShutterSpeed",  "ColorTempAsShot"))
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


def process_dcraw_opt(val, img, callexif=True, avg=False):
    try:
        return int(val)
    except ValueError:
        pass
    try:
        vals = [int(i) for i in val.strip().split()]
    except ValueError:
        pass
    else:
        if avg:
            return str(round(np.average(vals)))
        else:
            vals = " ".join([str(i) for i in vals])
            return f"'{vals}'"
    if callexif:
        rawinfo = Popen(shlex.split(f"exiftool -{val} -s {img}"),
                        stdout=PIPE).communicate()[0].decode(sys.stdin.encoding).split(":", 1)[1].strip()
        return process_dcraw_opt(rawinfo, img, False, avg)
    raise ValueError(f"Bad option given '{val}' could not be processed as a exiftool parameter")


def get_raw_frame(img, correct=True, overwrite=False, listonly=False, crop=None, bad_pixels=None,
                  black="PerChannelBlackLevel", white="LinearityUpperMargin", fo=None, ppm=None):
    black = process_dcraw_opt(black, img, avg=True)
    white = process_dcraw_opt(white, img, avg=True)
    if ppm is None:
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
    rawinfo = info_from_exif(img, correct, fo=fo)
    return ppm, *rawinfo


def process_colorspace_option(colorspace):
    cso = colorspace
    try:
        colorspace = [float(i) for i in colorspace.split()]
    except ValueError:
        if colorspace in ['rad', 'srgb']:
            return colorspace
    else:
        if len(colorspace) == 8:
            return colorspace[0:6], colorspace[6:]
    raise ValueError(f"bad value for colorspace: {cso}")

def pw_mtx(p, w):
    p = np.asarray(p).reshape(3, 2)
    wxyz = np.array([w[0], w[1], (1 - np.sum(w))])/w[1]
    z = (1 - np.sum(p, axis=1))[:, None]
    pxyz = np.hstack((p, z)).T
    c = np.einsum('ij,j->i', np.linalg.inv(pxyz), wxyz)
    return np.einsum('ij,j->ij', pxyz, c)


def cam_color_mtx(img, cs='rad'):
    # xyz->camRGB from adobeDNG/libraw/dcraw
    xyz_cam = get_xyz_cam(img)
    if cs == 'rad':
        cs = ((0.640, 0.330, 0.290, 0.600, 0.150, 0.060), (0.3333, 0.3333))
    elif cs == 'srgb':
        cs = ((0.64,  0.33,  0.3,  0.6,  0.15,  0.06), (0.3127, 0.329))
    rgb_xyz = pw_mtx(*cs)

    # rgb->camRGB
    rgb_cam = xyz_cam @ rgb_xyz
    # normalize
    n_cam_rgb = rgb_cam / np.sum(rgb_cam, axis=1, keepdims=True)
    # invert to camRGB->rgb
    cam_rgb = np.linalg.inv(n_cam_rgb)
    cam_rgbs = " ".join([f"{i:.05f}" for i in cam_rgb.ravel()])
    ps = " ".join([f"{i:.05f}" for i in cs[0]])
    ws = " ".join([f"{i:.05f}" for i in cs[1]])
    ls = " ".join([f"{i:.05f}" for i in rgb_xyz[1]])
    header = [f"# Camera2RGB= {cam_rgbs}",
              f"# TargetPrimaries= {ps}",
              f"# TargetWhitePoint= {ws}",
              f"# LuminanceRGB= {ls}"]
    return cam_rgb, header

def calibrate_frame(img, u, l, w, h, opts, bad_pixels, black="PerChannelBlackLevel",
                    white="LinearityUpperMargin", colorspace='rad', fo=None, scale=1.0, saturation=0.01, r=0.01):
    cam_rgb, header = cam_color_mtx(img, colorspace)
    ppm = img + "_calibrate.ppm"

    ppm, sh, ap, iso, _ = get_raw_frame(img, crop=(u,l, w, h), bad_pixels=bad_pixels, black=black, white=white, ppm=ppm, fo=fo)
    iso = iso/scale
    txt = img + "_calibrate.txt"
    f = open(txt, 'w')
    for h in header:
        print(h, file=f)
    print(f"{ppm} {iso} {ap:.03f} {1/sh:.08f}", file=f)
    f.close()
    opts += f" -o {saturation} -r {r}"
    f = Popen(shlex.split(f"linearhdr {opts} --exact --tsv {txt}"), stdout=PIPE)
    vals = Popen(shlex.split("total -m"), stdin=f.stdout, stdout=PIPE).communicate()[0].split()
    vals = [float(i) for i in vals] + [float(vals[-2]) + float(vals[-1])]
    os.remove(ppm)
    os.remove(txt)
    return ppm, sh, ap, iso, vals


def report(ppms, s=False, l=False, scale=1, sat_w=0.99, sat_b=.01, outf=None):
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


def report_calibrate(ppms, sort='shutter'):
    avg = 0
    div = 0
    minv = 1e9
    maxv = 0
    sorti = dict(image=0, shutter=1, aperture=2)
    si = sorti[sort]
    print("image\tiso\taperture\texposure_time\tsat_red\tsat_green\tsat_blue\tred\tgreen\tblue\tlum(assumes_D65)\tfrac_above\tfrac_below\tfrac_oor")

    ppms = sorted(ppms, key=lambda x: x[1])
    for ppm, sh, ap, iso, rgb in sorted(ppms, key=lambda x: x[si]):
        print(f"{ppm}\t{iso}\t{ap:.02f}\t{1/sh:.10f}\t" + "\t".join([f"{i:.04f}" for i in rgb]))
        if rgb[-2] == 0 and rgb[-1] == 0:
            minv = min(minv, rgb[-4])
            maxv = max(maxv, rgb[-4])
            avg += rgb[-4]
            div += 1
    print("Average value:", avg/div, file=sys.stderr)
    print("min-max for in range exposures:", minv, maxv, file=sys.stderr)
    print("adjust -o and -r settings of linearhdr to change viable range", file=sys.stderr)