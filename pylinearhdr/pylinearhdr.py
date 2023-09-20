# Copyright (c) 2023 Stephen Wasilewski, EPFL
# =======================================================================
# This program is free software: you can redistribute it and/or
# modify it under the terms of theGNU Lesser General Public License
# as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
# =======================================================================

from subprocess import Popen, PIPE
import shlex
import sys
import os
import re
from math import log2
from raytools import io
from raytools.mapper import ViewMapper
import numpy as np
from scipy.interpolate import UnivariateSpline


PREDEFINED_COLORS = dict(rad=((0.640, 0.330, 0.290, 0.600, 0.150, 0.060), (0.3333, 0.3333)),
                         srgb=((0.64,  0.33,  0.3,  0.6,  0.15,  0.06), (0.3127, 0.329)))


def camera_raw_values(img):
    tiff = img + "_raw.tiff"
    hdr = img + "_raw.hdr"
    Popen(shlex.split(f"rawconvert -disinterp -Z {tiff} {img}")).communicate()
    f = Popen(shlex.split(f"pfsin {tiff}"), stdout=PIPE)
    Popen(shlex.split(f"pfsout {hdr}"), stdin=f.stdout, stderr=PIPE).communicate()
    imrgb = io.hdr2carray(hdr)
    imar = np.max(imrgb, axis=0)
    os.remove(tiff)
    os.remove(hdr)
    return imar


def get_xyz_cam(img):
    rawinfo = Popen(shlex.split(f"rawconvert -identify {img}"),
                    stdout=PIPE).communicate()[0].decode(sys.stdin.encoding)
    a = rawinfo.rsplit("XYZ->CamRGB:", 1)[1].strip().split()[0:9]
    xyz_rgb = np.array([float(i) for i in a]).reshape(3, 3)
    return xyz_rgb


def get_d65_cam(img):
    rawinfo = Popen(shlex.split(f"rawconvert -identify {img}"),
                    stdout=PIPE).communicate()[0].decode(sys.stdin.encoding)
    a = rawinfo.rsplit("D65_multips:", 1)[1].strip().split()[0:3]
    xyz_rgb = np.array([float(i) for i in a])
    return xyz_rgb


def info_from_exif(img, correct, times=True, fo=None, shutterc=None):
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
        if shutterc is not None:
            shutter = shutter * np.exp(shutterc*shutter)
        apmat = np.isclose(aperture, anom)
        if np.any(apmat):
            aperture = np.array(aexc)[apmat][0]
        else:
            aperture = 2**(round(log2(aperture**2)*3,0)/6)
    if np.isinf(aperture):
        aperture = 1.0
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


def name_by_exif(img, prefix=None, aperture=True, iso=True, shutter=True):
    suff = img.rsplit(".", 1)[1]
    hinfo = []
    if iso:
        hinfo.append("ISO")
    if aperture:
        hinfo.append("Aperture")
    if shutter:
        hinfo.append("ShutterSpeed")
    if len(hinfo) > 0:
        hinfo = get_info(img, hinfo)
        hinfo = re.sub(r"\s+:\s+", " ", hinfo)
    else:
        hinfo = ""
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
            if v.lower() != "inf":
                v = "inf"
            else:
                v = str(int(round(float(v))))
        elif k == "ShutterSpeed":
            if "/" in v:
                v = v.split("/")[-1]
            k = "S-"
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


def get_raw_frame(img, correct=True, overwrite=False, listonly=False, crop=None, bad_pixels=None, bayer=False,
                  black="PerChannelBlackLevel", white="LinearityUpperMargin", fo=None, shutterc=None, tiff=None):
    correct = correct or fo or shutterc
    black = process_dcraw_opt(black, img, avg=True)
    white = process_dcraw_opt(white, img, avg=True)
    if tiff is None:
        tiff = img + ".tiff"
    if listonly:
        tiff = img
    elif overwrite or not os.path.isfile(tiff):
        cs = ""
        if crop is not None:
            cs = "-B {} {} {} {}".format(*crop)
        if bad_pixels is not None:
            cs += f" -P {bad_pixels}"
        if bayer:
            cs += " -disinterp"
        else:
            cs += " -q 11"
        Popen(shlex.split(f"rawconvert {cs} -Z {tiff} -k {black} -S {white} {img}")).communicate()
    rawinfo = info_from_exif(img, correct, fo=fo, shutterc=shutterc)
    return tiff, *rawinfo


def process_colorspace_option(colorspace):
    cso = colorspace
    try:
        colorspace = [float(i) for i in colorspace.split()]
    except ValueError:
        if colorspace in ['rad', 'srgb', 'raw']:
            return colorspace
    else:
        if len(colorspace) == 8:
            return colorspace[0:6], colorspace[6:]
    raise ValueError(f"bad value for colorspace: {cso}")


def pw_mtx(p, w):
    """calculate rgb->xyz matrix from primaries and whitepoint"""
    p = np.asarray(p).reshape(3, 2)
    wxyz = np.array([w[0], w[1], (1 - np.sum(w))])/w[1]
    z = (1 - np.sum(p, axis=1))[:, None]
    pxyz = np.hstack((p, z)).T
    c = np.einsum('ij,j->i', np.linalg.inv(pxyz), wxyz)
    return np.einsum('ij,j->ij', pxyz, c)


def mtx_pw(mtx):
    """calculate primaries and whitepoint from rgb->xyz matrix"""
    mtx = np.asarray(mtx).reshape(3,3)
    primaries = (mtx[0:2]/np.sum(mtx, axis=0)).T.ravel()
    whitepoint = np.sum(mtx, axis=1)[0:2]/np.sum(mtx)
    return primaries, whitepoint


def cam_color_mtx(xyzcam, cs='rad', cscale=None, normalize=True):
    """calculate camRGB->RGB from camera xyz->cam (from raw-identify or custom) and rgb primaries/whitepoint """
    # xyz->camRGB from adobeDNG/libraw/dcraw
    if cs == 'raw':
        return np.eye(3), [f"# Camera2RGB= 1 0 0 0 1 0 0 0 1"]
    if cs == 'rad':
        cs = PREDEFINED_COLORS['rad']
    elif cs == 'srgb':
        cs = PREDEFINED_COLORS['srgb']
    else:
        cs = (np.asarray(cs[0]).ravel(), np.asarray(cs[1]).ravel())
    rgb_xyz = pw_mtx(*cs)

    # rgb->camRGB
    rgb_cam = np.asarray(xyzcam).reshape(3, 3) @ rgb_xyz
    # normalize
    # if normalize:
    #     rgb_cam = rgb_cam / np.sum(rgb_cam, axis=1, keepdims=True)
    # invert to camRGB->rgb
    cam_rgb = np.linalg.inv(rgb_cam)
    cam_rgbs = " ".join([f"{i:.08f}" for i in cam_rgb.ravel()])

    ps = " ".join([f"{i:.04f}" for i in cs[0]])
    ws = " ".join([f"{i:.04f}" for i in cs[1]])
    ls = " ".join([f"{i:.08f}" for i in rgb_xyz[1]])
    header = [f"# Camera2RGB= {cam_rgbs}",
              f"# TargetPrimaries= {ps}",
              f"# TargetWhitePoint= {ws}",
              f"# LuminanceRGB= {ls}"]
    if cscale is not None:
        ccal = " ".join([str(i) for i in cscale])
        header.append(f"# RGBcalibration= {ccal}")
    return cam_rgb, header

# def cam_color_mtx(xyzcam, cs='rad', cscale=None):
#     """calculate camRGB->RGB from camera xyz->cam (from raw-identify or custom) and rgb primaries/whitepoint """
#     # xyz->camRGB from adobeDNG/libraw/dcraw
#     if cs == 'raw':
#         return np.eye(3), [f"# Camera2RGB= 1 0 0 0 1 0 0 0 1"]
#     if cs == 'rad':
#         cs = PREDEFINED_COLORS['rad']
#     elif cs == 'srgb':
#         cs = PREDEFINED_COLORS['srgb']
#     else:
#         cs = (np.asarray(cs[0]).ravel(), np.asarray(cs[1]).ravel())
#     rgb_xyz = pw_mtx(*cs)
#
#     # rgb->camRGB
#     rgb_cam = np.asarray(xyzcam).reshape(3, 3) @ rgb_xyz
#     cam_rgb = np.linalg.inv(np.asarray(xyzcam).reshape(3, 3)) @ np.linalg.inv(rgb_xyz)
#     print(xyzcam, rgb_xyz, cam_rgb, file=sys.stderr)
#     rgb_cam = rgb_xyz @ np.asarray(xyzcam).reshape(3, 3)
#     # normalize
#     # rgb_cam = rgb_cam / np.sum(rgb_cam, axis=1, keepdims=True)
#     # invert to camRGB->rgb
#     cam_rgb = np.linalg.inv(rgb_cam)
#     cam_rgbs = " ".join([f"{i:.08f}" for i in cam_rgb.ravel()])
#
#     ps = " ".join([f"{i:.04f}" for i in cs[0]])
#     ws = " ".join([f"{i:.04f}" for i in cs[1]])
#     ls = " ".join([f"{i:.08f}" for i in rgb_xyz[1]])
#     header = [f"# Camera2RGB= {cam_rgbs}",
#               f"# TargetPrimaries= {ps}",
#               f"# TargetWhitePoint= {ws}",
#               f"# LuminanceRGB= {ls}"]
#     if cscale is not None:
#         ccal = " ".join([str(i) for i in cscale])
#         header.append(f"# RGBcalibration= {ccal}")
#     return cam_rgb, header


def calibrate_frame(img, u, l, w, h, opts, bad_pixels, black="PerChannelBlackLevel", xyzcam=None, cscale=None, shutterc=None,
                    white="LinearityUpperMargin", colorspace='rad', fo=None, scale=1.0, saturation=0.01, r=0.01, bayer=False):
    if xyzcam is None:
        xyzcam = get_xyz_cam(img)
    cam_rgb, header = cam_color_mtx(xyzcam, colorspace, cscale=cscale)
    tiff = img + "_calibrate.tiff"
    tiff, sh, ap, iso, _ = get_raw_frame(img, crop=(u,l, w, h), bad_pixels=bad_pixels, black=black, white=white, tiff=tiff, fo=fo, shutterc=shutterc, bayer=bayer)
    iso = iso/scale
    txt = img + "_calibrate.txt"
    f = open(txt, 'w')
    for h in header:
        print(h, file=f)
    print(f"{tiff} {iso} {ap:.03f} {1/sh:.08f}", file=f)
    f.close()
    opts += f" -o {saturation} -r {r}"
    f = Popen(shlex.split(f"linearhdr {opts} --exact --tsv {txt}"), stdout=PIPE)
    vals = Popen(shlex.split("total -m"), stdin=f.stdout, stdout=PIPE).communicate()[0].split()
    vals = [float(i) for i in vals] + [float(vals[-2]) + float(vals[-1])]
    os.remove(tiff)
    os.remove(txt)
    return tiff, sh, ap, iso, vals


def report(tiffs, s=False, l=False, scale=1, sat_w=0.99, sat_b=.01, outf=None):
    if outf is None:
        outf = sys.stdout
    if l:
        print(f"Name Date ISO aperture etime shutter luminance range", file=sys.stderr)
    else:
        tiffs = sorted(tiffs, key=lambda x: x[1])
        if not s:
            capdate = sorted([tiffs[0][-1].strip(), tiffs[-1][-1].strip()])
            stop = capdate[-1].rsplit(" ", 1)[1]
            print(f"# CAPDATE= {capdate[0]}-{stop}", file=outf)
    for tiff, sh, ap, iso, time in tiffs:
        if l:
            fmax = scale * 100 * ap * ap * sh / iso
            print(f"{tiff} {time} {iso} {ap:.03f} {1/sh:.08f} {sh:.02f} = {sat_b*fmax:.02f} to {sat_w*fmax:.02f}", file=sys.stderr)
        elif s:
            print(f"pfsin {tiff} | pfstag --set 'ISO={iso/scale}' --set 'aperture={ap:.03f}' --set 'exposure_time={1/sh:.08f}'", file=outf)
        else:
            print(f"{tiff} {iso/scale} {ap:.03f} {1/sh:.08f}", file=outf)


def report_calibrate(tiffs, sort='shutter', target=None, header=True):
    avg = 0
    div = 0
    minv = 1e9
    maxv = 0
    sorti = dict(image=0, shutter=1, aperture=2)
    si = sorti[sort]
    if header:
        if target:
            print(f"image\tiso\taperture\texposure_time\tsat_red\tsat_green\tsat_blue\tred\tgreen\tblue\tlum\tfrac_above\tfrac_below\tfrac_oor\ttarget:{target}")
        else:
            print("image\tiso\taperture\texposure_time\tsat_red\tsat_green\tsat_blue\tred\tgreen\tblue\tlum\tfrac_above\tfrac_below\tfrac_oor")
    tiffs = sorted(tiffs, key=lambda x: x[1])
    for tiff, sh, ap, iso, rgb in sorted(tiffs, key=lambda x: x[si]):
        if target:
            print(f"{tiff}\t{iso}\t{ap:.02f}\t{1/sh:.10f}\t" + "\t".join([f"{i:.04g}" for i in rgb]) + f"\t{rgb[6]/target:.04g}")
        else:
            print(f"{tiff}\t{iso}\t{ap:.02f}\t{1/sh:.10f}\t" + "\t".join([f"{i:.04g}" for i in rgb]))
        if rgb[-2] == 0 and rgb[-1] == 0:
            minv = min(minv, rgb[-4])
            maxv = max(maxv, rgb[-4])
            avg += rgb[-4]
            div += 1
    print("Average value:", avg/div, file=sys.stderr)
    print("min-max for in range exposures:", minv, maxv, file=sys.stderr)
    print("adjust -o and -r settings of linearhdr to change viable range", file=sys.stderr)


def apply_vignetting_correction(img, vg, viewangle=180):
    vm = ViewMapper(viewangle=viewangle)
    imgv, vecs, _, _, _ = vm.init_img(img.shape[-1], features=img.shape[0])
    ang = vm.degrees(vecs)
    if vg.shape[1] == 4:
        for i in range(3):
            p = UnivariateSpline(vg[:, 0], vg[:, i+1], k=1, s=0)
            vf = p(ang).reshape(img.shape[1:])
            imgv[i] = img[i] * vf
    else:
        p = UnivariateSpline(vg[:, 0], vg[:, 1], k=1, s=0)
        vf = p(ang).reshape(img.shape[1:])
        imgv = img * vf
    return imgv


def str_primaries_2_mtx(inp, mtx=None):
    if mtx:
        rgb2xyz = np.linalg.inv(np.asarray([float(i) for i in mtx.strip().split()]).reshape(3,3))
        ps, ws = mtx_pw(rgb2xyz)
    elif inp in PREDEFINED_COLORS.keys():
        rgb2xyz = pw_mtx(*PREDEFINED_COLORS[inp])
        ps, ws = PREDEFINED_COLORS[inp]
    else:
        inp = np.fromstring(inp)
        rgb2xyz = pw_mtx(inp[0:6], inp[6:8])
        ps, ws = (inp[0:6], inp[6:8])
    return rgb2xyz, ps, ws
