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


import os
import shutil
import sys

import numpy as np
from clasp import click
import clasp.click_ext as clk
from clasp.script_tools import try_mkdir, clean_tmp, pipeline, sglob

import pylinearhdr
import raytools
from raytools import io, imagetools
from raytools.utility import pool_call
from pylinearhdr import shadowband as sb
from pylinearhdr import pylinearhdr as pl


def get_profiles():
    d = os.path.dirname(pylinearhdr.__file__)
    try_mkdir(f"{d}/pylinearhdr_profiles")
    profiles = sglob(f"{d}/pylinearhdr_profiles/*.cfg")
    pnames = [os.path.basename(p)[:-4] for p in profiles]
    return pnames


def get_profile(profile, forceexist=True):
    d = os.path.dirname(pylinearhdr.__file__)
    f = f"{d}/pylinearhdr_profiles/{profile}.cfg"
    if os.path.isfile(f) or not forceexist:
        return f
    else:
        print(f"{profile} profile not found in {d}/pylinearhdr_profiles/. choose from {global_profiles} or use -saveprofile", file=sys.stderr)
        raise click.Abort()


global_profiles = get_profiles()

@clk.pretty_name("HDR, TSV, FLOATS,FLOATS")
def hdr_data_load(ctx, param, s):
    """read np array from command line

    trys np.load (numpy binary), then np.loadtxt (space seperated txt file)
    then split row by spaces and columns by commas.
    """
    if s is None:
        return s
    if s == '-':
        s = clk.tmp_stdin(ctx)
    if os.path.exists(s):
        try:
            ar = io.load_txt(s)
        except ValueError:
            ar = s
        else:
            if len(ar.shape) == 1:
                ar = ar.reshape(-1, 3)
        return ar
    else:
        return np.array([[float(i) for i in j.split(',')] for j in s.split()]).reshape(-1, 3)


@click.group(invoke_without_command=True)
@click.option('-config', '-c', type=click.Path(exists=True),
              help="path of config file to load, if given, ignores profile")
@click.option('-profile', help=f"name of saved profile to load, options: {global_profiles}")
@click.option('-save',
              help="name of profile to save -config to.")
@click.option('-n', default=None, type=int,
              help='sets the environment variable RAYTOOLS_PROC_CAP set to'
                   ' 0 to clear (parallel processes will use cpu_limit)')
@click.option('--opts', '-opts', is_flag=True,
              help="check parsed options")
@click.option('--debug', is_flag=True,
              help="show traceback on exceptions")
@click.version_option(version=pylinearhdr.__version__)
@click.pass_context
def main(ctx, config=None, profile=None, save=None, n=None,  **kwargs):
    """the pylinearhdr executable is a command line interface to utility commands
    as part of the linearhdr package, some scripts available without any dependencies via cmake install. Run without
    a subcommand to check/set config and profiles.
    """
    if save is not None:
        if config is None:
            raise ValueError("-save requires -config")
        if save in global_profiles:
            if not click.confirm("Overwrite existing profile?"):
                while save in global_profiles:
                    save = click.prompt("Please enter a new name for profile:")
        nprofile = get_profile(save, False)
        shutil.copy(config, nprofile)
    if config is None and profile is not None:
        config = get_profile(profile)
    if ctx.invoked_subcommand is not None:
        os.environ['CLASP_PIPE'] = '1'
        raytools.io.set_nproc(n)
        ctx.info_name = 'pylinearhdr'
        ctx.obj = dict(temps=[])
        clk.get_config(ctx, config, None, None, None)
    else:
        print(f"config file loaded from: {config}", file=sys.stderr)
        print(open(config).read())


def is_file_or_tup_int(ctx, param, s):
    try:
        pix = clk.tup_int(ctx, param, s)
    except Exception:
        return clk.is_file(ctx, param, s)
    else:
        if pix is not None:
            path = clean_tmp(ctx)
            f = open(path, 'w')
            for p in pix:
                print(f"{p[0]}\t{p[1]}\t0", file=f)
            f.close()
            return path
    return None


vignetting_path = os.path.dirname(pylinearhdr.__file__) + f"/pylinearhdr_profiles"
vignetting_files = [i.rsplit("/", 1)[-1] for i in sglob(f"{vignetting_path}/*.txt")]

def is_vignette_file(ctx, param, s):
    """checks input file string with recursive prompt

    use os.environ['CLASP_PIPE'] = '1' in parent script
    or set CLASP_PIPE=1
    to disable prompt and avoid hanging process
    """
    if s in [None, 'None', 'none']:
        return None
    try:
        if os.path.exists(s):
            return s
        else:
            raise ValueError(s)
    except ValueError as e:
        s2 = f"{vignetting_path}/{s}"
        if os.path.exists(s2):
            return s2
        else:
            raise ValueError(f"{s} is not a file in current directory or {vignetting_path}")

shared_run_opts = [
    click.option("-badpixels", callback=is_file_or_tup_int,
                                help="file of bad pixels (rows: xpix ypix 0) where xpix is from left and ypix is "
                                     "from top or list of bad pixes x1,y1 x2,y2 etc."),
    click.option("-black", default="PerChannelBlackLevel",
                help="rawconvert darkness level. either a number(s) or exiftool key, it is critical this is kept "
                     "consistent between calibration and make_list/run. Possible options: AverageBlackLevel, "
                     "PerChannelBlackLevel, 2049 '2049 2049 2049 2049'"),
    click.option("-white", default="LinearityUpperMargin",
                help="rawconvert saturation level. either a number or exiftool key, it is critical this is kept "
                     " consistent between calibration and make_list/run. Possible options: NormalWhiteLevel, "
                     "SpecularWhiteLevel, LinearityUpperMargin, 10000"),
    click.option("-colorspace", default='rad',
                 help="by default outputs radiance RGB (defined by primaries/whitepoint:"
                      " ((0.640, 0.330, 0.290, 0.600, 0.150, 0.060), (0.3333, 0.3333)). other options are:"
                      " raw: do not convert colors from cam_raw"
                      " srgb: ((0.64,  0.33,  0.3,  0.6,  0.15,  0.06), (0.3127, 0.329) or a list of 8 values "
                      "(primaries + wp) the cam2rgb matrix and output primaries are written into the make_list header and "
                      " added to the output hdr header by linearhdr"),
    click.option("-xyzcam", callback=clk.split_float, help="custom xyz->cam matrix, if not given uses raw-identify"),
    click.option("-fo", "-f-overrides", callback=clk.tup_float,
                 help="if given, sets --correct to True. pairs of nominal/exact aperture values to correct. any aperture"
                      " not given will use the standard correction. for example give 11,11.4 22,23 to override standard"
                      " corrections for F11 and F22 (would be 11.31 and 22.63)"),
    click.option("-shutterc", "-sc", type=float,
                 help="if given, sets --correct to True. for correcting shutter speed. takes a single value 'A' "
                      "interpreted as: shutter=x*exp(A*x),  where x is the corrected (true power of 2) shutter speed. "
                      "this curve can be derived from a  sequence of hdr images using single files taken with "
                      "different shutter speeds. fit the function y=exp(A*x) where y = lum_s0/lum_st and lum_s0 is "
                      "the luminance of the  target at the longest exposure time. this can be done with an "
                      "exponential trendline in excel (make sure set intercept=1), or with "
                      "scipy.optimize.curve_fit(lambda t,a: np.exp(a*t),  x,  y, p0=(-1e-8,))"),
    click.option("-scale", default=1.0, help="calibration scale factor (applies to ISO, so do not use -s with linearhdr)"),
    click.option("-cscale", callback=clk.split_float, help="color calibration scale factor (applies via header, same as linearhdr -k)"),
    click.option("-nd", default=0.0, help="additional ND filter (applies to ISO, so do not use -s with linearhdr)"),
    click.option("-saturation", "-saturation-offset", "-s", default=0.01, help="saturation offset, if white is not LinearityUpperMargin, this must be changed"),
    click.option("-range", "-r", default=0.01, help="lower range of single raw exposure"),
    click.option("--verbose/--no-verbose", default=False, help="passed to linearhdr"),
    click.option("--interpfirst/--interpsecond", default=True, help="interpolate with raw convert (uses linear) or interpolate after merge (uses DHT)"),
    click.option("--bayer/--no-bayer", default=False, help="do not interpolate raw channels. forces -colorspace to 'raw' and ignores --interpfirst")
]


@main.command()
@click.argument("imgs", callback=clk.are_files)
@click.argument("crop", callback=clk.split_int)
@click.option("-hdropts", default="", help="options to pass to linearhdr (put in qoutes)")
@click.option("-sort", default="shutter", help="cane be image, aperture, shutter")
@click.option("-target", type=float, help="reference value")
@click.option("--header/--no-header", default=True, help="print column labels")
@clk.shared_decs(shared_run_opts)
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def calibrate(ctx, imgs, crop, badpixels=None, scale=1.0, nd=0.0, saturation=0.01, range=0.01, hdropts="",
              black="AverageBlackLevel", white="AverageBlackLevel", cscale=None, shutterc=None,
              colorspace='rad', fo=None, sort='shutter', target=None, header=True, xyzcam=None, verbose=False,
              bayer=False, interpfirst=True, **kwargs):
    """calibration routine, see README.rst

    imgs: list of images
    crop: help="<left> <upper> <width> <height>"
    """
    if verbose:
        hdropts += " --verbose"
    if bayer:
        colorspace = 'raw'
        hdropts += " -B"
    elif not interpfirst:
        hdropts += " -D"
    colorspace = pl.process_colorspace_option(colorspace)
    tiffs = pool_call(pl.calibrate_frame, imgs, *crop[0:4], hdropts, bad_pixels=badpixels, expandarg=False, black=black, pbar=False, shutterc=shutterc,
                     white=white, colorspace=colorspace, fo=fo, scale=scale * 10**nd, cscale=cscale, saturation=saturation, r=range, xyzcam=xyzcam,
                      bayer=bayer or (not interpfirst))
    pl.report_calibrate(tiffs, sort=sort, target=target, header=header)


make_list_opts = [
     click.argument("imgs", callback=clk.are_files),
     click.option("--shell/--no-shell", default=False,
                   help="output shell file for use with stdin of linearhdr: bash output.sh | linearhdr"),
    click.option("--fisheye/--no-fisheye", default=False,
                 help="apply fisheye_corr to 180 degree image (must be properly cropped and equiangular). "
                      "requires pcomb and RAYPATH"),
     click.option("--overwrite/--no-overwrite", default=False,
                   help="run rawconvert even if output file exists"),
     click.option("-header-line", '-hl', multiple=True,
                   help="lines to append to HDR header, e.g. LOCATION= 46.522833,6.580500"),
     click.option("--correct/--no-correct", default=True,
                   help="apply correction to nominal aperture and shutter speed values, use with linearhdr --exact"),
     click.option("--listonly/--no-listonly", default=False,
                   help="skip execution and just print metadata"),
     click.option("-hdropts", default="", help="additional options to linearhdr (with callhdr, overrides -r -s)"),
     click.option("-crop", callback=clk.split_int,
                   help="crop tiff (left upper W H)"),
]


def makelist_run(ctx, imgs, shell=False, overwrite=False, correct=False, listonly=False, scale=1.0, nd=0.0, saturation=0.01, range=0.01,
                 crop=None, badpixels=None, callhdr=False, hdropts="", fo=None, fisheye=False, xyzcam=None, cscale=None, shutterc=None,
                 black="AverageBlackLevel", white="AverageBlackLevel", colorspace='rad', clean=False, vfile=None, verbose=False, bayer=False,
                 interpfirst=True, **kwargs):
    """make list routine, use to generate input to linearhdr"""
    if verbose:
        hdropts += " --verbose"
    if bayer:
        colorspace = 'raw'
        hdropts += " -B"
    elif not interpfirst:
        hdropts += " -D"
    outf = sys.stdout
    outfn = "<makelist.txt>"
    if listonly:
        shell = False
        overwrite = False
    elif callhdr:
        outfn = clean_tmp(ctx)
        outf = open(outfn, 'w')
        if not correct:
            hdropts += " --nominal"
    tiffs = pool_call(pl.get_raw_frame, imgs, correct=correct, overwrite=overwrite, black=black, white=white, fo=fo, bayer=bayer or (not interpfirst),
                     shutterc=shutterc, listonly=listonly, crop=crop, bad_pixels=badpixels, expandarg=False, pbar=False)
    if xyzcam is None:
        xyzcam = pl.get_xyz_cam(imgs[0])
    cam_rgb, header = pl.cam_color_mtx(xyzcam, colorspace, cscale=cscale)
    for h in header:
        print(h, file=outf)
    pl.report(tiffs, shell, listonly, scale=scale * 10**nd, sat_w=1-saturation, sat_b=range, outf=outf)

    command = [f"linearhdr -r {range} -o {saturation} {hdropts} {outfn}"]
    if fisheye:
        command += ["pcomb -f fisheye_corr.cal -o - ", "getinfo -a 'VIEW= -vta -vh 180 -vv 180'"]
    if vfile is not None:
        ocmd = sys.argv
        if "run" in ocmd:
            ocmd = ocmd[:ocmd.index("run")]
        else:
            ocmd = ocmd[:ocmd.index("makelist")]
        ocmd = " ".join(ocmd) + f" vignette - -vfile {vfile}"
        command += [ocmd]

    if callhdr:
        outf.close()
        pipeline(command, outfile=sys.stdout)
        if clean:
            for tiff in tiffs:
                os.remove(tiff[0])
        clk.tmp_clean(ctx)
    else:
        print(" | ".join(command), file=sys.stderr)


@main.command()
@clk.shared_decs(make_list_opts + shared_run_opts)
@click.option("--callhdr/--no-callhdr", default=False, help="directly call linearhdr")
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def makelist(ctx, imgs, **kwargs):
    """make list routine, use to generate input to linearhdr"""
    makelist_run(ctx, imgs, **kwargs)


@main.command()
@clk.shared_decs(make_list_opts + shared_run_opts)
@click.option("--clean/--no-clean", default=True, help="delete tiff files after linearhdr")
@click.option("-vfile", callback=is_vignette_file, help="vignetting file, rows should be angle(degrees) factor(s) either one column"
                                                        "for luminance or 3 for RGB, apply in destination RGB space after lens "
                                                        " corrections. system stored files :\n" + "\n".join(vignetting_files))
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def run(ctx, imgs, **kwargs):
    """make list routine, use to generate input to linearhdr"""
    makelist_run(ctx, imgs, callhdr=True, **kwargs)


@main.command()
@click.argument("imgh", callback=clk.is_file)
@click.argument("imgv", callback=clk.is_file)
@click.argument("imgn", callback=clk.is_file)
@click.option("-outf", default="blended.hdr",
              help="output destination")
@click.option("-roh", default=0.0,
              help="rotation correction (degrees ccw) for horizontal band")
@click.option("-rov", default=0.0,
              help="rotation correction (degrees ccw) for vertical band")
@click.option("-sfov", default=2.0,
              help="field of view around shaded source that is valid in imgn (in case of partial ND filter)")
@click.option("-srcsize", default=6.7967e-05,
              help="solid angle of shaded source (steradians)")
@click.option("-bw", default=2.0,
              help="size (in degrees) of shadow band. Note: this value is doubled to account for error in centering on "
                   "source")
@click.option("-sunloc", callback=clk.split_int,
              help="pixel location of sun used to set band location. If not given auto-locates based on max pixels in"
                   "imgn.")
@click.option("--flip/--no-flip", default=False,
              help="by default the imgh will be used in the UL and LR quadrants, flip=True will use imgh UR and LL")
@click.option("-envmap",
              help="additionally save an hdr file suitable for use as an environment map, source information is in the header")
@click.option("-check",
              help="write out masks and other intermediate data (give as prefix)")
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def shadowband(ctx, imgh, imgv, imgn, outf="blended.hdr", roh=0.0, rov=0.0, sfov=2.0, srcsize=6.7967e-05, bw=2.0, flip=False, envmap=None, sunloc=None, check=None,
               **kwargs):
    """merge set of three images with horizontal, vertical and no shadow band all images are assumed to be 180 degree
    angular fisheye.

    imgh: hdr with horizontal aligned shadowband
    imgv: hdr with veritcal aligned shadowband
    imgn: hdr with no shadowband
    """
    if flip:
        rov = -rov
        roh = -roh
    hdata, hh = io.hdr2carray(imgh, header=True)
    vdata, hv = io.hdr2carray(imgv, header=True)
    sdata, hs = io.hdr2carray(imgn, header=True)
    blended, skyonly, source = sb.shadowband(hdata, vdata, sdata, roh=roh, rov=rov, sfov=sfov, srcsize=srcsize, bw=bw, flip=flip,
                                             envmap=envmap, sunloc=sunloc, check=check)
    vm = imagetools.hdr2vm(imgh)
    header = []
    if vm is not None:
        header = [vm.header()]
    header += hh + hv + hs
    io.carray2hdr(blended, outf, header)
    if skyonly is not None:
        io.carray2hdr(skyonly, envmap, header)
    # if envmap:
    #     radout = outf.rsplit(".", 1)[0] + ".rad"
    #     f = open(radout, 'w')
    #     srcrad = (srcsize/np.pi)**.5 * 180/np.pi * 2
    #     f.write("void light sun 0 0 3 {} {} {}\n".format(*blended[1][-1]))
    #     f.write("sun source solar 0 0 4 {} {} {} {}\n\n".format(*blended[1][0:3], srcrad))
    #     f.write(f"void colorpict imgfunc\n7 red green blue {outf} fisheye.cal fish_u fish_v\n0 0\n")
    #     f.write("imgfunc glow imgglow 0 0 4 1 1 1 0\n")
    #     f.write("imgglow source sky 0 0 4 0 1 0 180\n")
    #     io.carray2hdr(blended[0], outf, header)
    # else:
    #     io.carray2hdr(blended, outf, header)


@main.command()
@click.argument("imgs", callback=clk.are_files)
@click.option("-out", default="img", help="directory base name, becomes: img_XXX, use 'common'"
                                          " to auto-generate from common parts (seperated by _) or"
                                          "lcommon to shorten common parts by across group differences")
@click.option("-starti", default=0, help="start index")
@click.option("--ascend/--descend", default=True,
              help="shot order exposure time. ascending=shortest->longest")
@click.option("--r/--no-r", "--rename/--no-rename", default=False,
              help="also call rename")
@click.option("--preview/--no-preview", default=False,
              help="print outcome without moving")
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def sort(ctx, imgs, out="img", starti=0, ascend=True, preview=False, r=False, **kwargs):
    """sort sequence of images into folders based on changes in exposure time"""
    if r:
        rename.callback(imgs, **kwargs)
        imgs = ctx.obj['imgs']
    infos = pool_call(pl.get_raw_frame, imgs, listonly=True, correct=True, expandarg=False, pbar=False)
    i = starti
    if ascend:
        expt = 1e9
    else:
        expt = 0
    groups = [[]]
    for info in infos:
        if info[1] == expt:
            print("Duplicate!")
            continue
        if (info[1] < expt) != ascend:
            groups.append([])
            i += 1
        expt = info[1]
        groups[-1].append(info[0])
    outns = []
    for j, group in enumerate(groups):
        if out[-6:] == 'common':
            pieces = [g.rsplit(".", 1)[0].split("_") for g in group]
            pcount = min(len(g) for g in pieces)
            pieces = [g.rsplit(".", 1)[0].split("_", pcount) for g in group]
            match = [pieces[0][i] for i in range(pcount) if len(set([p[i] for p in pieces]))==1]
            outns.append("_".join(match))
        else:
            outns.append(out)
    if out[0] in ['l', 'n']:
        pieces = [g.rsplit(".", 1)[0].split("_") for g in outns]
        pcount = min(len(g) for g in pieces)
        pieces = [g.rsplit(".", 1)[0].split("_", pcount) for g in outns]
        match = [len(set([p[i] for p in pieces]))>1 for i in range(pcount)]
        outns = ["_".join([g[j] for j, m in enumerate(match) if m]) for g in pieces]
    for j, (group, outn) in enumerate(zip(groups, outns)):
        if preview:
            print(f"\n{outn}_{j:03d}/")
        else:
            try_mkdir(f"{outn}_{j:03d}")
        for g in group:
            if preview:
                print(f"\t{g}")
            else:
                shutil.move(g, f"{outn}_{j:03d}")
    # print([len(i) for i in groups])


@main.command()
@click.argument("imgs", callback=clk.are_files)
@click.option("-out", help="base name, becomes: img_XXX, if None, appends info to current name")
@click.option("--copy/--move", default=True, help="copy or move")
@click.option("--aperture/--no-aperture", default=True, help="name by aperture")
@click.option("--iso/--no-iso", default=False, help="name by iso")
@click.option("--shutter/--no-shutter", default=True, help="name by shutter")
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def rename(ctx, imgs, out=None, copy=True, aperture=True, shutter=True, iso=False, **kwargs):
    """rename raw files based on ISO_APERTURE_SHUTTER_CCT"""
    infos = pool_call(pl.name_by_exif, imgs, expandarg=False, prefix=out, pbar=False, aperture=aperture, shutter=shutter, iso=iso)
    copied = []
    for orig, dest in zip(imgs, infos):
        dest2 = dest
        i = 0
        while dest2 in copied:
            dest2 = dest.rsplit(".", 1)
            dest2 = "{}_{:02d}.{}".format(dest2[0], i, dest2[1])
            i += 1
        if copy:
            shutil.copy(orig, dest2)
        else:
            shutil.move(orig, dest2)
        copied.append(dest2)
    ctx.obj['imgs'] = copied


@main.command()
@click.argument("imgs", callback=clk.are_files)
@click.option("-threshold", default=0.01, help="pixels with raw signal greater than this value across all imgs are "
                                               "considered bad.")
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def badpixels(ctx, imgs, threshold=0.01, **kwargs):
    """find bad pixels and generate input for calibrate/makelist

    imgs: 2-6 raw format files taken with a range of shutter speeds with the lens cap on.
    """
    hdrs = pool_call(pl.camera_raw_values, imgs, expandarg=False)
    mask = np.all(np.stack(hdrs, axis=0) > threshold, axis=0)
    yres = mask.shape[1]
    badpixi = np.stack(np.unravel_index(np.arange(mask.size)[mask.ravel()], mask.shape)).T
    badpixi[:, 1] = yres - badpixi[:, 1] - 1
    for bp in badpixi:
        print("{}\t{}\t0".format(*bp))


@main.command()
@click.argument("img", callback=clk.is_file)
@click.option("-vfile", callback=is_vignette_file, help="vignetting file, rows should be angle(degrees) factor(s) either one column"
                                           "for luminance or 3 for RGB, apply in destination RGB space after lens "
                                           " corrections. system stored files :\n" + "\n".join(vignetting_files))
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def vignette(ctx, img, vfile=None, **kwargs):
    """apply vignetting correction file
    """
    img, header = io.hdr2carray(img, header=True)
    if vfile is None:
        print("Warning, no vignetting applied!", file=sys.stderr)
        imgv = img
    else:
        vg = np.loadtxt(vfile)
        imgv = pl.apply_vignetting_correction(img, vg)
    io.array2hdr(imgv, None, header=[f"VIGNETTING_CORRECTION= {vfile}"] + header)


@main.command()
@click.argument("img", callback=hdr_data_load)
@click.option("-inp", default='rad', help="input image primaries and wp. either 8 values, predefined aliases, or colorformat. Aliases:\n"
                                          "\t rad: (0.640, 0.330, 0.290, 0.600, 0.150, 0.060, 0.3333, 0.3333)\n"
                                          "\t srgb: (0.640,  0.330,  0.300,  0.600,  0.150,  0.060, 0.3127, 0.329)\ncolorformat:\n"
                                          "\t xyz yxy")
@click.option("-outp", default='srgb', help="output image primaries and wp. either 8 values or predefined aliases (see -inp)")
@click.option("-xyzrgb", default=None, help="alternative input as xyz->rgb matrix (overrides -inp)")
@click.option("-oxyzrgb", default=None, help="alternative output as xyz->rgb matrix (overrides -outp)")
@click.option("-rgbrgb", default=None, help="alternative input as rgb->rgb matrix (overrides -inp and -outp)")
@click.option("-vfile", callback=is_vignette_file, help="vignetting file, rows should be angle(degrees) factor(s) either one column"
                                                        "for luminance or 3 for RGB, apply in destination RGB space after lens "
                                                        " corrections. system stored files :\n" + "\n".join(vignetting_files))
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def color(ctx, img, inp='rad', outp='srgb', xyzrgb=None, oxyzrgb=None, rgbrgb=None, **kwargs):
    """apply color primary conversion
    """
    newheader = []
    if inp in ["xyz", "yxy"]:
        inp2xyz = np.eye(3)
    else:
        inp2xyz, _, _ = pl.str_primaries_2_mtx(inp, xyzrgb)

    if outp in ["xyz", "yxy"]:
        rgb2rgb = inp2xyz
        rgb2rgbs = " ".join([f"{i:.08f}" for i in rgb2rgb.ravel()])
        newheader.append(f"RGB2XYZ= {rgb2rgbs}")
    elif rgbrgb:
        rgb2rgb = np.fromstring(rgbrgb).reshape(3, 3)
        rgb2rgbs = " ".join([f"{i:.08f}" for i in rgb2rgb.ravel()])
        newheader.append(f"RGB2RGB= {rgb2rgbs}")
    else:
        orgb2xyz, ps, ws = pl.str_primaries_2_mtx(outp, oxyzrgb)
        rgb2rgb = np.linalg.inv(orgb2xyz) @ inp2xyz
        ps = " ".join([f"{i:.04f}" for i in ps])
        ws = " ".join([f"{i:.04f}" for i in ws])
        ls = " ".join([f"{i:.08f}" for i in orgb2xyz[1]])
        newheader += [f"TargetPrimaries= {ps}",
                      f"TargetWhitePoint= {ws}",
                      f"LuminanceRGB= {ls}"]
        rgb2rgbs = " ".join([f"{i:.08f}" for i in rgb2rgb.ravel()])
        newheader.append(f"RGB2RGB= {rgb2rgbs}")
    if type(img) == str:
        imgd, header = io.hdr2carray(img, header=True)
        if inp == "yxy":
            dy = imgd[0]
            dx = imgd[0]*imgd[1]/imgd[2]
            dz = (1-imgd[1]-imgd[2])*imgd[0]/imgd[2]
            imgd = np.stack((dx, dy, dz))
        rgb = np.einsum('ij,jkl->ikl', rgb2rgb, imgd)
        if outp == "yxy":
            dY = rgb[1]
            sxyz = np.sum(rgb, axis=0)
            dx = rgb[0]/sxyz
            dy = rgb[1]/sxyz
            rgb = np.stack((dY, dx, dy))
        try:
            header += [imagetools.hdr2vm(img).header()]
        except AttributeError:
            pass
        io.array2hdr(rgb, None, header=newheader + header)
    else:
        if inp == "yxy":
            dy = img[:, 0]
            dx = img[:, 0]*img[:, 1]/img[:, 2]
            dz = (1-img[:,1]-img[:,2])*img[:,0]/img[:,2]
            img = np.stack((dx,dy,dz)).T
        rgb = np.einsum('ij,kj->ki', rgb2rgb, img)
        if outp == "yxy":
            dY = rgb[:, 1]
            sxyz = np.sum(rgb, axis=1)
            dx = rgb[:, 0]/sxyz
            dy = rgb[:, 1]/sxyz
            rgb = np.stack((dY, dx, dy)).T
        for r in rgb:
            print(*r)

@main.result_callback()
@click.pass_context
def printconfig(ctx, returnvalue, **kwargs):
    """callback to clean up any temp files"""
    try:
        clk.tmp_clean(ctx)
    except Exception:
        pass


if __name__ == '__main__':
    main()