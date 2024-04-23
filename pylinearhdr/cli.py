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
import datetime
import itertools
import os
import shutil
import sys
from itertools import groupby

import numpy as np
from scipy.optimize import curve_fit
from clasp import click
import clasp.click_ext as clk
from clasp.script_tools import try_mkdir, clean_tmp, pipeline, sglob
from raytools.mapper import ViewMapper

import pylinearhdr
import raytools
from raytools import io, imagetools
from raytools.utility import pool_call
from pylinearhdr import shadowband as sb
from pylinearhdr import pylinearhdr as pl
from pylinearhdr import calibrate as cl


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
@click.option('-profile', '-p', help=f"name of saved profile to load, options: {global_profiles}")
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
                      "different shutter speeds. see the command 'pylinearhdr shutter'"),
    click.option("-scale", default=1.0,
                 help="calibration scale factor (applies to ISO, so do not use -s with linearhdr)"),
    click.option("-cscale", callback=clk.split_float,
                 help="color calibration scale factor (applies via header, same as linearhdr -k)"),
    click.option("-nd", default=0.0, help="additional ND filter (applies to ISO, so do not use -s with linearhdr)"),
    click.option("-saturation", "-saturation-offset", "-s", default=0.01,
                 help="saturation offset, if white is not LinearityUpperMargin, this must be changed"),
    click.option("-range", "-r", default=0.01,
                 help="lower range of single raw exposure"),
    click.option("-rawmultipliers", default=None, callback=clk.split_float,
                 help="override default premultipliers (1, 1, 1)"),
    click.option("--verbose/--no-verbose", default=False,
                 help="passed to linearhdr"),
    click.option("--interpfirst/--interpsecond", default=False,
                 help="default (false) interpolate after merge (uses DHT), if true,"
                      "interpolate with rawconvert (uses -interpq unless --half) or"),
    click.option("-interpq", default="DHT",
                 type=click.Choice(["linear", "VNG", "PPG", "AHD", "DCB", "DHT", "AAHD"], case_sensitive=False),
                 help="demosaicing algorithm, only used if --interpfirst, otherwise DHT"),
    click.option("--half/--no-half", default=False,
                 help="use half-scale output from rawconvert, disables --interpsecond and --rawgrid"),
    click.option("--rawgrid/--no-rawgrid", default=False,
                 help="do not interpolate raw channels. forces -colorspace to 'raw' and ignores --interpfirst"),
    click.argument("imgs", callback=clk.are_files),
    click.option("--fisheye/--no-fisheye", default=False,
                 help="apply fisheye_corr to 180 degree image (must be properly cropped and equiangular). "
                      "requires pcomb and RAYPATH"),
    click.option("--overwrite/--no-overwrite", default=False,
                 help="run rawconvert even if output file exists"),
    click.option("-header-line", '-hl', multiple=True,
                 help="lines to append to HDR header, e.g. LOCATION= 46.522833,6.580500"),
    click.option("-executable", default="linearhdr",
                 help="linearhdr executable"),
    click.option("--correct/--no-correct", default=True,
                 help="apply correction to nominal aperture and shutter speed values, use with linearhdr --exact"),
    click.option("--listonly/--no-listonly", default=False,
                 help="skip execution and just print metadata"),
    click.option("-hdropts", default="", help="additional options to linearhdr (with callhdr, overrides -r -s)"),
    click.option("-rawcopts", default="", help="additional options to rawconvert, for example: ' -C 0.999 9.999'"),
    click.option("-crop", callback=clk.split_int,
                 help="crop tiff (left upper W H)"),
]


def makelist_run(ctx, imgs, overwrite=False, correct=False, listonly=False, scale=1.0, nd=0.0, saturation=0.01, range=0.01,
                 crop=None, badpixels=None, callhdr=False, rawcopts="", hdropts="", fo=None, fisheye=False, xyzcam=None, cscale=None, shutterc=None,
                 black="AverageBlackLevel", white="AverageBlackLevel", colorspace='rad', clean=False, vfile=None, verbose=False, rawgrid=False,
                 interpfirst=False, header_line=None, half=False, rawmultipliers=None, interpq="DHT", executable="linearhdr", **kwargs):
    """make list routine, use to generate input to linearhdr"""
    if executable == "linearhdr":
        needs_full = "-wf" in hdropts or "-f" in hdropts
        needs_ext = "--we" in hdropts or ("-w" in hdropts and not needs_full)
        if needs_ext:
            executable = "linearhdr_full"
        elif not needs_full:
            executable = "linearhdr_slim"
    if half:
        interpfirst = True
        rawgrid = False
    if header_line is None:
        header_line = []
    else:
        header_line = list(header_line)
    multipliers = np.array([1, 1, 1])
    if verbose:
        hdropts += " --verbose"
    if rawgrid:
        colorspace = 'raw'
        if executable != "linearhdr_slim":
            hdropts += " -B"
    elif not interpfirst:
        hdropts += " -D"
    if executable == "linearhdr_slim" and interpfirst and not rawgrid:
        hdropts += " -B"
    outf = sys.stdout
    outfn = "<makelist.txt>"
    if listonly:
        overwrite = False
    elif callhdr:
        outfn = clean_tmp(ctx).rsplit("/", 1)[-1]
        outf = open(outfn, 'w')
        if not correct:
            hdropts += " --nominal"
    if xyzcam is None:
        xyzcam = pl.get_xyz_cam(imgs[0])
    xyzcam = np.asarray(xyzcam).reshape(3, 3)
    if rawmultipliers is not None:
        multipliers = np.array(rawmultipliers)[0:3]
    rawcopts += ' -r ' + " ".join([f"{i:.06f}" for i in multipliers]) + f" {multipliers[1]:.06f}"
    xyzcam = xyzcam * multipliers[:, None]
    rawconvertcom = pl.rawconvert_opts(imgs[0], crop=crop, bad_pixels=badpixels, rawgrid=rawgrid or (not interpfirst),
                                       black=black, white=white, rawcopts=rawcopts, half=half, interpq=interpq)
    tiffs = pool_call(pl.get_raw_frame, imgs, correct=correct, overwrite=overwrite, rawconvertcom=rawconvertcom, fo=fo,
                     shutterc=shutterc, listonly=listonly, expandarg=False, pbar=False)
    cam_rgb, header = pl.cam_color_mtx(xyzcam, colorspace, cscale=cscale)
    print(f"# {rawconvertcom}", file=outf)
    print(f"# pylinearhdr_VERSION= {pylinearhdr.__version__}", file=outf)
    if shutterc is not None:
        print(f"# SHUTTER_CORRECTION= {shutterc}", file=outf)
    if fo is not None:
        fostr = " ".join([f"{i[0]},{i[1]}" for i in fo])
        print(f"# APERTURE_OVERRIDES= {fostr}", file=outf)
    print(f"# XYZCAM= " + " ".join([f"{i:.08f}" for i in xyzcam.ravel()]), file=outf)
    print(f"# CAM_PREMULTIPLIERS= " + " ".join([f"{i:.08f}" for i in multipliers.ravel()]), file=outf)
    print("# CAPDATE= {}".format(datetime.datetime.now().strftime("%Y:%m:%d %H:%M:%S")), file=outf)
    rawstring = " ".join(imgs[0:1] + [i.rsplit("/", 1)[-1] for i in imgs[1:]])
    print(f"# RAW_IMAGES= {rawstring}", file=outf)
    for h in header:
        print(h, file=outf)
    for h in header_line:
        print(f"# {h}", file=outf)
    pl.report(tiffs, listonly, scale=scale * 10**nd, sat_w=1-saturation, sat_b=range, outf=outf)
    command = [f"{executable} -r {range} -o {saturation} {hdropts} {outfn}"]
    if fisheye:
        command += ["raytools solid2ang -"]
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
@clk.shared_decs(shared_run_opts)
@click.option("--callhdr/--no-callhdr", default=False, help="directly call linearhdr")
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def makelist(ctx, imgs, **kwargs):
    """make list routine, use to generate input to linearhdr"""
    makelist_run(ctx, imgs, **kwargs)


@main.command()
@clk.shared_decs(shared_run_opts)
@click.option("--clean/--no-clean", default=True, help="delete tiff files after linearhdr")
@click.option("-vfile", callback=is_vignette_file, help="vignetting file, rows should be angle(degrees) factor(s) either one column"
                                                        "for luminance or 3 for RGB, apply in destination RGB space after lens "
                                                        " corrections. system stored files :\n" + "\n".join(vignetting_files))
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def run(ctx, imgs, **kwargs):
    """make list routine, use to generate input to linearhdr"""
    makelist_run(ctx, imgs, callhdr=True, **kwargs)


@main.command()
@click.argument("imgh", callback=clk.are_files)
@click.argument("imgv", callback=clk.are_files)
@click.argument("imgn", callback=clk.are_files)
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
@click.option("-envmap",
              help="additionally save an hdr file suitable for use as an environment map, source information is in the header")
@click.option("-margin", default=20,
              help="pixel margin beyond 180 degrees in all input images, will crop and reproject final output")
@click.option("-check",
              help="write out masks and other intermediate data (give as prefix)")
@click.option("--align/--no-align", default=True,
              help="align H/V images and crop by margin (aligns H->V)")
@click.option("--fisheye/--no-fisheye", default=True,
              help="after alignment/crop reproject all inputs to angular fisheye")
@click.option("-rawext", default="CR2",
              help="when inputs are directories of raw files, the file extension to look for (do not include .)")
@click.option("-bandcfg", default="D70-1SB",
              help="config file for band image hdr merge (using pylinearhdr run) if argument includes a suffix, "
                   "load as config, otherwise as profile")
@click.option("-ndcfg", default="D70-1ND3",
              help="config file for band image hdr merge (using pylinearhdr run) if argument includes a suffix, "
                   "load as config, otherwise as profile")
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def shadowband(ctx, imgh, imgv, imgn, outf="blended.hdr", roh=0.0, rov=0.0, sfov=2.0, srcsize=6.7967e-05, bw=2.0,
               envmap=None, sunloc=None, check=None, margin=20, align=True, fisheye=True, rawext="CR2",
               bandcfg="D70-1SB", ndcfg="D70-1ND3", **kwargs):
    """merge set of three images with horizontal, vertical and no shadow band all images are assumed
    to be equisolid fisheye with an extra 'margin' around 180
    to allow for H/V image alignment (unless align/fisheye are false.

    imgh: hdr with horizontal aligned shadowband
    imgv: hdr with veritcal aligned shadowband
    imgn: hdr with no shadowband

    to recover environment map (outputs into southern hemisphere:

        getinfo testmap.hdr | grep SOURCE= | sed -E 's/[[:blank:]]*[A-Z]+SOURCE=[[:blank:]]*//g' | xform > test.rad

    or to rotate to west (for example):

        getinfo testmap.hdr | grep SOURCE= | sed -E 's/[[:blank:]]*[A-Z]+SOURCE=[[:blank:]]*//g' | xform -rz 90 > test.rad

    """
    sbobt = (f"SHADOWBAND= roh:{roh:.03f} rov:{rov:.03f} sfov:{sfov:.03f} srcsize:{srcsize:.04f} bw:{bw:.04f} "
             f"align:{align} fisheye:{fisheye} sunloc:{sunloc} margin:{margin}")

    def _prerun0(outf0, imgs0, cfg):
        if "." in cfg:
            runcom = f"pylinearhdr -c {cfg}"
        else:
            runcom = f"pylinearhdr -p {cfg}"
        runcom += f" run '{imgs0}'"
        pipeline([runcom], outf0, writemode='wb')
        return io.hdr2carray(outf0, header=True)

    def _prerun(imgd, cfg):
        imgd = imgd.rstrip("/")
        outf = imgd + ".hdr"
        imgs = f"{imgd}/*.{rawext}"
        return _prerun0(outf, imgs, cfg)

    try:
        if len(imgh) > 1:
            hdata, hh = _prerun0(outf.rsplit(".", 1)[0] + "_H.hdr", " ".join(imgh), bandcfg)
        else:
            hdata, hh = io.hdr2carray(imgh[0], header=True)
    except ValueError:
        hdata, hh = _prerun(imgh[0], bandcfg)
    try:
        if len(imgh) > 1:
            vdata, vv = _prerun0(outf.rsplit(".", 1)[0] + "_V.hdr", " ".join(imgv), bandcfg)
        else:
            vdata, hv = io.hdr2carray(imgv[0], header=True)
    except ValueError:
        vdata, hv = _prerun(imgv[0], bandcfg)
    try:
        if len(imgn) > 1:
            sdata, hs = _prerun0(outf.rsplit(".", 1)[0] + "_ND.hdr", " ".join(imgn), ndcfg)
        else:
            sdata, hs = io.hdr2carray(imgn[0], header=True)
    except ValueError:
        sdata, hs = _prerun(imgn[0], ndcfg)
    if align:
        if margin < 1:
            t = slice(None, None)
            margin = 0
        else:
            t = slice(margin, -margin)
        xo, yo = sb.align_images(hdata[0, t, t], vdata[0, t, t])
        sdata = sdata[:, t, t]
        vdata = vdata[:, t, t]
        if np.sum(np.abs((xo, yo))) > 0:
            xt = slice(max(0, margin + xo), min(hdata.shape[1], xo + hdata.shape[1] - margin))
            yt = slice(max(0, margin + yo), min(hdata.shape[2], yo + hdata.shape[2] - margin))
            xt2 = slice(0, xt.stop-xt.start)
            yt2 = slice(0, yt.stop-yt.start)
            hdata2 = np.zeros_like(sdata)
            hdata2[:, xt2, yt2] = hdata[:, xt, yt]
            hdata = hdata2
            hh.append(f"\tSHADOWBAND_IMAGE_ALIGN= {xo} {yo}")
        else:
            hdata = hdata[:, t, t]
        margin = 0
    if fisheye:
        if margin > 0:
            t = slice(margin, -margin)
            sdata = sdata[:, t, t]
            vdata = vdata[:, t, t]
            hdata = hdata[:, t, t]
        hdata, vdata, sdata = pool_call(imagetools.array_solid2ang, [hdata, vdata, sdata], expandarg=False, returnvm=False, pbar=False)
    blended, skyonly, source = sb.shadowband(hdata, vdata, sdata, roh=roh, rov=rov, sfov=sfov, srcsize=srcsize, bw=bw,
                                             envmap=envmap, sunloc=sunloc, check=check)
    vm = ViewMapper((0.0, -1.0, 0.0), viewangle=180)
    csp = "CENTRAL_SOURCE_PIXEL= {} {}".format(*vm.ray2pixel(source[0:3], hdata.shape[-1], False).ravel())
    header = [sbobt, csp] + hh + hv + hs + [vm.header()]

    io.carray2hdr(blended, outf, header, clean=True)
    if skyonly is not None:
        srcsize = (source[3] / np.pi) ** .5 * 360/np.pi
        header.append("SOLARSOURCE= " + "void light sun 0 0 3 {} {} {} sun source solar 0 0 4 {} {} {} {}".format(*source[4], *source[0:3], srcsize))
        header.append("SKYSOURCE= " + "void colorpict imgfunc 9 red green blue {} fisheye.cal fish_u fish_v -rz 180 0 0 imgfunc glow imgglow 0 0 4 1 1 1 0 imgglow source sky 0 0 4 0 -1 0 180".format(envmap))
        io.carray2hdr(skyonly, envmap, header, clean=True)


@main.command()
@click.argument("imgs", callback=clk.are_files)
@click.option("-out", default="img", help="directory base name, becomes: img_XXX, use 'common'"
                                          " to auto-generate from common parts (seperated by _) or"
                                          "lcommon to shorten common parts by across group differences")
@click.option("-starti", default=0, help="start index")
@click.option("-count", type=int, help="if given, overrides auto detection and groups into folders of count size"
                                       " based on order given (uses sorted wildcard expansion.")
@click.option("--ascend/--descend", default=True,
              help="shot order exposure time. ascending=shortest->longest")
@click.option("--r/--no-r", "--rename/--no-rename", default=False,
              help="also call rename")
@click.option("--preview/--no-preview", default=False,
              help="print outcome without moving")
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def sort(ctx, imgs, out="img", starti=0, ascend=True, preview=False, r=False, count=None, **kwargs):
    """sort sequence of images into folders based on changes in exposure time"""
    if r:
        rename.callback(imgs, **kwargs)
        imgs = ctx.obj['imgs']
    if count is None:
        infos = pool_call(pl.get_raw_frame, imgs, listonly=True, correct=True, expandarg=False, pbar=False)
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
            expt = info[1]
            groups[-1].append(info[0])
    else:
        def batched(iterable, n):
            it = iter(iterable)
            while batch := tuple(itertools.islice(it, n)):
                yield batch

        groups = list(batched(imgs, count))
    outns = []
    for j, group in enumerate(groups):
        if out[-6:] == 'common':
            pieces = [g.rsplit(".", 1)[0].split("_") for g in group]
            pcount = min(len(g) for g in pieces)
            pieces = [g.rsplit(".", 1)[0].split("_", pcount) for g in group]
            match = [pieces[0][i] for i in range(pcount) if len(set([p[i] for p in pieces])) == 1]
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
            print(f"\n{outn}_{j+starti:03d}/")
        else:
            try_mkdir(f"{outn}_{j+starti:03d}")
        for g in group:
            if preview:
                print(f"\t{g}")
            else:
                shutil.move(g, f"{outn}_{j+starti:03d}")
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
    """rename raw files based on ISO_APERTURE_SHUTTER"""
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
    io.array2hdr(imgv, None, header=[f"VIGNETTING_CORRECTION= {vfile}"] + header, clean=True)


@main.command()
@click.option("-seq", callback=clk.are_files_iter, multiple=True,
              help="sequence of images with different shutter speeds that all have same region (see crop)"
                   " in range. give multiple sequences with different lighting or crop regions to cover a wider "
                   "range of shutter speeds")
@click.option("-crop", multiple=True,
              help="crop region for -seq (give multiple in same order). give as: left upper width height, note this is "
                   "not the same as pcompos!). Use pfsview, photoshop or other raw image viewer to determine. region"
                   "should be consistently lit, color neutral and in range for all images in -seq. will duplicate last"
                   "-crop when there are fewer than -seq. if not given at all, uses full image (not recommended)")
@click.option("-runopts", default="", help="passed to pylinearhdr run")
@click.option("-channel", default="g", help="which color channel to use")
@click.option("-dataout",
              help="if given, save full data to this file")
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def shutter(ctx, seq=None, crop=None, runopts="", dataout=None, channel='g', **kwargs):
    """do relative shutter speed calibration

    When multiple sequences are given, they are sorted by the
    slowest in range shutter speed in each, then the relative
    values are aligned via the matching shutter speeds between
    sequences. This means that adjacent sequences must contain
    at least one overlapping shutterspeed, and that having
    multiple overlaps will yield a more robust alignment.

    Duplicate shutter values (within and across sequences) are
    averaged together prior to fitting, so that all shutter
    speeds included are equally weighted in the fit. Only images
    with 100% of pixels in range are included, so it is best to
    carefully crop the sequences to highly homogenous regions
    to maximize the number of in range frames.
    """
    if len(seq) == 0:
        print("No images given", file=sys.stderr)
        raise click.Abort
    if len(crop) == 0:
        crop = [None]
    results = []
    alldata = []
    # for each sequence calculate relative factors
    for i, sq in enumerate(seq):
        j = min(len(crop)-1, i)
        if crop[j] is None:
            croparg = ""
        else:
            croparg = f"-crop '{crop[j]}'"
        result = pool_call(cl.average_channel, sq, croparg, runopts=runopts, expandarg=False, channel=channel)
        # filter out of range
        result = np.array([[j[0], j[2]] for j in result if j[3] > 0.99])
        print(f"sequence {i}: {len(result)} out of {len(sq)} frames in range", file=sys.stderr)
        # sort by shutter speed
        result = result[np.argsort(result[:, 0])]
        alldata.append(result)
        # group and average by shutter speed
        result = [[x, np.mean([yi[1] for yi in y])]
                  for x, y in groupby(result, key=lambda x: f"{x[0]:.04f}")]
        results.append(result)
    # sort sequences by slowest shutter speed
    rsort = np.argsort([float(i[0][0]) for i in results])
    results = [results[i] for i in rsort]
    alldata = [alldata[i] for i in rsort]
    # these are the final lists
    shutters = [i[0] for i in results[0]]
    vals = [i[1] for i in results[0]]
    alldata_s = [alldata[0]]
    # scale subsequent sequences to first basis
    for i in range(1, len(results)):
        sf = []
        # find overlap to determine average scale factor
        for s, v in results[i]:
            try:
                si = shutters.index(s)
            except ValueError:
                pass
            else:
                sf.append(vals[si]/v)
        if len(sf) == 0:
            raise ValueError("sequences do not connect, at least one shutter "
                             "speed in each sequence must appear in a"
                             "previous (slower) sequence")
        sf = np.average(sf)
        shutters += [r[0] for r in results[i]]
        vals += [r[1] * sf for r in results[i]]
        s_data = np.copy(alldata[i])
        s_data[:, 1] *= sf
        alldata_s.append(s_data)
    alldata = np.vstack(alldata_s)
    alldata = alldata[np.argsort(alldata[:, 0])]
    results = np.stack(([float(i) for i in shutters], vals)).T
    # sort by shutter speed
    results_p = results[np.argsort(results[:, 0])]
    # group and average by shutter speed
    results = np.array([[float(x), np.mean([yi[1] for yi in y])]
                        for x, y in groupby(results_p, key=lambda x: f"{x[0]:.04f}")])
    bm = results[0, 1]
    results[:, 1] = results[:, 1] / bm
    coef, _ = curve_fit(lambda t,a,b: b*np.exp(a*t),  results[:, 0],  1/results[:, 1], p0=(-1e-8, 1))
    if dataout is not None:
        f = open(dataout, 'w')
        print("shutter\tfit\tavg\tsamples", file=f)
        for r, (x, y) in zip(results, groupby(alldata, key=lambda x: f"{x[0]:.04f}")):
            samples = [bm/yi[1] for yi in y]
            avg = 1/r[1]
            efit = coef[1]*np.exp(coef[0]*float(x))
            print(f"{x}\t{efit}\t{avg}\t" + "\t".join([f"{yi}" for yi in samples]), file=f)
        f.close()
    print("# calculated average relative scaling factors")
    for r in results:
        print(*r)
    print("# coefficients for exponential fit correction factor: B*e^(A*x)")
    print("# give 'A' as argument -shutterc to 'pylinearhdr run' or linearhdr. "
          "The 'B' coefficient can be ignored except for plotting/checking as "
          "this is only a relative correction.", file=sys.stderr)
    print(f"A: {coef[0]}", file=sys.stderr)
    print(f"B: {coef[1]}", file=sys.stderr)
    print(f"for configuration file:", file=sys.stderr)
    print(f"shutterc = {coef[0]}", file=sys.stderr)


@main.command()
@click.option("-seq", callback=clk.are_files_iter, multiple=True,
              help="sequence of images. give once for each aperture series.")
@click.option("-shutterc", "-sc", default=float,
              help="exponential shutter speed correction (results will not be valid unless given. give as 0.0"
                   " if camera has perfect shutter speed calibration.")
@click.option("-crop",
              help="crop region for -seq (give multiple in same order). give as: left upper width height, note this is "
                   "not the same as pcompos!). Use pfsview, photoshop or other raw image viewer to determine. region"
                   "should be consistently lit, color neutral and in range for all images in -seq. will duplicate last"
                   "-crop when there are fewer than -seq. if not given at all, uses full image (not recommended)")
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def aperture(ctx, seq=None, crop=None, shutterc=None, **kwargs):
    """do relative aperture calibration

    perform only after shutter correction. first sequence is taken as the relative reference (ie is assumed to be the
    "correct" aperture). As each sequence is assumed to be of the same exact lighting condition, crop can only be given
    once and applies to all sequences.
    """
    if len(seq) < 2:
        print("cannot calculate aperture correction without at least two sequences", file=sys.stderr)
        raise click.Abort
    if shutterc is None:
        print("-shutterc is a required option. use pylinearhdr shutter to compute", file=sys.stderr)
        raise click.Abort
    results = []
    if crop is None:
        croparg = ""
    else:
        croparg = f"-crop '{crop}'"
    # for each sequence calculate average value
    for sq in seq:
        result = cl.average_channel(" ".join(sq), croparg)
        if result[3] < 1:
            print(f"Warning! sequence '{sq}' does not yield an in range HDR image, {1-result[3]} is out of range",
                  file=sys.stderr)
        results.append([pl.info_from_exif(sq[0], False)[1], result[1], result[2]])
    refval = results[0][2]
    corrections = []
    print(f"Aperture corrections based on nominal aperture F-{results[0][0]} (exact: {results[0][1]})")
    print(f"Nominal:   Exact: Corrected:")
    for r in results[1:]:
        cf = np.sqrt(r[1]*r[1]*refval/r[2])
        print(f"{r[0]:>8.03f} {r[1]:>8.03f} {cf:>10.03f}")
        corrections.append(f"{r[0]},{round(cf,6)}")

    print(f"aperture correction string (option -fo)")
    print(f"fo = {' '.join(corrections)}")


@main.command()
@click.argument("reference", callback=clk.is_file)
@click.argument("test", callback=clk.is_file)
@click.option("-rc", callback=clk.is_file,
              help='file of pixel locations corresponding to target areas in the reference. each row has four numbers '
                   'x y w h, where x, y are the lower left corner. if not given script will interactively'
                   'generate. The first point should always be the "white-point"')
@click.option("-tc", callback=clk.is_file,
              help='file of pixel locations corresponding to target areas in the test. each row has four numbers '
                   'x y w h, where x, y are the lower left corner. if not given script will interactively'
                   'generate. The first point should always be the "white-point"')
@click.option("-refcol", default='rad',
              help="reference image primaries and wp. either 8 values, predefined aliases, or colorformat. Aliases:\n"
                   "\t rad: (0.640, 0.330, 0.290, 0.600, 0.150, 0.060, 0.3333, 0.3333)\n"
                   "\t srgb: (0.640,  0.330,  0.300,  0.600,  0.150,  0.060, 0.3127, 0.329)\ncolorformat:\n"
                   "\t xyz")
@click.option("-xyzcam", default=None, callback=clk.split_float,
              help="optional xyzcam for input to override/provide XYZCAM= header information. If not given, and not found,"
                   " refcol is assumed.")
@click.option("-refdataout", default=None,
              help="if given, save reference data to this path in XYZ")
@click.option("-testdataout", default=None,
              help="if given, save test data to this path (in input colorspace (raw))")
@click.option("--alternate/--no-alternate", default=False,
              help='if true and neither -rc or -tc are provided, opens both images simultaneously and expects '
                   'alternating inputs (cell 1 ref, cell 1 test, cell 2 ref, cell 2 test, etc.')
@click.option("-minimizer", type=click.Choice(['luv', 'lab'], case_sensitive=False), default='luv')
@click.option("--verbose/--no-verbose", default=False, help="print detailed results")
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def calibrate(ctx, reference, test, rc=None, tc=None, refcol='rad', xyzcam=None, alternate=False, refdataout=None, testdataout=None,
              verbose=False, minimizer='luv', **kwargs):
    """run color luminance calibration using reference data (or image) and test hdr
    (shutter and aperture corrected, but raw color)

    Arguments
    ---------

    reference: reference color HDR or data. if data, should be tab-seperated, with rows corresponding to tc. can
        can be either one column (luminance) or three column. If HDR or 3-column, colorspace must be
        given with -refcol. Note that the RGB channel scaling minimization will be done in the reference colorspace.
    test: merged hdr image of the calibration scene, should be self-calibrated for shutter and aperture, but output
        in raw (-colorspace raw in pylinearhdr run) if rgb=True or reference data is 3-channel. otherwise
        output in target color space (RGB or sRGB). If test was made with another program and does not have "XYZCAM"
        header line, then test is assumed to be in the same color space as reference (given by -refcol) unless xyzcam is given.
        Note that even if reference is greyscale, this will still apply (just not to the reference).

    """
    def _is_hdr(imgf):
        f = open(imgf, 'rb')
        ishdr = f.read(10) == b'#?RADIANCE'
        f.close()
        return ishdr

    refimg = _is_hdr(reference)
    testimg = _is_hdr(test)
    if not (refimg and testimg) and alternate:
        print("Warning reference or test is not an HDR, so setting alternate to False", file=sys.stderr)
        alternate = False
    if xyzcam:
        xyzcam = np.array(xyzcam).reshape(3,3)
    result, A, B = cl.calibrate(reference, test, rc, tc, alternate, refimg, refcol, xyz_cam=xyzcam, testimg=testimg, minimizer=minimizer)
    if refdataout:
        np.savetxt(refdataout, A.T)
    if testdataout:
        np.savetxt(testdataout, B.T)
    if verbose:
        k, v = list(result.items())[-1]
        print("final optimization results:")
        if minimizer == 'luv':
            print(" #        Dlum         Duv")
            for c, (i, j) in enumerate(v['duv_full'].T):
                print(f"{c:02d}  {i: >10.02%}  {j: >10.04f}")
        else:
            print(" #        De")
            for c, i in enumerate(v['dlab_full'].T):
                print(f"{c:02d}  {i: >10.02}")
    print("xyzcam matrix:")
    for k, v in result.items():
        print(f"{v['label']+':':<40} " + " ".join([f"{i:.06f}" for i in v['xyzcam'].ravel()]))


@main.command()
@click.argument("test", callback=clk.is_file)
@click.option("-tc", callback=clk.is_file,
              help='file of pixel locations corresponding to target areas in the test. each row has four numbers '
                   'x y w h, where x, y are the lower left corner. if not given script will interactively'
                   'generate.')
@click.option("-outf",
              help='data file to write. if not given, defaults to input + ".txt"')
@click.option("-checkimg",
              help='if given generate a masked image showing areas (values in overlapping areas get doubled)')
@click.option("-scale", default=179.0,
              help='scale data by this value (radiance compatible 179)')
@click.option("--stdout/--no-stdout", default=False,
              help='overrides outf, print to stdout')
@click.option("-zeromode", default=0, type=click.IntRange(0, 2),
              help='mode 0: treat zeros as normal part of data, include in average, '
                   'mode 1: take average on non-zero data, '
                   'mode 2: if area includes zeros, return zero for whole area, ')
@click.option("--mean/--no-mean", default=True,
              help='calculate mean (if true always first 3/4 columns')
@click.option("--stdev/--no-stdev", default=False,
              help='calculate stddev (if true always last 3/4 columns')
@click.option("--lum/--no-lum", default=False,
              help='include luminance column before RGB. if LuminanceRGB is not in header, assumes radiance colorspace')
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def getimgdata(ctx, test, tc=None, outf=None, stdout=False, zeromode=0, lum=False, checkimg=None, mean=True, stdev=False, scale=179.0, **kwargs):
    """get average squares of data from hdr"""
    if outf is None:
        outf = test + ".txt"
    tcells = cl.load_test_cells(test, tc)
    testd = cl.load_data(test, tcells, zero=zeromode, lum=lum, checkimg=checkimg, mean=mean, stdev=stdev, scale=scale)
    if mean or stdev:
        if stdout:
            for d in testd.T:
                print("\t".join([f"{i:.08g}" for i in d]))
        else:
            np.savetxt(outf, testd.T)
    elif not checkimg:
        print("Note: no data requested!", file=sys.stderr)


@main.command()
@click.argument("img", callback=hdr_data_load)
@click.option("-inp", default='rad', help="input image primaries and wp. either 8 values, predefined aliases, or colorformat. Aliases:\n"
                                          "\t rad: (0.640, 0.330, 0.290, 0.600, 0.150, 0.060, 0.3333, 0.3333)\n"
                                          "\t srgb: (0.640,  0.330,  0.300,  0.600,  0.150,  0.060, 0.3127, 0.329)\ncolorformat:\n"
                                          "\t xyz yxy yuv")
@click.option("-outp", default='srgb', help="output image primaries and wp. either 8 values or predefined aliases (see -inp)")
@click.option("-lab", default=None, callback=clk.split_float, help="override -outp and output La*b* color, give XYZ whitepoint")
@click.option("-xyzrgb", default=None, help="alternative input as xyz->rgb matrix (overrides -inp)")
@click.option("-oxyzrgb", default=None, help="alternative output as xyz->rgb matrix (overrides -outp)")
@click.option("-rgbrgb", default=None, help="alternative input as rgb->rgb matrix (overrides -inp and -outp)")
@click.option("--verbose/--no-verbose", default=False, help="print color transforn matrices to stderr")
@click.option("-maxv", type=float, help="cap max val, only applies to RGB or XYZ output")
@click.option("-minv", type=float, help="cap minv val, only applies to RGB or XYZ output")
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def color(ctx, img, inp='rad', outp='srgb', xyzrgb=None, oxyzrgb=None, rgbrgb=None, lab=None, verbose=True, maxv=None, minv=None, **kwargs):
    """apply color primary conversion
    """
    if lab is not None:
        lab = np.array(lab).ravel()[0:3]
        outp = 'xyz'
    if type(img) == str:
        imgd, header = io.hdr2carray(img, header=True)
        try:
            header += [imagetools.hdr2vm(img).header()]
        except AttributeError:
            pass
        rgb, outheader = pl.color_convert_img(imgd, header, inp, outp, xyzrgb=xyzrgb, rgbrgb=rgbrgb, oxyzrgb=oxyzrgb, verbose=verbose)
        if lab is not None:
            shape = rgb.shape
            rgb = cl.xyz_2_lab(rgb.reshape(3, -1), lab)
            rgb = rgb.reshape(shape)
        elif outp not in ["yxy", "yuv"]:
            if maxv is not None:
                rgb = np.minimum(maxv, rgb)
            if minv is not None:
                rgb = np.maximum(minv, rgb)
        io.array2hdr(rgb, None, header=outheader, clean=True)
    else:
        rgb2rgb, _ = pl.prep_color_transform(inp, outp, xyzrgb=xyzrgb, rgbrgb=rgbrgb, oxyzrgb=oxyzrgb, verbose=verbose)
        if inp == "yxy":
            dy = img[:, 0]
            dx = img[:, 0]*img[:, 1]/img[:, 2]
            dz = (1-img[:,1]-img[:,2])*img[:,0]/img[:,2]
            img = np.stack((dx,dy,dz)).T
        elif inp == "yuv":
            dy = img[:, 0]
            d = 9 * dy / img[:, 2]
            dx = img[:, 1] * 9 * dy / (4 * img[:, 2])
            dz = (d - dx - 15 * dy) / 3
            img = np.stack((dx, dy, dz)).T
        rgb = np.einsum('ij,kj->ki', rgb2rgb, img)
        if lab is not None:
            rgb = cl.xyz_2_lab(rgb.T, lab).T
        elif outp == "yxy":
            dY = rgb[:, 1]
            sxyz = np.sum(rgb, axis=1)
            sxyz[sxyz==0] = 1
            dx = rgb[:, 0]/sxyz
            dy = rgb[:, 1]/sxyz
            rgb = np.stack((dY, dx, dy)).T
        elif outp == 'yuv':
            d = rgb[:, 0] + 15 * rgb[:, 1] + 3 * rgb[:, 2]
            d[d == 0] = 1
            u = 4 * rgb[:, 0] / d
            v = 9 * rgb[:, 1] / d
            rgb = np.stack((rgb[:, 1], u, v)).T
        else:
            if maxv is not None:
                rgb = np.minimum(maxv, rgb)
            if minv is not None:
                rgb = np.maximum(minv, rgb)
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