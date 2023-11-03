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
import os
import shutil
import sys

import numpy as np
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
    click.option("-scale", default=1.0,
                 help="calibration scale factor (applies to ISO, so do not use -s with linearhdr)"),
    click.option("-cscale", callback=clk.split_float,
                 help="color calibration scale factor (applies via header, same as linearhdr -k)"),
    click.option("-nd", default=0.0, help="additional ND filter (applies to ISO, so do not use -s with linearhdr)"),
    click.option("-saturation", "-saturation-offset", "-s", default=0.01,
                 help="saturation offset, if white is not LinearityUpperMargin, this must be changed"),
    click.option("-range", "-r", default=0.01,
                 help="lower range of single raw exposure"),
    click.option("--verbose/--no-verbose", default=False,
                 help="passed to linearhdr"),
    click.option("--interpfirst/--interpsecond", default=True,
                 help="interpolate with rawconvert (uses linear) or interpolate after merge (uses DHT)"),
    click.option("--rawgrid/--no-rawgrid", default=False,
                 help="do not interpolate raw channels. forces -colorspace to 'raw' and ignores --interpfirst"),
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
    click.option("--premult/--no-premult", default=False,
                 help="normalize xyzcam matrix and premultiply raw channels to compensate. This will correct any"
                      " colorshift near out of bounds values at the expense of offsetting the valid raw camera range. "
                      "If the HDR sequence is safely in range, set this to false, otherwise true may yield better results"),
    click.option("-hdropts", default="", help="additional options to linearhdr (with callhdr, overrides -r -s)"),
    click.option("-crop", callback=clk.split_int,
                 help="crop tiff (left upper W H)"),
]


def makelist_run(ctx, imgs, shell=False, overwrite=False, correct=False, listonly=False, scale=1.0, nd=0.0, saturation=0.01, range=0.01,
                 crop=None, badpixels=None, callhdr=False, hdropts="", fo=None, fisheye=False, xyzcam=None, cscale=None, shutterc=None,
                 black="AverageBlackLevel", white="AverageBlackLevel", colorspace='rad', clean=False, vfile=None, verbose=False, rawgrid=False,
                 interpfirst=True, premult=False, header_line=None, **kwargs):
    """make list routine, use to generate input to linearhdr"""
    if header_line is None:
        header_line = []
    else:
        header_line = list(header_line)
    multipliers = np.array([1, 1, 1])
    if verbose:
        hdropts += " --verbose"
    if rawgrid:
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
    if xyzcam is None:
        xyzcam = pl.get_xyz_cam(imgs[0])
    xyzcam = np.asarray(xyzcam).reshape(3, 3)
    if premult:
        rowsum = np.sum(xyzcam, axis=1)
        multipliers = np.max(rowsum)/rowsum
    rawcopts = '-r ' + " ".join([f"{i:.06f}" for i in multipliers]) + f" {multipliers[1]:.06f}"
    xyzcam = xyzcam * multipliers[:, None]
    rawconvertcom = pl.rawconvert_opts(imgs[0], crop=crop, bad_pixels=badpixels, rawgrid=rawgrid or (not interpfirst),
                                       black=black, white=white, rawcopts=rawcopts)
    tiffs = pool_call(pl.get_raw_frame, imgs, correct=correct, overwrite=overwrite, rawconvertcom=rawconvertcom, fo=fo,
                     shutterc=shutterc, listonly=listonly, expandarg=False, pbar=False)
    cam_rgb, header = pl.cam_color_mtx(xyzcam, colorspace, cscale=cscale)
    # print(f"# pylinearhdr " + " ".join(sys.argv[1:]), file=outf)
    print(f"# pylinearhdr_VERSION= {pylinearhdr.__version__}", file=outf)
    print(f"# {rawconvertcom}", file=outf)
    print(f"# XYZCAM= " + " ".join([f"{i:.08f}" for i in xyzcam.ravel()]), file=outf)
    print(f"# CAM_PREMULTIPLIERS= " + " ".join([f"{i:.08f}" for i in multipliers.ravel()]), file=outf)
    print("# CAPDATE= {}".format(datetime.datetime.now().strftime("%Y:%m:%d %H:%M:%S")), file=outf)
    rawstring = " ".join(imgs[0:1] + [i.rsplit("/", 1)[-1] for i in imgs[1:]])
    print(f"# RAW_IMAGES= {rawstring}", file=outf)
    for h in header:
        print(h, file=outf)
    for h in header_line:
        print(f"# {h}", file=outf)
    pl.report(tiffs, shell, listonly, scale=scale * 10**nd, sat_w=1-saturation, sat_b=range, outf=outf)

    command = [f"linearhdr -r {range} -o {saturation} {hdropts} {outfn}"]
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
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def shadowband(ctx, imgh, imgv, imgn, outf="blended.hdr", roh=0.0, rov=0.0, sfov=2.0, srcsize=6.7967e-05, bw=2.0,
               envmap=None, sunloc=None, check=None, margin=20, align=True, fisheye=True, **kwargs):
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
             f"align:{align} fisheye:{fisheye}")
    hdata, hh = io.hdr2carray(imgh, header=True)
    vdata, hv = io.hdr2carray(imgv, header=True)
    sdata, hs = io.hdr2carray(imgn, header=True)
    if align:
        if margin < 1:
            t = slice(None, None)
            margin = 0
        else:
            t = slice(margin, -margin)
        xo, yo = sb.align_images(hdata[0, t, t], vdata[0, t, t])
        if np.sum(np.abs((xo, yo))) > 0:
            sdata = sdata[:, t, t]
            vdata = vdata[:, t, t]
            xt = slice(max(0, margin + xo), min(hdata.shape[1], xo + hdata.shape[1] - margin))
            yt = slice(max(0, margin + yo), min(hdata.shape[2], yo + hdata.shape[2] - margin))
            xt2 = slice(0, xt.stop-xt.start)
            yt2 = slice(0, yt.stop-yt.start)
            hdata2 = np.zeros_like(sdata)
            hdata2[:, xt2, yt2] = hdata[:, xt, yt]
            hdata = hdata2
            hh.append(f"\tSHADOWBAND_IMAGE_ALIGN= {xo} {yo}")
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
    header = [sbobt] + hh + hv + hs + [vm.header()]
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
    io.array2hdr(imgv, None, header=[f"VIGNETTING_CORRECTION= {vfile}"] + header, clean=True)


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
    if type(img) == str:
        imgd, header = io.hdr2carray(img, header=True)
        try:
            header += [imagetools.hdr2vm(img).header()]
        except AttributeError:
            pass
        rgb, outheader = pl.color_convert_img(imgd, header, inp, outp, xyzrgb=xyzrgb, rgbrgb=rgbrgb, oxyzrgb=oxyzrgb)
        io.array2hdr(rgb, None, header=outheader, clean=True)
    else:
        rgb2rgb, _ = pl.prep_color_transform(inp, outp, xyzrgb=xyzrgb, rgbrgb=rgbrgb, oxyzrgb=oxyzrgb)
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