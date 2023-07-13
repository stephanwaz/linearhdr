import os
import shutil
import sys

import numpy as np
from clasp import click
import clasp.click_ext as clk
from clasp.script_tools import try_mkdir

import pylinearhdr
import raytools
from raytools import io, imagetools
from raytools.utility import pool_call
from pylinearhdr import calibrate as ca
from pylinearhdr import make_list as ml
from pylinearhdr import shadowband as sb
from pylinearhdr import pylinearhdr as pl

@click.group()
@click.option('-config', '-c', type=click.Path(exists=True),
              help="path of config file to load")
@click.option('-n', default=None, type=int,
              help='sets the environment variable RAYTOOLS_PROC_CAP set to'
                   ' 0 to clear (parallel processes will use cpu_limit)')
@click.option('--opts', '-opts', is_flag=True,
              help="check parsed options")
@click.option('--debug', is_flag=True,
              help="show traceback on exceptions")
@click.version_option(version=pylinearhdr.__version__)
@click.pass_context
def main(ctx, config=None, n=None,  **kwargs):
    """the pylinearhdr executable is a command line interface to utility commands
    as part of the linearhdr package, some scripts available without any dependencies via cmake install.
    """
    os.environ['CLASP_PIPE'] = '1'
    raytools.io.set_nproc(n)
    ctx.info_name = 'pylinearhdr'
    clk.get_config(ctx, config, None, None, None)


@main.command()
@click.argument("imgs", callback=clk.are_files)
@click.argument("crop", callback=clk.split_int)
@click.option("-hdropts", default="", help="options to pass to linearhdr (put in qoutes)")
@click.option("-badpixels", callback=clk.is_file,
              help="file of bad pixels (rows: xpix ypix 0) where xpix is from left and ypix is from top")
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def calibrate(ctx, imgs, crop, badpixels=None, hdropts="", **kwargs):
    """calibration routine, see README.rst

    imgs: list of images
    crop: help="<upper_left> <upper_right> <width> <height>"
    """
    ppms = pool_call(ca.get_raw_frame, imgs, *crop[0:4], hdropts, bad_pixels=badpixels, expandarg=False)
    ca.report(ppms)


@main.command()
@click.argument("imgs", callback=clk.are_files)
@click.option("--shell/--no-shell", default=False,
              help="output shell file for use with stdin of linearhdr: bash output.sh | linearhdr")
@click.option("--overwrite/--no-overwrite", default=False,
              help="run dcraw_emu even if output file exists")
@click.option("-header-line", '-hl', multiple=True,
              help="lines to append to HDR header, e.g. LOCATION= 46.522833,6.580500")
@click.option("--correct/--no-correct", default=False,
              help="apply correction to nominal aperture and shutter speed values, use with linearhdr --exact")
@click.option("--listonly/--no-listonly", default=False,
              help="skip execution and just print metadata")
@click.option("-scale", default=1.0, help="calibration scale factor (applies to ISO, so do not use -s with linearhdr)")
@click.option("-nd", default=0.0, help="additional ND filter (applies to ISO, so do not use -s with linearhdr)")
@click.option("-saturation", "-saturation-offset", "-s", default=0.2, help="saturation offset (only needed by --listonly)")
@click.option("-r", "-range", default=0.01, help="lower range of single raw exposure (only needed by --listonly)")
@click.option("-crop", callback=clk.split_int,
              help="crop ppm (left upper W H)")
@click.option("-badpixels", callback=clk.is_file,
              help="file of bad pixels (rows: xpix ypix 0) where xpix is from left and ypix is from top")
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def makelist(ctx, imgs, shell=False, overwrite=False, correct=False, listonly=False, scale=1.0, nd=0.0, saturation=0.2, r=0.01,
             crop=None, badpixels=None, header_line=None, **kwargs):
    """make list routine, use to generate input to linearhdr"""
    if listonly:
        shell = False
        correct = True
        overwrite = False
    ppms = pool_call(ml.get_raw_frame, imgs, correct=correct, overwrite=overwrite,
                     listonly=listonly, crop=crop, bad_pixels=badpixels, expandarg=False)
    print(ml.header_info(imgs[0]))
    for li in header_line:
        print(f"# {li}")
    ml.report(ppms, shell, listonly, scale=scale * 10**nd, sat_w=1-saturation, sat_b=r)


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
@click.option("-sfov", default=180.0,
              help="field of view around shaded source that is valid in imgn (in case of partial ND filter)")
@click.option("-srcsize", default=6.7967e-05,
              help="solid angle of shaded source (steradians)")
@click.option("-bw", default=2.0,
              help="size (in degrees) of shadow band. Note: this value is doubled to account for error in centering on "
                   "source")
@click.option("--flip/--no-flip", default=False,
              help="by default the imgh will be used in the UL and LR quadrants, flip=True will use imgh UR and LL")
@click.option("--envmap/--no-envmap", default=False,
              help="do not add source to image, instead, return as radiance source description")
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def shadowband(ctx, imgh, imgv, imgn, outf="blended.hdr", roh=0.0, rov=0.0, sfov=4.0, srcsize=6.7967e-05, bw=2.0, flip=False, envmap=False,
               **kwargs):
    """merge set of three images with horizontal, vertical and no shadow band all images are assumed to be 180 degree
    angular fisheye.

    imgh: hdr with horizontal aligned shadowband
    imgv: hdr with veritcal aligned shadowband
    imgn: hdr with no shadowband
    """
    hdata = io.hdr2carray(imgh)
    vdata = io.hdr2carray(imgv)
    sdata = io.hdr2carray(imgn)
    blended = sb.shadowband(hdata, vdata, sdata, roh=roh, rov=rov, sfov=sfov, srcsize=srcsize, bw=bw, flip=flip,
                            envmap=envmap)
    header = imagetools.hdr2vm(imgh).header()
    if envmap:
        radout = outf.rsplit(".", 1)[0] + ".rad"
        f = open(radout, 'w')
        srcrad = (srcsize/np.pi)**.5 * 180/np.pi * 2
        f.write("void light sun 0 0 3 {} {} {}\n".format(*blended[1][-1]))
        f.write("sun source solar 0 0 4 {} {} {} {}\n\n".format(*blended[1][0:3], srcrad))
        f.write(f"void colorpict imgfunc\n7 red green blue {outf} fisheye.cal fish_u fish_v\n0 0\n")
        f.write("imgfunc glow imgglow 0 0 4 1 1 1 0\n")
        f.write("imgglow source sky 0 0 4 0 1 0 180\n")
        io.carray2hdr(blended[0], outf, [header])
    else:
        io.carray2hdr(blended, outf, [header])


@main.command()
@click.argument("imgs", callback=clk.are_files)
@click.option("-out", default="img", help="directory base name, becomes: img_XXX")
@click.option("-starti", default=0, help="start index")
@click.option("--ascend/--descend", default=True,
              help="shot order exposure time. ascending=shortest->longest")
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def sort(ctx, imgs, out="img", starti=0, ascend=True, **kwargs):
    """sort sequence of images into folders based on changes in exposure time"""
    infos = pool_call(ml.get_raw_frame, imgs, listonly=True, correct=True, expandarg=False)
    i = starti
    if ascend:
        expt = 1e9
    else:
        expt = 0
    for info in infos:
        if info[1] == expt:
            print("Duplicate!")
            continue
        if (info[1] < expt) != ascend:
            i += 1
        expt = info[1]
        try_mkdir(f"{out}_{i:03d}")
        shutil.move(info[0], f"{out}_{i:03d}")

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