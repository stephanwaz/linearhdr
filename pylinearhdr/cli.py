import sys

import numpy as np

from clasp import click
import clasp.click_ext as clk
from clasp.script_tools import pipeline

import pylinearhdr
import raytools
from raytools import io, imagetools
from raytools.utility import pool_call
from pylinearhdr import calibrate as ca
from pylinearhdr import make_list as ml
from pylinearhdr import shadowband as sb

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
    raytools.io.set_nproc(n)
    ctx.info_name = 'raytools'
    clk.get_config(ctx, config, None, None, None)


@main.command()
@click.argument("imgs", callback=clk.are_files)
@click.argument("crop", callback=clk.split_int)
@click.option("-hdropts", default="", help="options to pass to linearhdr (put in qoutes)")
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def calibrate(ctx, imgs, crop, hdropts="", **kwargs):
    """calibration routine, see README.rst

    imgs: list of images
    crop: help="<upper_left> <upper_right> <width> <height>"
    """
    ppms = pool_call(ca.get_raw_frame, imgs, *crop[0:4], hdropts, expandarg=False)
    ca.report(ppms)


@main.command()
@click.argument("imgs", callback=clk.are_files)
@click.option("--shell/--no-shell", default=False,
              help="output shell file for use with stdin of linearhdr: bash output.sh | linearhdr")
@click.option("--overwrite/--no-overwrite", default=False,
              help="run dcraw_emu even if output file exists")
@click.option("--correct/--no-correct", default=False,
              help="apply correction to nominal aperture and shutter speed values, use with linearhdr --exact")
@click.option("--listonly/--no-listonly", default=False,
              help="skip execution and just print metadata")
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def makelist(ctx, imgs, shell=False, overwrite=False, correct=False, listonly=False, **kwargs):
    """make list routine, use to generate input to linearhdr"""
    if listonly:
        shell = False
        correct = False
        overwrite = False
    ppms = pool_call(ml.get_raw_frame, imgs, correct=correct, overwrite=overwrite, listonly=listonly, expandarg=False)
    ml.report(ppms, shell, listonly)

@main.command()
@click.argument("imgh", callback=clk.is_file)
@click.argument("imgv", callback=clk.is_file)
@click.argument("imgn", callback=clk.is_file)
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
@clk.shared_decs(clk.command_decs(pylinearhdr.__version__, wrap=True))
def shadowband(ctx, imgh, imgv, imgn, roh=0.0, rov=0.0, sfov=4.0, srcsize=6.7967e-05, bw=2.0, flip=False, **kwargs):
    """merge set of three images with horizontal, vertical and no shadow band all images are assumed to be 180 degree
    angular fisheye.

    imgh: hdr with horizontal aligned shadowband
    imgv: hdr with veritcal aligned shadowband
    imgn: hdr with no shadowband
    """
    hdata = io.hdr2carray(imgh)
    vdata = io.hdr2carray(imgv)
    sdata = io.hdr2carray(imgn)
    blended = sb.shadowband(hdata, vdata, sdata, roh=roh, rov=rov, sfov=sfov, srcsize=srcsize, bw=bw, flip=flip)
    header = imagetools.hdr2vm(imgh).header()
    io.carray2hdr(blended, None, [header])


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