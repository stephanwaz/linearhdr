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
import sys
from subprocess import Popen, PIPE

import numpy as np
from clasp import click
from clasp.script_tools import pipeline

from raytools import io
from scipy import optimize
from pylinearhdr import pylinearhdr as pl


def hdr_get_primaries(img, defp=(0.640, 0.330, 0.290, 0.600, 0.150, 0.060),
                      defwp=(0.3333, 0.3333)):
    pri = io.hdr_header(img, items=["TargetPrimaries", "TargetWhitePoint"])
    if len(pri[0]) == 0:
        pri[0] = defp
    else:
        pri[0] = [float(i) for i in pri[0].split()]
    if len(pri[1]) == 0:
        pri[1] = defwp
    else:
        pri[1] = [float(i) for i in pri[1].split()]
    return pri


def calibrate(ref, test, rc, tc, alternate=False, refimg=True, refcol='rad'):
    if alternate and rc is None and tc is None:
        rcells, tcells = load_test_cells_simul(ref, test)
    else:
        if refimg:
            rcells = load_test_cells(ref, rc)
        else:
            rcells = None
        tcells = load_test_cells(test, tc)
    refd = load_data(ref, rcells)
    testd = load_data(test, tcells)

    # get reference colorspace
    rgb_xyz, _ = pl.prep_color_transform(refcol, 'xyz')
    xyz_rgb = np.linalg.inv(rgb_xyz)


    cst = hdr_get_primaries(test, None, None)
    xyz_cam = io.hdr_header(test, items=["XYZCAM"])[0]
    # check if test is in a colorspace
    if cst[0] is not None:
        cam_xyz = pl.pw_mtx(*cst)
    # check in test has initial transform matrix
    elif len(xyz_cam) > 0:
        xyz_cam = np.array([float(i) for i in xyz_cam.split()]).reshape(3, 3)
        cam_xyz = np.linalg.inv(xyz_cam)
    # assume test is in same space as reference
    else:
        cam_xyz = rgb_xyz

    popts = np.get_printoptions()
    np.set_printoptions(5, suppress=True)

    # test in raw RGB
    B = testd

    if len(refd.shape) == 1:
        # luminance scaling:
        B_y = np.einsum('ij,jk->ik', cam_xyz, B)[1]
        lum_scale = np.linalg.lstsq(B_y[:, None], refd, rcond=None)[0][0]

        result = dict(start=dict(xyzcam=np.linalg.inv(cam_xyz)),
                      lum=dict(xyzcam=np.linalg.inv(cam_xyz * lum_scale)))

        print("************************* RESULTS *************************")
        for k, v in result.items():
            B_xyz = np.einsum('ij,jk->ik', np.linalg.inv(v['xyzcam']), B)
            v['dlum'] = rel_lum_diff(refd, B_xyz[1])
            v['dlum_full'] = rel_lum_diff(refd, B_xyz, False)
            print(f"lum diff for {k}:\t{v['dlum']:.05f}")
        print(f"\nluminance scale factor: {lum_scale:.05f}")
        print("***********************************************************")
    else:
        # reference data in XYZ
        A = np.einsum('ij,jk->ik', rgb_xyz, refd)
        # luminance scaling:
        B_xyz = np.einsum('ij,jk->ik', cam_xyz, B)
        lum_scale = np.linalg.lstsq(B_xyz[1][:, None], A[1], rcond=None)[0][0]

        # reference channel scaling:
        B_rgb = np.einsum('ij,jk->ik', xyz_rgb, B_xyz)
        rf = np.linalg.lstsq(B_rgb[0][:, None], refd[0], rcond=None)[0]
        gf = np.linalg.lstsq(B_rgb[1][:, None], refd[1], rcond=None)[0]
        bf = np.linalg.lstsq(B_rgb[2][:, None], refd[2], rcond=None)[0]
        rgb_scale = np.array([rf, gf, bf]).ravel()

        # linear color matrix optimization:
        cam_xyz_opt = np.linalg.solve((B@B.T).T, (A@B.T).T).T
        # non-linear color matrix optimization:
        cam_xyz_opt2 = optimize.minimize(color_diff_mtx, cam_xyz_opt.ravel(), args=(A, B),
                                         method='SLSQP', options=dict(maxiter=100)).x.reshape(3, 3)

        result = dict(start=dict(xyzcam=np.linalg.inv(cam_xyz)),
                      lum=dict(xyzcam=np.linalg.inv(cam_xyz * lum_scale)),
                      rgb=dict(xyzcam=np.linalg.inv(cam_xyz * rgb_scale[:, None])),
                      lopt=dict(xyzcam=np.linalg.inv(cam_xyz_opt)),
                      nopt=dict(xyzcam=np.linalg.inv(cam_xyz_opt2)))

        print("************************* RESULTS *************************")
        for k, v in result.items():
            B_xyz = np.einsum('ij,jk->ik', np.linalg.inv(v['xyzcam']), B)
            v['dlab'] = color_diff_lab(A, B_xyz)
            v['duv'] = color_diff_uv(A, B_xyz)
            v['dxy'] = color_diff_xy(A, B_xyz)
            v['dlab_full'] = color_diff_lab(A, B_xyz, False)
            v['duv_full'] = color_diff_uv(A, B_xyz, False)
            v['dxy_full'] = color_diff_xy(A, B_xyz, False)
            print(f"diff for {k}:\tlab:{v['dlab']:.05f}\tuv:{v['duv']:.05f}\txy:{v['dxy']:.05f}")
        print("***********************************************************")

    # print(result['lopt']['dlab_full'])
    np.set_printoptions(**popts)
    return result, refd, B


def load_data(img, cells, zero=False):
    if cells is None:
        return np.loadtxt(img).T
    img = io.hdr2carray(img) * 179
    data = []
    for cell in cells:
        x1, y1, x2, y2 = cell
        if zero and np.any(img[:, x1:x1+x2, y1:y1+y2] == 0):
            data.append(np.zeros(3))
        else:
            rgb = np.average(img[:, x1:x1+x2, y1:y1+y2], axis=(1, 2))
            data.append(rgb)
    return np.asarray(data).T


def xyz_2_xy(xyz):
    d = np.sum(xyz, axis=0)
    return xyz[0:2]/d


def xyz_2_yuv(xyz):
    d = xyz[0] + 15 * xyz[1] + 3 * xyz[2]
    u = 4 * xyz[0] / d
    v = 9 * xyz[1] / d
    return np.stack((xyz[1], u, v))


def xyz_2_lab(A, wp):
    Ar = A/wp[:, None]
    e = 216/24389
    k = 24389/27
    fA = (k * Ar + 16)/116
    fA.flat[Ar.ravel() > e] = np.cbrt(Ar)
    L = 116*fA[1] - 16
    a = 500 * (fA[0] - fA[1])
    b = 200 * (fA[1] - fA[2])
    return np.stack((L, a, b))


def color_diff_uv(A, A2, total=True):
    Auv = xyz_2_yuv(A)
    A2uv = xyz_2_yuv(A2)

    duv = A2uv[1:] - Auv[1:]
    result = np.sqrt(np.sum(np.square(duv), axis=0))
    if total:
        return np.sqrt(np.average(np.square(result)))
    else:
        return result


def color_diff_yuv(A, A2, total=True):
    Auv = xyz_2_yuv(A)
    A2uv = xyz_2_yuv(A2)

    duv = A2uv[1:] - Auv[1:]
    duv2 = np.sqrt(np.sum(np.square(duv), axis=0))
    dl = (A2uv[0] - Auv[0])/Auv[0]
    result = np.stack((dl, duv2))
    if total:
        return np.sqrt(np.average(np.square(result)))
    else:
        return result

def rel_lum_diff(A, A2, total=True):
    result = (A2 - A)/A
    if total:
        return np.sqrt(np.average(np.square(result)))
    else:
        return result

def color_diff_lab(A, A2, total=True):
    if A[1, 0] > A2[1, 0]:
        wp = A[:, 0] * 1.2
    else:
        wp = A2[:, 0] * 1.2
    la = xyz_2_lab(A, wp)
    la2 = xyz_2_lab(A2, wp)
    jnd = 2.3
    result = np.sqrt(np.sum(np.square(la2-la), axis=0))/jnd
    if total:
        return np.sqrt(np.average(np.square(result)))
    else:
        return result


def color_diff_xy(A, A2, total=True):
    xya = xyz_2_xy(A)
    xya2 = xyz_2_xy(A2)
    result = xya2 - xya
    if total:
        return np.sqrt(np.average(np.square(result)))
    else:
        return result


def color_diff_mtx(mtx, A, B, total=True):
    mtx = np.asarray(mtx).reshape(3,3)
    A2 = np.einsum('ij,jk->ik', mtx, B)
    return color_diff_lab(A, A2, total)


def color_diff_mtx_uv(mtx, A, B, total=True):
    mtx = np.asarray(mtx).reshape(3,3)
    A2 = np.einsum('ij,jk->ik', mtx, B)
    return color_diff_yuv(A, A2, total)


def color_diff_mtx_xy(mtx, A, B, total=True):
    mtx = np.asarray(mtx).reshape(3,3)
    A2 = np.einsum('ij,jk->ik', mtx, B)
    return color_diff_xy(A, A2, total)


def get_test_cells(img, outf=None):
    if outf is None:
        outf = img.rsplit(".", 1)[0] + ".txt"
        click.confirm(f'No cells given for {img}, generate?', default=True, abort=True)
        outf = click.prompt(f"name for cell location file", default=outf)
    f = open(outf, 'w')
    g = Popen(['ximage', '-op', img], stdout=PIPE, stderr=PIPE)

    print("In the same order as any other images (test, reference)\n"
          "click on (and press 't') the lower left and then upper right corner of each\n"
          "area to use for calibration. when finished, press 'q'")
    i = 0
    x, y, w, h = (0, 0, 0, 0)
    j = 0
    while True:
        pt = [int(k) for k in g.stdout.readline().decode('UTF-8').strip().split()]
        if len(pt) != 2:
            break
        if i == 0:
            x, y = pt
            i += 1
            print(f"cell {j}:\t{x}\t{y}", end="")
            j += 1
        else:
            w, h = (pt[0] - x, pt[1] - y)
            print(f"\t{w}\t{h}")
            i = 0
            f.write(f"{x}\t{y}\t{w}\t{h}\n")
        if g.poll() is not None:
            break
    f.close()
    return outf


def get_test_cells_simul(img1, img2, outf1=None, outf2=None):
    if outf1 is None:
        outf1 = img1.rsplit(".", 1)[0] + ".txt"
        outf1 = click.prompt(f"cell location file for {img1}", default=outf1)
    if outf2 is None:
        outf2 = img2.rsplit(".", 1)[0] + ".txt"
        outf2 = click.prompt(f"cell location file for {img2}", default=outf2)
    f1 = open(outf1, 'w')
    f2 = open(outf2, 'w')
    g = Popen(['ximage', '-op', img2, img1], stdout=PIPE, stderr=PIPE)

    print("Alternating between each image (starting with the reference,\n"
          "click on (and press 't') the lower left and then upper right corner of each\n"
          "area to use for calibration. when finished, press 'q'")
    i = 0
    x, y, w, h = (0, 0, 0, 0)
    j = 0
    cf = True
    while True:
        pt = [int(k) for k in g.stdout.readline().decode('UTF-8').strip().split()]
        if len(pt) != 2:
            break
        if cf:
            f = f1
            img = img1
        else:
            f = f2
            img = img2
        if i == 0:
            x, y = pt
            i += 1
            print(f"cell {j} for {img}:\t{x}\t{y}", end="")
            sys.stdout.flush()
            j += not cf
        else:
            w, h = (pt[0] - x, pt[1] - y)
            print(f"\t{w}\t{h}")
            i = 0
            cf = not cf
            f.write(f"{x}\t{y}\t{w}\t{h}\n")
        if g.poll() is not None:
            break
    f1.close()
    f2.close()
    return outf1, outf2


def load_test_cells(img, outf):
    confirm = False
    if outf is None or not os.path.isfile(outf):
        outf = get_test_cells(img, outf)
        confirm = True
    try:
        cells = np.loadtxt(outf, dtype=int)
        if cells.shape[1] != 4:
            raise ValueError
    except ValueError:
        click.confirm(f'Bad cell file, regenerate (this will erase {outf})?', default=True, abort=True)
        os.remove(outf)
        return load_test_cells(img, outf)
    cellarea = cells[:, 2] * cells[:, 3]
    print(f"statistics for {len(cells)} cells in {outf}:", file=sys.stderr)
    print(f"width (min, med, max): {np.percentile(cells[:, 2], (0, 50, 100))}", file=sys.stderr)
    print(f"height (min, med, max): {np.percentile(cells[:, 3], (0, 50, 100))}", file=sys.stderr)
    print(f"area (min, med, max): {np.percentile(cellarea, (0, 50, 100))}", file=sys.stderr)
    if confirm:
        click.confirm(f'proceed?', default=True, abort=True)
    return cells


def load_test_cells_simul(img1, img2):
    confirm = False
    outf1, outf2 = get_test_cells_simul(img1, img2)
    try:
        cells1 = np.loadtxt(outf1, dtype=int)
        cells2 = np.loadtxt(outf2, dtype=int)
        if cells1.shape[1] != 4 or cells2.shape[1] != 4 or len(cells1) != len(cells2):
            raise ValueError
    except ValueError:
        print(f'Bad cell file(s), or files do not match')
        raise click.Abort
    cellarea = cells1[:, 2] * cells1[:, 3]
    print(f"statistics for {len(cells1)} cells in {outf1}:")
    print(f"width (min, med, max): {np.percentile(cells1[:, 2], (0, 50, 100))}")
    print(f"height (min, med, max): {np.percentile(cells1[:, 3], (0, 50, 100))}")
    print(f"area (min, med, max): {np.percentile(cellarea, (0, 50, 100))}")
    cellarea = cells2[:, 2] * cells2[:, 3]
    print(f"statistics for {len(cells2)} cells in {outf2}:")
    print(f"width (min, med, max): {np.percentile(cells2[:, 2], (0, 50, 100))}")
    print(f"height (min, med, max): {np.percentile(cells2[:, 3], (0, 50, 100))}")
    print(f"area (min, med, max): {np.percentile(cellarea, (0, 50, 100))}")
    if confirm:
        click.confirm(f'proceed?', default=True, abort=True)
    return cells1, cells2


def average_green(img, croparg):
    """return average green value and fraction in range for a single frame and crop area"""
    runcom = [f"pylinearhdr run -colorspace raw -hdropts ' -m 0 -x 0 --tsv' --half {croparg} '{img}'",
              "getinfo -d -", "rcalc -e '$1=$2;cond=$2;$2=1'", 'total -m']
    info = list(pl.info_from_exif(img.split()[0], True, False))
    result = [float(i) for i in pipeline(runcom).strip().split()]
    if len(result) != 2:
        result = [0.0, 0.0]
    return info[0:2] + result
