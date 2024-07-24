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
import sys

import numpy as np
from raytools import io, translate
from raytools import imagetools as im
from raytools.mapper import ViewMapper
from scipy import ndimage, stats
from raytraverse.lightpoint import SrcViewPoint, LightPointKD

def profile_angles(a, rh=0.0, rv=0.0):
    a = np.atleast_2d(a)
    if rh == 0:
        ha = a
    else:
        ha = translate.rotate_elem(np.atleast_2d(a), -rh, 1)
    if rv == 0:
        va = a
    else:
        va = translate.rotate_elem(np.atleast_2d(a), -rv, 1)
    return np.stack((np.arctan2(va[:, 0], va[:, 1]), np.arctan2(ha[:, 2], ha[:, 1]))).T  * 180/np.pi


def blend_bands(x, b, c=None):
    if c is None:
        c = b
    return np.where(x > b, np.where(x < b+c, np.cos((x-b)*np.pi/c)/2+.5, 0), 1)


def generate_masks(v, sb_cpxyz, roh=0, rov=0, res=1000, bw=2.0):
    # blend along midpoint between two shadowbands
    rotation = 45+np.average((roh, rov))
    # profile angles of sun
    sangles45 = profile_angles(sb_cpxyz, rotation, rotation).reshape(1, 1, 2)
    # profile angles of pixels
    pangles45 = profile_angles(v, rotation, rotation).reshape(res, res, 2)
    pdiff45 = sangles45 - pangles45
    mask = np.zeros((res, res))

    # generate masks for each quadrant
    ul = np.logical_and(pdiff45[..., 0] > 0, pdiff45[..., 1] <= 0)
    lr = np.logical_and(pdiff45[..., 0] <= 0, pdiff45[..., 1] > 0)
    ur = np.logical_and(pdiff45[..., 0] <= 0, pdiff45[..., 1] <= 0)
    ll = np.logical_and(pdiff45[..., 0] > 0, pdiff45[..., 1] > 0)

    mask[ul] = 0
    mask[ur] = 1
    mask[ll] = 1
    mask[lr] = 0

    # profile angles of sun (now along shadowband orientation)
    sangles = profile_angles(sb_cpxyz, roh, rov).reshape(1, 1, 2)
    # profile angles of pixels
    pangles = profile_angles(v, roh, rov).reshape(res, res, 2)

    pdiff = np.abs(sangles - pangles)
    band = bw/180*res
    # don't blur near intersection too much, blur more once its safe
    init_blend_mask = np.any(pdiff < bw * 1.5, axis=2)
    # add a buffer of band# pixels
    init_blend_mask = ndimage.maximum_filter(init_blend_mask, band)
    mask = ndimage.uniform_filter(mask, band)

    vangle = sangles[..., 0]/90
    sv = np.abs(np.sin(sangles * np.pi/180)).ravel()
    angrat = sv[0]/np.sum(sv)
    # offset away from more glancing and keep some blur
    mask2 = ndimage.uniform_filter((mask > angrat).astype(float), band/4)
    mask[init_blend_mask] = mask2[init_blend_mask]

    vms = ViewMapper(dxyz=sb_cpxyz, viewangle=180)
    # for primary interpolation make a diamond shape mask around sun to match band intersection
    pri_src_mask = np.logical_and(pdiff[..., 0] < bw * 0.75,
                                  pdiff[..., 1] < bw * 0.75 * (1 + abs(vangle)))
    pri_src_mask = ndimage.maximum_filter(pri_src_mask, band)
    # also include perimeter around sun to avoid overflow and marginal flare
    sec_src_data = vms.degrees(v).reshape(mask.shape)
    sec_src_mask = sec_src_data < bw * 1.5
    pri_src_mask = np.logical_or(pri_src_mask, sec_src_mask)
    pri_src_mask = ndimage.uniform_filter(pri_src_mask.astype(float), band / 2)
    return mask, pri_src_mask, sec_src_data


def shadowband(hdata, vdata, sdata, roh=0.0, rov=0.0, sfov=2.0, srcsize=6.7967e-05, bw=2.0,
               envmap=False, sunloc=None, check=None):
    header = []
    vm = ViewMapper(viewangle=180)
    res = hdata.shape[-1]
    band = bw/vm.viewangle*res
    v = vm.pixelrays(res).reshape(-1, 3)
    viewmask = vm.in_view(v, indices=False).reshape(res, res)
    omega = vm.pixel2omega(vm.pixels(res), res).ravel()
    lum = np.squeeze(sdata.reshape(-1, res**2))
    if lum.shape[0] == 3:
        lum = io.rgb2rad(lum.T)
    # find peak in sdata (un shaded image)
    if sunloc is not None:
        sb_cpxyz = vm.pixel2ray(np.array(sunloc[0:2])[None], res)
        up = translate.degrees(sb_cpxyz.ravel(), v) < 3
    else:
        up = v[:, 2] > 0
        sb_cpxyz = None
    # only use peak if ~> 1000000 cd/m^2
    pklum = lum[up] > 6000
    hassource = np.sum(pklum) > 0
    pxyz = None
    if hassource:
        try:
            pxyz = im.find_peak(v[up][pklum], omega[up][pklum], lum[up][pklum], peaka=srcsize)[0][0]
        except TypeError:
            print("Warning Zero Sun Energy! check -sunloc parameter", file=sys.stderr)
            hassource = False
    else:
        print("Warning Zero Sun Energy! check -sunloc parameter", file=sys.stderr)
    if sb_cpxyz is None:
        sb_cpxyz = pxyz
    if sb_cpxyz is None:
        raise ValueError("No Peak found in ND image and no solar location give. Cannot merge.")
    mask, pri_src_mask, sec_src_mask = generate_masks(v, sb_cpxyz, roh, rov, res, bw=bw)

    pri_src_data = np.logical_and(np.logical_and(pri_src_mask < .5, pri_src_mask > (band - 4)/(2 * band)), viewmask)
    sec_src_data = np.logical_and(np.logical_and(sec_src_mask > bw * 1.5, sec_src_mask < bw * 2), viewmask)
    # interpolate everything in outer regions to allow for blending
    src_mask = np.logical_and(pri_src_mask > .01, viewmask)

    # first make initial estimate (quadrants for "safe" values, maximum for near values)
    blend = (vdata * mask + hdata * (1-mask)) * (1 - pri_src_mask) + pri_src_mask * np.maximum(vdata, hdata)

    if check is not None:
        io.array2hdr(mask, f"{check}_mask.hdr")
        io.array2hdr(pri_src_mask, f"{check}_srcmask.hdr")
        io.carray2hdr(hdata, f"{check}_H.hdr")
        io.carray2hdr(vdata, f"{check}_V.hdr")
        io.carray2hdr(sdata, f"{check}_ND.hdr")

    lum_u = blend[:, pri_src_data].reshape(3, -1).T
    bw2 = min(300, int(np.sum(pri_src_data)*.2))
    lp_u = LightPointKD(None, v[pri_src_data.ravel()], lum_u, vm=vm, features=3, calcomega=False, write=False)

    i_u, w_u = lp_u.interp(v[src_mask.ravel()], bw2, lum=False, angle=False, dither=False, offset=True)
    clum = lp_u.apply_interp(i_u, lum_u, w_u)

    # this is now a usable result, but in some cases we should assume a higher peak in the middle
    blend[:, src_mask] = blend[:, src_mask] * (1 - pri_src_mask[src_mask][None]) + pri_src_mask[src_mask][None] * clum.T

    # so... interpolation step 2: fit linear model for brightening closer to sun
    # interpolate between edge of interpolation frame and extrapolated center

    lum = blend[:, sec_src_data].reshape(3, -1).T
    d = sec_src_mask[sec_src_data].ravel()
    _, ir, rr, _, _ = stats.linregress(d, lum[:, 0])
    _, ig, rg, _, _ = stats.linregress(d, lum[:, 1])
    _, ib, rb, _, _ = stats.linregress(d, lum[:, 2])
    # only do something if correlated in right direction (negative means peak in the middle)
    if np.min((rr, rg, rb)) < -0.25:
        src_mask = np.logical_and(sec_src_mask < bw * 1.5, viewmask)
        dn = bw * 1.5
        # distance from center to pixel
        d3 = sec_src_mask[src_mask].ravel()
        # distance from center to pixel scaled by inner mask region
        df = np.clip(d3/dn, 0, 1)[None]
        # this is the extrapolated value at peak
        clum = np.array([ir, ig, ib])[:, None]
        # make sure we don't reduce the value, so take max of existing, and blend + (clum - blend) (1-df) (so we only add slope but not base)
        # which is rearranged to: blend * df + clum * (1 - df) so that the new part can be blended in with outer_src_mask
        avg = clum.ravel() - np.average(blend[:, src_mask], axis=1)
        blend[:, src_mask] = np.maximum(blend[:, src_mask], blend[:, src_mask] * df + clum * (1 - df))
        header.append("EXTRAPOLATED_CIRCUMSOLAR= {:.05G},{:.05G},{:.05G} / {:.02f} degrees".format(*avg, bw * 1.5))
    if envmap is not None:
        skyonly = np.copy(blend)
    else:
        skyonly = None

    if hassource:
        # isolate area around source to gather lens flare energy
        vm_valid = ViewMapper(pxyz, sfov).in_view(v, indices=False)
        flare = sdata.reshape(3, -1)[:, vm_valid]
        skyflare = blend.reshape(3, -1)[:, vm_valid]
        sol_lumrgb = np.sum(flare * omega[None, vm_valid], axis=1) / srcsize
        sol_lum = io.rgb2rad(sol_lumrgb)
        sky_lum = io.rgb2rad(np.sum(skyflare * omega[None, vm_valid], axis=1) / srcsize)
        cf = (sol_lum - sky_lum) / sol_lum
        sol_lumrgb *= cf
        if check is not None:
            flare2 = np.copy(sdata).reshape(3, -1) * cf
            flare2[:, np.logical_not(vm_valid)] = 0
            io.array2hdr(flare2.reshape(sdata.shape), f"{check}_src.hdr")

        opxyz = translate.rotate_elem(pxyz[None], 180).ravel()
        source = (*opxyz, srcsize, sol_lumrgb)
        # draw source on image
        mask = vm.in_view(v)
        src = SrcViewPoint(None, np.asarray(pxyz).reshape(-1, 3), sol_lumrgb, res=1)
        src.add_to_img(blend, v, mask, vm=vm)
        vm = ViewMapper((0.0, -1.0, 0.0), viewangle=180)
        header.append("CENTRAL_SOURCE_PIXEL= {} {}".format(*vm.ray2pixel(source[0:3], hdata.shape[-1], False).ravel()))
    else:
        source = None
    return blend, skyonly, source, header
