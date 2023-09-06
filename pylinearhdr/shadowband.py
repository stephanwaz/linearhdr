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
    return np.stack((np.arctan2(va[:, 0], va[:, 1]), np.arctan2(ha[:, 2], ha[:, 1]))).T


def blend_bands(x, b, c=None):
    if c is None:
        c = b
    return np.where(x > b, np.where(x < b+c, np.cos((x-b)*np.pi/c)/2+.5, 0), 1)


def shadowband(hdata, vdata, sdata, roh=0.0, rov=0.0, sfov=2.0, srcsize=6.7967e-05, bw=2.0, flip=False,
               envmap=False, sunloc=None, check=None):
    vm = ViewMapper(viewangle=180)
    res = hdata.shape[-1]
    band = bw/vm.viewangle*res

    v = vm.pixelrays(res).reshape(-1, 3)
    omega = vm.pixel2omega(vm.pixels(res), res).ravel()
    lum = np.squeeze(sdata.reshape(-1, res**2))
    if lum.shape[0] == 3:
        lum = io.rgb2rad(lum.T)
    # find peak in sdata (un shaded image)
    if sunloc is not None:
        pxyz = vm.pixel2ray(np.array(sunloc[0:2])[None], res)
    else:
        pxyz = im.find_peak(v, omega, lum, peaka=srcsize)[0][0]

    # blend along midpoint between two shadowbands
    rotation = 45+np.average((roh, rov))
    # profile angles of sun
    sangles = profile_angles(pxyz, rotation, rotation).reshape(1, 1, 2)
    # profile angles of pixels
    pangles = profile_angles(v, rotation, rotation).reshape(res, res, 2)
    pdiff = sangles - pangles
    mask = np.zeros((res, res))

    # generate masks for each quadrant
    ul = np.logical_and(pdiff[..., 0] > 0, pdiff[..., 1] <= 0)
    lr = np.logical_and(pdiff[..., 0] <= 0, pdiff[..., 1] > 0)
    ur = np.logical_and(pdiff[..., 0] <= 0, pdiff[..., 1] <= 0)
    ll = np.logical_and(pdiff[..., 0] > 0, pdiff[..., 1] > 0)

    mask[ul] = 0
    mask[ur] = 1
    mask[ll] = 1
    mask[lr] = 0
    mask = ndimage.uniform_filter(mask, band/4)

    # profile angles of sun (now along shadowband orientation)
    sangles = profile_angles(pxyz, roh, rov).reshape(1, 1, 2)
    # profile angles of pixels
    pangles = profile_angles(v, roh, rov).reshape(res, res, 2)
    pdiff = np.abs(sangles - pangles) * 180/np.pi
    # anything within 3 degrees of both shadowbands
    outer_src_mask = np.all(pdiff < bw * 1.5, axis=2)
    # blur to create tolerance regions
    outer_src_mask = ndimage.uniform_filter(outer_src_mask.astype(float), band/4)
    # use this area for interpolation
    inter_src_mask = np.logical_and(outer_src_mask < .5, outer_src_mask > .01)

    # first make initial estimate (quadrants for "safe" values, maximum for near values)
    blend = (vdata * mask + hdata * (1-mask)) * (1 - outer_src_mask) + outer_src_mask * np.maximum(vdata, hdata)

    # interpolate everything in outer regions to allow for blending
    src_mask = outer_src_mask > .01

    if check is not None:
        mask3 = np.copy(np.broadcast_to(mask, (3, *mask.shape)))
        mask3[0, src_mask] = 0
        mask3[1, src_mask] = 1 - mask3[1, src_mask]
        mask3[2, src_mask] = mask3[2, src_mask]
        mask3[1:, outer_src_mask >= .5] = 0
        mask3[0, outer_src_mask >= .5] = 1
        # srcmaskimg = np.copy(blend)
        # srcmaskimg *= (1-src_mask + .2)/1.2
        # srcmaskimg[0, inter_src_mask] *= 20

        io.array2hdr(mask3, f"{check}_mask.hdr")
        # io.array2hdr(srcmaskimg, f"{check}_srcmask.hdr")

    vms = ViewMapper(dxyz=pxyz, viewangle=45)
    lum = blend[:, inter_src_mask].reshape(3, -1).T

    lp = LightPointKD(None, v[inter_src_mask.ravel()], lum, vm=vm, features=3, calcomega=False, write=False)
    i, w = lp.interp(v[src_mask.ravel()], 20, lum=False, angle=False)

    clum = lp.apply_interp(i, lum, w)

    # interpolate between edge of interpolation frame and extrapolated center
    # distance from frame to center
    i, d = lp.query_ray(v[src_mask.ravel()])
    d2 = vms.radians(lp.vec[i])
    # distance from center to pixel
    d3 = vms.radians(v[src_mask.ravel()])
    df = np.clip(d3/d2, 0, 1)[:, None]

    # extrapoloate center
    d = vms.degrees(v[inter_src_mask.ravel()])
    sr, ir, _, _, _ = stats.linregress(d, lum[:, 0])
    sg, ig, _, _, _ = stats.linregress(d, lum[:, 1])
    sb, ib, _, _, _ = stats.linregress(d, lum[:, 2])
    clum = clum * df + np.array([(ir, ig, ib)]) * (1-df)
    # blend in combined interpolation
    blend[:, src_mask] = blend[:, src_mask] * (1 - outer_src_mask[src_mask][None]) + outer_src_mask[src_mask][None] * clum.T

    # isolate area around source to gather lens flare energy
    vm_valid = ViewMapper(pxyz, sfov).in_view(v, indices=False)
    flare = sdata.reshape(3, -1)[:, vm_valid] - blend.reshape(3, -1)[:, vm_valid]
    # discard under exposed pixels
    flare[flare < 0] = 0
    sol_luminance = np.sum(flare * omega[None, vm_valid], axis=1) / srcsize
    if check is not None:
        flare2 = np.maximum(sdata - blend, 0).reshape(3,-1)
        flare2[:, np.logical_not(vm_valid)] = 0
        io.array2hdr(flare2.reshape(sdata.shape), f"{check}_src.hdr")

    if envmap is not None:
        skyonly = np.copy(blend)
        source = (*pxyz, srcsize, sol_luminance)
    else:
        skyonly = source = None
    # draw source on image
    mask = vm.in_view(v)
    src = SrcViewPoint(None, np.asarray(pxyz).reshape(-1, 3), sol_luminance, res=1)
    src.add_to_img(blend, v, mask, vm=vm)

    return blend, skyonly, source
