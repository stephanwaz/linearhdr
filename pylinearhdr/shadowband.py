import sys

import numpy as np
from raytools import io, translate
from raytools import imagetools as im
from raytools.mapper import ViewMapper
from scipy import ndimage
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
               envmap=False, sunloc=None):
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

    # calculate profile angles relative to shadow bands
    sangles = profile_angles(pxyz, roh, rov).reshape(1, 1, 2)
    pangles = profile_angles(v, roh, rov).reshape(res, res, 2)
    pdiff = sangles - pangles
    mask = np.zeros((res, res))

    # generate masks for each quadrant
    ul = np.logical_and(pdiff[..., 0] > 0, pdiff[..., 1] <= 0)
    lr = np.logical_and(pdiff[..., 0] <= 0, pdiff[..., 1] > 0)
    ur = np.logical_and(pdiff[..., 0] <= 0, pdiff[..., 1] <= 0)
    ll = np.logical_and(pdiff[..., 0] > 0, pdiff[..., 1] > 0)

    # find location of band at each row/column
    x_max = np.sum(pdiff[..., 0] > 0, axis=0, keepdims=True)
    y_max = np.sum(pdiff[..., 1] > 0, axis=1, keepdims=True)
    idx = np.arange(res)

    # offset band width into quadrants based on distance from band
    # 0 is H 1 is V
    if flip:
        mask[ul] = 1 - blend_bands((x_max - idx[:, None])[ul], band)
        mask[lr] = 1 - blend_bands((idx[:, None] - x_max)[lr], band)
        mask[ur] = blend_bands((idx[None] - y_max)[ur], band)
        mask[ll] = blend_bands((y_max - idx[None])[ll], band)
    else:
        mask[ul] = blend_bands((idx[None] - y_max)[ul], band)
        mask[lr] = blend_bands((y_max - idx[None])[lr], band)
        mask[ur] = 1 - blend_bands((idx[:, None] - x_max)[ur], band)
        mask[ll] = 1 - blend_bands((x_max - idx[:, None])[ll], band)

    # outer peak area (for max replacement)
    peakr = 2
    # area around peak to replace with interpolation
    rpeakr = 1.5
    # inner peak area (for donut interpolation)
    cpeakr = 1

    # isolate pixels near peak where blending interacts
    vms = ViewMapper(dxyz=pxyz, viewangle=bw*peakr*2)
    angle = vms.degrees(v)
    outer_src_mask = angle < bw * peakr
    source_mask = ndimage.uniform_filter(outer_src_mask.reshape(mask.shape).astype(float), band/2)

    # blend on mask, falling back to maximum value near peak (which gets interpolated over anyways)
    blend = (vdata * mask + hdata * (1-mask)) * (1 - source_mask) + source_mask * np.maximum(vdata, hdata)

    # identify inner source area for interpolation
    replace_mask = angle < bw * rpeakr
    # isolate donut area used for interpolation
    src_mask = np.logical_and(angle > bw * cpeakr, outer_src_mask)
    lum = blend.reshape(3, -1).T[src_mask]

    # make light point from this donut of values for interpolation
    lp = LightPointKD(None, v[src_mask], lum, vm=vms, features=3, calcomega=False, write=False)
    # interpolate middle
    i, w = lp.interp(v[replace_mask], 40, lum=False, angle=False)

    # calculate slope towards center
    mp = (peakr + cpeakr) * bw / 2
    dp = (peakr - cpeakr) * bw / 2
    # interpolate between within src_mask
    ib = np.average(lum[angle[src_mask] < mp], axis=0)
    ob = np.average(lum[angle[src_mask] >= mp], axis=0)
    slope = (ib - ob) / dp

    maxangle = bw * cpeakr
    angle2 = np.minimum(angle[replace_mask][:, None], maxangle)
    # combine donut interpolation with slope towards middle (plus noise for blending
    clum = ((lp.apply_interp(i, lum, w) + slope * (maxangle - angle2)) *
            (1 + np.random.default_rng().normal(scale=.01, size=angle2.shape)))

    # replace center values with new interpolation
    blend[:, replace_mask.reshape(res, res)] = clum.T
    # isolate area around source to gather lens flare energy
    vm_valid = ViewMapper(pxyz, sfov).in_view(v, indices=False)
    flare = sdata.reshape(3, -1)[:, vm_valid] - blend.reshape(3, -1)[:, vm_valid]
    # discard under exposed pixels
    flare[flare < 0] = 0
    sol_luminance = np.sum(flare * omega[None, vm_valid], axis=1) / srcsize

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
