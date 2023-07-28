from subprocess import Popen, PIPE
import shlex
import os
from raytools import io
import numpy as np


def camera_raw_values(img):
    ppm = img + "_raw.ppm"
    hdr = img + "_raw.hdr"
    Popen(shlex.split(f"dcraw_emu -disinterp -4 -o 0 -w -Z {ppm} {img}")).communicate()
    f = Popen(shlex.split(f"pfsin {ppm}"), stdout=PIPE)
    Popen(shlex.split(f"pfsout {hdr}"), stdin=f.stdout, stderr=PIPE).communicate()
    imrgb = io.hdr2carray(hdr)
    imar = np.max(imrgb, axis=0)
    os.remove(ppm)
    os.remove(hdr)
    return imar