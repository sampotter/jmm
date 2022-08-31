#!/usr/bin/env python

import numpy as np
import sys

from PIL import Image

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f'usage: {sys.argv[0]} <BIN_PATH> <IMG_PATH>')
        exit(1)
    shape = (256, 256, 4)
    arr = np.fromfile(sys.argv[1]).reshape(shape)
    img = Image.fromarray(np.uint8(255*(arr/arr.max())))
    img.save(sys.argv[2])
