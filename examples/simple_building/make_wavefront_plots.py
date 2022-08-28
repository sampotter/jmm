# NOTE: some of this stuff is hardcoded!

import matplotlib.pyplot as plt
import numpy as np

plt.ion()

from pathlib import Path

m, n = 256, 256

path = Path('../../build/examples/simple_building')

T = [
    4.37777,
    5.57171,
    6.76565,
    7.95959,
    9.15352,
    10.3475,
    11.5414,
    12.7353,
    13.9293,
    15.1232,
    16.3172,
    17.5111,
]

plt.figure(figsize=(8, 8))
img = np.fromfile(path/f'image0003.bin').reshape(m, n, 4)
plt.imshow(img, interpolation='antialiased')
plt.axis('off')
plt.gca().set_aspect('equal')
plt.tight_layout()
plt.show()

plt.figure(figsize=(8, 12))
for i in range(12):
    img = np.fromfile(path/f'image{i:04d}.bin').reshape(m, n, 4)
    plt.subplot(3, 4, i + 1)
    plt.imshow(img, interpolation='antialiased')
    plt.axis('off')
    # plt.tick_params(left=False,
    #                 bottom=False,
    #                 labelleft=False,
    #                 labelbottom=False)
    plt.gca().set_aspect('equal')
    plt.title(rf'$T = {T[i]/340:0.3g}$s')
plt.tight_layout()
plt.show()
