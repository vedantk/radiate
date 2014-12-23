#!/usr/bin/python2

import argparse
import random
import itertools
from math import e, log, sqrt

import cv2
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay

desc = '''Generate a 3d model out of a face'''

def zlog(x):
    return log(2*e + -(0.5 - x)**2)

def make_simulacrum(img, outpath, plot):
    h, w = map(float, img.shape[:2])

    R = 0.35               # Target face radius (normalized).
    pr_keep_face = 0.005   # Probability of keeping face vertices.
    pr_keep_bkgd = 0.001   # Probability of keeping background vertices.

    xs, ys, zs = [], [], []
    within_face = lambda r, c: (r/h - 0.5)**2 + (c/w - 0.5)**2 < R**2
    depth_noise = lambda: (random.random() - 0.5) / 750.0
    on_border = lambda r, c: (r == 0 or r == h-1) or (c == 0 or c == w-1)
    z = lambda r, c: sqrt(zlog(r/h) + zlog(c/w)) + depth_noise()
    z0 = z(0, 0)

    for (r, c) in itertools.product(range(int(h)), range(int(w))):
        pr = random.random()
        if within_face(r, c):
            if pr > pr_keep_face:
                continue
        else:
            if pr > pr_keep_bkgd and not on_border(r, c):
                continue

        xs.append(r)
        ys.append(c)

        if on_border(r, c):
            zs.append(z0)
        else:
            zs.append(z(r, c))

    D = Delaunay(zip(xs, ys))

    xs = xs + xs
    ys = ys + ys
    zs = zs + map(lambda x: z0 - (x - z0), zs)

    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(ys, xs, zs, c='g', marker='+')
        plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('img')
    parser.add_argument('out')
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()

    img = cv2.imread(args.img)
    make_simulacrum(img, args.out, args.plot)
