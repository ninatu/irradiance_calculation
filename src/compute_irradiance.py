# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import argparse
import os


def get_coords(side, count, angle=0):
    coords = np.empty(count)

    d_side = side / count
    for i in range(count):
        coords[i] = (- side / 2 + d_side / 2 + i * d_side) * np.cos(angle)
    return coords


def create_scene_figure(points_plat, points_light, max_X):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x, y, z = [points_plat[:, :, i].flatten() for i in range(3)]
    ax.scatter(x, y, z, c='r', marker='o', s=1)

    x, y, z = [points_light[:, :, i].flatten() for i in range(3)]
    ax.scatter(x, y, z, c='b', marker='o', s=1)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlim(-max_X, max_X)
    ax.set_ylim(-max_X, max_X)
    return fig


def create_irradiance_figure(irradiance, max_X):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Y')
    ax.set_ylabel('X')
    ax.set_xlim(-max_X, max_X)
    ax.set_ylim(-max_X, max_X)
    ax.imshow(irradiance, cmap='gray', extent=(-max_X, max_X, -max_X, max_X))
    return fig


def print_irradiance_to_file(irradiance, path):
    with open(path, 'w') as fout:
        for row in irradiance:
            for elem in row:
                fout.write("{:.2f} ".format(elem))
            fout.write("\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Calculation of the irradiance for a flat lambert light")
    parser.add_argument('--F', type=float, help='Total radiant flux (Watt)')
    parser.add_argument('--side_plat', type=float, help='Size of platform side (meter)')
    parser.add_argument('--side_light', type=float, help='Size of light source side (meter)')
    parser.add_argument('--angle_light', type=float, help='Angle of inclination of the light source (degree)')
    parser.add_argument('--n_plat', type=int, help='Number of cells of platform')
    parser.add_argument('--n_light', type=int, help='Number of cells of light source')
    parser.add_argument('--path_V_lambda', type=str, help='Path to V(lambda) curve')
    parser.add_argument('--path_RSPD', type=str, help='Path to relative SPD')
    parser.add_argument('--output_dir', type=str, help='Path directory to save results')
    args = parser.parse_args()

    lamb_min = 390
    lamb_max = 831
    d_lamb = 5

    side_light = args.side_light
    angle_light = args.angle_light

    side_plat = args.side_plat
    n_plat = args.n_plat
    n_light = args.n_light

    F = args.F

    path_V_lambda = args.path_V_lambda
    path_RSPD = args.path_RSPD
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    V_labmda = pd.read_excel(path_V_lambda, header=None, names=["x"], index_col=0)['x']
    rel_SPD = pd.read_excel(path_RSPD, header=None, names=["x"], index_col=0)["x"]

    V_labmda = np.array([V_labmda[i] for i in range(lamb_min, lamb_max, d_lamb)])
    rel_SPD = np.array([rel_SPD[i] for i in range(lamb_min, lamb_max, d_lamb)])

    # calculation of cell coordinates
    xs_plat = get_coords(side_plat, n_plat)
    ys_plat = get_coords(side_plat, n_plat)
    normal_plat = np.array([0, 0, 1])

    xs_light = get_coords(side_light, n_light)
    ys_light = get_coords(side_light, n_light, np.pi * (angle_light / 180.0))
    zs_light = side_plat / 2 + get_coords(side_light, n_light, np.pi * ((90 - angle_light) / 180.0))
    normal_light = np.array([0, np.sin(np.pi * angle_light / 180.0), -np.cos(np.sin(np.pi * angle_light / 180.0))])

    points_plat = np.zeros((n_plat, n_plat, 3))
    points_plat[:, :, 0] = xs_plat[..., np.newaxis]
    points_plat[:, :, 1] = ys_plat[np.newaxis, ...]

    points_light = np.zeros((n_light, n_light, 3))
    points_light[:, :, 0] = xs_light[..., np.newaxis]
    points_light[:, :, 1] = ys_light[np.newaxis, ...]
    points_light[:, :, 2] = zs_light[np.newaxis, ...]

    # saving a scene image
    max_half_side = max(side_plat, side_light) / 2
    fig = create_scene_figure(points_plat, points_light, max_half_side)
    fig.savefig(os.path.join(output_dir, "scene.png"))

    # computation of radiant flux
    F_cell = F / (n_light * n_light)
    norm_constant = rel_SPD.sum() * d_lamb / F_cell
    F_lambda = rel_SPD / norm_constant
    # computation I_0 of lambert source
    I_o_lambda = F_lambda / (4 * np.pi)

    # computation vector from platform point to light point
    points_plat = points_plat.reshape((n_plat * n_plat, 3))
    points_light = points_light.reshape((n_light * n_light, 3))
    plat_to_light = points_light[np.newaxis, :, :] - points_plat[:, np.newaxis, :]
    del points_plat, points_light

    # computation distance between platform point and light point
    r = np.linalg.norm(plat_to_light, axis=2)

    # cosine of angle between normal of surface of light source and direction to platform
    cos_phi = (- plat_to_light * normal_light[np.newaxis, np.newaxis, :]).sum(axis=2) / r

    # cosine of angle between normal of platform and direction to light
    cos_alpha = (plat_to_light * normal_plat[np.newaxis, np.newaxis, :]).sum(axis=2) / r

    del plat_to_light

    # going over all wavelengths
    E_lamb = I_o_lambda[np.newaxis, np.newaxis, :] * (cos_phi * cos_alpha / r ** 2)[..., np.newaxis]
    del cos_phi, cos_alpha, r
    irradiance = 683 * (E_lamb * V_labmda[np.newaxis, np.newaxis, :] * d_lamb).sum(axis=(1, 2))
    irradiance = irradiance.reshape((n_plat, n_plat))

    fig_irradiance = create_irradiance_figure(irradiance, max_half_side)
    fig_irradiance.savefig(os.path.join(output_dir, "irradiance.png"))
    print_irradiance_to_file(irradiance, os.path.join(output_dir, "irradiance.txt"))
