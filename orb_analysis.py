#!/usr/bin/env python
# coding: utf-8
import json
import matplotlib.pyplot as plt
import os.path
import numpy as np
import itertools
from cycler import cycler

from orb_func import compute_energies, generate_density_labels, generate_orbital_labels
from orb_func import read_density_plots
from orb_func import read_orbital_plots
import scipy.optimize as opt

from orbital_asymptotics import orbital_scrubbing, energies_scrubbing, compute_alpha_nl, compute_b_nl, \
    compute_asymptotic, plot_ground_asymptotics, plot_response_asymptotics, compute_fits
from read_madness_json import read_molrespone_json, tensor_to_numpy

response_dir = os.path.curdir
moldft_dir = response_dir + "/.."
mol_name = os.path.realpath(moldft_dir).split("/")[-1]
moldft_json = moldft_dir + "/calc_info.json"
with open(moldft_json, "r") as json_file:
    moldft_j = json.load(json_file)

response_json = "response_base.json"

with open(response_json, "r") as json_file:
    response_j = json.load(json_file)
    r_params, proto_data = read_molrespone_json(response_j)

num_states = r_params["states"]
num_orbitals = r_params["num_orbitals"]
num_orders = 2

orb_energies = tensor_to_numpy(moldft_j["scf_eigenvalues_a"])

# [num_states, num_orbs] = get_num_states(response_dir)
x_offset = 24.0
rho = read_density_plots(response_dir, num_orders, num_states, x_offset)
phi = read_orbital_plots(response_dir, num_orders, num_states, num_orbitals, x_offset)

[c_labels, r_labels] = generate_orbital_labels(2, num_states, num_orbitals)

plot_path = response_dir + "/plots"

l_angular = [0, 0, 1, 0, 1, ]
orb_labels = {1: '1s', 2: '2s', 5: '2p', 6: '3s', 9: '3p'}
shell_sizes = {'1s': 1, '2s': 2, '2p': 3, '3s': 4, '3p': 5, '4s': 6}
size = shell_sizes[orb_labels[num_orbitals]]

num_orb_i = [1, 1, 3, 1, 3, 1, 5]
[x_df, y_df, z_df] = orbital_scrubbing(phi, c_labels, num_orbitals)
energies = energies_scrubbing(size, orb_energies)
ground_data = x_df["000"]
r = ground_data[:, 0]
orb_data = ground_data[:, 1:]
[alpha_nl, beta_nl] = compute_alpha_nl(size, energies)
rin = 5
rout = 7
radius = 1
[b_nl, b_nl_opt] = compute_b_nl(r, orb_data, alpha_nl, beta_nl, size, rin, rout, radius)

plot_lo = 1
plot_hi = 8
plot_ground_asymptotics(r, orb_data, size, b_nl, alpha_nl, beta_nl, plot_lo, plot_hi, mol_name)

x0_data = x_df['x10'][:, 1:]
omega = response_j['response_parameters']['omega']
rin = 3.0
rout = 8
radius = 1
[alpha_nl_new, beta_nl] = compute_alpha_nl(size, energies + omega, shift=1)
[b_nl, b_nl_opt] = compute_b_nl(r, x0_data, alpha_nl, beta_nl, size, rin + 3, rout + 10, radius)

#plot_response_asymptotics(r, x0_data, size, b_nl, alpha_nl, beta_nl, plot_lo, 20, mol_name + '-x0')

fit_p = [2.5, 4, 5, 8, 11]
rad_p = [3.0, 1, 1, 6, 3]
fits = compute_fits(r, orb_data, alpha_nl, beta_nl, size, fit_p, rad_p)

plot_ground_asymptotics(r, orb_data, size, fits[:, 0], fits[:, 1], beta_nl, 0.75, 30, mol_name + '-fit')
print("_______________response_____________________")
fit_r = [1.5, 4, 8, 14, 18]
rad_r = [5.0, 3, 3, 4, 4]
r_fits = compute_fits(r, x0_data, alpha_nl, beta_nl, size, fit_r, rad_r)
plot_response_asymptotics(r, x0_data, size, r_fits[:, 0], r_fits[:, 1], beta_nl, 0.75, 30, mol_name + '-fit')

