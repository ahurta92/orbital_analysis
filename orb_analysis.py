#!/usr/bin/env python
# coding: utf-8
import json
from orb_func import compute_orbital_energies
from orb_func import read_density_plots
from orb_func import read_orbital_plots
from read_madness_json import read_molrespone_json

root_dir = "/home/adrianhurtado/projects/madness-test-suite/tests_response/orbital_analysis"
molecule = "/10_Be"

working_dir = root_dir + molecule

response_dir = working_dir + "/excited_state"
response_json = response_dir + "/response.json"

with open(response_json, "r") as json_file:
    response_j = json.load(json_file)
    r_params, proto_data = read_molrespone_json(response_j)

num_states = r_params["states"]
num_orbitals = r_params["num_orbitals"]
num_orders = 2

# [num_states, num_orbs] = get_num_states(response_dir)
x_offset = 24.0
rho = read_density_plots(response_dir, num_orders, num_states, x_offset)
phi = read_orbital_plots(response_dir, num_orders, num_states, num_orbitals, x_offset)
plot_path = response_dir + "/plots"

lo_val = 1e-6
hi_val = 8e-6
p_lo = -24
p_hi = 0
center_p = 0

#                                   0, plot_path)
e_orbs = compute_orbital_energies(phi, num_orders, num_states, num_orbitals, lo_val, hi_val,
                                  center_p, plot_path)
