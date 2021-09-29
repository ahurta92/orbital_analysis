#!/usr/bin/env python
# coding: utf-8
from orb_func import read_orb_plots
from orb_func import read_ground_orb_plots
from orb_func import get_num_states
from orb_func import compute_energies
from read_output import harvest_outfile_pass
import glob

molecule = "molecules/10_Be"
print(molecule)
# test = "dynamic_p1"
# test = "dynamic_p2"o
# o_0000=0.000
# o_0025=0.025
frequency = "o_0000"
response_dir = molecule + "/" + frequency
tests = glob.glob(molecule + "/o_*")

response_file = response_dir + "/response.out"
print(response_file)

psivar = harvest_outfile_pass(response_file)
o_e = psivar["ORBITAL ENERGIES"]

# ground plots
[num_states, num_orbs] = get_num_states(response_dir)
# g_x = read_ground_orb_plots(num_orbs, response_dir, "x")
# g_y = read_ground_orb_plots(num_orbs, response_dir, "y")
# g_z = read_ground_orb_plots(num_orbs, response_dir, "z")

# direction x (1 state m oritals
# x_x = read_orb_plots("x", response_dir, "x", num_states, num_orbs)
# y_x = read_orb_plots("y", response_dir, "x", num_states, num_orbs)

# direction y
# x_y = read_orb_plots("x", response_dir, "y", num_states, num_orbs)
# y_y = read_orb_plots("y", response_dir, "y", num_states, num_orbs)
# direction z
# y_z = read_orb_plots("y", response_dir, "z", num_states, num_orbs)
# x_z = read_orb_plots("x", response_dir, "z", num_states, num_orbs)

# plot_path = response_dir + "/plots"

lo_val = 0.02
hi_val = 0.025
p_lo = 0
p_hi = 24

# e_xx_df = compute_orbital_energies(x_x, lo, hi, plot_path)
# e_yx_df = compute_orbital_energies(y_x, lo, hi, plot_path)
# e_gx_df = compute_energies(g_x, lo_val, hi_val, p_hi, plot_path, o_e)
# e_xx_df = compute_energies(x_x, lo_val, hi_val, p_hi, plot_path, o_e)
# e_xy_df = compute_energies(x_y, lo_val, hi_val, p_hi, plot_path, o_e)
# e_xz_df = compute_energies(x_z, lo_val, hi_val, p_hi, plot_path, o_e)
