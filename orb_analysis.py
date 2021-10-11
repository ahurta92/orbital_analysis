#!/usr/bin/env python
# coding: utf-8
# from orb_func import read_orb_plots
from orb_func import read_density_plots
from orb_func import compute_energies_density
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
num_states = int(psivar["NUM STATES"])
num_orbitals = int(psivar["NUM ORBITALS"])
num_orders=2

# ground plots
# [num_states, num_orbs] = get_num_states(response_dir)


rho_x = read_density_plots(response_dir,num_orders , num_states)
plot_path = response_dir + "/plots"

lo_val = 0.02
hi_val = 0.025
p_lo = 0
p_hi = 24

e_rho_x = compute_energies_density(rho_x, num_orders,num_states,lo_val, hi_val, p_hi, plot_path)
# e_xx_df = compute_orbital_energies(x_x, lo, hi, plot_path)
# e_yx_df = compute_orbital_energies(y_x, lo, hi, plot_path)
# e_gx_df = compute_energies(g_x, lo_val, hi_val, p_hi, plot_path, o_e)
# e_xx_df = compute_energies(x_x, lo_val, hi_val, p_hi, plot_path, o_e)
# e_xy_df = compute_energies(x_y, lo_val, hi_val, p_hi, plot_path, o_e)
# e_xz_df = compute_energies(x_z, lo_val, hi_val, p_hi, plot_path, o_e)
