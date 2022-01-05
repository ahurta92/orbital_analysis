#!/usr/bin/env python
# coding: utf-8
# from orb_func import read_orb_plots
from orb_func import read_density_plots
from orb_func import plot_single_orbital
from orb_func import read_orbital_plots
from orb_func import compute_energies_density
from orb_func import compute_orbital_energies
from read_output import harvest_outfile_pass
import glob
import numpy as np
import scipy.optimize as opt
import math
import matplotlib.pyplot as plt
import glob
import pandas as pd
# molecule = "molecules/10_Be"
molecule = "molecules/11_Ne"
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
num_orders = 2

# ground plots
# [num_states, num_orbs] = get_num_states(response_dir)
xoffset = 24.0
rho = read_density_plots(response_dir, num_orders, num_states, xoffset)
phi = read_orbital_plots(response_dir, num_orders, num_states, num_orbitals, xoffset)
plot_path = response_dir + "/plots"

lo_val = 1e-6
hi_val = 8e-6
p_lo = -24
p_hi = 0
center_p = 0

#                                   0, plot_path)
e_orbs = compute_orbital_energies(phi, num_orders, num_states, num_orbitals, lo_val, hi_val,
                                  center_p, plot_path)
