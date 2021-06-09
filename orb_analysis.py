#!/usr/bin/env python
# coding: utf-8
from orb_func import read_orb_plots
from orb_func import read_ground_orb_plots
from orb_func import get_num_states
from orb_func import compute_energies
from orb_func import make_plots
import matplotlib.pyplot as plt

test = "static"
#test = "dynamic_p1"
#test = "dynamic_p2"

# ground plots
[num_states, num_orbs] = get_num_states(test)
g_x = read_ground_orb_plots(num_orbs, test, "x")

# direction x
x_x = read_orb_plots("x", test, "x", num_states, num_orbs)
y_x = read_orb_plots("y", test, "x", num_states, num_orbs)

# direction y
x_y = read_orb_plots("x", test, "y", num_states, num_orbs)
y_y = read_orb_plots("y", test, "y", num_states, num_orbs)
# direction z
y_z = read_orb_plots("y", test, "z", num_states, num_orbs)
x_z = read_orb_plots("x", test, "z", num_states, num_orbs)

plot_path = "plots/" + test + "/"
lo_val = 0.0005
hi_val = 0.0009
p_lo=0
p_hi=24

#plot_x_x = make_plots(x_x, plt.plot, p_lo, p_hi, plot_path)
#plot_y_x = make_plots(y_x, plt.plot, p_lo, p_hi,plot_path)

# e_xx_df = compute_orbital_energies(x_x, lo, hi, plot_path)
# e_yx_df = compute_orbital_energies(y_x, lo, hi, plot_path)
e_gx_df = compute_energies(g_x, lo_val, hi_val, p_hi, plot_path)
e_xx_df = compute_energies(x_z, lo_val, hi_val,  p_hi, plot_path)
g_x0=g_x.loc['ground_x_0'][0]