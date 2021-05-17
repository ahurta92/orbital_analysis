#!/usr/bin/env python
# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import pandas as pd
from orb_func import read_orb_plots
from orb_func import make_plots

molecule = "N2"

x_x = read_orb_plots("x", "Ne", "x")
x_y = read_orb_plots("x", "Ne", "y")
x_z = read_orb_plots("x", "Ne", "z")

y_x = read_orb_plots("y", "Ne", "x")
y_y = read_orb_plots("y", "Ne", "y")
y_z = read_orb_plots("y", "Ne", "z")

plot_path = "plots/" + molecule + "/"
plot_x_x = make_plots(x_x, plt.plot, 25, 40, plot_path)
plot_y_x = make_plots(y_x, plt.plot, 25, 40, plot_path)
