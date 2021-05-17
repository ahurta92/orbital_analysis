import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import pandas as pd


def read_orb_plots(xy, molecule, di):
    mol_path = molecule + "/"
    x_states = glob.glob(mol_path + "x*_x*_*")
    num_states = 0
    num_orbitals = 0
    # find number for res states and orbitals
    for x_s in x_states:
        x_split = x_s.split("_")
        num_states = max(int(x_split[4]), num_states)
        num_orbitals = max(int(x_split[6]), num_orbitals)
    print("num states", num_states)
    print("num orbitals", num_orbitals)
    cols = []
    for j in range(num_orbitals):
        cols.append(str(j))
    o_dict = {}
    for i in range(num_states):
        p_1 = xy + "_state_" + str(i)
        o_dict[p_1] = []
        for j in range(num_orbitals):
            p_2 = (
                mol_path + xy + "_direction_" + di + "_res_" + str(i) + "_orb_" + str(j)
            )
            o_dict[p_1].append(np.genfromtxt(p_2, dtype=float))
    return pd.DataFrame.from_dict(o_dict, orient="index", columns=cols)


def make_plots(df, plot_func, low_range, high_range, plot_path):

    columns = list(df)
    rows = df.index.copy()
    plot_df = pd.DataFrame(index=rows, columns=columns)

    for idx in rows:
        print("idx", idx)
        for col in columns:
            data = df.loc[idx][col]
            r = data[:, 0]
            phi = data[:, 1]
            r_index = np.argwhere((r > low_range) & (r < high_range))
            r = r[r_index]
            phi = phi[r_index]
            plot_df.loc[idx][col] = plot_func(r, phi)
            title = idx + "_" + str(col)
            plt.title(title)
            plt.savefig(plot_path + title + ".png")
            plt.show()

    return plot_df
