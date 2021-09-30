import numpy as np
import matplotlib.pyplot as plt
import glob
import pandas as pd
from sklearn.linear_model import LinearRegression


def read_ground_orb_plots(num_orbitals, test, di):
    plot_path = test + "/plots/ground/"
    o_dict = {}

    cols = []
    for j in range(num_orbitals):
        cols.append(str(j))
    for i in range(1):
        p_1 = "ground_" + di + "_" + str(i)
        o_dict[p_1] = []
        for j in range(num_orbitals):
            p_2 = "ground_" + di + "_" + str(j)
            p_3 = plot_path + p_2 + ".plt"
            o_dict[p_1].append(np.genfromtxt(p_3, dtype=float))
    return pd.DataFrame.from_dict(o_dict, orient="index", columns=cols)


def read_ground_density_plots(test_dir):
    # test directory + /plots +/ground
    plot_path = test_dir + "/plots/ground/"
    o_dict = {}

    cols = []
    directions = "xyz"
    for d in directions:
        p_1 = "ground_density_" + d
        o_dict[p_1] = []
        p_3 = plot_path + p_1 + ".plt"
        o_dict[p_1].append(np.genfromtxt(p_3, dtype=float))
    return pd.DataFrame.from_dict(o_dict, orient="index")


def read_transition_density_plot(test_dir, num_states, direction):
    # test directory + /plots +/ground
    plot_path = test_dir + "/plots/xy/"
    o_dict = {}

    cols = []
    for i in range(num_states):
        p_1 = "rho1_direction_" + direction + "_res_" + str(i)
        o_dict[p_1] = []
        p_3 = plot_path + p_1 + ".plt"
        o_dict[p_1].append(np.genfromtxt(p_3, dtype=float))
    return pd.DataFrame.from_dict(o_dict, orient="index")


def get_num_states(test):
    mol_path = test + "/plots/xy/"
    x_states = glob.glob(mol_path + "x*_x*_*")
    # print(x_states)
    num_states = 0
    num_orbitals = 0
    # find number for res states and orbitals
    for x_s in x_states:
        x_split = x_s.split("_")
        num_xs = [s for s in x_split if s.isdigit()]
        # print(num_xs)
        num_states = max(int(num_xs[0]), num_states)
        num_orbitals = max(int(num_xs[1]), num_orbitals)
    num_orbitals = num_orbitals + 1
    num_states = num_states + 1
    print("num states", num_states)
    print("num orbitals", num_orbitals)
    return [num_states, num_orbitals]


def read_orb_plots(xy, test, di, num_states, num_orbitals):
    mol_path = test + "/plots/xy/"
    # find number for res states and orbitals
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
            plot_df.loc[idx][col] = plot_func(r, np.abs(phi))
            title = idx + "_" + str(col)
            plt.title(title)
            plt.savefig(plot_path + title + ".png")
            plt.show()

    return plot_df


# TODO figure out fit to exp(-alpha r) alpha=sqrt(-2*e_i)
# exp(-alpha*r)/r^2
def compute_energies(df, lo_val, hi_val, center_p, plot_path, o_e) -> object:
    # these are the orbitals
    columns = list(df)
    # these are the states
    rows = df.index.copy()
    eng_df = pd.DataFrame(index=rows, columns=columns)

    lo_val = np.log(np.abs(lo_val))
    print(lo_val)
    hi_val = np.log(np.abs(hi_val))
    num_orbs = df.columns.size
    labels = ["log_phi", "log_phi_predict"]
    mark_size = 0.25
    for idx in rows:
        plt_count = 1
        for col in columns:
            plt.subplot(num_orbs, 1, plt_count)
            title = idx
            y_label = "orbital " + str(int(col) + 1)
            if plt_count == 1:
                plt.title(title)
                plt.legend(labels)
            data = df.loc[idx][col]
            r = data[:, 0]
            phi = data[:, 1]
            log_phi = np.log(np.abs(phi))
            # check if log(abs(phi)) is > eps
            if log_phi.max() > hi_val:
                # first point where value is greater than low
                lo_idx = np.argwhere(log_phi > lo_val)[0][0]
                hi_idx = np.argwhere(log_phi > hi_val)[0][0]
                # we get the values in between
                r_far = r[lo_idx:hi_idx]
                r_far = r_far.reshape(-1, 1)
                log_phi_i = log_phi[lo_idx:hi_idx]
                log_phi_i = log_phi_i.reshape(-1, 1)
                reg = LinearRegression().fit(r_far, log_phi_i)

                r_index = np.argwhere((r > r[lo_idx]) & (r < center_p))
                # get regression score and energy
                score = reg.score(r_far, log_phi_i)
                momentum = reg.coef_[0]
                Energy = -(momentum[0] ** 2) / 2

                print(
                    "Regression Score | orbital energy " + title + "_" + str(col),
                    " | ",
                    str(round(score, 2)),
                    " ",
                    str(round(Energy, 3)),
                )
                # get predicted values
                log_phi_predict = reg.predict(r_far)
                # plot log_phi
                plt.plot(r[r_index], log_phi[r_index], "ro", markersize=mark_size)
                plt.plot(r_far, log_phi_predict, "bo", markersize=mark_size)

                plt.ylabel(y_label)
                # plt.text(r[lo_idx], hi_val + 1, "Score = " + str(round(score, 2)))
                plt.text(r[lo_idx], lo_val + 1, "E = " + str(round(Energy, 3)))
                # plt.plot(r, 1/r,'go')
                plt.savefig(plot_path + title + ".png")
                eng_df.loc[idx][col] = reg
            else:
                eng_df.loc[idx][col] = 0
                plt.plot(r, log_phi, "ro", markersize=mark_size)
                plt.ylabel(y_label)
                plt.savefig(plot_path + title + ".png")

            plt_count += 1
        plt.show()
    return eng_df
