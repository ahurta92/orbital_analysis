import numpy as np
import matplotlib.pyplot as plt
import glob
import pandas as pd
from sklearn.linear_model import LinearRegression


# return dataframe of densities for given direction
# order by states
# zero order has a single state
# all other order has num_states
#
def read_density_plots(response_dir, num_orders, num_states):
    plot_path = response_dir + "/plots/densities/"

    cols = ["x", "y", "z"]
    idx = ["00"]
    for o in range(1, num_orders):
        for states in range(num_states):
            idx.append(str(o) + str(states))

    print(idx)
    densities_df = pd.DataFrame(index=idx, columns=cols)

    for o in idx:
        for di in cols:
            order = o[0]
            state = o[1]
            p_path = "rho" + order + "_" + di + "_" + state + ".plt"
            p_3 = plot_path + p_path
            print(p_3)
            densities_df.loc[o, di] = np.genfromtxt(p_3, dtype=float)
    return densities_df


# TODO figure out fit to exp(-alpha r) alpha=sqrt(-2*e_i)
# exp(-alpha*r)/r^2
def compute_energies_density(df, num_orders, num_states, lo_val, hi_val,
                             center_p, plot_path) -> object:

    cols = ["x", "y", "z"]
    idx = ["00"]
    for o in range(1, num_orders):
        for states in range(num_states):
            idx.append(str(o) + str(states))

    num_idx = idx.__len__()
    num_cols = cols.__len__()
    print(num_orders)
    print(num_states)
    print(idx)

    eng_df = pd.DataFrame(index=idx, columns=cols)

    lo_val = np.log(np.abs(lo_val))
    hi_val = np.log(np.abs(hi_val))
    labels = ["log_phi", "log_phi_predict"]
    mark_size = 0.25
    plt_count = 1
    plt.subplot(111)
    plt.title("density plots")

    for o in idx:
        order = o[0]
        state = o[1]
        for di in cols:

            plt.subplot(num_idx, num_cols, plt_count)
            title = "rho" + state + "_" + di + "_" + state
            y_label = o
            x_label = di
            if plt_count == 1:
                plt.legend(labels)
            data = df.loc[o, di]
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
                Energy = -(momentum[0]**2) / 2

                print(
                    "Regression Score | density energy " + title + "_" +
                    str(state),
                    " | ",
                    str(round(score, 2)),
                    " ",
                    str(round(Energy, 3)),
                )
                # get predicted values
                log_phi_predict = reg.predict(r_far)
                # plot log_phi
                plt.plot(r[r_index],
                         log_phi[r_index],
                         "ro",
                         markersize=mark_size)
                plt.plot(r_far, log_phi_predict, "bo", markersize=mark_size)
                plt.ylabel(y_label)
                plt.xlabel(x_label)
                # plt.text(r[lo_idx], hi_val + 1, "Score = " + str(round(score, 2)))
                plt.text(r[lo_idx], lo_val + 1, "E = " + str(round(Energy, 3)))
                # plt.plot(r, 1/r,'go')
                plt.savefig(plot_path + title + ".png")
                eng_df.loc[o, di] = reg
            else:
                eng_df.loc[o, di] = 0
                plt.plot(r, log_phi, "ro", markersize=mark_size)
                plt.ylabel(y_label)
                plt.xlabel(x_label)

            plt_count += 1
    plt.savefig(plot_path + "density_plots" + ".png")
    plt.show()
    return eng_df


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
            p_2 = (mol_path + xy + "_direction_" + di + "_res_" + str(i) +
                   "_orb_" + str(j))
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
