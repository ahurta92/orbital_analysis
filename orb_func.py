import numpy as np
import scipy.optimize as opt
import math
import matplotlib.pyplot as plt
import glob
import pandas as pd
from sklearn.linear_model import LinearRegression


def plot_single_orbital(df, row, col):
    data = df.at[row, col]
    r = data[:, 0]
    phi = data[:, 1]
    log_phi = np.log(np.abs(phi))

    r_index = np.argwhere((r > 16) & (r < 40))
    # we get the values in between
    fig = plt.plot(
        r[r_index],
        log_phi[r_index],
        "ro",
        markersize=0.04,
    )
    plt.show()
    return r, phi, fig


# return dataframe of densities for given direction
# order by states
# zero order has a single state
# all other order has num_states
def generate_density_labels(num_orders, num_states):
    cols = ["x", "y", "z"]
    rows = ["00"]
    for o in range(1, num_orders):
        for states in range(num_states):
            rows.append(str(o) + str(states))

    print(rows)
    return [rows, cols]


def read_density_plots(response_dir, num_orders, num_states, xoffset):
    plot_path = response_dir + "/plots/densities/"
    [idx, cols] = generate_density_labels(num_orders, num_states)
    densities_df = pd.DataFrame(index=idx, columns=cols)

    for o in idx:
        for di in cols:
            order = o[0]
            state = o[1]
            p_path = "rho" + order + "_" + di + "_" + state + ".plt"
            p_3 = plot_path + p_path
            print(p_3)
            densities_df.loc[o, di] = np.loadtxt(p_3, dtype=float)
            densities_df.loc[o,
                             di][:,
                                 0] = densities_df.loc[o, di][:, 0] - xoffset
    return densities_df


def generate_orbital_labels(num_orders, num_states, num_orbitals):
    cols = []
    directions = ["x", "y", "z"]
    xstate = ["x", "y"]
    for d in directions:
        for o in range(num_orbitals):
            cols.append(d + str(o))

    rows = ["000"]
    for o in range(1, num_orders):
        for states in range(num_states):
            for s in xstate:
                rows.append(s + str(o) + str(states))
    print(rows)
    print(cols)
    return [rows, cols]


def read_orbital_plots(response_dir, num_orders, num_states, num_orbitals,
                       xoffset):
    plot_path = response_dir + "/plots/orbitals/"
    [idx, cols] = generate_orbital_labels(num_orders, num_states, num_orbitals)
    orb_df = pd.DataFrame(index=idx, columns=cols)
    print(idx)
    print(cols)
    print(orb_df)

    for row in idx:
        for col in cols:
            orb_type = row[0]  # type 0 x y
            order = row[1]  # order
            state = row[2:]  # state 0 x y
            di = col[0]  # direction xyz
            orb = col[1:]  # orb number
            p_path = "phi" + orb_type + "_" + di + "_" + state + "_" + orb + ".plt"
            p_3 = plot_path + p_path
            print(p_3)
            orb_df.at[row, col] = np.loadtxt(p_3, dtype=float)
            orb_df.at[row, col][:, 0] = orb_df.at[row, col][:, 0] - xoffset

    return orb_df


def f(r, coeff, alpha, beta):
    return coeff * r**alpha * np.exp(-beta * r)
def logf(r,a,b,c):
    return a+b*np.log(r)-c*r



# TODO figure out fit to exp(-alpha r) alpha=sqrt(-2*e_i)
# exp(-alpha*r)/r^2
def compute_energies_density(df, num_orders, num_states, lo_val, hi_val,
                             center_p, plot_path) -> object:
    [idx, cols] = generate_density_labels(num_orders, num_states)
    num_idx = idx.__len__()
    num_cols = cols.__len__()
    print(num_orders)
    print(num_states)
    eng_df = pd.DataFrame(index=idx, columns=cols)
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
            data = df.at[o, di]
            r = data[:, 0]
            phi = data[:, 1]
            absphi = np.abs(phi)
            log_phi = np.log10(np.abs(phi))
            # check if log(abs(phi)) is > eps
            if absphi.max() > hi_val:
                # first point where value is greater than low
                lo_idx = np.argwhere(lo_val < absphi)[0][0]
                hi_idx = np.argwhere(hi_val < absphi)[0][0]
                # we get the values in between
                r_far = np.array(r[lo_idx:hi_idx])
                print("r_lo ", r_far[0], "r_hi", r_far[-1], "rfar.size() ",
                      r_far.size)
                # r_far = r_far.reshape(-1, 1)
                phi_far = np.array(phi[lo_idx:hi_idx])
                log_phi_far=np.log(np.abs(phi_far))
                # phi_far = phi_far.reshape(-1, 1)
                popt, pcov = opt.curve_fit(logf, r_far, log_phi_far, [0.0, 0.0, 0.0])
                #popt, pcov = opt.curve_fit(f, r_far, phi_far, [1.0, 1.0, 1.0])

                coeff, alpha, beta = popt
                coeff=np.exp(-coeff)

                print("coeff", coeff, "alpha", alpha, "beta", beta)
                e_i = -(beta**2) / 2

                r_index = np.argwhere((r > r[lo_idx]) & (r < center_p))
                # get predicted values
                fit = f(r_far, coeff, alpha, beta)
                logfit = np.log10(np.abs(fit))
                # plot log_phi
                plt.plot(r[r_index],
                         log_phi[r_index],
                         "ro",
                         markersize=mark_size)
                plt.plot(r_far, logfit, "bo", markersize=mark_size)
                plt.ylabel(y_label)
                plt.xlabel(x_label)
                # plt.text(r[lo_idx], hi_val + 1, "Score = " + str(round(score, 2)))
                plt.text(r[lo_idx], lo_val + 1, "E = " + str(round(e_i, 3)))
                # plt.plot(r, 1/r,'go')
                eng_df.loc[o, di] = e_i
            else:
                eng_df.loc[o, di] = 0
                plt.plot(r, log_phi, "ro", markersize=mark_size)
                plt.ylabel(y_label)
                plt.xlabel(x_label)

            plt_count += 1
    plt.savefig(plot_path + "density_plots" + ".png")
    plt.show()
    return eng_df


def generate_state_orbital_labels(num_orders, state, num_orbitals):
    cols = []
    directions = ["x", "y", "z"]
    for d in directions:
        for o in range(num_orbitals):
            cols.append(d + str(o))

    rows = ["000"]
    xstate = ["x", "y"]
    for o in range(1, num_orders):
        for s in xstate:
            rows.append(s + str(o) + str(state))
    print(rows)
    print(cols)
    return [rows, cols]


# TODO figure out fit to exp(-alpha r) alpha=sqrt(-2*e_i)
# exp(-alpha*r)/r^2
def compute_orbital_energies(df, num_orders, num_states, num_orbitals, lo_val,
                             hi_val, center_p, plot_path) -> object:
    labels = ["log_phi", "log_phi_predict"]
    mark_size = 0.15
    plt_count = 1
    o_e = []
    for s in range(num_states):
        [rows, cols] = generate_state_orbital_labels(num_orders, s,
                                                     num_orbitals)
        fig, axs = plt.subplots(nrows=3,
                                ncols=3 * num_orbitals,
                                sharex=False,
                                sharey=False)
        eng_df = pd.DataFrame(index=rows, columns=cols)
        i = 0
        for row in rows:
            j = 0
            otype = row[0]
            order = row[1]
            state = row[2:]
            for c in cols:
                direc = c[0]
                orb_num = c[1:]
                x_label = c
                y_label = row
                data = df.at[row, c]

                r = data[:, 0]
                phi = data[:, 1]
                absphi = np.abs(phi)
                log_phi = np.log(np.abs(phi))
                log10phi = np.log10(np.abs(phi))
                # check if log(abs(phi)) is > eps
                print(absphi.max())
                if absphi.max() >.00005 :
                    # first point where value is greater than low
                    hi_idx = np.argwhere(absphi > lo_val)[-1][0]
                    lo_idx = np.argwhere(absphi > hi_val)[-1][0]
                    r_far = r[lo_idx:hi_idx]
                    # r_far = r_far.reshape(-1, 1)
                    phi_far = phi[lo_idx:hi_idx]
                    log_phi_far =np.log(np.abs(phi_far))
                    # phi_far = phi_far.reshape(-1, 1)
                    popt, pcov = opt.curve_fit(logf, r_far, log_phi_far,
                                               [1.0, 0.0, 0.0])

                    coeff, alpha, beta = popt
                    coeff=np.exp(coeff)

                    print("coeff", coeff, "alpha", alpha, "beta", beta)
                    e_i = -(beta**2) / 2
                    print("row", row, "col", c, "e_i", e_i)

                    # we get the values in between
                    r_index = np.argwhere((r > center_p) & (r < r[hi_idx]))
                    fit = f(r_far, coeff, alpha, beta)
                    logfit = np.log10(np.abs(fit))
                    # get regression score and energy
                    # plot log_phi
                    axs[i, j].plot(r[r_index],
                                   log10phi[r_index],
                                   "ro",
                                   markersize=mark_size)
                    axs[i, j].plot(r_far, logfit, "bo", markersize=mark_size)
                    # plt.text(r[lo_idx], hi_val + 1, "Score = " + str(round(score, 2)))
                    estring = str(round(e_i, 3))
                    xpos = 0
                    ypos = np.log10(hi_val)
                    axs[i, j].text(xpos,
                                   ypos,
                                   estring,
                                   weight="bold",
                                   color="midnightblue")
                    eng_df.loc[row,c] = e_i
                else:
                    lo_idx = np.argwhere(r > 20)[0][0]
                    r_index = np.argwhere((r > 20) & (r < center_p))
                    # we get the values in between
                    axs[i, j].plot(
                        r[r_index],
                        np.log10(0.5 * np.ones(r[r_index].size)),
                        "ro",
                        markersize=0.04,
                    )
                j += 1
            i += 1
        fig.suptitle("State " + str(s) + " Orbitals")
        fig.savefig("state_" + str(s) + "orbs.svg", )
        fig.show()
        o_e.append(eng_df)
    return o_e

    # fig.supxlabel(cols)
    # fig.supylabel(idx)
