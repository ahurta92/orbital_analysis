import json
import matplotlib.pyplot as plt
import os.path
import pandas as pd
import numpy as np
import scipy.optimize as opt

# We are going to sphericaly average the orbitals with l>0
# For this reason we will only include a single averaged value each orbital type
# 1s 2s 2p 3s 3p 4s 3d 3p...
# 0  0  1  0  1  0  2  1
#
# Only capable of working with Ar
l_angular = [0, 0, 1, 0, 1, ]
r_orbs = ['1s', '2s', '2p', '3s', '3p', '4s', '3d', '4p']
orb_labels = {1: '1s', 2: '2s', 5: '2p', 6: '3s', 9: '3p'}
shell_sizes = {'1s': 1, '2s': 2, '2p': 3, '3s': 4, '3p': 5, '4s': 6}

num_orb_i = [1, 1, 3, 1, 3, 1, 5]


# take df and average out the
def orbital_scrubbing(psi_df: pd.DataFrame, c_labels: list[str], num_orbitals: int):
    size = shell_sizes[orb_labels[num_orbitals]]
    x_dict = {}
    y_dict = {}
    z_dict = {}
    for c in c_labels:
        # grab the x y and z data
        x_series = psi_df[c][0:num_orbitals]
        y_series = psi_df[c][num_orbitals:(2 * num_orbitals)]
        z_series = psi_df[c][2 * num_orbitals:(3 * num_orbitals)]
        # make the size Npts x num_orbitals+1
        npts = x_series[0].shape[0]
        x_array = np.zeros((npts, size + 1))
        y_array = np.zeros((npts, size + 1))
        z_array = np.zeros((npts, size + 1))
        r = x_series[0][:, 0]
        # get the r value
        x_array[:, 0] = r
        y_array[:, 0] = r
        z_array[:, 0] = r
        # get the values for the array
        j = 1
        m = 0

        while j < size + 1:
            start_index = m
            print("start index", start_index)
            for jj in range(num_orb_i[j - 1]):
                i = jj + start_index
                print("i", i)
                x_array[:, j] += x_series[i][:, 1] ** 2
                y_array[:, j] += y_series[i][:, 1] ** 2
                z_array[:, j] += z_series[i][:, 1] ** 2
                m += 1
            x_array[:, j] = np.sqrt(x_array[:, j])
            y_array[:, j] = np.sqrt(y_array[:, j])
            z_array[:, j] = np.sqrt(z_array[:, j])
            j += 1
        x_dict[c] = x_array
        y_dict[c] = y_array
        z_dict[c] = z_array
    return [x_dict, y_dict, z_dict]


def energies_scrubbing(size, energies):
    j = 0
    m = 0
    averaged_energies = np.zeros(size)
    while j < size:
        start_index = m
        for jj in range(num_orb_i[j]):
            i = jj + start_index
            averaged_energies[j] += energies[i]
            m += 1
        averaged_energies[j] /= num_orb_i[j]

        j += 1
    return averaged_energies


# assumes that p orbitals have been averaged
# only working with closed shell atoms for this reason.
def compute_alpha_nl(num_orbitals, energies):
    # First step is to compute alpha_h
    beta_h = np.sqrt(-2 * energies[num_orbitals - 1])
    alpha_h = 1 / beta_h - 1
    # now we need to compute alpha_nl based on the value of l
    l_homo = l_angular[num_orbitals - 1]
    homo_label = r_orbs[num_orbitals - 1]
    homo_type = homo_label[1]
    homo_level = homo_label[0]  # get the level of the valence
    alpha_nl = []
    beta_nl = []
    # if atom only consists of s orbitals then alpha_nl = 1/beta_nl
    # beta_nl=sqrt(-2 e_nl)
    if homo_label == '1s' or homo_label == '2s':
        for oj in range(num_orbitals):
            beta_nli = np.sqrt(-2 * energies[oj])
            alpha_i = (1 / beta_nli) - 1
            alpha_nl.append(alpha_i)
            beta_nl.append(beta_nli)
    else:
        for oj in range(num_orbitals - 1):
            l_orbital = l_angular[oj]
            if l_orbital != l_homo:
                alpha_nl.append(alpha_h - np.abs(l_orbital - l_homo) - 1)
            # if the l_homo is 0
            else:
                if l_homo == 0:
                    alpha_nl.append(
                        alpha_h - 2 * (1 + 1))  # 1 should be lmin the lowest nonzer orbital amoung occupied orbitals
                else:
                    alpha_nl.append(alpha_h - 3)
            # each orbital should have the same value of beta
            beta_nl.append(beta_h)
        alpha_nl.append(alpha_h)
        beta_nl.append(beta_h)
    return alpha_nl, beta_nl


def compute_b_nl(r, chi, alpha_nl, beta_nl, size, r_inner, r_outer, radius=0.05):
    print(chi.shape)
    homo_label = r_orbs[size - 1]  # get the label of the homo
    homo_type = homo_label[1]
    homo_level = homo_label[0]  # get the level of the valence
    b_nl = []
    b_nl_opt = []

    alpha = 0
    beta = 0

    def logchi(r, c):
        return c + alpha * np.log(r) - beta * r

    for k in range(size):
        orb_label = r_orbs[k]
        orb_type = orb_label[0]

        if orb_type == homo_type:
            rp = r_outer
        else:
            rp = r_inner
        r_lo = rp - radius
        r_hi = rp + radius
        r_index = np.argwhere((r > r_lo) & (r < r_hi))
        r_long = r[r_index].flatten()
        log_chi_far = np.log(np.abs(chi[r_index, k])).flatten()

        alpha = alpha_nl[k]
        beta = beta_nl[k]
        popt, pcov = opt.curve_fit(logchi, r_long, log_chi_far, [1.0])
        print(popt)

        b = chi[r_index, k] / r[r_index] ** alpha_nl[k] / np.exp(-beta_nl[k] * r_long)
        print(b)
        b_nl.append(b[0][0])
        b_nl_opt.append(popt[0])

    return b_nl, b_nl_opt


def compute_asymptotic(r, alpha_nl, beta_nl, b_nl, size):
    x = np.empty((r.shape[0], size))
    for i in range(size):
        x[:, i] = b_nl[i] * r ** alpha_nl[i] * np.exp(-beta_nl[i] * r)
    return x


def plot_ground_asymptotics(r, phi, size, alpha_nl, beta_nl, b_nl, plot_lo, plot_hi, mol_name):
    p_index = np.argwhere((r > plot_lo) & (r < plot_hi))
    r_plot = r[p_index].flatten()
    phi_plot = phi[p_index][:]
    sshape = (phi_plot.shape[0], phi_plot.shape[2])
    phi_plot = phi_plot.reshape(sshape)
    shell_markers = [".", "o", 'v', '*', '^', '1', '2', '3', '4'][0:size]
    colors = "rbygco"[0:size]
    x_asymp = compute_asymptotic(r_plot, alpha_nl, beta_nl, b_nl, size)
    logphi = np.log(np.abs(phi_plot))
    logxasmp = np.log(np.abs(x_asymp))

    orb_names = list(shell_sizes.keys())[0:size]
    for k in range(size):
        y = np.empty((logphi.shape[0], 2))
        y[:, 0] = logphi[:, k]
        y[:, 1] = logxasmp[:, k]
        print(y.shape)
        plt.plot(r_plot, y[:, 0],
                 label=orb_names[k], color=colors[k], marker=shell_markers[k], markevery=20, linestyle='None')
        plt.plot(r_plot, y[:, 1],
                 label=orb_names[k], color=colors[k], linestyle="--")

    plt.title("Asymptotic Region of Ground Orbitals ")
    plt.xlabel("r(a.u)")
    plt.ylabel("$\ln|\chi_{nl}|$")
    plt.legend()
    plt.savefig(mol_name + ".svg")
    plt.show()
