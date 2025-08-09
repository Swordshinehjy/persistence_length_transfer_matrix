from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import eigvals
from scipy.integrate import cumulative_trapezoid, quad
from scipy.interpolate import interp1d
from scipy.spatial.transform import Rotation as R


def read_data(file_name):
    data = np.loadtxt(file_name)
    data = np.reshape(data, (-1, 2))
    if data[:, 0].max() - data[:, 0].min() != 360:
        # symmetric, only calculate half of the dihedral potential
        if data[:, 0].min() == 0 or data[:, 0].min(
        ) == 180:  # from 0 to 180 or 180 to 360
            mirrored = np.column_stack((-data[:, 0] + 360, data[:, 1]))
            combined = np.vstack((data, mirrored))
        elif data[:, 0].min() == -180:  # from -180 to 0
            mirrored = np.column_stack((-data[:, 0], data[:, 1]))
            combined = np.vstack((np.column_stack(
                (mirrored, data[:, 0] + 360, data[:, 1]))))
        combined = np.unique(combined, axis=0)
        return combined[np.argsort(combined[:, 0])]
    else:  # full dihedral potential
        return data[np.argsort(data[:, 0])]


def setup_interpolation(data_label, kTval):
    data = read_data(Path(f"{data_label['label']}.txt"))
    fitf = interp1d(data[:, 0],
                    data[:, 1],
                    kind='cubic',
                    fill_value="extrapolate")
    norm_val, _ = quad(lambda x: np.exp(-fitf(x) / kTval), 0, 360)
    x_values = np.linspace(0, 360, 1000)
    prob_vals = np.exp(-fitf(x_values) / kTval) / norm_val
    cum_dist = cumulative_trapezoid(prob_vals, x_values, initial=0)
    norm_cdf = cum_dist / cum_dist[-1]
    unique_cdf_vals, unique_indices = np.unique(norm_cdf, return_index=True)
    corresponding_x_vals = x_values[unique_indices]
    inv_cdf = interp1d(unique_cdf_vals,
                       corresponding_x_vals,
                       kind='cubic',
                       fill_value="extrapolate")
    return {
        'data': data,
        'fitf': fitf,
        'prob_vals': prob_vals,
        'x_values': x_values,
        'cum_dist': cum_dist
    }


def plot_dihedral_potentials(interp_data):
    """Plot dihedral potentials and their probability distributions."""
    plt.figure(figsize=(18, 5))
    plt.subplot(1, 3, 1)
    for key, data in interp_data.items():
        plt.plot(data['data'][:, 0],
                 data['data'][:, 1],
                 f"{data['color']}o",
                 label=data['label'])
        plt.plot(data['x_values'],
                 data['fitf'](data['x_values']),
                 f"{data['color']}--")
    format_subplot("Dihedral Angle [Deg.]", "Dihedral Potential (kJ/mol)",
                   "Dihedral Potentials")
    plt.subplot(1, 3, 2)
    for key, data in interp_data.items():
        plt.plot(data['x_values'],
                 data['prob_vals'],
                 f"{data['color']}-",
                 label=data['label'])
    format_subplot("Angle [deg.]", "Probability", "Probability Distributions")
    plt.subplot(1, 3, 3)
    for key, data in interp_data.items():
        plt.plot(data['cum_dist'] / data['cum_dist'][-1],
                 data['x_values'],
                 f"{data['color']}-",
                 label=data['label'])
    format_subplot("Probability", "Dihedral Angle [deg.]",
                   "Cumulative Probability Distributions")

    plt.tight_layout()
    plt.show()


def format_subplot(xlabel, ylabel, title):
    """Format subplot with consistent styling."""
    plt.xlabel(xlabel, fontsize=16, fontfamily="Helvetica")
    plt.ylabel(ylabel, fontsize=16, fontfamily="Helvetica")
    plt.xticks(fontsize=14, fontfamily="Helvetica")
    plt.yticks(fontsize=14, fontfamily="Helvetica")
    plt.legend(fontsize=14, prop={'family': 'Helvetica'})
    plt.grid(True)
    plt.minorticks_on()
    plt.title(title, fontsize=18, fontfamily="Helvetica")


def compute_persistence_in_repeats(all_data, l_list, Angle_rad, rotation_types,
                                   kTval):
    """
    Calculates the persistence length in repeat units

    Parameters:
    all_data: dict, the all_data you constructed previously (including fitf)
    l_list: list or array, the length of each segment in the primitive (actually irrelevant here)
    angle_rad: array, the bond angle of each segment in the primitive, in radians (e.g., your angle)
    rotation_types: array, the rotation id corresponding to each segment in the primitive (0 means no rotation)

    Returns:
    (lp_in_repeats, lambda_max, Mmat)
    """
    M = len(l_list)
    # L_rep = float(np.sum(l_list))
    A_list = []
    for i in range(M):
        rot_id = int(rotation_types[i])
        theta = float(Angle_rad[i])  # rad
        if rot_id == 0:
            # fixed conformation：phi = 0（m=1,s=0）
            m_i = 1.0
            s_i = 0.0
        else:
            fitf = all_data[rot_id]['fitf']
            # weight w(phi) = exp(-V(phi)/kT)
            Z, _ = quad(lambda phi_deg: np.exp(-fitf(phi_deg) / kTval), 0, 360)

            # calculate <cos(phi)> and <sin(phi)> using quad
            m_i, _ = quad(
                lambda phi_deg: np.cos(np.deg2rad(phi_deg)) * np.exp(-fitf(
                    phi_deg) / kTval), 0, 360)
            m_i /= Z
            s_i, _ = quad(
                lambda phi_deg: np.sin(np.deg2rad(phi_deg)) * np.exp(-fitf(
                    phi_deg) / kTval), 0, 360)
            s_i /= Z

        # create S_i and R_y(theta) to obtain A_i = S_i @ R_y(theta)
        S = np.array([[m_i, -s_i, 0.0], [s_i, m_i, 0.0], [0.0, 0.0, 1.0]])
        c = np.cos(theta)
        s = np.sin(theta)
        R_y = np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])
        A_i = S @ R_y
        A_list.append(A_i)

    # transfer matrix Mmat = A_{M-1} ... A_1 A_0
    Mmat = np.eye(3)
    for A in A_list:
        Mmat = A @ Mmat

    # eigen values
    eigs = eigvals(Mmat)
    lambda_max = float(np.max(np.abs(eigs)))

    # if lambda_max is close to 1, numerical stability may be an issue; if >=1 (numerical error), clip to 1 - eps
    if lambda_max >= 1.0:
        eps = 1e-12
        if lambda_max > 1.0 + 1e-8:
            print("Warning: lambda_max > 1 (numerical error)")
        lambda_max = min(lambda_max, 1.0 - eps)

    # persistence length measured in number of repeat units:
    lp_in_repeats = -1.0 / np.log(lambda_max)
    # persistence length in physical length (same units as l_list)
    # lp_in_length = lp_in_repeats * L_rep
    return lp_in_repeats, lambda_max
