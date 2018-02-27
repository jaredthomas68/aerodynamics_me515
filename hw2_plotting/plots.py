import numpy as np
from matplotlib import pylab as plt

def plot_cl_cmac():
    taf = np.loadtxt("thin_airfoil_cl_cmac.txt")
    taf_alpha = taf[:, 0]
    taf_cl = taf[:, 1]
    taf_cmac = taf[:, 2]
    panel = np.loadtxt("panel_cl_cmac.txt")
    p_alpha = panel[:, 0]
    p_cl= panel[:, 1]
    p_cmac= panel[:, 2]

    plt.plot(taf_alpha, taf_cl, 'r', label="$C_L$ (thin air foil theory)")
    plt.plot(taf_alpha, taf_cmac, 'b', label="$C_{mac}$ (thin air foil theory)")
    plt.plot(p_alpha, p_cl, '--r', label="$C_L$ (panel method)")
    plt.plot(p_alpha, p_cmac, '--b', label="$C_{mac}$ (panel method)")

    plt.legend(loc=2)
    plt.xlabel("Alpha (radians)")
    plt.ylabel("Coefficient value")

    plt.savefig('cl_cmac.pdf')
    plt.show()

def plot_convergence():

    data = np.loadtxt("panel_convergence_cl_cmac_cd_alpha0.000.txt")

    panels = data[:, 0]
    cl = data[:, 1]
    cmac = data[:, 2]
    cd = data[:, 3]

    plt.plot(panels, cl, label="$C_L$")
    plt.plot(panels, cd, label="$C_D$")
    plt.plot(panels, cmac, label="$C_{mac}$")

    plt.legend()
    plt.xlabel("Panels")
    plt.ylabel("Coefficient value")

    plt.savefig('convergence.pdf')
    plt.show()

def plot_Cp():

    panel_ideal = np.loadtxt("panel_cp_alphac0.009.txt")
    panel_cl1 = np.loadtxt("panel_cp_alphac0.087.txt")
    taf_ideal = np.loadtxt("thin_airfoil_cp_alphac0.009.txt")
    taf_cl1 = np.loadtxt("thin_airfoil_cp_alphac0.087.txt")
    cfd_ideal = np.loadtxt("cp_cfd_a00.txt")
    cfd_cl1 = np.loadtxt("cp_cfd_cl1.txt")

    x_pi = panel_ideal[:, 0]
    cp_pi = panel_ideal[:, 1]

    x_pcl1 = panel_cl1[:, 0]
    cp_pcl1 = panel_cl1[:, 1]

    x_ti = taf_ideal[:, 0]
    cp_ti_b = taf_ideal[:, 1]
    cp_ti_t = taf_ideal[:, 2]

    x_tcl1 = taf_cl1[:, 0]
    cp_tcl1_b = taf_cl1[:, 1]
    cp_tcl1_t = taf_cl1[:, 2]

    x_cfd_ideal = cfd_ideal[:, 0]
    cp_cfd_ideal = cfd_ideal[:, 1]

    x_cfd_cl1 = cfd_cl1[:, 0]
    cp_cfd_cl1 = cfd_cl1[:, 1]


    plt.plot(x_pi, cp_pi, 'r', label="$C_P$ (panel)")
    plt.plot(x_ti, -cp_ti_b, 'b--', label="$C_P$ (thin airfoil)")
    plt.plot(x_ti, cp_ti_t, 'b--')
    plt.scatter(x_cfd_ideal, cp_cfd_ideal, c='g', label="$C_P$ (CFD)", edgecolors='g')

    plt.legend()
    plt.xlabel("X postion")
    plt.ylabel("$C_P$")
    # plt.show()
    plt.ylim([2,-2])


    plt.savefig("cp_compare_ideal.pdf")

    plt.figure()

    plt.plot(x_pcl1, cp_pcl1, 'r', label="$C_P$ (panel)")
    plt.plot(x_tcl1, -cp_tcl1_b, 'b--', label="$C_P$ (thin airfoil)")
    plt.plot(x_tcl1, cp_tcl1_t, 'b--')
    plt.scatter(x_cfd_cl1, cp_cfd_cl1, c='g', label="$C_P$ (CFD)", edgecolors='g')
    plt.legend()
    plt.xlabel("X postion")
    plt.ylabel("$C_P$")

    plt.ylim([2,-2])

    plt.savefig("cp_compare_cl1.pdf")
    plt.show()

if __name__ == "__main__":

    plot_cl_cmac()
    # plot_convergence()
    # plot_Cp()