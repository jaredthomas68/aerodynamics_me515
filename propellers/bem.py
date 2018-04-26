import numpy as np
from scipy.interpolate import UnivariateSpline, interp1d

def get_clcd_splines():

    clt = np.loadtxt("cl_data_tangler.csv", delimiter=',', encoding='utf-8-sig')
    clx = np.loadtxt("cl_data_xfoil.csv", delimiter=',', encoding='utf-8-sig')
    clo = np.loadtxt("cl_data_ostowari.csv", delimiter=',', encoding='utf-8-sig')

    allcl_alpha = np.hstack([clt[:, 0]*np.pi/180., clx[:, 0]*np.pi/180., clo[:, 0]*np.pi/180.])
    allcl_cl = np.hstack([clt[:, 1], clx[:, 1], clo[:, 1]])
    # allcl_spline = UnivariateSpline(allcl_alpha, allcl_cl, k=5)
    allcl_spline = interp1d(allcl_alpha, allcl_cl, bounds_error=False, fill_value='extrapolate')

    cdt = np.loadtxt("cd_data_tangler.csv", delimiter=',', encoding='utf-8-sig')
    cdx = np.loadtxt("cd_data_xfoil.csv", delimiter=',', encoding='utf-8-sig')
    cdo = np.loadtxt("cd_data_ostowari.csv", delimiter=',', encoding='utf-8-sig')

    allcd_alpha = np.hstack([cdt[:, 0]*np.pi/180., cdx[:, 0]*np.pi/180., cdo[:, 0]*np.pi/180.])
    allcd_cd = np.hstack([cdt[:, 1], cdx[:, 1], cdo[:, 1]])
    # allcd_spline = UnivariateSpline(allcd_alpha, allcd_cd)
    allcd_spline = interp1d(allcd_alpha, allcd_cd, bounds_error=False, fill_value='extrapolate')

    return allcl_spline, allcd_spline

def get_induction():

    return

def get_thrust():

    return

def get_torque():

    return

def phi_resid(phi, theta, M, Re, c, r, Rb, B, omega, Vinf, cl_spline, cd_spline, return_all=False):
    """

    :param phi: inflow angle (rad)
    :param theta: twist (rad)
    :param M: mach number
    :param Re: reynolds number
    :param c: chord length (m)
    :param r: local radius (m)
    :param Rb: propeller radius (m)
    :param B: number of blades
    :param omega: rotation rate (rad/s)
    :param Vinf:
    :return: residual
    """

    # angle of attack
    alpha = theta - phi
    alphad =alpha*180./np.pi
    cl = cl_spline(alpha)
    cd = cd_spline(alpha)

    cl = cl/np.sqrt(1.-M**2)

    cn = cl*np.cos(phi)-cd*np.sin(phi)
    ct = cl*np.sin(phi)+cd*np.cos(phi)

    # tip losses
    insidef = np.sin(phi)
    f = (B/2.)*(Rb-r)/(r*np.abs(np.sin(phi)))
    inside = np.exp(-f)
    F = (2./np.pi)*np.arccos(inside)

    # local solidity
    sigma_p = B * c / (2. * np.pi * r)

    # induction factor
    a = 1./(4.*F*np.sin(phi)**2/(sigma_p*cn)-1.)

    # tangential induction factor
    a_p = 1./(4.*F*np.sin(phi)*np.cos(phi)/(sigma_p*ct)+1.)

    # solution terms
    lambdar = omega*r/Vinf
    Rs = np.sin(phi)/(1.+a) - np.cos(phi)/(lambdar*(1.-a_p))

    if return_all:
        return cn, ct, a, a_p
    else:
        return Rs

def integrand(r, phi0, theta, M, Re, rho, c, Rb, B, omega, Vinf, cl_spline, cd_spline, False, output='thrust'):

    from scipy.optimize import root

    res = root(phi_resid, phi0, (theta, M, Re, c, r, Rb, B, omega, Vinf, cl_spline, cd_spline, False))

    phi_star = res['x']

    cn, ct, a, a_p = phi_resid(phi_star, theta, M, Re, c, r, Rb, B, omega, Vinf, cl_spline, cd_spline, True)

    w_squared = (Vinf * (1. + a)) ** 2 + (omega * r * (1 - a_p)) ** 2

    if output is 'thrust':
        return B * cn * 0.5 * rho * w_squared * c
    elif output is 'torque':
        return B * ct * 0.5 * rho * w_squared * c * r
    else:
        raise ValueError*("incorrect output type")

def bem(J, r, c, beta, Rb=1.):
    from scipy.integrate import quad
    cl_spline, cd_spline = get_clcd_splines()

    B = 2.
    Db = Rb*2.
    speed_of_sound = 343. #m/s
    phi0 = 2.*np.pi/180.
    rho = 1.23
    mu = 1.79E-5

    omega_m = 5000.
    omega = omega_m*(1./60.)*2.*np.pi

    n = omega/(2.*np.pi)

    rh = np.min(r)
    rt = Rb

    CT = np.zeros_like(J)
    CP = np.zeros_like(J)
    CQ = np.zeros_like(J)
    eta = np.zeros_like(J)

    for i in np.arange(0, J.size):

        theta = beta[i]*np.pi/180.

        Vinf = J[i] * n * Db

        M = Vinf / speed_of_sound

        Re = rho*Vinf*c[i]/mu

        # get thrust
        res = quad(integrand, rh, rt, (phi0, theta, M, Re, rho, c[i], Rb, B, omega, Vinf, cl_spline, cd_spline, False, 'thrust'))
        T = res[0]

        # get torque
        res = quad(integrand, rh, rt,
                   (phi0, theta, M, Re, rho, c[i], Rb, B, omega, Vinf, cl_spline, cd_spline, False, 'torque'))
        Q = res[0]

        # get power
        P = Q*omega

        # get coefficients
        CT[i] = T/(rho*n**2*Db**4)
        CQ[i] = Q/(rho*n**2*Db**5)
        CP[i] = P/(rho*n**3*Db**5)

        # get efficiency
        eta[i] = J[i]*CT[i]/CP[i]

    return CT, CQ, CP, eta

def plot_clcd():
    clt = np.loadtxt("cl_data_tangler.csv", delimiter=',', encoding='utf-8-sig')
    clx = np.loadtxt("cl_data_xfoil.csv", delimiter=',', encoding='utf-8-sig')
    clo = np.loadtxt("cl_data_ostowari.csv", delimiter=',', encoding='utf-8-sig')

    allcl_alpha = np.hstack([clt[:, 0], clx[:, 0], clo[:, 0]])
    allcl_cl = np.hstack([clt[:, 1], clx[:, 1], clo[:, 1]])
    # allcl_spline = UnivariateSpline(allcl_alpha, allcl_cl, k=5)
    allcl_spline = interp1d(allcl_alpha, allcl_cl)

    cdt = np.loadtxt("cd_data_tangler.csv", delimiter=',', encoding='utf-8-sig')
    cdx = np.loadtxt("cd_data_xfoil.csv", delimiter=',', encoding='utf-8-sig')
    cdo = np.loadtxt("cd_data_ostowari.csv", delimiter=',', encoding='utf-8-sig')

    allcd_alpha = np.hstack([cdt[:, 0], cdx[:, 0], cdo[:, 0]])
    allcd_cd = np.hstack([cdt[:, 1], cdx[:, 1], cdo[:, 1]])
    # allcd_spline = UnivariateSpline(allcd_alpha, allcd_cd)
    allcd_spline = interp1d(allcd_alpha, allcd_cd)

    import matplotlib.pyplot as plt

    labels = ['Tangler 2005', 'XFOIL', 'Ostowari 1985', 'Interpolation']

    plt.scatter(clt[:, 0], clt[:, 1], color='r', label=labels[0])
    plt.scatter(clx[:, 0], clx[:, 1], color='b', label=labels[1])
    plt.scatter(clo[:, 0], clo[:, 1], color='y', label=labels[2])

    clx = np.linspace(-28, 29.5, 100)
    cly = allcl_spline(clx)
    plt.plot(clx, cly, label=labels[3])
    plt.xlabel('Angle of Attack (deg)')
    plt.ylabel('Lift Coefficient')

    plt.legend(frameon=False, loc=4)
    plt.savefig('cl.pdf')

    plt.figure()

    plt.semilogy(cdt[:, 0], cdt[:, 1], 'or', label=labels[0])
    plt.semilogy(cdx[:, 0], cdx[:, 1], 'sb', label=labels[1])
    plt.semilogy(cdo[:, 0], cdo[:, 1], '*y', label=labels[2])

    cdx = np.linspace(-29, 30, 100)
    cdy = allcd_spline(cdx)
    plt.plot(cdx, cdy, label=labels[3])
    plt.xlabel('Angle of Attack (deg)')
    plt.ylabel('Drag Coefficient')
    plt.legend(frameon=False, loc=4)

    plt.savefig('cd.pdf')

    plt.show()

    return 0

def plot_results():

    data = np.loadtxt('apce_10x5_kt0821_5000.txt')
    data_g = np.loadtxt('apce_10x5_geom.txt')

    r = data_g[:, 0]
    c = data_g[:, 1]
    beta = data_g[:, 2]

    J = data[:, 0]
    CTd = data[:, 1]
    CPd = data[:, 2]
    etad = data[:, 3]

    Jr = J
    CT, CQ, CP, eta = bem(Jr, r, c, beta)

    import matplotlib.pyplot as plt

    labels = ['BEM', 'UIUC']
    plt.plot(Jr, CT, label=labels[0])
    plt.scatter(J, CTd, label=labels[1])
    plt.ylabel('CT')
    plt.xlabel('J')
    plt.legend(frameon=False)

    plt.savefig('CT.pdf')

    plt.figure()
    plt.plot(Jr, CP, label=labels[0])
    plt.scatter(J, CPd, label=labels[1])
    plt.ylabel('CP')
    plt.xlabel('J')
    plt.legend(frameon=False)

    plt.savefig('CP.pdf')

    plt.figure()
    plt.plot(Jr, eta, label=labels[0])
    plt.scatter(J, etad, label=labels[1])
    plt.ylabel('$\eta$')
    plt.xlabel('J')
    plt.ylim([0, 1])
    plt.legend(frameon=False)

    plt.savefig('eta.pdf')

    plt.show()

    return


if __name__ == "__main__":

    # plot_clcd()

    plot_results()
