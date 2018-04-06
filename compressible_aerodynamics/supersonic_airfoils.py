import numpy as np
from scipy.optimize import root

def cot(angle):

    f = 1./np.tan(angle)

    return f

def get_theta(tmax, c, alpha, top=True, front=True):
    theta_foil = np.arctan2(tmax, c)
    if front:
        if top:
            theta = np.arctan2(tmax, c) - alpha
            # theta = np.arctan(tmax/c)-alpha
        else:
            theta = np.arctan2(tmax, c) + alpha
            # theta = np.arctan(tmax/c) + alpha
    else:
        theta = -2.*np.arctan2(tmax, c)
        # theta = -2.*np.arctan(tmax/c)

    return theta

def theta_beta_m(beta, M1, gamma):
    tan_theta = 2. * cot(beta) * (M1 **2 * np.sin(beta) ** 2 - 1.) / (M1 ** 2 * (gamma + np.cos(2. * beta)) + 2.)
    return tan_theta

def get_beta(theta, M1, gamma, weak_shock=True):

    def root_func(beta, theta, M1, gamma):
        theta_f = np.arctan(theta_beta_m(beta,M1,gamma))
        f = np.tan(theta) - theta_f

        return f

    # decide which root is of interest and initialize beta0
    if weak_shock:
        beta0 = np.pi/12.
    else:
        beta0 = np.pi/2.

    # solve for beta
    result = root(root_func, beta0, args=(theta, M1, gamma))
    beta = result['x'][0]

    return beta

def get_m_shock(theta, beta, gamma, m1):

    mn1 = m1*np.sin(beta)

    mn2 = np.sqrt((2.+(gamma-1.)*mn1**2)/(2.*gamma*mn1**2-(gamma-1.)))

    m2 = mn2/np.sin(beta-theta)

    return m2

def get_nu(m, gamma):
    nu = np.sqrt((gamma + 1.) / (gamma - 1.)) * np.arctan2(np.sqrt((m ** 2 - 1.) * (gamma - 1.) / (gamma + 1.)),1.) - \
         np.arctan2(np.sqrt(m ** 2 - 1.),1)

    return nu

def get_m_fan(theta, m1, gamma):

    def root_func_fan(m2, m1, theta, gamma):
        nu1 = get_nu(m1, gamma)
        nu2 = get_nu(m2, gamma)
        f = -theta + nu2 - nu1
        return f

    result = root(root_func_fan, m1, args=(m1, abs(theta), gamma))
    m2 = result['x'][0]

    return m2

def get_pressure_oblique_shock(p1, m1, beta, gamma):

    # normal component of upstream mach number
    mn1 = m1*np.sin(beta)

    # pressure behind the oblique shock (see eq. 9.16, pg. 569 in Anderson)
    p2 = p1*(1.+(2.*gamma/(gamma+1.))*(mn1**2-1.))

    return p2

def get_pressure_expansion_fan(p1, m1, m2, gamma):

    # pressure behind the expansion fan shock (see eq. 9.45, pg. 595 in Anderson)
    p2 = p1 * ((1.+((gamma-1.)/2.)*m1**2)/(1.+((gamma-1.)/2.)*m2**2))**(gamma/(gamma-1.))

    return p2

def get_properties(theta, m1, p1, gamma):

    if theta > 0.0:
        beta = get_beta(theta, m1, gamma)
        m2 = get_m_shock(theta, beta, gamma, m1)
        p2 = get_pressure_oblique_shock(p1, m1, beta, gamma)
    else:
        m2 = get_m_fan(theta, m1, gamma)
        p2 = get_pressure_expansion_fan(p1, m1, m2, gamma)

    return m2, p2

def get_dynamic_pressure(gamma, P, M):

    q = (gamma/2.)*P*M**2

    return q

def get_coefficients_shock_theory(d, p1, m1, gamma, tmax, c, alpha):

    theta_front_top = get_theta(tmax, c, alpha, front=True, top=True)
    m2top, p2top = get_properties(theta_front_top, m1, p1, gamma)
    q2top = get_dynamic_pressure(gamma, p2top, m2top)

    theta_back_top = get_theta(tmax, c, alpha, front=False, top=True)
    m3top, p3top = get_properties(theta_back_top, m2top, p2top, gamma)
    q3top = get_dynamic_pressure(gamma, p3top, m3top)

    theta_front_bot = get_theta(tmax, c, alpha, front=True, top=False)
    m2bot, p2bot = get_properties(theta_front_bot, m1, p1, gamma)
    q2bot = get_dynamic_pressure(gamma, p2bot, m2bot)

    theta_back_bot = get_theta(tmax, c, alpha, front=False, top=False)
    m3bot, p3bot = get_properties(theta_back_bot, m2bot, p2bot, gamma)
    q3bot = get_dynamic_pressure(gamma, p3bot, m3bot)

    # print theta_front_top, theta_back_top, theta_front_bot, theta_back_bot
    # print m2top, m3top, m2bot, m3bot
    # print p2top, p3top, p2bot, p3bot

    theta_foil = 0.5*tmax/(c/2.)
    xi1 = theta_foil-alpha
    xi2 = theta_foil-alpha-2.*theta_foil
    xi3 = theta_foil+alpha
    xi4 = theta_foil+alpha-2.*theta_foil

    D1 = p2top*np.sin(xi1)
    D2 = p3top*np.sin(xi2)
    D3 = p2bot*np.sin(xi3)
    D4 = p3bot*np.sin(xi4)
    D = d*(D1+D2+D3+D4)
    # D = d*(-p2top*np.sin(theta_foil-alpha)+p3top*np.sin(-theta_foil-alpha)-p2bot*np.sin(-theta_foil-alpha)+p3bot*np.sin(theta_foil-alpha))

    L1 = -p2top*np.cos(xi1)
    L2 = -p3top*np.cos(xi2)
    L3 = p2bot*np.cos(xi3)
    L4 = p3bot*np.cos(xi4)
    L = d*(L1+L2+L3+L4)

    # L = d*(p2top*np.cos(abs(theta_foil-alpha))-p3top*np.cos(abs(-theta_foil-alpha))+p2bot*np.cos(abs(-theta_foil-alpha))-p3bot*np.cos(abs(theta_foil-alpha)))

    q = get_dynamic_pressure(gamma, p1, m1)

    x1 = (d / 2.) * np.cos(theta_foil - alpha)
    dx1 = -1.
    y1 = (d / 2.) * np.sin(theta_foil - alpha)
    dy1 = 0.0

    x2 = x1 + (c / 2.) * np.cos(alpha)
    dx2 = dx1
    y2 = y1 - (c / 2.) * np.sin(alpha)
    dy2 = 0.0

    x3 = (d / 2.) * np.cos(-theta_foil - alpha)
    dx3 = -1.
    y3 = (d / 2.) * np.sin(-theta_foil - alpha)
    dy3 = 0.0

    x4 = x3 + (c / 2.) * np.cos(alpha)
    dx4 = dx3
    y4 = y3 - (c / 2.) * np.sin(alpha)

    xcp = (x1 * L1 - y1 * D1 + x2 * L2 + y2 * D2 + x3 * L3 + y3 * D3 + x4 * L4 - y4 * D4)/(L1+L2+L3+L4)

    cl = L/(q*c)
    cd = D/(q*c)

    return cl, cd, xcp

def get_coefficients_thin_airfoil_theory(m1, alpha, tmax,c):

    # lift coefficient
    cl = 4*alpha/np.sqrt(m1**2-1.)

    # drag coefficient
    dycdx_bar = 0. # no camber
    dytdx_bar = 4.*tmax**2/c**2 # average airfoil thickness
    cd = (4./np.sqrt(m1**2-1.))*(alpha**2+dycdx_bar**2+dytdx_bar**2)

    cmle = -2.*alpha/np.sqrt(m1**2-1.)

    def root_func_moment(x, cmle, cn, c):
        f = cmle + (x/c)*cn
        return f

    result = root(root_func_moment, 1., (cmle, cl, c))
    xcp = result['x'][0]

    return cl, cd, xcp

def plot_theta_beta_m():
    gamma = 1.4

    import matplotlib.pyplot as plt

    betav = np.linspace(0, np.pi/2, 200)

    for M1 in np.array([1.5, 2., 3., 9E16]):
        thetav = np.arctan(theta_beta_m(betav, M1, gamma))

        plt.plot(thetav * 180. / np.pi, betav * 180. / np.pi, label='M=%0.1f' % M1)
    plt.xlim([0, 50])
    plt.xlabel('deflection angle')
    plt.ylabel('shock angle')
    plt.legend(loc=4, frameon=False)
    plt.show()

def plot_cl_cd(c, tmax):
    import matplotlib.pyplot as plt

    # tmax = 0.05*c
    d = np.sqrt((c / 2.) ** 2 + (tmax / 2.) ** 2)

    M1 = 2.0
    gamma = 1.4
    P1 = 1.01325E5  # N/m**2

    num = 500
    alpha = np.linspace(-15., 15., num) * np.pi / 180.
    cl_ta = np.zeros(num)
    cd_ta = np.zeros(num)
    xcp_ta = np.zeros(num)
    cl_s = np.zeros(num)
    cd_s = np.zeros(num)
    xcp_s = np.zeros(num)

    for i in np.arange(0, num):
        cl_ta[i], cd_ta[i], xcp_ta[i] = get_coefficients_thin_airfoil_theory(M1, alpha[i], tmax, c)
        cl_s[i], cd_s[i], xcp_s[i] = get_coefficients_shock_theory(d, P1, M1, gamma, tmax, c, alpha[i])

    colors = ['r', 'b']
    labels = ['$c_l$ Thin Airfoil Theory', '$c_d$ Thin Airfoil Theory', '$c_l$ Shock Theory', '$c_d$ Shock Theory']
    # alpha = alpha*180./np.pi
    plt.plot(alpha, cl_ta, '--'+colors[0])
    # plt.plot(alpha, cd_ta, colors[0])
    plt.plot(alpha, cl_s, '--'+colors[1])
    # plt.plot(alpha, cd_s, colors[1])

    plt.xlabel('Angle of Attack (deg)')
    plt.ylabel('Coefficient Value')

    plt.legend(labels, loc=4, frameon=False)
    plt.savefig("cl-tmax%.2f.pdf" % tmax, transparent=True)

    # drag polar
    plt.figure()
    plt.plot(cd_ta, cl_ta, colors[0], label="Thin Airfoil Theory")
    plt.plot(cd_s, cl_s, colors[1], label="Shock Theory")
    plt.ylabel("$c_l$")
    plt.xlabel("$c_d$")
    plt.legend(loc=3, frameon=False)
    plt.savefig("drag-polar-tmax%.2f.pdf" % tmax, transparent=True)

    # center of pressure
    plt.figure()
    plt.plot(alpha, xcp_ta, label="Thin Airfoil Theory")
    plt.plot(alpha, xcp_s, label="Shock Theory")
    plt.xlabel('Angle of Attack (deg)')
    plt.ylabel('$X_{cp}$ (percent of chord)')
    plt.ylim([0., 1.])
    plt.legend(frameon=False)
    plt.savefig("center-of-pressure-tmax%.2f.pdf" % tmax, transparent=True)

    plt.show()


if __name__ == "__main__":

    c = 1.0
    tmax = 0.01 * c
    # tmax = 0.05 * c

    plot_cl_cd(c, tmax)

    # plot_theta_beta_m()

    print 'theta_foil: ', np.tan(tmax/c)
    d = np.sqrt((c / 2.) ** 2 + (tmax / 2.) ** 2)

    M1 = 2.0
    gamma = 1.4
    alpha = 0. * np.pi / 180.

    P1 = 1.01325E5  # N/m**2
    # P1 = 1.0E5  # N/m**2

    cl_s, cd_s, xcp_s = get_coefficients_shock_theory(d, P1, M1, gamma, tmax, c, alpha)
    cl_ta, cd_ta, xcp_ta = get_coefficients_thin_airfoil_theory(M1, alpha, tmax, c)
    print cl_s, cd_s
    print cl_ta, cd_ta