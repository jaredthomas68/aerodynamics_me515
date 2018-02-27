import numpy as np
from scipy.integrate import quad
from matplotlib import pylab as plt
from scipy.interpolate import interp1d



def camber_func(x, epsilon, p):

    # find the value of dy/dx
    if x < 0.0:
        raise ValueError('x must be on the interval [0, 1]')
    if x <= p:
        y = (epsilon/p**2)*(2.*p*x-x**2)
        dydx = (epsilon/p**2)*(2.*p-2.*x)
    elif x <= 1:
        y = (epsilon/(1.-p)**2)*(1.-2.*p+2.*p*x-x**2)
        dydx = (epsilon/(1-p)**2)*(2.*p-2.*x)
    else:
        raise ValueError('x must be on the interval [0, 1]')

    return y, dydx

def thickness_func(x, tau, c):

    # calculate the location for a unit chord
    loc = x/c

    # calculate the thickness of the airfoil at x
    thickness = 10.*tau*(0.2969*np.sqrt(loc)-0.1260*loc-0.3516*loc**2 + 0.2843*loc**3-0.1015*loc**4)

    # calculate the derivative of the thickness distribution with respect to x
    if loc > 0:
        dtdx = 10.*tau*(0.5*0.2969*loc**(-0.5)-0.1260-2.*0.3516*loc + 3.*0.2843*loc**2-4.*0.1015*loc**3)
    else:
        dtdx = 0.0

    return thickness, dtdx

def q_func(x, tau, c, v0):

    _, dtdx = thickness_func(x, tau, c)

    return v0 * dtdx

def q_func_integrand(s, args):
    x = args[0]
    tau = args[1]
    c = args[2]
    v0 = args[3]
    _, dtdx = thickness_func(s, tau, c)

    return v0*dtdx/(x-s)

def lift_func(A0, A1, V0, c, density):

    # calculate the lift per unit span
    Lp = np.pi*density*(V0**2)*c*(A0+0.5*A1)

    # calculate the lift coefficient
    cl = 2.*np.pi*(A0 + 0.5*A1)

    return Lp, cl

def moment_func(A0, A1, A2, x, density, V0, c):

    # calculate the moment at x
    M = np.pi*density*(V0**2)*c*(A0*(x-c/4.)+A1*(x/2.-c/4.)+A2*c/8.)

    # calculate the coefficient
    cm = 2.*np.pi*(A0*(x-c/4.)+A1*(x/2.-c/4.)+A2*c/8.)

    return M, cm

def cp_func(xx, alpha, gamma, A, epsilon, p, tau=0.12, c=1., V0=10., top=True):

    u0 = V0*np.cos(alpha)
    v0 = V0*np.sin(alpha)

    if np.size(xx) > 1:
        cp = np.zeros_like(xx)
        for i, x in zip(np.arange(0, xx.size), xx):
            # print "HERE", x
            uc, vc, ut, vt = velocity_components_func(x, A, tau, c, V0, top)

            # Riegel correction
            _, dydx = camber_func(x, epsilon, p)
            _, dTdx = thickness_func(x, tau, c)
            Vmod = ((u0+ut+uc)**2+(v0+vt+vc)**2)/np.sqrt(1.+(dydx+0.5*dTdx)**2)
            # calculate the pressure at a given point
            # cp[i] = 1.-((u0+ut+uc)**2+(v0+vt+vc)**2)/(u0**2+v0**2)
            cp[i] = 1.-Vmod/V0**2

    else:
        uc, vc, ut, vt = velocity_components_func(xx, A, tau, c, V0, top)

        # calculate the pressure at a given point
        # cp = 1. - ((u0 + ut + uc) ** 2 + (v0 + vt + vc) ** 2) / V0**2

        _, dydx = camber_func(xx, epsilon, p)
        _, dTdx = thickness_func(xx, tau, c)
        Vmod = ((u0 + ut + uc) ** 2 + (v0 + vt + vc) ** 2) / np.sqrt(1. + (dydx + 0.5 * dTdx) ** 2)

        cp = 1. - Vmod / V0**2

    return cp

def velocity_components_func(x, A, tau=0.12, c=1., v0=10., top=True):
    tol = 1E-6

    # solve for velocity related to camber
    if top:
        uc = gamma_pre_computed_function(x, c)/2.
    else:
        uc = -gamma_pre_computed_function(x, c)/2.

    # print x
    vc = (-1./(2.*np.pi))*quad(gamma_integrand, tol, c, args=[A, x, c, v0], limit=1000, points=[x])[0]

    # solve for velocity related to thickness
    ut = (1./(2.*np.pi))*quad(q_func_integrand, tol, c, args=[x, tau, c, v0], limit=1000, points=[x])[0]

    if top:
        vt = q_func(x, tau, c, v0)/2.
    else:
        vt = -q_func(x, tau, c, v0)/2.


    return uc, vc, ut, vt

# def boundary_velocity_func(V0, alpha, dydx, surface):
#
#     if surface == 0:
#         vt = 0.5*V0*dTdx
#         uc =
#     elif surface == 1:
#         vt = -0.5 * V0 * dTdx
#
#     vc = -V0*(alpha-dydx)

def gamma_integrand(s, args):
    A = args[0]
    x = args[1]
    c = args[2]
    v0 = args[3]
    gamma_int = gamma_pre_computed_function(s, c)/(x-s)
    if abs(x-s) < 1E-12:
        print "in gamma integrand, %0.10f, %0.19f" %(x, s)
    return gamma_int


def gamma_func(x, A, v0, c=1.):

    # find angle
    phi = np.arccos(1.-2*x/c)
    # cosinex = 1. - 2.*x/c
    # sinex = 2.*np.sqrt((x/c)*(1-x/c))

    # calculate gamma
    if np.size(phi) > 1:
        gamma = np.zeros_like(phi)
        for p, i in zip(phi, np.arange(0, phi.size)):
            print np.arange(1, A.size)
            gamma[i] = 2.*v0*(A[0]*(1.+np.cos(p))+np.sum(A[1:]*np.sin(np.arange(1,A.size)*p)*np.sin(p)))/np.sin(p)
    else:
        print np.arange(1, A.size)
        gamma = 2.*v0 * (A[0] * (1. + np.cos(phi)) + np.sum(A[1:] * np.sin(np.arange(1, A.size) * phi) * np.sin(phi))) / np.sin(phi)
    # gamma = 2.*v0*(A[0]*(1+cosinex)+np.sum(A[1:]*np.sin(np.arange(1,A.size)*phi)*np.sin(phi)))/np.sin(phi)

    return gamma

def gamma_pre_computed_function(x, c=1.):
    global gamma_interp_func
    # find angle
    phi = np.arccos(1. - 2. * x / c)

    # get gamma value
    gamma = gamma_interp_func(phi)

    return gamma

def pre_compute_gamma(phi, A, v0):

    GG = np.zeros_like(phi)

    for p, i in zip(phi, np.arange(0, phi.size)):

        if p < 1E-10:
            GG[i] = 10000.
        else:
            GG[i] = 2. * v0 * (
                    A[0] * (1. + np.cos(p)) + np.sum(A[1:] * np.sin(np.arange(1, A.size) * p) * np.sin(p))) / np.sin(p)
        # if GG[i] < 0:
        #     GG[i] = 0.0
    gamma_data_func = interp1d(phi, GG, kind='cubic')

    return gamma_data_func

def integrand_func(theta, args):

    # parse input args
    epsilon = args[0]
    p = args[1]
    c = args[2]
    N = args[3]

    # calculate x from theta
    x = (c/2.)*(1.-np.cos(theta))

    _, b = camber_func(x, epsilon, p)

    # adjust the integrand according to which fourier coefficient is being calculated
    if N == 0:
        integrand = b
    else:
        integrand = b*np.cos(N*theta)

    return integrand

def fourier_coefficient_func(alpha, btheta, btheta_args, Ncoefficients):

    # initalize fourier coefficient array
    fourier_coefficients = np.zeros(Ncoefficients)

    # calc first fourier coefficient
    btheta_args[-1] = 0
    fourier_coefficients[0] = alpha - (1./np.pi)*(quad(btheta, 0.0, np.pi, btheta_args)[0])

    # calc the second to N fourier coefficients
    for i in np.arange(1, Ncoefficients):
        btheta_args[-1] = i
        fourier_coefficients[i] = (2./np.pi)*(quad(btheta, 0.0, np.pi, btheta_args)[0])

    return fourier_coefficients

def find_ideal_angle_of_attack():

    # get fourier coefficients
    alpha = 0.0866455743735  # angle of attack where Cl=1.0
    # alpha = 0.0
    epsilon = 0.04
    p = 0.4
    c = 1.0
    tau = 0.12
    N = 2

    V0 = 10.
    density = 1.2

    btheta = integrand_func
    btheta_args = [epsilon, p, c, 0]

    alphas = np.linspace(-0.1, 2., 100)
    fc = np.zeros([alphas.size,N])

    for i, a in zip(np.arange(0, alphas.size), alphas):
        fc[i, :] = fourier_coefficient_func(a, btheta, btheta_args, N)

    from scipy.interpolate import UnivariateSpline

    spl = UnivariateSpline(fc[:,0], alphas)



    print spl(0)


def plot_cp(save_data):

    # get fourier coefficients
    # alpha = 0.0866455743735 # angle of attack where Cl=1.0
    # alpha = 0.00898577292547 # ideal angle of attack, where A0 = 0.0
    alpha = 0.0
    epsilon = 0.04
    p = 0.4
    c = 1.0
    tau = 0.12
    N = 100

    V0 = 10.
    density = 1.2

    btheta = integrand_func
    btheta_args = [epsilon, p, c, 0]

    fc = fourier_coefficient_func(alpha, btheta, btheta_args, N)
    print "after FC"
    print fc
    quit()
    # pre-compute gamma
    samples = 10000
    phi0 = np.linspace(0.0,np.pi,samples)
    global gamma_interp_func
    gamma_interp_func = pre_compute_gamma(phi0, fc, V0)

    print "gamma computed"

    # x1 = np.linspace(0, 1, 100)
    # gamma = gamma_func(x1, fc, V0, c)

    # x2 = np.linspace(0., 1, 100)
    # phi2 = np.arccos(1.-2*x2/c)
    # phi = np.arccos(1. - 2 * x1 / c)
    # plt.scatter(phi, gamma)
    # plt.plot(phi0, gamma_interp_func(phi0))
    # plt.show()
    # quit()
    tol = 1E-6
    # x2 = np.linspace(tol, c, 200)
    # cp_top = cp_func(x2, alpha, gamma_integrand, fc, tau, c, V0, top=True)
    # cp_bottom = cp_func(x2, alpha, gamma_integrand, fc,tau,  c, V0, top=False)
    #
    # plt.plot(x2, cp_top, '-')
    # plt.plot(x2, cp_bottom, '-')

    x2 = np.linspace(tol, c, 30)
    cp_top = cp_func(x2, alpha, gamma_integrand, fc, epsilon, p, tau, c, V0, top=True)
    cp_bottom = cp_func(x2, alpha, gamma_integrand, fc, epsilon, p, tau, c, V0, top=False)

    if save_data:
        np.savetxt("thin_airfoil_cp_alphac%.3f.txt" % alpha, np.c_[x2, cp_bottom, cp_top], header="x, cp_bottom, cp_top")

    plt.plot(x2, cp_top, '--')
    plt.plot(x2, cp_bottom, '--')
    plt.ylim([2, -2])

    plt.show()


def plot_cl_and_cmac(save_data):
    alpha = 0.0866455743735
    alpha = 0.
    epsilon = 0.04
    p = 0.4
    c = 1.0
    N = 100

    V0 = 10.
    density = 1.2

    btheta = integrand_func
    btheta_args = [epsilon, p, c, 0]

    fc = fourier_coefficient_func(alpha, btheta, btheta_args, N)
    fc
    print fc
    quit()

    samples = 20
    cl = np.zeros(samples)
    cm = np.zeros(samples)
    alpha = np.linspace(-np.pi / 8., np.pi / 8., samples)
    for i, a in zip(np.arange(0, samples), alpha):
        fc = fourier_coefficient_func(a, btheta, btheta_args, N)
        _, cl[i] = lift_func(fc[0], fc[1], V0, c, density)
        _, cm[i] = moment_func(fc[0], fc[1], fc[2], c / 4., density, V0, c)

    from scipy.interpolate import UnivariateSpline
    # spl =
    # spl = UnivariateSpline(cl, alpha)
    # alpha_cl1 = spl(1.)
    #
    # spl1 = UnivariateSpline(alpha, cl)
    # print spl1(alpha_cl1), alpha_cl1

    if save_data:
        np.savetxt("thin_airfoil_cl_cmac.txt" % alpha, np.c_[alpha, cl, cm], header="alpha, cl, cm")
    print cl, cm
    plt.plot(alpha, cl, label='$C_L$')
    plt.plot(alpha, cm, label='$C_{mac}$')
    plt.legend(loc=2)
    plt.xlabel('Angle of Attack (radians)')
    plt.ylabel('Coefficient Value')
    plt.show()

if __name__ == "__main__":
    save_data = False
    # find_ideal_angle_of_attack()
    # plot_cl_and_cmac(save_data)
    plot_cp(save_data)