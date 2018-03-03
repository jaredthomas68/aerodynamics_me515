import numpy as np
from matplotlib import pylab as plt
from scipy.integrate import quad, odeint, trapz
import cmath as cm

def plot_results(xdata, ydata, labels, xlabel='', ylabel='', title='', linestyles=None, logy=False, logx=False, legend_loc=0, ylim=None, save_plot=False, filename=''):

    for x, y, linestyle in zip(xdata, ydata, linestyles):
        if logy and not logx:

            plt.semilogy(np.transpose(xdata), np.transpose(ydata))

        elif logx and not logy:

            plt.semilogx(np.transpose(xdata), np.transpose(ydata))

        elif logy and logx:

            plt.loglog(np.transpose(xdata), np.transpose(ydata))

        else:

            plt.plot(x, y, linestyle=linestyle)

    plt.legend(labels=labels, loc=legend_loc, frameon=False)

    if ylim is None:
        ylim = [np.min(ydata), np.max(ydata)]
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.ylim(ylim)

    if save_plot:
        plt.savefig(filename, transparent=True)

    plt.show()

    return


def reynolds_number(V, L, nu=0, mu=0, rho=0):

    if nu == 0:
        R = rho*V*L/mu
    else:
        R = V*L/nu

    return R

def shape_factor_tlambda(tlambda):

    if tlambda >= 0.0 and tlambda <= 0.1:
        L = 0.22+1.57*tlambda-1.8*tlambda**2
        H = 2.61 - 3.75*tlambda + 5.24*tlambda*2
    elif tlambda >= -0.1 and tlambda < 0.0:
        L = 0.22 + 1.402*tlambda + 0.018*tlambda/(0.107+tlambda)
        H = 2.088 + 0.0731/(0.14+tlambda)
    else:
        raise ValueError("tlambda must be in the interval [-0.1 and 0.1]")

    return H, L

def shape_factor_dstar_theta(delta_star, theta):

    H = delta_star/theta

    return H

def shape_factor_of_h1(H1):

    if H1 >= 5.3:
        H = 1.1 + 0.86*(H1-3.3)**-0.777

    else:
        H = 0.6778 + 1.1538*cm.exp(-0.326*cm.log(H1-3.3))

    return np.real(H)

def fric_coefficient_of_h(H, Ret):

    Cf = 0.246*10**(-0.678*H)*Ret**-0.268

    return Cf

def shape_factor_1(H):

    if H <= 1.6:
        H1 = 0.8234*(H-1.1)**-1.287+3.3
    else:
        H1 = 1.5501*(H-0.6778)**-3.064 + 3.3

    return H1

def edge_velocity_integrand(x, args):

    # check if flat plate
    if args[0]:
        ve = args[1] # if flat plate, ve is freestreem velocity

    return ve

def thwaite_boundary(x, L, v0, ve, ve0, nu, flat_plate=True, dve_dx0=0.):

    x_star = x/L
    ve_star = ve/v0
    ve0_star = ve0/v0
    v0_star = v0/v0

    # Reynold's number based on length
    Rel = reynolds_number(v0, L, nu)

    if flat_plate:
        theta0 = 0.
    else:
        if dve_dx0 == 0.0:
            raise ValueError('dve_dx0 must be specified if not solving a flat plate')

        theta0 = np.sqrt(0.075*nu/(dve_dx0))

    #initialize disp thickness
    theta = np.zeros_like(x_star)
    delta_star = np.zeros_like(x_star)
    for i in np.arange(0, x_star.size):
        result = quad(edge_velocity_integrand, 0., x_star[i], np.array([flat_plate, v0_star]))
        ve_integral = result[0]
        # displacement thickness
        theta[i] = L*np.sqrt((0.45)/(Rel*ve_star**6)*ve_integral+((theta0/L)**2)*(ve0_star/ve_star)**6)

        # Reynold's number based on displacement thickness
        Ret = reynolds_number(v0, L, nu=nu)

        tlambda = (theta[i]**2/nu)*dve_dx0

        shape_factor, _ = shape_factor_tlambda(tlambda)

        delta_star[i] = shape_factor*theta[i]

    #TODO fix cf
    Rex = reynolds_number(v0,x,nu)
    cf = 0.664/np.sqrt(Rex)

    CF = trapz(cf, x, x[1]-x[0])

    return delta_star, theta, np.ones_like(x)*shape_factor, cf, CF

def head_ode_func(y, x, v0, nu, ve, dve_dx):
    theta =y[0]
    H1 = y[1]
    H = shape_factor_of_h1(H1)

    Ret = reynolds_number(v0, theta, nu)
    cf = fric_coefficient_of_h(H,Ret)

    dtheta_dx = cf/2. - (theta/ve)*(H+2.)*dve_dx

    dH1_dx = (1./(ve*theta))*(-ve*H1*dtheta_dx-theta*H1*dve_dx) + ve*0.0396*cm.exp(-0.6169*cm.log(H1-3.))

    return np.array([dtheta_dx, np.real(dH1_dx)])

def head_boundary(x, v0, ve, dve_dx, nu, theta0, H10):

    result = odeint(head_ode_func, np.array([theta0, H10]), x, args=(v0, nu, ve, dve_dx))
    theta = result[:, 0]
    H1 = result[:, 1]

    H = np.zeros_like(H1)
    for i in np.arange(0, H1.size):
        H[i] = shape_factor_of_h1(H1[i])

    Ret = reynolds_number(v0, theta, nu)
    cf = fric_coefficient_of_h(H, Ret)

    delta_star = theta*H

    CF = trapz(cf, x, x[1]-x[0])

    return delta_star, theta, H, cf, CF

def blasius_boundary(x, Rex, Rel):

    # boundary layer thickness
    delta = 5.*x/np.sqrt(Rex)

    # displacement thickness
    delta_star = 1.72*x/np.sqrt(Rex)

    # momentum thickness
    theta = 0.664*x/np.sqrt(Rex)

    # local friction coefficient
    cf = 0.664/np.sqrt(Rex)

    # skin friction drag coefficient
    CF = 1.328/np.sqrt(Rel)

    shape_factor = shape_factor_dstar_theta(delta_star, theta)

    return delta_star, theta, shape_factor, cf, CF

def schlichting_boundary(x, Rex, Rel):

    # boundary layer thickness
    delta = 0.37*x/Rex**0.2

    # displacement thickness
    delta_star = 0.046*x/Rex**0.2

    # momentum thickness
    theta = 0.036*x/Rex**0.2

    # local friction coefficient
    cf = 0.0592/Rex**0.2

    # skin friction drag coefficient
    CF = 0.074/Rel**0.2

    shape_factor = shape_factor_dstar_theta(delta_star, theta)

    return delta_star, theta, shape_factor, cf, CF