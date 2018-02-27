import numpy as np
import matplotlib.pylab as plt

def cart_to_polar(x, y):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan(y/x)

    return r, theta

def uniform_flow(r, theta, v0, kappa=1.):

    psi = v0*r*np.sin(theta)

    return psi

def doublet(r, theta, kappa=1.):

    psi = -kappa*np.sin(theta)/(2.*np.pi*r)

    return psi


def cylinder_flow(x, y, v0, R):

    kappa = 2.*np.pi*v0*R**2

    r, theta = cart_to_polar(x, y)

    psi_uniform = uniform_flow(r, theta, v0, kappa)

    psi_doublet = doublet(r, theta, kappa)

    psi_total = psi_uniform + psi_doublet

    Vr = (1./r)*(v0*r*np.cos(theta))

    Vt = -(1.+R**2/r**2)*v0*np.sin(theta)

    V = Vt*np.cos(theta) + Vr*np.cos(theta)


    return psi_total, V


def stagnation_point():


    return

if __name__ == "__main__":

    x = np.arange(-2., 2., 0.01)
    y = np.arange(-2., 2., 0.01)

    xx, yy = np.meshgrid(x, y)

    psi, V = cylinder_flow(xx, yy, v0=1., R=1.)

    psi[xx**2+yy**2<=1.] = 0

    plt.contour(xx, yy, V)

    plt.show()

