import numpy as np

def cart_to_polar(x, y):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(x,y)

    return r, theta

def uniform_flow(x, y, v0, kappa=1.):

    r, theta = cart_to_polar(x, y)

    psi = v0*r*np.sin(theta)

    return psi

def doublet_flow(x, y, kappa=1.):

    r, theta = cart_to_polar(x, y)

    psi = -kappa*np.sin(theta)/(2.*np.pi*r)

    return psi

def vortex_flow(xx, yy, x, y, gamma=1.):


    # convert to polar
    r, theta = cart_to_polar((xx-x), (yy-y))

    # stream function

    psi = gamma*np.log(r)/(2.*np.pi)

    Vt = -gamma/(2.*np.pi*r)

    Vx = -Vt * np.cos(theta)

    Vy = Vt * np.sin(theta)

    return psi, Vx, Vy
