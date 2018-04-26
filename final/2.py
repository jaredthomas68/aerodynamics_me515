import numpy as np

def root_function(m1, gamma, m2, beta):

    a = (1./m2)*np.sqrt((2.+(gamma-1.)*(m1*np.sin(beta))**2)/(2.*gamma*(m1*np.sin(beta))**2-(gamma-1.)))
    b = 2.*(1./np.tan(beta))*(m1**2*np.sin(beta)**2-1.)/(m1**2*(gamma+np.cos(2.*beta))+2.)

    res = np.tan(beta-np.arcsin(a)) - b

    return res

if __name__ == "__main__":

    gamma = 1.4
    m2 = 1.5
    beta = 50.*np.pi/180. #shock angle in radians

    from scipy.optimize import root

    res = root(root_function, 2., (gamma, m2, beta))

    print res

