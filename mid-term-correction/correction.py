import numpy as np


if __name__ == "__main__":

    alpha = 0.1
    v = 2.
    c = 2.5

    from scipy.integrate import quad

    def integrand(x, c, v, alpha):
        g = 2.*v*alpha*np.sqrt(1.-x/c)/np.sqrt(x/c)

        return g

    res = quad(integrand, 0, c, (c, v, alpha))

    print res[0]/(2.*c)

