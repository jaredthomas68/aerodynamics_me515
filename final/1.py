import numpy as np

def xle(y, c0):

    x = -0.25*c0*np.sqrt(1.-(2.*y/b)**2)

    return x

def xte(y, c0):

    x = 0.75 * c0 * np.sqrt(1. - (2. * y / b) ** 2)

    return x

def chord(y, c0):

    te = xte(y, c0)
    le = xle(y, c0)
    # print le, te
    c = te - le

    # c = xte(y, c0) - xle(y, c0)

    return c

if __name__ == "__main__":


    b = 10.     # wing span (m)
    S = 20.     # wind area (m**2)
    P = 1E5     # free stream pressue (Pa)
    M = 0.4     # mach number
    alpha = 2.  # degrees
    AR = b**2/S # aspect ratio (m)
    c0 = 4.*S/(np.pi*b)
    gamma = 1.4

    from scipy.integrate import quad

    hSc = quad(chord, 0., b/2., c0)

    Sc = 2.*hSc[0]

    print S, Sc


    alpha *= np.pi/180. # convert angle of attack to radians

    m = 2.*np.pi

    cl0 = m*alpha*(1.-S*m/(b**2*np.pi+S*m))

    print "cl0: ", cl0

    cl = cl0/np.sqrt(1.-M**2)

    print "cl: ", cl

    L = 0.5*cl*(M**2)*gamma*P*S

    print "L: ", L

    print np.sin(alpha), alpha
    print np.cos(alpha)

    y = 4.999
    ch = chord(y, c0)
    print "c0: ",  xle(y, c0)/ch, xte(y, c0)/ch

    Mo = L*0.5
    print "Mo: ", Mo

    from scipy.integrate import quad
    integral = quad(chord, 0, b/2., c0)
    mcl = integral[0]/(b/2.)
    print "mcl, c0:, ",mcl, c0

    cM = Mo/(0.5*M**2*gamma*P*S*mcl)

    print "cM: ", cM


