import numpy as np

if __name__ == "__main__":

    # experimental values
    v1e = 2. #ft/s
    d1e = 2. #in
    d2e = 2.5 #in
    d3e = 4.5 #in

    v2e = v1e*d1e**2/(d3e**2-d2e**2) # m/s

    print v2e

    rhow = 2. # slugs/ft**3
    muw = 2.2E-5 #lbf*s/ft**2
    Ree = rhow*v2e*d2e*(1./12.)/muw
    print "Ree: ", Ree
    # real values
    d2r = 10. #in

    mua = 3.737E-7 # lbf*s/ft**2
    rhoa = 2.3769E-3 # slugs/ft**3

    Rer = Ree
    v2r = Ree*mua/(rhoa*d2r*(1./12.)) # ft/s

    print "v2r: ", v2r

    d1r = d1e*d2r/d2e
    print d1r

    d3r = d3e*d2r/d2e
    print d3r

    v1r = v2r*(d3r**2-d2r**2)/d1r**2
    print v1r

    p1r = 14.696 # lb/in**2
    # p2r = p1r + (0.5*rhoa*v1r**2-0.5*rhoa*v2r**2)*(1./12.)**2
    p2r = p1r + (-0.5*rhoa*v2r**2)*(1./12.)**2 #lb/in**2
    print "p2r: ", p2r

    A2 = (np.pi*d2r**2/4.) # in**2
    Dr = A2*(p2r-p1r)
    print "Dr: ", Dr