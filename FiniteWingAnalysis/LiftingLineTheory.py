import numpy as np

def plot_lift_and_cl(V, L, rho, cl_max, b, gamma):
    from matplotlib import pylab as plt

    def s_ellipse(V, L, rho, cl_max, b):
        s = 4. * L / (cl_max * rho * V ** 2 * b)
        return s
    def c_ellipse(y, b, V, L, rho, cl_max):
        c = 2.*L*np.sqrt(1.-(y/(b/2.))**2)/(cl_max*V**2*rho*b*np.pi/4.)
        return c
    def cl_ellipse(y, cl_max):
        cl = np.ones_like(y)*cl_max
        return cl
    def wing_area_ellipse(b, V, L, rho, cl_max):
        from scipy.integrate import quad
        results = quad(c_ellipse, 0, b/2., (b, V, L, rho, cl_max))
        return 2.*results[0]

    def lift(y, b, L):
        l = L*np.sqrt(1.-(y/(b/2.))**2)/(b*np.pi/4.)
        return l

    def c_rec(V, L, rho, b, cl_max):
        c = 2*L/(rho*V**2*cl_max*b*np.pi/4.)
        return c
    def s_rec(V, L, rho, cl_max):
        s = 2.*L/(rho*V**2*cl_max*np.pi/4.)
        return s
    def cl_rec(y, V, L, rho, b, c):
        cl = 2.*L*np.sqrt(1-(y/(b/2.))**2)/(rho*V**2*c*b*np.pi/4.)
        return cl
    def wing_area_rec(b, c):
        area = b*c
        return area


    # def cl_max_tap():
    #     2.*L*np.sqrt(1.-(y/(b/2.))**2)/(c*)
    def y_cl_max_tap(b, gamma):
        y_cl_max = 0.5*b*(1-gamma)
        return y_cl_max
    def cr_tap(yclmax, clmax, gamma, L, V, b, rho):
        gamma0 = L/(rho*V*b*np.pi/4.)
        cr = 2.*gamma0*np.sqrt(1.-(yclmax/(b/2.))**2)/((1.-(1.-gamma)*yclmax/(b/2.))*clmax*V)
        return cr
    def cl_tap(L, rho, V, b, y, cr):
        gamma0 = L / (rho * V * b * np.pi / 4.)
        cl = 2.*gamma0*np.sqrt(1.-(y/(b/2.))**2)/(V*cr*(1.-(1.-gamma)*y/(b/2.)))
        return cl
    def wing_area_tapered(b, cr, gamma):
        area = b*cr*(1.-(1.-gamma)/2.)
        return 2.*area
    def c_tap(y, b, gamma, cr):
        c = cr*(1.-(1.-gamma)*y/(b/2.))
        return c

    # for ellipse
    s_e = s_ellipse(V, L, rho, cl_max, b)
    y = np.linspace(0, b/2., 1000)
    c_e = c_ellipse(y, b, V, L, rho, cl_max)
    cl_e = cl_ellipse(y, cl_max)
    l_e = lift(y, b, L)
    a_e = wing_area_ellipse(b, V, L, rho, cl_max)



    # for rectangle
    c_r = c_rec(V, L, rho, b, cl_max)
    cl_r = cl_rec(y, V, L, rho, b, c_r)
    l_r = lift(y, b, L)
    s_r = s_rec(V, L, rho, cl_max)
    a_r = wing_area_rec(b, c_r)


    # for tapered
    y_cl_max_t = y_cl_max_tap(b, gamma)
    cr_t = cr_tap(y_cl_max_t, cl_max, gamma, L, V, b, rho)
    l_t = lift(y, b, L)
    cl_t = cl_tap(L, rho, V, b, y, cr_t)
    a_t = wing_area_tapered(b, cr_t, gamma)
    c_t = c_tap(y, b, gamma, cr_t)



    # plot cl

    fig, ax = plt.subplots()
    ax.plot(y, cl_e, label="Elliptic planform $c_l$")
    ax.plot(y, cl_r, label="Rectagular planform $c_l$")
    ax.plot(y, cl_t, label="Tapered planform $c_l$")
    plt.xlabel("Distance from root (m)")
    plt.ylabel("Section lift coefficient")
    ax.set_ylim([0, 0.75])
    plt.legend(frameon=False)
    plt.savefig('p1_section_cl.pdf', transparent=True)

    # plot normalized cl
    fig, ax = plt.subplots()
    ax.plot(y, cl_e*c_e/(a_e/b), label="Elliptic planform $c_l$")
    ax.plot(y, cl_r*c_r/(a_r/b), label="Rectagular planform $c_l$")
    ax.plot(y, cl_t*c_t/(a_t/b), label="Tapered planform $c_l$")
    plt.xlabel("Distance from root (m)")
    plt.ylabel("Normalized Section lift coefficient ($c_lc/\overline{c}$)")
    plt.legend(frameon=False)
    plt.savefig('p1_normalized_section_cl.pdf', transparent=True)

    # plot L prime
    fig, ax = plt.subplots()
    ax.plot(y, l_e, label="Elliptic planform $L^\prime$", linewidth=5)
    ax.plot(y, l_r, '--',  label="Rectangular planform $L^\prime$", linewidth=5)
    ax.plot(y, l_t, ':', label="Tapered planform $L^\prime$", linewidth=5)
    plt.xlabel("Distance from root (m)")
    plt.ylabel("Section lift")
    plt.legend(frameon=False)
    plt.savefig('p1_section_l.pdf', transparent=True)

    print "Ellipse Area: %0.3f" % a_e
    print "Rectangle Area: %0.3f" % a_r
    print "Tapered Area: %0.3f" % a_t

    plt.show()




if __name__ == "__main__":
    # design requirements
    # L: 25 kN
    # V: 60 m/s
    # b: 15 m
    # density: 1.05 kg/m**3
    # planforms: (1) elliptical, (2) rectangular, (3) tapered
    # taper ratio: 0.5
    # m: 2*pi (lift curve slope)
    # CL_max: 0.5

    # find:
    #   (a) L(y) and CL(y) that minimize induced drag for the given span
    #   (b) min wing area
    L = 25000.
    V = 60.
    b = 15.
    rho = 1.05
    gamma = 0.5
    m = 2.*np.pi
    cl_max = 0.5


    plot_lift_and_cl(V, L, rho, cl_max, b, gamma)


