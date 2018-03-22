import numpy as np


def get_panels(span, Npanels):

    nodes = np.linspace(0,span/2, Npanels+1)

    mid_points = np.zeros(nodes.size-1)

    for i in np.arange(0, mid_points.size):
        mid_points[i] = (nodes[i]+nodes[i+1])/2.

    return nodes, mid_points

def get_k(y, ybar):

    k_pos = np.zeros([y.size-1, y.size])
    k_neg = np.zeros([y.size-1, y.size])

    for i in np.arange(0, y.size-1):
        for j in np.arange(0, y.size):
            # print j
            k_pos[i,j] = (y[i+1]-y[i])/(ybar[i]-y[j])
            k_neg[i,j] =(y[i+1]-y[i])/(ybar[i]+y[j])

    return k_pos, k_neg

def get_dic(k_pos, k_neg, rho, Npanels):

    dic = np.zeros([Npanels, Npanels])
    # print k_neg.shape, k_pos.shape, dic.shape
    for i in np.arange(0, dic[0].size):
        for j in np.arange(0, dic[0].size):
            dic[i,j] = (rho/(2.*np.pi))*(k_pos[i,j]-k_pos[i, j+1] - k_neg[i,j] + k_neg[i, j+1])

    return dic

def get_lic(rho, v, y, Npanels):

    lic = np.zeros(Npanels)

    for i in np.arange(0, lic.size):
        lic[i] = 2.*rho*v*(y[i+1]-y[i])

    return np.array([lic])


def get_mic(rho, v, y, Npanels):

    mic = np.zeros(Npanels)

    for i in np.arange(0, mic.size):
        mic[i] = .5*rho*v*(y[i+1]**2-y[i]**2)

    return np.array([mic])

def get_gamma_ref(b, L, rho, V, ybar, span=15.):
    gamma0 = L / (rho * V * b * np.pi / 4.)
    gamma_ref = gamma0*np.sqrt(1.-(ybar/(span/2.))**2)
    return gamma_ref

def VLM(rho, v, L, span, Npanels, Mref=None):

    # get panel nodes and mid-points
    y, ybar = get_panels(span, Npanels)

    # calculate entries of K matrix
    k_pos, k_neg = get_k(y, ybar)

    # calculate values of induced drag matrix (DIC)
    dic = get_dic(k_pos, k_neg, rho, Npanels)

    # calculate values of lift matrix (LIC)
    lic = get_lic(rho, v, y, Npanels)

    # calculate values of root bending moment matrix (MIC)
    mic = get_mic(rho, v, y, Npanels)

    # get Mref if not provided
    if Mref is None:
        gamma_ref = get_gamma_ref(span, L, rho, v, ybar, span)
        Mref = np.zeros([1])
        Mref[0] = np.array([np.dot(mic[0], gamma_ref)])

    # get Lref
    Lref = np.ones([1])*L #np.matmul(lic[0], gamma_ref)

    # solve for entries of Gamma matrix
    ## build A
    A = np.concatenate((dic, lic.T, mic.T), axis=1)
    azeros = np.zeros([1,2])
    A = np.concatenate((A, np.concatenate((lic, azeros), axis=1)), axis=0)
    A = np.concatenate((A, np.concatenate((mic, azeros), axis=1)), axis=0)
    ## build b
    bzeros = np.zeros([1, Npanels])
    b = np.concatenate((bzeros[0], Lref, Mref), axis=0)
    ## solve the system
    x = np.linalg.solve(A, b)
    ## extract gamma
    gamma_solve = x[:Npanels]

    # calculate induced drag
    di = np.matmul(gamma_solve.T, np.matmul(dic, gamma_solve))

    # calculate lift
    Lprime = rho*v*gamma_solve

    # get span efficiency
    q = 0.5*rho*v**2
    einv = L**2/(q*np.pi*span**2*di)

    return di, Lprime, einv, ybar, Mref

def a_through_c():
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
    m = 2. * np.pi
    cl_max = 0.5

    Npanels = 100

    wingA = VLM(rho, V, L, b, Npanels)
    wingB = VLM(rho, V, L, b * 1.05, Npanels, Mref=wingA[4])

    dia = wingA[0]
    La = wingA[1]
    ea = wingA[2]
    ya = wingA[3]

    dib = wingB[0]
    Lb = wingB[1]
    eb = wingB[2]
    yb = wingB[3]

    print "Efficiency of wing A: ", ea
    print "Efficiency of wing B: ", eb

    print "Percent drag reduction: ", 100.*(dia-dib)/dia

    from matplotlib import pyplot as plt

    plt.plot(ya, La, label='Wing A')
    plt.plot(yb, Lb, label='Wing B')
    plt.xlabel('Distance from root (m)')
    plt.ylabel('$L\prime$ (N)')
    plt.legend(frameon=False)
    plt.savefig('p2_lift_distributions.pdf')
    plt.show()

def partd():
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
    be = 15.
    rho = 1.05
    gamma = 0.5
    m = 2. * np.pi
    cl_max = 0.5

    Npanels = 100

    wingA = VLM(rho, V, L, be, Npanels)
    die = wingA[0]

    Ntests = 100

    di = np.zeros(Ntests)
    spans = np.linspace(.85, 5.0, Ntests)

    for i in np.arange(0, Ntests):
        wingB = VLM(rho, V, L, be * spans[i], Npanels, Mref=wingA[4])
        di[i] = wingB[0]

    def jones(be, b):
        ratio = 8.*(be/b)**4 - 16.*(be/b)**3 + 9.*(be/b)**2
        return ratio

    di_ratio = jones(be, spans*be)

    from matplotlib import pyplot as plt

    plt.plot(1./spans, di/die, '--o', label='VLM', markerfacecolor='none')
    plt.plot(1./spans, di_ratio, label='R. T. Jones')
    plt.xlabel('Span Ratio ($b_e/b$)')
    plt.ylabel('Induced drag ratio ($D_{i}/D_{ie}$)')
    plt.legend(frameon=False, loc=2)
    plt.savefig('p2_induced_drag_ratio.pdf', transparent=True)
    plt.show()

if __name__ == "__main__":

    # a_through_c()
    partd()