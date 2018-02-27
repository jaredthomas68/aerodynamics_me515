import numpy as np
from matplotlib import pylab as plt

def polar_to_cartesian_x(theta, r=1):

    x = (r*np.cos(theta) + 1.0)/2.

    return x

def camber_func(x, epsilon, p):

    y = np.zeros_like(x)
    dydx = np.zeros_like(x)

    for i, xv in zip(np.arange(0, x.size), np.abs(x)):
        # find the value of dy/dx
        if xv < 0.0:
            raise ValueError('x must be on the interval [0, 1]')
        if xv <= p:
            y[i] = (epsilon/p**2)*(2.*p*xv-xv**2)
            dydx[i] = (epsilon/p**2)*(2.*p-2.*xv)
        elif xv <= 1:
            y[i] = (epsilon/(1.-p)**2)*(1.-2.*p+2.*p*xv-xv**2)
            dydx[i] = (epsilon/(1-p)**2)*(2.*p-2.*xv)
        else:
            raise ValueError('x must be on the interval [0, 1]')

    return y, dydx

def thickness_func(x, tau, c):

    # calculate the location for a unit chord
    locs = np.abs(x)/c

    thickness = np.zeros_like(x)
    dtdx = np.zeros_like(x)

    for i, loc in zip(np.arange(0, x.size), locs):

        # calculate the thickness of the airfoil at x
        thickness[i] = 10.*tau*(0.2969*np.sqrt(loc)-0.1260*loc-0.3516*loc**2 + 0.2843*loc**3-0.1015*loc**4)

        # calculate the derivative of the thickness distribution with respect to x
        if loc > 0:
            dtdx[i] = 10.*tau*(0.5*0.2969*loc**(-0.5)-0.1260-2.*0.3516*loc + 3.*0.2843*loc**2-4.*0.1015*loc**3)
        else:
            dtdx[i] = 0.0

    return thickness, dtdx

def get_theta(x, y):

    # get the sin and cos of theta
    # sintheta, costheta = get_sin_cos(x, y)

    theta = np.zeros(x.size-1)
    # calculate the values of theta
    for i in np.arange(0, x.size-1):
        theta[i] = np.arctan2(y[i+1]-y[i], x[i+1]-x[i])

    return theta

def get_nodes(Npanels, epsilon, p, tau, c):

    # determine the number of nodes
    Nnodes = Npanels + 1

    # set x values
    theta = np.linspace(0, np.pi, Npanels/2. + 1)
    # print theta
    x_one_side = polar_to_cartesian_x(theta)
    # plt.scatter(x_one_side, np.zeros(x_one_side.size))
    # plt.show()
    x = np.hstack([x_one_side, np.flip(x_one_side[:-1], 0)])
    sign = np.hstack([-np.ones_like(x_one_side), np.ones_like(x_one_side[:-1])])

    # find camber line
    y_camber, _ = camber_func(x, epsilon, p)

    # find thickness
    thickness, _ = thickness_func(x, tau, c)

    y = np.zeros(x.size)

    for i in np.arange(0, Npanels+1):

        y[i] = y_camber[i] + sign[i]*thickness[i]

    return np.abs(x), y

def get_xbar_ybar(x, y):

    Npanels = x.size - 1

    xbar = np.zeros(Npanels)
    ybar = np.zeros(Npanels)

    for i in np.arange(0, Npanels):
        xbar[i] = (x[i]+x[i+1])/2.
        ybar[i] = (y[i]+y[i+1])/2.

    return xbar, ybar

def get_panel_lengths(x, y):

    Npanels = x.size - 1

    length = np.zeros(Npanels)

    for i in np.arange(0, Npanels):
        length[i] = np.sqrt((y[i+1]-y[i])**2 + (x[i+1]-x[i])**2)

    return length

def get_sin_cos(x, y):

    Npanels = x.size - 1
    sintheta = np.zeros(Npanels)
    costheta = np.zeros(Npanels)

    length = get_panel_lengths(x, y)

    for i in np.arange(0, Npanels):
        sintheta[i] = (y[i+1]-y[i])/length[i]
        costheta[i] = (x[i+1]-x[i])/length[i]

    return sintheta, costheta

def get_tangent_normal_venctors(x, y):

    Npanels = x.size - 1

    tangent = np.zeros([Npanels, 2])
    normal = np.zeros([Npanels, 2])

    sintheta, costheta = get_sin_cos(x, y)

    for i in np.arange(0, Npanels):
        tangent[i] = np.array([costheta[i], sintheta[i]])
        normal[i] = np.array([-sintheta[i], costheta[i]])

    return tangent, normal

def get_beta(x, y, xbar, ybar):

    # get the number of panels
    Npanels = x.size - 1

    # initialize beta
    beta = np.zeros([Npanels, Npanels])

    # loop over rows of beta
    for i in np.arange(0, Npanels):
        # loop over columns of beta
        for j in np.arange(0, Npanels):

            # calculate beta
            if i == j:
                beta[i, j] = np.pi
            else:
                beta[i, j] = np.arctan2(((x[j]-xbar[i])*(y[j+1]-ybar[i])-(y[j]-ybar[i])*(x[j+1]-xbar[i])),
                                        ((x[j]-xbar[i])*(x[j+1]-xbar[i])+(y[j]-ybar[i])*(y[j+1]-ybar[i])))

    return beta

def get_r(x, y, xbar, ybar):

    # get the number of panels
    Npanels = x.size - 1

    # initialize r
    r = np.zeros([Npanels, Npanels+1])

    # loop over rows of r
    for i in np.arange(0, Npanels):
        # loop over columns of r
        for j in np.arange(0, Npanels+1):

            # calculate r
            r[i,j] = np.sqrt((xbar[i]-x[j])**2 + (ybar[i] - y[j])**2)

    return r

def panel_solve(V0 = 10., alpha=0.0, epsilon=0.04, p=0.4, c=1.0, tau=0.12, Npanels=100):

    x, y = get_nodes(Npanels, epsilon, p, tau, c)

    # x = np.array([1,   1.00000000000000e-05,    1.00000000000000e-05,    1])
    # y = np.array([-0.00125999999999998,    -0.000560572146286566 ,   0.000564572096286566 ,   0.00125999999999998])

    # find number of panels
    Npanels = x.size - 1

    # initialize coefficient matrix (A)
    A = np.zeros([Npanels+1, Npanels+1])

    # initialize result vector (b)
    b = np.zeros(Npanels+1)

    # calculate the location of the midpoint of each panel
    xbar, ybar = get_xbar_ybar(x, y)

    # solve for theta
    theta = get_theta(x, y)
    # print ybar
    # quit()
    # solve for r
    r = get_r(x, y, xbar, ybar)

    # solve for beta
    beta = get_beta(x, y, xbar, ybar)

    # loop over row in A and b
    for i in np.arange(0, Npanels):

        # loop over columns of A
        for j in np.arange(0, Npanels):

            # solve for the internal elements of A
            A[i, j] = np.sin(theta[i]-theta[j])*np.log(r[i,j+1]/r[i,j]) + np.cos(theta[i] - theta[j])*beta[i,j]

            # solve for last element in each internal row of A
            A[i, Npanels] += np.cos(theta[i] - theta[j]) * np.log(r[i, j + 1] / r[i, j]) - np.sin(theta[i] - theta[j]) * beta[i, j]

        # solve for elements of b
        b[i] = 2.*np.pi*V0*np.sin(theta[i]-alpha)

    # loop over columns of A
    for j in np.arange(0, Npanels):

        # solve for the last element in each internal column of A
        for k in np.array([0, Npanels-1]):
            A[Npanels, j] += np.sin(theta[k] - theta[j]) * beta[k, j] - np.cos(theta[k] - theta[j]) * np.log(r[k, j + 1] / r[k, j])

    # solve for A[N, N]
    for k in np.array([0, Npanels-1]):
        for j in np.arange(0, Npanels):
            A[Npanels, Npanels] += np.cos(theta[k] - theta[j]) * beta[k, j] + np.sin(theta[k] - theta[j]) * np.log(r[k, j + 1] / r[k, j])

    # solve for b[Npanels]
    b[Npanels] = - 2.*np.pi*V0*(np.cos(theta[0]-alpha) + np.cos(theta[Npanels-1] - alpha))

    # solve for q and gamma (A*[q,gamma] = b
    X = np.linalg.solve(A, b)

    # extract q
    q = X[0:Npanels]

    # extract gamma
    gamma = X[-1]

    # initialize Vt
    Vt = np.zeros(Npanels)

    # initialize Cp
    Cp = np.zeros(Npanels)

    # loop over the panels
    for i in np.arange(0, Npanels):

        # calculate the part of Vt that is only dependent on i
        Vt[i] = V0 * np.cos(theta[i] - alpha)

        # loop over the panels
        for j in np.arange(0, Npanels):

            # calculate Vt
            Vt[i] += (1./(2.*np.pi))*q[j]*(np.sin(theta[i]-theta[j])*beta[i,j]-np.cos(theta[i]-theta[j])*np.log(r[i,j+1]/r[i,j]))
            Vt[i] += (gamma/(2.*np.pi))*(np.cos(theta[i]-theta[j])*beta[i,j]+np.sin(theta[i]-theta[j])*np.log(r[i,j+1]/r[i,j]))

        # calculate Cp
        # calculate Cp
        Cp[i] = 1. - (Vt[i]/V0)**2

    # calculate cl
    cl = 0.0
    cmac = 0.0
    cd = 0.0
    length = get_panel_lengths(x, y)
    for i in np.arange(0, Npanels):
        cl -= Cp[i]*length[i]*np.cos(theta[i]-alpha)
        cd -= Cp[i]*length[i]*np.sin(theta[i]-alpha)
        cmac -= (0.25 - x[i])*Cp[i]*length[i]*np.cos(theta[i]-alpha) + y[i]*Cp[i]*length[i]*np.sin(theta[i]-alpha)


    return x, y, xbar, ybar, Cp, cl, cmac, cd

def plot_cl_cmac(save_data):
    Npanels = 100
    # alpha = 0.00898577292547  # ideal angle of attack, where A0 = 0.0
    alpha = 0.0866  # ideal angle of attack, where A0 = 0.0
    # alpha = 0.0
    epsilon = 0.04
    p = 0.4
    c = 1.0
    tau = 0.12
    V0 = 10.



    samples = 10
    cl = np.zeros(samples)
    cm = np.zeros(samples)
    alpha = np.linspace(-np.pi / 8., np.pi / 8., samples)
    for i, a in zip(np.arange(0, samples), alpha):
        x, y, xbar, ybar, Cp, cl[i], cm[i], cd = panel_solve(V0, a, epsilon, p, c, tau, Npanels)

    print (cl[1]-cl[0])/(alpha[1]-alpha[0])
    from scipy.interpolate import UnivariateSpline

    # spl = UnivariateSpline(cl, alpha)
    # alpha_cl1 = spl(1.)
    #
    # spl1 = UnivariateSpline(alpha, cl)
    # print spl1(alpha_cl1), alpha_cl1

    if save_data:
        np.savetxt("panel_cl_cmac.txt", np.c_[alpha, cl, cm], header="alpha, cl, cm")

    plt.plot(alpha, cl, label='$C_L$')
    plt.plot(alpha, cm, label='$C_{mac}$')
    plt.legend(loc=2)
    plt.xlabel('Angle of Attack (radians)')
    plt.show()

def run_panel(save_data):

    Npanels = 100
    # alpha = 0.00898577292547  # ideal angle of attack, where A0 = 0.0
    alpha = 0.0866  # ideal angle of attack, where A0 = 0.0
    # alpha = 0.0
    epsilon = 0.04
    p = 0.4
    c = 1.0
    tau = 0.12
    V0 = 10.

    x, y, xbar, ybar, Cp, cl, cmac, cd = panel_solve(V0, alpha, epsilon, p, c, tau, Npanels)

    plt.plot(xbar, Cp)

    plt.gca().invert_yaxis()
    plt.show()

    x, y = get_nodes(Npanels, epsilon, p, tau, c)

    # x = np.array([1, 1.00000000000000e-05, 1.00000000000000e-05, 1])
    # y = np.array([-0.00125999999999998, -0.000560572146286566, 0.000564572096286566, 0.00125999999999998])

    xbar, ybar = get_xbar_ybar(x, y)

    tangent, normal = get_tangent_normal_venctors(x, y)

    theta = get_theta(x, y)
    # print theta
    print Cp.shape, x.shape
    if save_data:
        np.savetxt("panel_cp_alphac%.3f.txt" % alpha, np.c_[xbar, Cp], header="x, Cp")
    plt.plot(x, y)
    plt.scatter(x, y)
    plt.axis('equal')
    plt.show()

def plot_convergence(save_data):
    Npanels = 100
    alpha = 0.0866455743735  # angle of attack where Cl=1.0
    # alpha = 0.00898577292547 # ideal angle of attack, where A0 = 0.0
    alpha = 0.0
    epsilon = 0.04
    p = 0.4
    c = 1.0
    tau = 0.12
    V0 = 10.

    samples = 100
    cl = np.zeros(samples)
    cmac = np.zeros(samples)
    cd = np.zeros(samples)
    Npanels = np.arange(0, samples)*10+2
    # print Npanels
    for i, N in zip(np.arange(0, samples), Npanels):
        print i, N
        x, y, xbar, ybar, Cp, cl[i], cmac[i], cd[i] = panel_solve(V0, alpha, epsilon, p, c, tau, N)

    from scipy.interpolate import UnivariateSpline

    # spl = UnivariateSpline(cl, alpha)
    # alpha_cl1 = spl(1.)
    #
    # spl1 = UnivariateSpline(alpha, cl)
    # print spl1(alpha_cl1), alpha_cl1

    if save_data:
        np.savetxt("panel_convergence_cl_cmac_cd_alpha%.3f.txt" % alpha, np.c_[Npanels, cl, cmac, cd], header="panels, cl, cmac, cd")

    plt.plot(Npanels, cl)
    plt.plot(Npanels, cd)
    plt.plot(Npanels, cmac)
    # plt.plot(alpha, cl, label='$C_L$')
    # plt.plot(alpha, cm, label='$C_{mac}$')
    # plt.legend(loc=2)
    # plt.xlabel('Angle of Attack (radians)')
    # plt.show()

def plot_airfoil(save_data):

    Npanels = 100
    alpha = 0.00898577292547  # ideal angle of attack, where A0 = 0.0
    # alpha = 0.0
    epsilon = 0.04
    p = 0.4
    c = 1.0
    tau = 0.12

    x, y = get_nodes(Npanels, epsilon, p, tau, c)
    # print x

    plt.plot(x, y, label="Panels")
    plt.scatter(x, y, label="Nodes")
    plt.axis('equal')
    plt.xlabel("x location")
    plt.ylabel("y location")
    plt.legend()
    plt.show()

if __name__ == "__main__":

    save_data = False

    run_panel(save_data)
    # plot_convergence(save_data)
    # plot_cl_cmac(save_data)
    # plot_airfoil(save_data)