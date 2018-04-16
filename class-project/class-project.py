import numpy as np

def make_csv():

    data = np.loadtxt('cp-xfoil.txt')
    data2 = np.loadtxt('x-foil-cp-2.txt')

    np.savetxt('cp-xfoil.csv', np.c_[data[:,0], data[:, 1]], delimiter=', ')

    import matplotlib.pyplot as plt
    plt.plot(data[:,0], data[:,1])
    plt.plot(data2[:,0], data2[:,1], '--')
    plt.show()

def plot_results():
    data_xfoil = np.loadtxt('cp-xfoil.csv', delimiter=", ")
    data_cfd = np.loadtxt('CFDAbsolutePressureWithMesh2.csv', delimiter=",")
    data_wt = np.loadtxt('wt-p4.csv', delimiter=",")

    Pg = 1.03
    V1wt, rhowt, Pgwt = wt_vel(Pg)

    qwt = 0.5*rhowt*V1wt**2

    xwt = data_wt[:, 0]
    ywt = data_wt[:, 2]

    pwt = data_wt[:, 1]

    cpwt = (pwt-1.03)/qwt



    P = 1.01325E5   # Pa
    rho = 1.18415   # kg/m^3
    Re = 200000.
    mu = 1.75E-5
    c = 0.1524
    v = Re*mu/(rho*c)

    # L = 15.675
    L = 15.892
    # D = 0.6721
    D_t = 0.6265
    D_p = 0.2862

    q = 0.5*rho*v**2

    x = data_cfd[:, 0]/c
    p = (data_cfd[:, 1] - P)/q

    from scipy.integrate import trapz

    Lwt = -np.sum(pwt*(xwt[2]-xwt[1]))*c*np.cos(5.*np.pi/180.)
    Lwt = trapz(pwt,xwt)*c*np.cos(5.*np.pi/180.)
    Dwt = trapz(pwt,ywt)*c*np.sin(5.*np.pi/180.)

    cl_wt = Lwt/c
    cd_wt = Dwt/c

    print Lwt

    cl_cfd = L/(q*c)
    cd_cfd = D_p/(q*c)
    print "CL CFD: ", cl_cfd
    print "CD CFD: ", cd_cfd

    print "CL wt: ", cl_wt
    print "CD wt: ", cd_wt


    import matplotlib.pyplot as plt
    plt.plot(data_xfoil[:,0], data_xfoil[:,1], label='xfoil')
    plt.scatter(x, p, label='CFD')
    plt.plot(xwt, pwt, 's', label='WT')
    # plt.plot(data_wt[:,0], data_wt[:,1], label='Wind Tunnel')
    plt.xlabel('$x/c$')
    plt.ylabel('$C_P$')
    plt.legend(frameon=False)
    plt.gca().invert_yaxis()

    plt.savefig('cp.png', transparent=True)
    plt.show()

    return 0

def wt_vel(Pg):

    T = 22.  # Ce
    T += 273.15  # K
    mmHG = 72.8  # pressure
    R = 287.058  # J/(kgK)
    # P = 9706.  # Pa
    P = 89700.  # Pa

    Pg *= 248.84

    rho = P / (R * T)

    V1 = np.sqrt(2. * (Pg) / rho)

    print V1

    return V1, rho, Pg

if __name__ == "__main__":

    Pg = 1.03
    # wt_vel(Pg)

    plot_results()
    # make_csv()
    #
    # Re = 200000.
    # co = 0.1524
    # cx = 1.
    # # mu = 1.75E-5
    # # rho =
    #
    # Re_xfoil = Re*cx/co
    #
    # print Re_xfoil

