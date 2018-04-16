import numpy as np

def get_pressure(T, R, Re, M, gamma, c, mu):

    p = np.sqrt(T)*np.sqrt(R)*Re*mu/(M*np.sqrt(gamma)*c)

    return p

def get_mu(T, T0, mu0, S):

    mu = mu0*(T/T0)**(3./2.)*(T0+S)/(T+S)

    return mu

def plot_airfoil():

    import matplotlib.pyplot as plt

    data = np.loadtxt('rae2822-geom2.csv', delimiter=',')
    x = data[:, 0]
    y = data[:, 1]

    plt.plot(x, y)
    plt.axis('equal')
    plt.show()

def get_blayer_thickness():
    from math import log
    T = 0.016
    y = 4.2927E-6
    s = 1.2
    n = log((T / y) * (s - 1.) + 1., 1.2)
    return n


def complete_points():

    data = np.loadtxt('rae2822-geom.txt')
    n = data.shape[0]
    x1 = data[:n/2, 0]
    x2 = data[n/2:, 0]
    x = np.hstack([x1[:-1],np.flip(x2, 0)])
    y1 = data[:n/2, 1]
    y2 = data[n/2:, 1]
    y = np.hstack([y1[:-1], np.flip(y2, 0)])
    z = np.zeros(n-1)

    np.savetxt('rae2822-geom2.csv', np.c_[x, y, z], delimiter=', ', header='x, y, z')

def plot_cp():
    n = get_blayer_thickness()
    print n
    c = 1.

    Re = 6.5E6
    R = 286.9
    gamma = 1.4
    rho = 1.1606E-1

    mu0 = 1.716E-5
    T0 = 273.15
    T = 216.66
    S = 110.4

    mu = get_mu(T, T0, mu0, S)
    nu = mu / rho

    a = np.sqrt(gamma * R * T)
    M = 0.729
    v = M * a

    p = get_pressure(T, R, Re, M, gamma, c, mu)

    print "Velocity: ", v
    print "Dynamic viscosity: ", mu
    print "Kinematic viscosity: ", nu
    print "Mach: ", M
    print "Velocity: ", v
    print "Speed of sound: ", a
    print "Temperature: ", T
    print "Pressure: ", p
    print "Density: ", rho
    print "Re: ", Re

    q = (gamma / 2.) * p * M ** 2
    # L = 7157.
    # L = 6700.
    L = 6698.4
    # D = 172.
    # D = 152.
    D = 142.3
    print "CL: ", L / (q * c)
    print 'CD: ', D / (q * c)

    data = np.loadtxt("2822_cp.txt")
    cfd = np.loadtxt("cfd_results_p.csv", delimiter=',')

    import matplotlib.pyplot as plt
    plt.scatter(cfd[:,0], cfd[:,1]/q, color='b', label='CFD')
    plt.scatter(data[:,0], -data[:,1], color='r', label='Experimental Data')
    plt.xlabel('Relative Chord Position')
    plt.ylabel('Pressure Coefficient')
    ax = plt.gca()
    ax.set_ylim(ax.get_ylim()[::-1])
    plt.legend(loc=4, frameon=False)
    plt.savefig('cfd_cp.pdf', transparent=True)
    plt.show()

if __name__ == "__main__":
    plot_cp()
    # altitude = 18300 m
    # complete_points()
    # plot_airfoil()
    n = get_blayer_thickness()
    print n
    c = 1.

    Re = 6.5E6
    R = 286.9
    gamma = 1.4
    rho = 1.1606E-1

    mu0 = 1.716E-5
    T0 = 273.15
    T = 216.66
    S = 110.4

    mu = get_mu(T, T0, mu0, S)
    nu = mu/rho

    a = np.sqrt(gamma*R*T)
    M = 0.729
    v = M*a

    p = get_pressure(T, R, Re, M, gamma, c, mu)

    print "Velocity: ", v
    print "Dynamic viscosity: ", mu
    print "Kinematic viscosity: ", nu
    print "Mach: ", M
    print "Velocity: ", v
    print "Speed of sound: ", a
    print "Temperature: ", T
    print "Pressure: ", p
    print "Density: ", rho
    print "Re: ", Re


    q = (gamma/2.)*p*M**2
    # L = 7157.
    # L = 6700.
    L = 6697.16
    # L = 7760
    # D = 172.
    # D = 152.
    D = 142.2
    print "CL: ", L/(q*c)
    print 'CD: ', D/(q*c)