import numpy as np

def get_pressure(T, R, Re, M, gamma, c, mu):

    p = np.sqrt(T)*np.sqrt(R)*Re*mu/(M*np.sqrt(gamma)*c)

    return p

def get_mu(T, T0, mu0, S):

    mu = mu0*(T/T0)**(3./2.)*(T0+S)/(T+S)

    return mu

def plot_airfoil():

    import matplotlib.pyplot as plt

    data = np.loadtxt('rae2822-geom.txt')
    x = data[:, 0]
    y = data[:, 1]

    plt.scatter(x, y)
    plt.axis('equal')
    plt.show()

def complete_points():

    data = np.loadtxt('rae2822-geom.txt')
    x = data[:, 0]
    y = data[:, 1]
    z = np.zeros(x.size)

    np.savetxt('rae2822-geom.txt', np.c_[x, y, z], header='x, y, z')


if __name__ == "__main__":

    # altitude = 18300 m
    plot_airfoil()

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

    a = np.sqrt(gamma*R*T)
    M = 0.729
    # M = 3.

    p = get_pressure(T, R, Re, M, gamma, c, mu)

    print M, p