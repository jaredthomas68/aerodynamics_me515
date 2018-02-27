import numpy as np
import matplotlib.pylab as plt
from matplotlib import animation
from pot_flow_func import *

def leap_frog():
    d = 1.
    start = np.array([10., 10.])

    time_step = 0.01
    steps = 5000

    #  Location of points
    #  1   2
    #  0   3

    loc = np.zeros([steps + 1, 4, 2])
    loc[0] = np.array(([start[0], start[1]],
                       [start[0], start[1] + d],
                       [start[0] + d, start[1] + d],
                       [start[0] + d, start[1]]))

    gamma = np.array([1., -1., -1., 1.])

    v = np.zeros_like(loc[0])

    for t in np.arange(0, steps):
        v[:] = 0.0
        for i in np.arange(0, 4):
            count = 0

            for j in np.arange(0, 4):
                if i == j:
                    continue
                else:
                    _, vx, vy = vortex_flow(loc[t, i, 0], loc[t, i, 1], loc[t, j, 0], loc[t, j, 1], gamma=gamma[j])
                    v[i, 0] += vx
                    v[i, 1] += vy

            # print v
        # print loc[t], v
        loc[t + 1] = loc[t] + time_step * v
        # print v, time_step*v

    print loc.shape
    plt.plot(loc[:, 0, 0], loc[:, 0, 1], 'r', label='Vortex Ring A')
    plt.plot(loc[:, 1, 0], loc[:, 1, 1], 'r')
    plt.plot(loc[:, 2, 0], loc[:, 2, 1], 'b', label='Vortex Ring B')
    plt.plot(loc[:, 3, 0], loc[:, 3, 1], 'b')
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    plt.axis('equal')
    plt.legend(frameon=False)
    plt.savefig('leap_frog.pdf', transparent=True)
    plt.show()

def visualize():

    x = np.arange(3,7.,0.15)
    y = np.arange(3,7.,0.5)

    xx, yy = np.meshgrid(x, y)

    d = 1.
    start = np.array([4.75, 4.75])

    #  Location of points
    #  1   2
    #  0   3

    loc = np.array(([start[0], start[1]],
                    [start[0], start[1] + d],
                    [start[0] + d, start[1] + d],
                    [start[0] + d, start[1]]))
    psi = np.zeros_like(xx)
    vxx = np.zeros_like(xx)
    vyy = np.zeros_like(xx)

    gamma = np.array([1., -1., -1., 1.])
    for i in np.arange(0,4):
        p, vx, vy = vortex_flow(xx, yy, loc[i, 0], loc[i, 1], gamma=gamma[i])

        psi += p
        vxx += vx
        vyy += vy

    # plt.contour(xx, yy, vxx, 100)
    # plt.contourf(xx, yy, vyy, 100)
    # plt.contourf(xx, yy, psi, 100)
    plt.contour(xx, yy, np.sqrt(vxx**2+vyy**2),100)
    plt.quiver(xx, yy, vxx, vyy)
    plt.scatter(loc[:,0], loc[:,1])
    plt.show()
    return


if __name__ == "__main__":

    leap_frog()
    # visualize()


    #
    # psi = np.zeros([4, y.size, x.size])
    # v = np.zeros([4, 2, y.size, x.size])
    #
    # # # print v[0,0,:].shape
    # # # print v[0,1,:].shape
    # # # print v[1,1,:].shape
    # # # print v[2,1,:].shape
    # # # print psi[0,:].shape
    # # #
    # # # quit()
    # #
    # # psi[0,:], v[0,0,:], v[0,1,:] = vortex_flow(xx, yy, loc[0,0], loc[0,1])
    # # psi[1], v[1,0], v[1,1] = vortex_flow(xx, yy, loc[1,0], loc[1,1])
    # # psi[2], v[2,0], v[2,1] = vortex_flow(xx, yy, loc[2,0], loc[2,1])
    # # psi[3], v[3,0], v[3,1] = vortex_flow(xx, yy, loc[3,0], loc[3,1])


    # v_tot = np.sqrt(vx**2+vy**2)

    # plt.contour(xx, yy, v)

    # plt.show()