from utilities import *
from matplotlib import pylab as plt


def get_results():

    Vel = 1.0                   # free stream velocity
    nu = 1.57E-4                # kinematic viscosity
    Npoints = 1000               # how many point to include along the plate
    epsilon = 1E-6              # start a little past the leading edge
    x_max = 1E6*nu/Vel
    x = np.linspace(epsilon, x_max, Npoints)    # locations of interest along the plate
    save_figs = True

    L = x_max

    Rel = reynolds_number(Vel, L, nu=nu)
    Rex = reynolds_number(Vel, x, nu=nu)

    # print Rel, Rex

    # run Thwaite's method
    delta_star_T, theta_T, H_T, cf_T, CF_T = thwaite_boundary(x, L, Vel, Vel, Vel, nu, flat_plate=True, dve_dx0=0.0)

    # run Balsious' method
    delta_star_B, theta_B, H_B, cf_B, CF_B = blasius_boundary(x, Rex, Rel)

    # run Head's method
    H0 = 1.4
    theta0 = theta_B[0]
    delta_star_H, theta_H, H_H, cf_H, CF_H = head_boundary(x, Vel, Vel, 0., nu, theta0, H0)
    # print delta_star_H.size, x.size
    # quit()

    # run Schlichting's method
    delta_star_S, theta_S, H_S, cf_S, CF_S = schlichting_boundary(x, Rex, Rel)

    # assign plotting labels
    labels = ["Thwaite", "Blausius", "Head", "Schlichting"]

    # define x data
    x_data = [x, x, x, x]

    # define data to be plotted together
    dt_data = [delta_star_T, delta_star_B, delta_star_H, delta_star_S]
    mt_data = [theta_T, theta_B, theta_H, theta_S]
    cf_data = [cf_T, cf_B, cf_H, cf_S]
    h_data = [H_T, H_B, H_H, H_S]

    # set axis labels
    xlabel = 'X location along the plate'
    linestyles = ['-', '--', '-', '--']

    ylabel = 'Displacement thickness ($\delta^*$)'
    filename = "displacement-thickness.pdf"
    plot_results(x_data, dt_data, labels, xlabel=xlabel, ylabel=ylabel, legend_loc=2, linestyles=['-', '--', '-', '--'],
                 filename=filename, save_plot=save_figs)

    ylabel = r'Momentum thickness ($\theta$)'
    filename = 'momentum-thickness.pdf'
    plot_results(x_data, mt_data, labels, xlabel=xlabel, ylabel=ylabel, legend_loc=2, linestyles=['-', '--', '-', '--'],
                 filename=filename, save_plot=save_figs)

    ylabel = 'Local friction coefficient ($c_f$)'
    filename = 'local-friction-coefficient.pdf'
    plot_results(x_data, cf_data, labels, xlabel=xlabel, ylabel=ylabel, legend_loc=1, ylim=[0, .02], linestyles=['-', '--', '-', '--'],
                 filename=filename, save_plot=save_figs)

    ylabel = 'Shape factor ($H$)'
    filename = 'shape-factor.pdf'
    plot_results(x_data, h_data, labels, xlabel=xlabel, ylabel=ylabel, legend_loc=2, linestyles=['-', '--', '-', '--'],
                 filename=filename, save_plot=save_figs, ylim = [0, 4])

    return


if __name__ == "__main__":

    get_results()