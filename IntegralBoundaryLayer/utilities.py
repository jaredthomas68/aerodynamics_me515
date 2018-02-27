import numpy as np

def reynolds_number(V, L, nu=0, mu=0, rho=0):

    if nu == 0:
        R = rho*V*L/mu
    else:
        R = V*L/nu

    return R

def shape_factor(tlambda):

    if tlambda >= 0.0 and tlambda <= 0.1:
        L = 0.22+1.57*tlambda-1.8*tlambda**2
        H = 2.61 - 3.75*tlambda + 5.24*tlambda*2
    elif tlambda >= -0.1 and tlambda < 0.0:
        L = 0.22 + 1.402*tlambda + 0.018*tlambda/(0.107+tlambda)
        H = 2.088 + 0.0731/(0.14+tlambda)

    return H, L

def thwaite_boundary(flow_condition='laminar'):


    return displacement_thickness, momentum_thickness, shape_factor, local_skin_friction_coefficient

def head_boundary():

    return 0

def blasius_boundary():

    return 0

def schlichting_boundary():