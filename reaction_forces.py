from force import *
from Input_variables import *
import numpy as np
from inertia import *
'''
Calculates the reaction forces of the statically indeterminate structure using equation of motion 
and displacement relations
'''
#inertia coefficients
Iyy, Izz = get_inertia(Booms)
Iyy, Izz, Izy = rotate_inertia(theta, Iyy, Izz)

#Izz = 2.67475559e-05 #1.31296155e-05 #0.00002674754 test values
#Iyy =  6.93287501e-05 #3.04460444e-05  #0.0000693287
#Izy = 3.15646082e-05 #1.28363308e-05#-0.00003156458
C1 = Izy/(E*(Izz*Iyy-Izy**2.))
C2 = Izz/(E*(Izz*Iyy-Izy**2.))
C3 = Iyy/(E*(Izz*Iyy-Izy**2.))
C4 = (Izz*Iyy-Izy**2.)
#force definition in y
x_R_1y = l_a-x_1
x_R_2y = l_a-x_2
x_R_3y = l_a-x_3
R_1y = force(x_R_1y, sign=-1.)
R_2y = force(x_R_2y, sign=-1.)
R_3y = force(x_R_3y, sign=-1.)
q_coeff = distributed_load(0, sign=1.)

#force definition in z
x_F_I = l_a-(x_2-x_a/2.)
x_P = l_a-(x_2+x_a/2.)
R_1z = force(x_R_1y, sign=-1.)
R_2z = force(x_R_2y, sign=-1.)
F_I = force(x_F_I, sign=-1.)
P_coeff = force(x_P, sign=-1.)

'''Reaction forces in z from EoM'''

def reaction_z(angle):
    '''Given angle outputs R1z, R2z, FI '''
    theta = np.pi*angle/180.
    condition_1z = np.array([-1., -1., -1.]), np.array([P])
    condition_2z = np.array([(x_2 - x_1), 0., x_a / 2.]), np.array([P * (x_a / 2.)])
    condition_3z = np.array([d_1, 0., h_a * (-np.sin(theta) + np.cos(theta)) / 2.]), \
                   - np.array([(P * h_a) / 2. * (-np.sin(theta) + np.cos(theta)) + q * l_a * (
                           0.25 * c_a - h_a / 2.) * np.cos(theta)])
    unknowns_z = np.vstack((np.vstack((condition_1z[0], condition_2z[0])), condition_3z[0]))
    coeff_z = np.vstack((np.vstack((condition_1z[1], condition_2z[1])), condition_3z[1]))
    RF_z = np.linalg.solve(unknowns_z, coeff_z)
    return RF_z

RF_z = reaction_z(theta)
forces_z = np.append(RF_z, P) #R1z, R2z, FI, P

'''Reaction forces in y'''

def coefficients_boundary(position, C1, C2):
    array_displacement_Mz = C1 * np.array([R_1y.displacement(position), R_2y.displacement(position),
                                      R_3y.displacement(position), q_coeff.displacement(position)]) #uses forces in y

    array_displacement_My = C2 * np.array([R_1z.displacement(position), R_2z.displacement(position),
                                      F_I.displacement(position), P_coeff.displacement(position)]) #uses forces in x
    coeff_A = position
    return array_displacement_Mz, array_displacement_My, coeff_A

def solve_boundary_y(boundary_0, boundary_1, boundary_2, distance):

    coeff = - (np.dot(boundary_1, forces_z) + boundary_0[3] * q) + distance
    unknowns = np.append(boundary_0[:3], np.append(boundary_2, 1.))
    return unknowns, coeff

def reaction_forcesy( C1, C2):
    #1st boundary_condition
    # y = D3 @ x = la - x3
    boundary_h3y = coefficients_boundary(x_R_3y, C1, C2)
    #2nd boundary_condition
    # y = 0 @ x = la - x2
    boundary_h2y = coefficients_boundary(x_R_2y, C1, C2)
    #3rd boundary_condition
    # y = D1 @ x = la - x1
    boundary_h1y = coefficients_boundary(x_R_1y, C1, C2)



# Linear system for y-displacement





    condition_1y = solve_boundary_y(boundary_h1y[0], boundary_h1y[1], boundary_h1y[2], d_1) # order R_1y, R_2y, R_3y
    condition_2y = solve_boundary_y(boundary_h2y[0], boundary_h2y[1], boundary_h2y[2], 0.)
    condition_3y = solve_boundary_y(boundary_h3y[0], boundary_h3y[1], boundary_h3y[2], d_3)

    # Conditions from EoM
    condition_4y = np.array([1., 1., 1., 0., 0.]), np.array([q*l_a])
    condition_5y = np.array([-(x_2-x_1), 0., -(-x_3+x_2), 0., 0.]), np.array([q*l_a*(l_a/2-x_2)])

    unknowns_y = np.vstack((np.vstack((condition_1y[0], condition_2y[0])), np.vstack((condition_3y[0],
                 np.vstack((condition_4y[0], condition_5y[0]))))))
    coeff_y = np.vstack((np.vstack((condition_1y[1], condition_2y[1])), np.vstack((condition_3y[1],
                 np.vstack((condition_4y[1], condition_5y[1]))))))

    RF_y = np.linalg.solve(unknowns_y, coeff_y) # R_1y, R_2y, R_3y, A, B
    return RF_y

RF_y = reaction_forcesy(-C3, C1)
forces_y = np.append(RF_y[:3], q)



'''Displacement in z at hinge 3'''
#need to solve for boundary conditions

#1st boundary_condition
# y = 0 @ x = la - x1
boundary_h1z = coefficients_boundary(x_R_1y, C1, -C2)
#2nd boundary_condition
# y = 0 @ x = la - x2
boundary_h2z = coefficients_boundary(x_R_2y, C1, -C2)




def solve_boundary_z(boundary_0, boundary_1, boundary_2):

    coeff = - (np.dot(boundary_0, forces_y) + np.dot(boundary_1, forces_z))
    unknowns = np.append(boundary_2, 1.)
    return unknowns, coeff


condition_1z = solve_boundary_z(boundary_h1z[0], boundary_h1z[1], boundary_h1z[2])  # order A, B
condition_2z = solve_boundary_z(boundary_h2z[0], boundary_h2z[1], boundary_h2z[2])
unknowns_z = np.vstack((condition_1z[0], condition_2z[0]))
coeff_z = np.vstack((condition_1z[1], condition_2z[1]))
int_constants_z = np.linalg.solve(unknowns_z, coeff_z)  # A, B
