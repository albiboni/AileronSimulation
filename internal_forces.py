from force import *
from Input_variables import *
import numpy as np
from inertia import *
from reaction_forces import *
'''
Calculates internal forces. functions have some unused arguments to facilitate plotting
'''

def get_displacement(position, C1, C2, int_constants, forces_y):
    '''given position on x axis, constants C1 and C2, and integrations constants.
    the output is displacement in y and z, in any x'''
    array_displacement_Mz = C1 * np.array([R_1y.displacement(position), R_2y.displacement(position),
                                           R_3y.displacement(position),
                                           q_coeff.displacement(position)])  # uses forces in y

    array_displacement_My = C2 * np.array([R_1z.displacement(position), R_2z.displacement(position),
                                           F_I.displacement(position),
                                           P_coeff.displacement(position)])  # uses forces in x

    displacement = np.dot(array_displacement_Mz, forces_y) + np.dot(array_displacement_My, forces_z) \
                   + int_constants[0]*position + int_constants[1]
    return displacement

hinge_3z_displacement = get_displacement(x_R_3y, C1, -C2, int_constants_z, forces_y) #displacement at hinge 3

'''Functions for moment, displacement, slope'''

def get_slope(position, C1, C2, int_constants):
    '''given position on x axis, constants C1 and C2, and integrations constants.
    the output is slope in y and z, in any x'''
    array_slope_Mz = C1 * np.array([R_1y.slope(position), R_2y.slope(position),
                                           R_3y.slope(position),
                                           q_coeff.slope(position)])  # uses forces in y

    array_slope_My = C2 * np.array([R_1z.slope(position), R_2z.slope(position),
                                           F_I.slope(position),
                                           P_coeff.slope(position)])  # uses forces in x

    slope = (np.dot(array_slope_Mz, forces_y) + np.dot(array_slope_My, forces_z)) + int_constants[0]
    return slope

def get_curvature(position, C1, C2, int_constants):
    '''given position on x axis, constants C1 and C2, and integrations constants.
    the output is curvature in y and z, in any x'''
    array_curvature_Mz = C1 * np.array([R_1y.moment(position), R_2y.moment(position),
                                           R_3y.moment(position),
                                           q_coeff.moment(position)])  # uses forces in y

    array_curvature_My = C2 * np.array([R_1z.moment(position), R_2z.moment(position),
                                           F_I.moment(position),
                                           P_coeff.moment(position)])  # uses forces in x

    curvature = (np.dot(array_curvature_Mz, forces_y) + np.dot(array_curvature_My, forces_z))
    return curvature

def get_moment_Mz(position, C1, C2, int_constants):
    array_Mz = np.array([R_1y.moment(position), R_2y.moment(position), R_3y.moment(position), q_coeff.moment(position)])
    moment = np.dot(array_Mz, forces_y)
    return moment

def get_moment_My(position, C1, C2, int_constants):
    array_My = np.array([R_1z.moment(position), R_2z.moment(position), F_I.moment(position), P_coeff.moment(position)])
    moment = np.dot(array_My, forces_z)
    return moment

def get_shear_Fy(position, C1, C2, int_constants):
    array_Fy = np.array([R_1y.coeff(position), R_2y.coeff(position), R_3y.coeff(position), q_coeff.coeff(position)])
    shear = -np.dot(array_Fy, forces_y)
    return shear

def get_shear_Fz(position, C1, C2, int_constants):
    array_Fz = np.array([R_1z.coeff(position), R_2z.coeff(position), F_I.coeff(position), P_coeff.coeff(position)])
    shear = -np.dot(array_Fz, forces_z)
    return shear

def get_moment_Mx(position, C1, C2, int_constants):
    '''
    shear center is assumed to be the hinge line because the actual rotation happens there
    negative moment means that it is a clockwise internal moment
    :param position:
    :param C1:
    :param C2:
    :param int_constants:
    :return:
    '''
    #shear_center = 0.20
    #hingetocenter= shear_center-h_a/2
    hinge = h_a/2
    theta_rad = theta*np.pi/180.

    M_q = -np.dot(rotation_matrix(theta),np.array([0., q]))[1]*position*(0.25*c_a-hinge)
    M_P = - P*(hinge*(np.cos(theta_rad)-np.sin(theta_rad)))*np.heaviside(position-x_P, 0.)
    #M_R3y = -np.dot(rotation_matrix(theta),np.array([0., forces_y[2]]))[1]*(hingetocenter)*np.heaviside(position-x_R_3y, 0.)
    M_FI = - forces_z[2]*(hinge*(np.cos(theta_rad)-np.sin(theta_rad)))*np.heaviside(position-x_F_I, 0.)
    #M_R2y = -np.dot(rotation_matrix(theta),np.array([0., forces_y[1]]))[1]*(hingetocenter)*np.heaviside(position-x_R_2y, 0.)
    #M_R1y = -np.dot(rotation_matrix(theta),np.array([0., forces_y[0]]))[1]*(hingetocenter)*np.heaviside(position-x_R_1y, 0.)
    #M_R2z = np.dot(rotation_matrix(theta), np.array([forces_z[1], 0.]))[1] * (hingetocenter) * np.heaviside(position - x_R_2y, 0.)
    #M_R1z = np.dot(rotation_matrix(theta), np.array([forces_z[0], 0.]))[1] * (hingetocenter) * np.heaviside(position - x_R_1y, 0.)
    M = M_q + M_P + M_FI  #+M_R3y + M_R2y + M_R1y#+ M_R1z+ M_R2z
    return M

def rotate_points(z,y):
    radians = theta * np.pi / 180.
    rotation_matrix = np.array([[np.cos(radians), -np.sin(radians)],
                                [np.sin(radians), np.cos(radians)]])
    matrix = np.dot(rotation_matrix, np.array([z,y]))
    z,y = matrix[0], matrix[1]
    return z,y


def get_normal_stress(position_x):
    Cz = centroid(Booms)
    position_y = Booms[:,1]
    position_z = Booms[:,0] - Cz
    moment_My = get_moment_My(position_x, 1., 1., 1.)
    moment_Mz = get_moment_Mz(position_x, 1., 1., 1.)
    normal_array = np.zeros((len(position_y),3))
    for i in range(len(position_y)):
        z,y = rotate_points(position_z[i], position_y[i])
        normal_stress = moment_Mz*(Iyy*y-Izy*z)/C4 + moment_My*(Izz*z- Izy*y)/C4
        normal_array[i] = np.array([z, y, normal_stress])
    return normal_array

'''Print max normal stress at the ribs
x_points= np.array([x_R_1y, x_F_I, x_P, x_R_2y])
n_points = 1000
#x_points = np.linspace(0, l_a, num=n_points)
y_points = np.zeros((len(x_points),4))
for i in range(np.size(x_points)):
    run = get_normal_stress(x_points[i])
    idx = np.argmax(run[:, 2])
    y_points[i] = np.append(x_points[i], run[idx])
print(y_points)

idx = np.argmax(y_points[:,3])
massimo = y_points[idx]
print(massimo)
for i in range(np.size(x_points)):
    run = get_normal_stress(x_points[i])
    idx = np.argmax(run[:, 2])
    while run[idx][2] >1.122e+09:
        run[idx][2]=0
        idx = np.argmax(run[:, 2])
        y_points[i] = np.append(x_points[i], run[idx])
    else:
        y_points[i] = np.append(x_points[i], run[idx])
idx = np.argmax(y_points[:,3])
massimo = y_points[idx]
print(massimo)
'''
'''Plot bending'''
import matplotlib.pyplot as plt
n_points = 1000
x_points = np.linspace(0, l_a, num=n_points)
y_points = np.zeros((n_points))
y_points1 = np.zeros((n_points))
y_points2 = np.zeros((n_points))

for i in range(np.size(x_points)):
    y_points[i] = get_displacement(x_points[i], -C3, C1, RF_y[3:],forces_y)
    y_points1[i] = get_displacement(x_points[i], -C3, C1, RF_y[3:],forces_y)-(h_a/2.*np.sin(theta*np.pi/180.))
    y_points2[i] = get_displacement(x_points[i], -C3, C1, RF_y[3:],forces_y) + ((c_a-h_a / 2.) * np.sin(theta * np.pi / 180.))

plt.plot(x_points,y_points2, label='trailing-edge')
plt.plot(x_points,y_points, label= 'hinge-line')
plt.plot(x_points,y_points1, label='leading-edge')
plt.title('Deflection of aileron in the y axis along the span')
plt.xlabel('Span-wise location[m], from hinge 3 to 1', fontsize=10)
plt.ylabel('Aileron deflection in y axis [m]', fontsize=10)
plt.legend()
plt.show()


