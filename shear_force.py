import numpy as np
from Boom import *
from inertia import *
from reaction_forces import *
from internal_forces import *
import matplotlib.pyplot as plt
Iyy, Izz = get_inertia(Booms) #coordinate frame crosssection
Cz = centroid(Booms)
'''
Calculates Shear flow in the structure
'''


def get_forces(position, theta):
    '''
    this function comuptes the sum of forces in y and z in the cross_section reference frame at given x position
    :param position: position in x direction in global reference frame from hinge3 --> hinge1
    :param theta: rotation angle of cross_section
    :return: sum of forces in y and z in the cross_section reference system at the specific x position
    '''
    rotate_y = np.dot(rotation_matrix(theta),np.array([0., get_shear_Fy(position, 1,1,1)]))
    rotate_z = np.dot(rotation_matrix(theta),np.array([get_shear_Fz(position, 1,1,1), 0.]))
    Vy = rotate_y[1]+ rotate_z[1]
    Vz = rotate_y[0] + rotate_z[0]
    return Vy, Vz


def get_qb(NSpS, Iyy, Izz, Vz1, Vy1, BrA):
    '''
    Get open section shear flow q_B with imaginary cut.
    :param NSpS: half number of booms
    :param Iyy: inertia in y
    :param Izz: inertia in z
    :param Vz1: sum of forces in z
    :param Vy1: sum of forces in y
    :param BrA: boom area
    :return: shear flow q_B
    '''

    Qb = np.zeros(NSpS+3)


    Qb[7] = -Vz1/Iyy*BrA[7]*Sty[7] - Vy1/Izz*BrA[7]*(StZ[7] - Cz)  + Qb[8]
    Qb[6] = -Vz1/Iyy*BrA[6]*Sty[6] - Vy1/Izz*BrA[6]*(StZ[6] - Cz)  + Qb[7]
    Qb[5] = -Vz1/Iyy*BrA[5]*Sty[5] - Vy1/Izz*BrA[5]*(StZ[5] - Cz)  + Qb[6]
    Qb[4] = -Vz1/Iyy*BrA[4]*Sty[4] - Vy1/Izz*BrA[4]*(StZ[4] - Cz)  + Qb[5]
    Qb[3] = -Vz1/Iyy*BrA[3]*Sty[3] - Vy1/Izz*BrA[3]*(StZ[3] - Cz)  + Qb[4]
    Qb[1] = (-Vz1/Iyy*BrA[1]*Sty[1] - Vy1/Izz*BrA[1]*(StZ[1] - Cz))  + Qb[0]
    Qb[2] = (-Vz1/Iyy*BrA[2]*Sty[2] - Vy1/Izz*BrA[2]*(StZ[2] - Cz))  + Qb[1]
    Qb[9] = (-Vz1/Iyy*BrA[8]*Sty[8] - Vy1/Izz*BrA[8]*(StZ[8] - Cz))  + Qb[1] - Qb[2]

    return Qb


def get_Qs01(ha, q, Sa, skt, spt, s): #closed section shear flow 1 for the LE of the Aileron
    '''
    closed section shear flow 1 for the semicircle part
    :param ha: aileron height [m]
    :param q: open section shear flow
    :param Sa: half of the surface distance of the Aileron
    :param skt: skin thickness
    :param spt: Spar thickness
    :param s: Stringer Location
    :return: closed section shear flow 1 for the semicircle part
    '''

    q1c1 = 2*(q[1]*(s[8] - s[1])) + q[9]*ha*skt/spt
    q1c2 = -(2*((s[8] - s[0])) + ha*skt/spt)

    q2c1 = 2*(q[3]*(s[3] - s[2]) + q[2]*(s[2] - s[8]) + q[4]*(s[4] - s[3]) + q[5]*(s[5] - s[4]) + q[6]*(s[6] - s[5]) + q[7]*(s[7] - s[6]) + q[8]*(Sa - s[7])) - q[9]*ha*skt/spt
    q2c2 = -(2*(Sa - s[8]) + ha*skt/spt)
    h3 = ha*skt/spt
    qs02 = (q2c1 - (q1c1/q1c2)*h3)/(q2c2 + h3**2/q1c2)
    qs01 = (q1c1 - qs02*h3)/q1c2
    return qs01


def get_Qs02(ha, q, Sa, skt, spt, s):
    '''
    closed section shear flow 2 for the triangular part
    :param ha: aileron height [m]
    :param q: open section shear flow
    :param Sa: half of the surface distance of the Aileron
    :param skt: skin thickness
    :param spt: Spar thickness
    :param s: Stringer Location
    :return: closed section shear flow 2 for the triangular part
    '''

    q1c1 = 2*(q[1]*(s[8] - s[1])) + q[9]*ha*skt/spt
    q1c2 = -(2*((s[8] - s[0])) + ha*skt/spt)

    q2c1 = 2*(q[3]*(s[3] - s[2]) + q[2]*(s[2] - s[8]) + q[4]*(s[4] - s[3]) + q[5]*(s[5] - s[4]) + q[6]*(s[6] - s[5]) + q[7]*(s[7] - s[6]) + q[8]*(Sa - s[7])) - q[9]*ha*skt/spt
    q2c2 = -(2*(Sa - s[8]) + ha*skt/spt)
    h3 = ha*skt/spt
    qs02 = (q2c1 - (q1c1/q1c2)*h3)/(q2c2 + h3**2/q1c2)
    qs01 = (q1c1 - qs02*h3)/q1c2
    return qs02

def get_q_shear(NSpS, Iyy, Izz, Vz1, Vy1, BrA, qs01, qs02):
    '''
    Get open section shear flow q_B with imaginary cut.
    :param NSpS: half number of booms
    :param Iyy: inertia in y
    :param Izz: inertia in z
    :param Vz1: sum of forces in z
    :param Vy1: sum of forces in y
    :param BrA: boom area
    :return: shear flow q_B
    '''

    Qb = np.zeros((NSpS+3))
    Qb[0] = qs01
    Qb[8] = qs02
    Qb[7] = -Vz1/Iyy*BrA[7]*Sty[7] - Vy1/Izz*BrA[7]*(StZ[7] - Cz)  + Qb[8]
    Qb[6] = -Vz1/Iyy*BrA[6]*Sty[6] - Vy1/Izz*BrA[6]*(StZ[6] - Cz)  + Qb[7]
    Qb[5] = -Vz1/Iyy*BrA[5]*Sty[5] - Vy1/Izz*BrA[5]*(StZ[5] - Cz)  + Qb[6]
    Qb[4] = -Vz1/Iyy*BrA[4]*Sty[4] - Vy1/Izz*BrA[4]*(StZ[4] - Cz)  + Qb[5]
    Qb[3] = -Vz1/Iyy*BrA[3]*Sty[3] - Vy1/Izz*BrA[3]*(StZ[3] - Cz)  + Qb[4]
    Qb[1] = (-Vz1/Iyy*BrA[1]*Sty[1] - Vy1/Izz*BrA[1]*(StZ[1] - Cz))  + Qb[0]
    Qb[2] = (-Vz1/Iyy*BrA[2]*Sty[2] - Vy1/Izz*BrA[2]*(StZ[2] - Cz))  + Qb[1]
    Qb[9] = (-Vz1/Iyy*BrA[8]*Sty[8] - Vy1/Izz*BrA[8]*(StZ[8] - Cz))  + Qb[1] - Qb[2]
    #rest= -Qb[:-1] #better positive
    #q_total = np.append(Qb, rest)

    return Qb


def Shear_moment_arm(NSpS, StZ, Sty, ha, Qsc):  # finding the moment arm for the moment equation for shear flows.
    q = Qsc

    Qsct = 0
    Qscy = np.zeros(NSpS + 3)
    Qscy[0] = 2 * Qsc[0] * (StZ[1] - StZ[0]) * Sty[0]
    Qscy[1] = 2 * Qsc[1] * np.abs(StZ[8] - StZ[1]) * Sty[1]
    Qscy[2] = 2 * Qsc[2] * np.abs(StZ[2] - StZ[8]) * Sty[8]
    Qscy[3] = 2 * Qsc[3] * np.abs(StZ[3] - StZ[2]) * Sty[2]
    Qscy[4] = 2 * Qsc[4] * np.abs(StZ[4] - StZ[3]) * Sty[3]
    Qscy[5] = 2 * Qsc[5] * np.abs(StZ[5] - StZ[4]) * Sty[4]
    Qscy[6] = 2 * Qsc[6] * np.abs(StZ[6] - StZ[5]) * Sty[5]
    Qscy[7] = 2 * Qsc[7] * np.abs(StZ[7] - StZ[6]) * Sty[6]
    Qscy[8] = 0
    Qscy[9] = 0

    Qscz = np.zeros(NSpS + 3)
    Qscz[0] = 2 * Qsc[0] * (Sty[1] - Sty[0]) * (StZ[0] - ha / 2)
    Qscz[1] = 2 * Qsc[1] * (Sty[8] - Sty[1]) * (StZ[1] - ha / 2)
    Qscz[2] = 2 * Qsc[2] * (Sty[2] - Sty[8]) * (StZ[8] - ha / 2)
    Qscz[3] = 2 * Qsc[3] * (Sty[3] - Sty[2]) * (StZ[2] - ha / 2)
    Qscz[4] = 2 * Qsc[4] * (Sty[4] - Sty[3]) * (StZ[3] - ha / 2)
    Qscz[5] = 2 * Qsc[5] * (Sty[5] - Sty[4]) * (StZ[4] - ha / 2)
    Qscz[6] = 2 * Qsc[6] * (Sty[6] - Sty[5]) * (StZ[5] - ha / 2)
    Qscz[7] = 2 * Qsc[7] * (Sty[7] - Sty[6]) * (StZ[6] - ha / 2)
    Qscz[8] = 0
    Qscz[9] = Qsc[9] * (ha) * (StZ[8] - ha / 2)

    for i in xrange(0, NSpS + 3, 1):
        Qsct += Qscy[i] + Qscz[i]

    return Qsct


def ShearCenterZ(Qsct, qs01, qs02, Ac, At, ha, Vy):
    VZsc = ((Qsct + 2*Ac*qs01 + 2*At*qs02))/Vy - ha/2
    #VYsc = ((Qsct + 2*Ac*qs01 + 2*At*qs02))/Vzsc
    #print "this should be 0 for shear center Y: ", VYsc/Vzsc
    return VZsc

def get_q_torque(position):
    Torque = -get_moment_Mx(position, 1., 1., 1.)
    triangle_term = (2*Lsl/t_sk+h_a/t_sp)/(2.*At*G)
    circle_term = (CirC/2*t_sk + h_a/t_sp)/(2*Ac*G)
    A_torque = np.array([[1., 1., 0, 0],
                  [-1., 0, 2. * Ac, 0],
                  [0, -1., 0, 2. * At],
                  [0, 0, circle_term, -triangle_term]])

    b_torque = np.array([Torque, 0., 0.,0.])
    sol_torq = np.linalg.solve(A_torque, b_torque) #Tc, Tt, qTc, qTt
    return sol_torq


def get_twist(position):
    q = get_q_torque(position)[3]
    circle_term = (CirC / 2 * t_sk + h_a / t_sp) / (2 * Ac * G)
    triangle_term = (2 * Lsl / t_sk + h_a / t_sp) / (2. * At * G)
    theta_distribution = q*triangle_term
    return theta_distribution*180./np.pi

'''Print shear center
Vy, Vz = get_forces(1.5, theta)
qb = get_qb(NSpS, Iyy, Izz, Vz, Vy, BrA)
qs01 = get_Qs01(h_a, qb, Sa, t_sk, t_sp, StrinLoc)
qs02 = get_Qs02(h_a, qb, Sa, t_sk, t_sp, StrinLoc)
moment_arm = Shear_moment_arm(NSpS, StZ, Sty, h_a, qb)
print(ShearCenterZ(moment_arm, qs01, qs02, Ac, At, h_a, Vy))
'''

'''Print Twist
n_points = 1000
x_points = np.linspace(0, l_a, num=n_points)
step = 2.6e-3
y_points = np.zeros((n_points))
start_angle = 28
for i in range(np.size(x_points)):

    start_angle += get_twist(x_points[i])*step
    y_points[i] += start_angle

plt.plot(x_points,y_points)
plt.title("Aileron's twist along the span")
plt.xlabel('Span-wise location[m], from hinge 3 to 1', fontsize=10)
plt.ylabel('Aileron twist [deg]', fontsize=10)
plt.legend()
plt.show()
'''

def compute_q(position):
    Vy, Vz = get_forces(position, theta) #first entry refers to the position
    qb = get_qb(NSpS, Iyy, Izz, Vz, Vy, BrA)
    qs01 = get_Qs01(h_a, qb, Sa, t_sk, t_sp, StrinLoc)
    qs02 = get_Qs02(h_a, qb, Sa, t_sk, t_sp, StrinLoc)
    q_shear = get_q_shear(NSpS, Iyy, Izz, Vz, Vy, BrA, qs01, qs02)
    q_T=get_q_torque(position)
    q_circle = q_shear[:2]+q_T[2],
    q_spar = q_shear[9]+q_T[2]-q_T[3]
    q_triangle = q_shear[2:9] + q_T[3]
    q_tot = np.append(q_circle,np.append(q_triangle,q_spar))
    return q_tot

'''Get max shear in ribs
n_points = 1000
x_points = np.linspace(0, l_a, num=n_points)
y_points = np.zeros((n_points,3))
for i in range(np.size(x_points)):
    run = compute_q(x_points[i])
    idx = np.argmax(run)

    y_points[i] = np.append(x_points[i], np.append(idx, run[int(idx)]))
idx = np.argmax(y_points[:,2])
y_points = y_points
massimo = y_points[idx]
print(massimo)
x_points= np.array([x_R_1y, x_F_I, x_P, x_R_2y])
n_points = 1000
#x_points = np.linspace(0, l_a, num=n_points)
y_points = np.zeros((len(x_points),3))
for i in range(np.size(x_points)):
    run = compute_q(x_points[i])
    idx = np.argmax(run)
    y_points[i] = np.append(x_points[i], np.append(idx, run[idx]))
y_points = y_points
print(y_points)
'''
