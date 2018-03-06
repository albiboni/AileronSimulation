import numpy as np
from Boom import *

'''
It calculates the inertia from teh booms and it rotates it between different coordinate frames, body frame and general
'''

#summary booms, position is given from leading edge
Booms = np.zeros((15,3))
Booms[0] = np.array([Sty[0], Sty[0], BrA[0]])
for i in xrange(1, len(StZ)-1):
    Booms[i] = np.array([StZ[i], Sty[i], BrA[i]])
    Booms[15-i] = np.array([StZ[i], -Sty[i], BrA[i]])
Booms= np.insert(Booms, 2, [StZ[8],Sty[8],BrA[8]], axis=0)
Booms= np.insert(Booms, 15, [StZ[8],-Sty[8],BrA[8]], axis=0)


def centroid(Booms): #centroid
    SumBY = 0
    SumY = 0
    for c in range(np.shape(Booms)[0]):
        SumBY += Booms[c][2]*Booms[c][0]
        SumY += Booms[c][2]
    Cz = SumBY/SumY
    return Cz


def get_inertia(Booms): #moment of inertia around centroid around y, Iyy, and around x, Ixx
    Cz = centroid(Booms)
    Iyy = 0
    Izz = 0
    for c in range(np.shape(Booms)[0]):
        Iyy += Booms[c][2]*(Booms[c][0] - Cz)**2
        Izz += Booms[c][2]*(Booms[c][1])**2
    return Iyy, Izz


def rotation_matrix (theta):
    'clockwise rotation theta in degrees'
    #http://mathworld.wolfram.com/RotationMatrix.html
    radians = theta*np.pi/180.
    rotation_matrix = np.array([[np.cos(radians), np.sin(radians)],
                                [-np.sin(radians), np.cos(radians)]])
    return rotation_matrix


def rotate_inertia (theta, Iyy, Izz, Izy=0., Iyz=0.):

    '''produces clockwise rotation of theta, considers inertia from centroid
    and neutral axis (Izz and Iyy). Izy is 0 because in this reference frame there is a symmetry plane
    output:we are interested only interested in 3 entries of this matrix '''
    #http: // www.kwon3d.com / theory / moi / iten.html
    #http://www.kwon3d.com/theory/moi/triten.html
    inertia_rotated = np.dot(np.dot(rotation_matrix(theta).T, np.array([[Iyy, Iyz], [Izy, Izz]])),
                             rotation_matrix(theta))
    return inertia_rotated[0,0], inertia_rotated[1,1], inertia_rotated[0,1]

#print(rotate_inertia(28, 3.7271242564*10**(-5),6.30441732535*10**(-6)))

#print(rotate_inertia(28, 8.611195*10**(-5),9.964356*10**(-6)))  Ale

#Iyy, Izz, Izy = rotate_inertia(theta, Iyy, Izz)
#print(get_inertia(Booms))
#http://calcresource.com/moment-of-inertia-rotation.html #for testing
