'''
INPUT VARIABLES given for the assignment Boeing 737
'''


#geometry properties aileron in general
c_a = 0.605                 #chord length aileron [m]
l_a = 2.661                 #span of the aileron [m]
x_1 = 0.172                 #x-location of hinge 1 [m]
x_2 = 1.211                 #x-location of hinge 2 [m]
x_3 = 2.591                 #x-location of hinge 3 [m]
x_a = 35.0 * 10**(-2)        #distance between actuator 1 and 2 [m]
h_a = 20.5 * 10**(-2)        #aileron height [m]
t_sk = 1.1 * 10**(-3)        #skin thickness [m]
E = 73.1 * 10**(9)            #E-modulus Al2024-T3 [Pa]
G = 28.e9

#geometry properties spar
t_sp = 2.8 * 10**(-3)        #spar thickness [m]

#geometry properties stiffener
t_st = 1.2 * 10**(-3)        #stiffener thickness [m]
h_st = 1.6 * 10**(-2)        #stiffener height [m]
w_st = 1.9 * 10**(-2)        #stiffener width [m]
n_st = 15                   #number of stiffeners (equally spaced) [-]

#displacement variables
d_1 = 11.54e-2 /2.54       #vertical displacement hinge 1 [m]
d_3 = 18.40e-2 /2.54      #vertical displacement hinge 3 [m]
theta = 28.                  #maximum upward deflection [deg]

#load variables
P = 97.4 * 10**(3)           #load in actuator 2 [N]
q = 5.54 * 10**(3)           #net aerodynamic load [N/m]


