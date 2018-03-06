import numpy as np
'''
It defines coefficients of the force for different uses, especially in the calculation of the internal forces
'''
class force(object):

    def __init__(self, position, sign = 1.):
        '''sign might be +1 or -1
        force position refers to the x position of the force in the reference frame specified in the report
        x line on hinge line, '''
        self.force_position = position
        self.force_sign = sign

    def coeff(self, current_position):
        '''output is coefficient of the force, givrn position considered'''
        return self.force_sign * np.heaviside(current_position-self.force_position, 0.)

    def moment(self, current_position):
        '''output is moment coefficient of the force, givrn position considered'''

        return self.force_sign*(current_position - self.force_position) \
               * np.heaviside(current_position-self.force_position, 0.)

    def slope(self, current_position):
        '''output is slope coefficient'''

        return self.force_sign*((current_position - self.force_position)**2.)/2.\
               * np.heaviside(current_position-self.force_position, 0.)

    def displacement(self, current_position):
        '''output is displacement coefficient'''

        return self.force_sign * ((current_position - self.force_position) ** 3.) / 6.\
               * np.heaviside(current_position-self.force_position, 0.)

class distributed_load(object):

    def __init__(self, position, sign = 1.):
        '''sign might be +1 or -1
        force position refers to the x position of the force in the reference frame specified in the report
        x line on hinge line, '''
        self.force_position = position
        self.force_sign = sign

    def coeff(self, current_position):
        '''output is coefficient of the force, givrn position considered'''

        return self.force_sign*(current_position - self.force_position)\
               * np.heaviside(current_position-self.force_position, 0.)

    def moment(self, current_position):
        '''output is moment coefficient of the force, givrn position considered'''

        return self.force_sign*((current_position - self.force_position)**2.)/2.\
               * np.heaviside(current_position-self.force_position, 0.)

    def slope(self, current_position):
        '''output is slope coefficient'''

        return self.force_sign * ((current_position - self.force_position) ** 3.) / 6.\
               * np.heaviside(current_position-self.force_position, 0.)

    def displacement(self, current_position):
        '''output is displacement coefficient'''

        return self.force_sign * ((current_position - self.force_position) ** 4.) / 24.\
               * np.heaviside(current_position-self.force_position, 0.)
