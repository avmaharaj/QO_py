import numpy as np


class T2g_SpinOrbit_Model(object):
    """
    A T2g_SpinOrbit_Model produces the diagonal and off diagonal blocks of a Hamiltonian
    corresponding to electrons hopping on a 3D lattice, in the dxz, dyz and dxy orbitals
    with spin orbit coupling includded.
    """


    def __init__(self, spin_orbit, theta, phi, mu, delta, t2):
        """
        Construct a new T2g_SpinOrbit_Model instance.

        Required arguments:
        - spin_orbit: The strength of spin orbit coupling
        - theta: the polar angle of the magnetic field
        - phi: the azimuthal angle of the magnetic field
        - mu: the chemical potential
        - delta: the value of the broadening
        - t2: the second neighbor hopping
        """


        #Initialize the external parameters
        self.spin_orbit = spin_orbit
        self.theta = theta
        self.phi = phi
        self.mu = mu
        self.broadening = delta
        self.t2 = t2

        #Initialize the bfield to 1 - the solver will update this
        #passed for now

        #Initialize the ky and kz values (not going to be used for now)
        self.kz = 0.0
        self.ky = 0.0

        #The size of a block
        self.dim = 6



    def get_dim(self):
    
        return self.dim
    

    def diagonal_block(self, xpos, bfield):
        """
        returns a diagonal block of the Hamiltonian, given the current position
        
        Inputs - the x position (integer)
        
        Outputs - a numpy array of size (self.dim x self.dim) which is a diagonal
                    block of the Hamiltonian
        """


        #TODO: figure out logic and exact nature of this this term
        block = np.zeros((self.dim, self.dim),dtype = np.complex128)

        block[0,0] = 2.0*np.cos(2.0*np.pi*bfield*xpos*np.cos(self.theta) - self.ky) \
                    + 2.0*np.cos(2.0*np.pi*bfield*xpos*np.sin(self.theta) - self.kz) \
                    + self.mu + self.broadening*1j
        block[3,3] = block[0,0]

        block[1,1] = 2.0*self.t2*np.cos(2.0*np.pi*bfield*xpos*np.cos(self.theta) - self.ky) \
                    +2.0*np.cos(2.0*np.pi*bfield*xpos*np.sin(self.theta) - self.kz) \
                    + self.mu + self.broadening*1j
        block[4,4] = block[1,1]

        block[2,2] = 2.0*np.cos(2.0*np.pi*bfield*xpos*np.cos(self.theta) - self.ky) \
                    +2.0*self.t2*np.cos(2.0*np.pi*bfield*xpos*np.sin(self.theta) - self.kz) \
                    + self.mu + self.broadening*1j
        block[5,5] = block[2,2]

        block[0,1] = self.spin_orbit*1j
        block[1,0] = -1.0*block[0,1]

        block[0,5] = -1.0*self.spin_orbit
        block[5,0] = block[0,5]

        block[1,5] = self.spin_orbit*1j
        block[5,1] = -1.0*block[1,5]

        block[2,3] = self.spin_orbit
        block[3,2] = block[2,3]

        block[2,4] = -1.0*self.spin_orbit*1j
        block[4,2] = -1.0*block[2,4]

        block[3,4] = -1.0*self.spin_orbit*1j
        block[4,3] = -1.0*block[3,4]



        return block


    def right_offdiag_block(self, xpos, bfield):
        """
        returns upper off diagonal block of the Hamiltonian, given the current position
        
        Inputs - the x position (integer)
        
        Outputs - a numpy array of size (self.dim X self.dim)
        """
        #TODO: figure out what this term looks like
        block = np.zeros((self.dim, self.dim),dtype = np.complex128)
        block[0,0] = self.t2
        block[3,3] = self.t2

        block[1,1] = 1.0
        block[4,4] = 1.0

        block[2,2] = 1.0
        block[5,5] = 1.0


    

        return block


    def left_offdiag_block(self, xpos, bfield):
        """
        returns upper off diagonal block of the Hamiltonian, given the current position
        
        Inputs - the x position (integer)
        
        Outputs - a numpy array of size (self.dim X self.dim)
        """
        #TODO: figure out what this term looks like
        block = np.zeros((self.dim, self.dim),dtype = np.complex128)
        block[0,0] = 1.0 + 0.0*1j
        block[1,1] = 1.0 + 0.0*1j
        

        block[0,0] = self.t2
        block[3,3] = self.t2

        block[1,1] = 1.0
        block[4,4] = 1.0

        block[2,2] = 1.0
        block[5,5] = 1.0
        
        return block







