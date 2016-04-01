import numpy as np
from scipy import linalg as slinalg

class RecursiveSolver(object):
    """
        A Solver encapsulates all the logic necessary for calculating
        the density of states recursively. 
        
        The solver accepts the range of magnetic fields, and optionally, number of b-field points
        
        The Solver accepts a model object (which consists of the diagonal and
        off diagonal blocks of the hamiltonian), and uses the recursive d.o.s.
        algorithm to calcualte the density of states.
        
        Calling run() method returns a numpy array with the density of states vs.
        field strength
        
        
        Example usage as follows:
        
        model = MyAwesomeModel(spin_orbit = , density = , theta = , phi = )
        solver = RecursiveSolver(model,system_size = , min_field = , max_field =  , mu = )
        solver.run()
        
        ##TODO: Parallelize the run loop implementation.
        

    """


    def __init__(self, model, system_size = 8192, min_invfield = 0.001, max_invfield = 1000,  **kwargs):
        """
        Construct a new RecursiveSolver instance.

        Required arguments:
        - model: A model object conforming to the API described above
        - system_size : the linear dimension of the sytem
        - min_field : the minimum field strength
        - max_field : the maximum field_strength
        - mu: the energy/chemical potential to be used
        OPTIONALLY:
        - num_steps = the number of bpts to include
        - print_every = the number of iterations at which to print
        - verbose = whether or not to print updates of progress
        """
        #Set up the solver's attributes
        self.model = model
        
        self.system_size = system_size
        self.min_invfield = min_invfield
        self.max_invfield = max_invfield
        
        
        # Unpack keyword arguments
        self.num_steps = kwargs.pop('num_steps', 1000)
        self.print_every = kwargs.pop('print_every', 100)
        self.verbose = kwargs.pop('verbose', True)
        
        self.dim = model.get_dim()

        # Throw an error if there are extra keyword arguments
        if len(kwargs) > 0:
          extra = ', '.join('"%s"' % k for k in kwargs.keys())
          raise ValueError('Unrecognized arguments %s' % extra)

        self.bfield = 0.0
        self.dLeft = np.zeros((self.system_size,self.dim,self.dim),dtype = np.complex128)
        self.cLeft = np.zeros((self.system_size,self.dim,self.dim), dtype = np.complex128)
        self.dRight = np.zeros((self.system_size,self.dim,self.dim), dtype = np.complex128)
        self.cRight = np.zeros((self.system_size,self.dim,self.dim), dtype = np.complex128)
          







    def run(self):
        """
        Runs the solver to evaluate d.o.s. as a function of field
        
        """
        
        #first set up the field values (could move this to init)
        lower_inv_field = self.min_invfield
        upper_inv_field = self.max_invfield
        dfield = (upper_inv_field - lower_inv_field) / self.num_steps
        
        output = np.zeros(( self.num_steps, 2),dtype = np.complex128)


        #now loop over values of the magnetic field
        #TODO: Parallelize this loop in some way
        for b in xrange(self.num_steps):

            dos = 0.0 + 0.0*1j
            self.bfield = 1.0*np.pi/(lower_inv_field + b * dfield)
            dos = self._get_dos()

            output[b,0], output[b,1] = 1.0/self.bfield, dos
        
        
            if self.verbose and (b+1)%self.print_every == 0 :
                print('Completed iteration %i of %i' %(b+1, self.num_steps))


        return output



    def _get_dos(self):
        """
        Performs the upwards and downwards sweep to get the density of states
        for a single value of the bfield
            
        """
        self._downsweep()
        self._upsweep()
        dos = self._calc_dos()

        return dos



    def _downsweep(self):
        """
        Sweeps down the hamiltonian, filling in self.dLeft, self.cLeft
        """
        
        model = self.model
        for i in xrange(self.system_size-1):
            if i == 0:
                self.dLeft[i] = model.diagonal_block( i, self.bfield )
            else:
                self.dLeft[i] = model.diagonal_block( i, self.bfield ) \
                        + np.dot(self.cLeft[i-1], model.right_offdiag_block(i-1,self.bfield))


            dinv = np.linalg.inv(self.dLeft[i])
            self.cLeft[i] = -1.0 * np.dot(model.left_offdiag_block(i+1, self.bfield),dinv)
    
        i = self.system_size -1
        self.dLeft[i] = model.diagonal_block( i, self.bfield ) \
            + np.dot(self.cLeft[i-1], model.right_offdiag_block(i-1,self.bfield))



    def _upsweep(self):
        """
        Sweeps up the hamiltonian, filling in self.dright
        """
        model = self.model
        for i in reversed(xrange(1,self.system_size)):
            if i == self.system_size - 1:
                self.dRight[i] = model.diagonal_block( i, self.bfield )
            else:
                self.dRight[i] = model.diagonal_block( i, self.bfield ) + \
                        np.dot(self.cRight[i+1], model.left_offdiag_block(i+1,self.bfield))
    
            dinv =  np.linalg.inv(self.dRight[i])
            self.cRight[i] = -1.0 * np.dot(model.right_offdiag_block(i-1, self.bfield),dinv)

        i = 0
        self.dRight[i] = model.diagonal_block( i, self.bfield ) + \
            np.dot(self.cRight[i+1], model.left_offdiag_block(i+1,self.bfield))




    def _calc_dos(self):
        """
        Finds trace of the Green's function
        """
        model = self.model
        dos = 0.0 + 0.0*1j
        
        
        for i in xrange(self.system_size):
            ginv = -1.0*model.diagonal_block( i, self.bfield ) + self.dLeft[i] + self.dRight[i]
            
            g = slinalg.inv(ginv)
        
            dos -= (np.trace(g)).imag
        
        dos /= np.pi * self.system_size

        return dos
