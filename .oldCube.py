from __future__ import print_function
import numpy as np
import sys

class Cube:
    ''' Class to read in a cube file, and split up real/imag and alpha/beta
        components. It can plot some less common properties, like GHF 
        magnetization or electrostatic potential gradients, as well as write the
        cubes out to file for use in GaussView or other cube-viewing programs, 
        such as PyMol
    '''
    def __init__(self, filename,vspin=False,cubeSelect=1):
        self.comment = []
        self.sign    = None # sign on number atoms determines if we are dealing
                            # with total density or a molecular orbital.
        self.natoms  = None
        self.ncubes  = None
        self.origin  = None 
        self.nx      = None 
        self.ny      = None 
        self.nz      = None 
        self.x       = None 
        self.y       = None 
        self.z       = None 
        self.atoms   = []
        # GHF only; do cubegen with vspin to extract magnetization
        self.vspin   = vspin
        with open(filename,'r') as f:
           # Save comments (first two lines)
           for i in xrange(2):
               self.comment.append(f.readline().strip('\n')) 
           # Grab number of atoms and the cubefile origin
           line = f.readline().split()
           self.natoms = abs(int(line[0]))
           self.sign   = np.sign(int(line[0]))
           self.origin = np.array([float(line[1]),
                                   float(line[2]),
                                   float(line[3])])
           # Grab number voxels and axis vector
           for vector in [['nx','x'],['ny','y'],['nz','z']]:
               line = f.readline().split()
               self.__dict__[vector[0]] = int(line[0])
               self.__dict__[vector[1]] = np.array([float(line[1]),
                                                    float(line[2]),
                                                    float(line[3])])
           # Grab atoms and coordinates 
           for atom in xrange(self.natoms):
               line = f.readline().split()
               self.atoms.append([line[0],line[1],line[2],line[3],line[4]])
           # How many cubes are in file?
           if self.sign < 0:
               line = f.readline().split()
               self.ncubes = int(line[0])
           if cubeSelect > self.ncubes:
               print("Warning, have selected a cube beyond end of file, setting\
                      to last")
               cubeSelect = self.ncubes
           if cubeSelect < 1:
               print("Warning, cubeSelect index starts at one; adjusting")
               cubeSelect = 1
           print("NCUBES: ", self.ncubes, " SELECTED: ", cubeSelect)

           # Grab volumetric data, which finishes out our file
           vals = np.asarray([float(v) for s in f for v in s.split()])
           # now we need to see what kind of data we have
           # (real/complex?,two-component)

           if(self.vspin):
               if len(vals) == 4*self.nx*self.ny*self.nz * self.ncubes:
                   print("Case1")
                   self.N = np.zeros((self.nx,self.ny,self.nz))
                   self.Mx = np.zeros((self.nx,self.ny,self.nz))
                   self.My = np.zeros((self.nx,self.ny,self.nz))
                   self.Mz = np.zeros((self.nx,self.ny,self.nz))

                   for idx,v in enumerate(vals.reshape(-1,self.ncubes*4) \
                                      [:,(cubeSelect-1)*4:cubeSelect*4].\
                                      reshape(-1)):
                       i = int(np.floor(idx/4))
                       if idx % 4 == 0:
                           self.N[i/(self.ny*self.nz),
                                     (i/self.nz)%self.ny,
                                      i%self.nz] = float(v)
                       elif idx % 4 == 1:
                           self.Mx[i/(self.ny*self.nz),
                                     (i/self.nz)%self.ny,
                                      i%self.nz] = float(v)
                       elif idx % 4 == 2:
                           self.My[i/(self.ny*self.nz),
                                     (i/self.nz)%self.ny,
                                      i%self.nz] = float(v)
                       elif idx % 4 == 3:
                           self.Mz[i/(self.ny*self.nz),
                                     (i/self.nz)%self.ny,
                                      i%self.nz] = float(v)
               
           elif len(vals) == self.nx*self.ny*self.nz * self.ncubes:
               print("Case2")
               self.volRA = np.zeros((self.nx,self.ny,self.nz))
               
               for idx,v in enumerate(vals.reshape(-1,self.ncubes) \
                                      [:,(cubeSelect-1):cubeSelect].\
                                      reshape(-1)):
                   self.volRA[idx/(self.ny*self.nz),
                           (idx/self.nz)%self.ny,
                            idx%self.nz] = float(v)
               self.volNorm = np.sqrt(self.volRA**2)
           elif len(vals) == 2*self.nx*self.ny*self.nz * self.ncubes:
               print("Case3")
               # This is the complex one component case
               self.volRA = np.zeros((self.nx,self.ny,self.nz))
               self.volIA = np.zeros((self.nx,self.ny,self.nz))
               
               for idx,v in enumerate(vals.reshape(-1,self.ncubes*2) \
                                      [:,(cubeSelect-1)*2:cubeSelect*2].\
                                      reshape(-1)):
                   i = int(np.floor(idx/2))
                   if idx % 2 == 0:
                       self.volRA[i/(self.ny*self.nz),
                                 (i/self.nz)%self.ny,
                                  i%self.nz] = float(v)
                   elif idx % 2 == 1:
                       self.volIA[i/(self.ny*self.nz),
                                 (i/self.nz)%self.ny,
                                  i%self.nz] = float(v)
               self.volNorm = np.sqrt(self.volRA**2 +
                                      self.volIA**2 )
           elif len(vals) == 4*self.nx*self.ny*self.nz * self.ncubes:
               print("Case4")
               # This is the complex 2-component case, so we have to split the
               # volumetric data into four parts. The density is stored real
               # alpha, imaginary alpha, real beta, imaginary beta 
               self.volRA = np.zeros((self.nx,self.ny,self.nz))
               self.volIA = np.zeros((self.nx,self.ny,self.nz))
               self.volRB = np.zeros((self.nx,self.ny,self.nz))
               self.volIB = np.zeros((self.nx,self.ny,self.nz))

               #In order to get 2nd cube need to reshape to ncube*4 
               #(4 mo1 4 mo2) and select appropriate range:

               for idx,v in enumerate(vals.reshape(-1,self.ncubes*4) \
                                      [:,(cubeSelect-1)*4:cubeSelect*4].\
                                      reshape(-1)):
                   i = int(np.floor(idx/4))
                   if idx % 4 == 0:
                       self.volRA[i/(self.ny*self.nz),
                                 (i/self.nz)%self.ny,
                                  i%self.nz] = float(v)
                   elif idx % 4 == 1:
                       self.volIA[i/(self.ny*self.nz),
                                 (i/self.nz)%self.ny,
                                  i%self.nz] = float(v)
                   elif idx % 4 == 2:
                       self.volRB[i/(self.ny*self.nz),
                                 (i/self.nz)%self.ny,
                                  i%self.nz] = float(v)
                   elif idx % 4 == 3:
                       self.volIB[i/(self.ny*self.nz),
                                 (i/self.nz)%self.ny,
                                  i%self.nz] = float(v)
               self.volNorm = np.sqrt(self.volRA**2 +
                                      self.volRB**2 +
                                      self.volIA**2 +
                                      self.volIB**2) 
       
           else:
               print("NX,Y,Z,vals = ", self.nx, " ", self.ny, " ", self.nz, " ", len(vals))
               raise NameError, "cube file not valid"

    def plot_property(self,prop=None):
        ''' This function will plot vector fields for you to enjoy. Right now,
            if you feed it an electrostatic potential (ESP), it will compute and
            plot the gradient. If you have a noncollinear solution and run VSPIN
            with cubegen, it will plot the magnetization (MAG)
            You need Mayavi to plot in 3D.
        '''
        try:
            import mayavi.mlab as mlab 
        except ImportError:
            print("[ImportError] You need Mayavi \n \
                   (http://docs.enthought.com/mayavi/mayavi/)") 
        if not prop:
            sys.exit('No property chosen!')
        if prop == 'ESP':
            volume = np.gradient(self.volRA,edge_order=2)
            volumeX = volume[0]
            volumeY = volume[1]
            volumeZ = volume[2]
            import mayavi.mlab as mlab
            vol = mlab.flow(volumeX,volumeY,volumeZ,\
                seed_resolution=25,\
                seed_visible=True,\
                seedtype='plane',\
                seed_scale=2,\
                colormap='bone',\
                integration_direction='both',\
                line_width=1.0)
            #mlab.quiver3d(volumeX,volumeY,volumeZ)
            mlab.show()
        elif (self.vspin) and (prop == 'MAG'):
            import mayavi.mlab as mlab
            import matplotlib.pyplot as plt 
            plt.quiver(self.Mx[:,self.ny/2,:],self.Mz[:,self.ny/2,:],scale=10.0)
            plt.show()
            vol = mlab.flow(self.Mx,self.My,self.Mz,\
                seed_resolution=25,\
                seed_visible=True,\
                seedtype='plane',\
                seed_scale=2,\
                colormap='bone',\
                integration_direction='both',\
                line_width=1.0)
            #mlab.quiver3d(self.Mx,self.My,self.Mz,mask_points=50)
            #mlab.show()

    def write_out(self,filename,data='RA'):
        '''Generate new cube file
           String variable 'data' specifies which component you want
           RA = real alpha
           IA = imag alpha
           RB = real beta 
           IB = imag beta
        '''
        if data == 'RA':
            volume = self.volRA
        elif data == 'IA':
            volume = self.volIA
        elif data == 'RB':
            volume = self.volRB
        elif data == 'IB':
            volume = self.volIB
        elif data == 'NORM':
            volume = self.volNorm 
        # Below are some options to print the Magnitude (Mag) and Argument (Arg)
        # of complex cubes, e.g. if z = a + i*b -> |z|*exp(i*theta) where |z| is 
        # the Magnitude (Mag) and theta is the Argument (Arg)
        elif data == 'MagA':
            volume = np.sqrt(self.volRA**2 + self.volIA**2)
        elif data == 'ArgA':
            volume = (np.arctan2(self.volRA,self.volIA))
        elif data == 'MagB':
            volume = np.sqrt(self.volRB**2 + self.volIB**2)
        elif data == 'ArgB':
            volume = (np.arctan2(self.volRB,self.volIB))
        # this the magnitude and angle between alpha and beta components for GHF
        elif data == 'MagABr':
            volume = np.sqrt(self.volRA**2 + self.volRB**2)
        elif data == 'ArgABr':
            volume = (np.arctan2(self.volRA,self.volRB))
        # if GHF and VSPIN, you can write out Magnetization densities a la eq 18 
        # of Bulik, et al, PHYSICAL REVIEW B 87, 035117 (2013)
        # self.N  = P_{aa} + P_{bb}
        # self.Mx = 2*Re(P_{ab})
        # self.My = 2*Im(P_{ab})
        # self.Mz = P_{aa} - P_{bb}
         
        if self.vspin:
            if data == 'N':
                volume = self.N
            elif data == 'Mx':
                volume = self.Mx
            elif data == 'My':
                volume = self.My
            elif data == 'Mz':
                volume = self.Mz

        with open(filename, 'w') as f:
            for i in self.comment:
                print(str(i),file=f)
            print(" %4d %.6f %.6f %.6f" % (self.sign*self.natoms, 
                self.origin[0], self.origin[1],self.origin[2]),file=f)
            print(" %4d %.6f %.6f %.6f" % (self.nx, self.x[0], self.x[1],
                self.x[2]), file=f)
            print(" %4d %.6f %.6f %.6f" % (self.ny, self.y[0], self.y[1],
                self.y[2]), file=f)
            print(" %4d %.6f %.6f %.6f" % (self.nz, self.z[0], self.z[1],
                self.z[2]), file=f)
            for atom in self.atoms:
                print(" %s %s %s %s %s" % (atom[0], atom[1], atom[2], atom[3],
                    atom[4]), file=f)
            if self.sign < 0:
                print("    1    "+str(self.natoms),file=f)
            for ix in xrange(self.nx):
                for iy in xrange(self.ny):
                    for iz in xrange(self.nz):
                        print(" %.5e " % volume[ix,iy,iz],end="",file=f)
                        if (iz % 6 == 5):
                            print('',file=f)
                    print('',file=f)
       

if __name__ == '__main__':
    # Basic usage example: dump complex GHF MO data to separate real and # imaginary, alpha and beta, cubes for visualizing each separately.
    filename=sys.argv[-1]
    #for filename in listcubes:
    for selection in range(1,3):
      #atom = Cube(filename,cubeSelect=selection)
      #atom.write_out(filename+"twoc_ra_%i.cube" %selection,data='RA')
      #atom.write_out(filename+'twoc_ia.cube',data='IA')
      #atom.write_out(filename+"twoc_rb_%i.cube" %selection,data='RB')
      #atom.write_out(filename+'twoc_ib.cube',data='IB')
      #atom.write_out(filename+"twoc_maga_%i.cube" %selection,data='MagA')
      #atom.write_out(filename+"twoc_magb_%i.cube" %selection,data='MagB')
      #atom.write_out(filename+"twoc_arga_%i.cube" %selection,data='ArgA')
      #atom.write_out(filename+"twoc_argb_%i.cube" %selection,data='ArgB')
      #atom.write_out(filename+"twoc_argabr_%i.cube" %selection,data='ArgABr')
      #atom.write_out(filename+"twoc_magabr_%i.cube" %selection,data='MagABr')



    
    ## GHF cubegen w/ VSPIN example to plot and dump magnetization densities
    magnetization = Cube('ghf_Sing.cube',vspin=True)
    #magnetization.plot_property('MAG')
    magnetization.write_out('Mz.cube',data='N')
    # electrostatic potential example
    #ESP = Cube('esp-example.cube')
    #ESP.plot_property('ESP')
    

