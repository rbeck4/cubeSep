#!/usr/env python
#Gonna need to pray to St. Isidore for this one

import numpy as np
import sys
import time

class Cube:
  '''Class to read in cube information.  It can split up the real/imag and
     alpha/beta components.  It can plot less common properties, such as GHF
     magnetization or electrostatic potential gradients, as well as write the
     resulting cube files out to be opened in programs such as GaussView or
     PyMol.
  '''
  def __init__(self, filename)
    self.comment = []
    self.sign   = None #Num. Atms sign det. if total den or molec. orb.
    self.natoms = None
    self.ncubes = None
    self.nMOs   = None
    self.origin = None
    self.nx     = None
    self.ny     = None
    self.nz     = None
    self.x      = None
    self.y      = None
    self.z      = None
    self.nC     = None
    self.atoms  = []

    with open(filename, 'r') as f:
      tick = time.time()
      #Save Comments:
      for i in range(2):
        self.comment.append(f.readline().strip('\n'))
      
      #Number atoms and cubefile origin
      line = f.readline().split()
      self.natoms = abs(int(line[0]))
      self.sign   = np.sign(int(line[0]))
      self.origin = np.array([float(line[x]) for x in range(1,4)])
      #Number of cubes in file:
      self.ncubes = int(line[4])
      
      #Number voxels and axis vect
      for vector in [['nx','x'],['ny','y'],['nz','z']]:
        line = f.readline().split()
        self.__dict__[vector[0]] = int(line[0])
        self.__dict__[vector[1]] = np.array([float(line[x]) for x in range(1,4)])
  
      #Atoms and Coords
      for atom in range(self.natoms):
        line = f.readline().split()
        self.atoms.append(line)
  
      #How many cubes are in file?
      line = f.readline().split()
      if self.sign < 0:
        #Mult MOs
        self.nMOs = int(line[0])
      
      #Get the volumetric data (to the end of the file)
      #Probably not fastest way of doing things, but w/o foreknowledge abt the
      #type of job we can't pre-allocate the array very well...
      vals = np.asarray([float(v) for s in f for v in s.split()])

      #Determine 
      self.nC = int(len(vals) / (self.nx * self.ny * self.nz * self.ncubes))

      #Cube File Valid?
      if self.nC != 1 or self.nC != 2 or self.nC !=4:
        raise NameError("Num. Comp. of %i in cubfile not (yet) valid" %self.nC)

      tock = time.time()
      print("Initialized in %.2f s" %(tock-tick))

  def get_numCubes(self):
    '''
    Getter fxn for number of cubes
    '''
    return self.numCubes

  def do_vspin(self, cubeSelect=1):
    '''
    For Mag. information
    '''
    
    #Only def. in o.g. for complex. two-comp case
    if self.nC != 4:
      raise NameError("Can't do for num Comp. = %i." %self.nC)
    
    #Select cube of interest, may not be most mem efficient, but can read so...
    data = vals.reshape(-1,self.ncubes*self.nC)\
                        [:,(cubeSelect-1)*self.nC:cubeSelect*self.nC]\
                        .reshape(-1)
    #Return N, Mx, My, Mz:
    return [data[0::self.nC], data[1::self.nC], \
            data[2::self.nC], data[3::self.nC]]

if __name__ == '__main__':
  filename = sys.argv[-1]

  atom = Cube(filename)
  
