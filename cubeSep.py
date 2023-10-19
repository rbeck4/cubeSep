#!/usr/env python
#Gonna need to pray to St. Isidore for this one

import numpy as np
import os.path
import psutil
import struct
import sys
import time

class Cube:
  '''Class to read in cube information.  It can split up the real/imag and
     alpha/beta components.  It can plot less common properties, such as GHF
     magnetization or electrostatic potential gradients, as well as write the
     resulting cube files out to be opened in programs such as GaussView or
     PyMol.
  '''
  def __init__(self, filename):
    '''
    Initialize, read in and set up cube file.  Need to call specific method with
    this version of the script (between real/complex1c/complex2c).
    '''
    self.comment = []
    self.sign   = None #Num. Atms sign det. if total den or molec. orb.
    self.natoms = None
    self.ncubes = None
    self.nMOs   = None
    self.MOnames= []
    self.origin = None
    self.nx     = None
    self.ny     = None
    self.nz     = None
    self.x      = None
    self.y      = None
    self.z      = None
    self.nC     = None
    self.atoms  = []
    self.vals   = None
    self.tmpfl  = None

    #Ensure supplied file exists and isn't script:
    if filename == 'cubeSep.py' or not os.path.exists(filename):
      raise NameError("File %s either does not exist or is self." %filename)

    #In case of RAM limitations:
    if 2 * os.path.getsize(filename) >= (0.6) * psutil.virtual_memory()[1]:
      print("Large cube, using memmap")
      self.tmpfl = 'cubeSep_tempFile.dat'

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
  
      #MO labels (if MO), can stretch multiple lines
      fullMO = 0
      while fullMO < 1: 
        line = f.readline().split()
        #Get number of elements (should be eq. to ncubes but...)
        if self.sign < 0 and not self.nMOs:
          #Mult MOs
          self.nMOs = int(line.pop(0))
        for x in line:
          self.MOnames.append(x)
        if len(self.MOnames) >= self.nMOs:
          fullMO = 10
      
      #Get the volumetric data (to the end of the file)
      #We will need to have either into ram or into file depending on file
      #size, and computer memory... if you have a better way of doing this
      #let me know!
      if not self.tmpfl:
        self.vals = np.asarray([float(v) for s in f for v in s.split()])
      else:
        tmpLength = 0
        with open(self.tmpfl, 'wb') as binFl:
          for line in f:
            load = line.split()
            tmpLength += len(load)
            for i in range(len(load)):
              binFl.write(struct.pack("f", float(load[i])))
        self.vals = np.memmap(self.tmpfl, dtype=np.float32, mode='r')

      #Determine 
      self.nC = int(len(self.vals) / \
                    (self.nx * self.ny * self.nz * self.ncubes))

      #Cube File Valid?
      if self.nC not in [1, 2, 4]:
        raise NameError("Num. Comp. of %i in cubfile not (yet) valid" %self.nC)

      tock = time.time()
      print("Initialized in %.2f s." %(tock-tick))

  
  def __enter__(self):
    '''
    Needed to use with?
    '''
    return self


  def __exit__(self, exc_type, exc_value, traceback):
    '''
    If used memmap remove the file so doesn't clog up file system
    '''
    if self.tmpfl:
      os.unlink(self.tmpfl)


  def get_numCubes(self):
    '''
    Getter fxn for number of cubes
    '''
    return self.ncubes

  def cube_select(self, cubeSelect=1):
    '''
    Return the data array for the selected cube of interest)
    May not be the most memory efficient, but is readable so....
    '''
    #Ensure cubeSelect is within range:
    if cubeSelect < 1:
      print("Warning, the selected cube (%i) is out of range, setting to 1" \
            %cubeSelect)
      cubeSelect = 1
    
    if cubeSelect > self.ncubes:
      print("Warning, the selected cube (%i) is out of range, setting to max" \
            %cubeSelect)
      cubeSelect = self.ncubes
    
    return self.vals.reshape(-1,self.ncubes*self.nC)\
                        [:,(cubeSelect-1)*self.nC:cubeSelect*self.nC]\
                        .reshape(-1)
  
  def get_vspin(self, cubeSelect=1):
    '''
    For Mag. information
    '''
    #Only def. in o.g. for complex. two-comp case
    if self.nC != 4:
      raise NameError("Can't do for num Comp. = %i." %self.nC)
    
    data = self.cube_select(cubeSelect)
    #Return N, Mx, My, Mz:
    return [data[0::self.nC], data[1::self.nC], \
            data[2::self.nC], data[3::self.nC]]

  
  def get_volRA(self, cubeSelect=1):
    '''
    General volRA (and normRA) return)
    '''
    data = self.cube_select(cubeSelect)
    #RA is always first, separation would be nC, but know that too:
    return data[0::self.nC]


  def get_volRB(self, cubeSelect=1):
    '''
    General volRB (and normRB) return)
    '''
    if self.nC < 4:
      raise NameError("B orbitals not included for numComponent = %i" %self.nC)

    data = self.cube_select(cubeSelect)
    #RA is always first, separation would be nC, but know that too:
    return data[2::self.nC]


  def get_volIA(self, cubeSelect=1):
    '''
    General volIA (and normIA) return)
    '''
    if self.nC < 2:
      raise NameError("IA orbitals not included for numComponent = %i" %self.nC)

    data = self.cube_select(cubeSelect)
    #RA is always first, separation would be nC, but know that too:
    return data[1::self.nC]


  def get_volIB(self, cubeSelect=1):
    '''
    General volIB (and normIB) return)
    '''
    if self.nC < 4:
      raise NameError("IB orbitals not included for numComponent = %i" %self.nC)

    data = self.cube_select(cubeSelect)
    #RA is always first, separation would be nC, but know that too:
    return data[3::self.nC]


  def write_to_file(self, cubeData, fileNameOut="cub.cube", cubeSelect=1):
    '''
    Write single cube file to fileNameOut, need selected cube if you want the
    title to carry over from the combined cube file
    '''
    tick = time.time()
    with open(fileNameOut, 'w') as f:
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
          print("    1    "+self.MOnames[cubeSelect-1],file=f)
      lineidx = 0
      for i in range(len(cubeData)):
        lineidx += 1
        print( "%.5e " % cubeData[i], file=f, end='')
        if (lineidx == 6):
          lineidx = 0
          print('', file=f)
        if (i % self.nz == self.nz-1):
          lineidx = 0
          print('', file=f)

    tock = time.time()
    print("Finished writing %s in %.2f s." %(fileNameOut, tock-tick))


if __name__ == '__main__':
  '''
  Please view https://doi.org/10.1002/jcc.26196.
  There are getter functions for the select volumes real/imag. A/B. Printing is
  general, and just needs a data array passed in, it will format as a Gaussian
  cube file.  This script is able to handle cube files that contain multiple
  cube files (e.g. cubegen 1 MO=1-3 test.fchk test.cube).  The example uses a
  'with' statement so that in the case of large cube files, the generated 
  temporary files can be cleaned from the file system.
  
  Potential cubes of interest:
    RA : Real alpha part of cube (get_volRA()[0])
    IA : Imaginary alpha part of cube (get_volIA()[0]) (num. component >= 2)
    RB : Real beta part of cube (get_volRB()[0])       (num. component >= 4)
    IB : Imaginary beta part of cube (get_volIB()[0])  (num. component >= 4)
    ArgA : Argument alpha (np.arctan2(get_volRA()[0], get_volIA()[0])
                                                       (num. component >= 2)
    ArgB : Argument beta (np.arctan2(get_volRB()[0], get_volIB()[0])) 
                                                       (num. component >= 4)
    MagA : Magnitude alpha sqrt((get_volRA()**2 + get_volIA()**2))                                                      
                                                       (num. component >= 2)
    MagB : Magnitude beta (sqrt(get_volRB()**2 + get_volIB()**2))
                                                       (num. component >= 4)
    MagAB : Magnitude between alpha and beta in GHF (sqrt(RA**2 + IA**2 + 
             RB**2 + IB**2))                           (num. component >= 4)
    ArgAB : Argument between alpha and beta in GHF (np.arctan2(sqrt(MagA), 
             sqrt(MagB)))                              (num. component >= 4)
    '''                                                       
  filename = sys.argv[-1]
  basename = filename.replace(".cube",'')
  with Cube(filename) as atom:
    ncubs = atom.get_numCubes()
    for i in range(1,ncubs+1):
      atom.write_to_file(atom.get_volRA(cubeSelect=i), \
                         basename+"_splitRa_%i.cube" %i, cubeSelect=i)
      atom.write_to_file(atom.get_volRB(cubeSelect=i), 
                         basename+"_splitRb_%i.cube" %i, cubeSelect=i)
      atom.write_to_file(np.sqrt(np.square(atom.get_volRA(cubeSelect=i)) + \
                         np.square(atom.get_volIA(cubeSelect=i))), \
                         basename+"_splitMaga_%i.cube" %i, cubeSelect=i)
      atom.write_to_file(np.sqrt(np.square(atom.get_volRB(cubeSelect=i)) + \
                         np.square(atom.get_volIB(cubeSelect=i))), \
                         basename+"_splitMagb_%i.cube" %i, cubeSelect=i)
      atom.write_to_file(np.arctan2(atom.get_volIA(cubeSelect=i), \
                         atom.get_volRA(cubeSelect=i)), \
                         basename+"_splitArga_%i.cube" %i, cubeSelect=i)
      atom.write_to_file(np.arctan2(atom.get_volIB(cubeSelect=i), \
                         atom.get_volRB(cubeSelect=i)), \
                         basename+"_splitArgb_%i.cube" %i, cubeSelect=i)
      atom.write_to_file(np.sqrt(np.square(atom.get_volRA(cubeSelect=i)) + \
                         np.square(atom.get_volIA(cubeSelect=i)) + \
                         np.square(atom.get_volRB(cubeSelect=i)) + \
                         np.square(atom.get_volIB(cubeSelect=i))), \
                         basename+"_splitMagab_%i.cube" %i, cubeSelect=i)
      atom.write_to_file(np.arctan2(np.sqrt( \
                         np.square(atom.get_volRA(cubeSelect=i) + \
                         np.square(atom.get_volIA(cubeSelect=i)))), \
                         np.sqrt(np.square(atom.get_volRB(cubeSelect=i)) + \
                         np.square(atom.get_volIB(cubeSelect=i)))), \
                         basename+"_splitArgab_%i.cube" %i, cubeSelect=i)
