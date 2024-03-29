#!/usr/env python
#Gonna need to pray to St. Isidore for this one

import numpy as np
import os.path
import psutil
import sys

def readCub(flNM):
  headder = []
  cube    = None
  val     = None

  with open(flNM, 'r') as f:
    #Comments:
    for i in range(2):
      headder.append(f.readline().strip('\n'))

    #Number atoms, orgin, n pts/vox?
    line = f.readline().strip('\n').split()
    if len(line) < 5:
      line.append('4')
    if 1 == int(line[4]):
      line[4] = 1
    natms = int(line[0])
    if natms < 1:
      natms = natms * -1
    else:
      line[0] = natms * -1
      print("WARNING, for use with GaussView we're changing this back to a", \
            "negative sign.  Sometimes this means that the units are", \
            "changed, which I don't do.  Just be careful!")
    headder.append(line)
    if not 4 == int(line[4]):
      raise ValueError("Currently files with more than 1 cube not supported")

    npts = 1
    for i in range(3):
      #Get nx, ny, nz
      line = f.readline().strip('\n').split()
      headder.append(line)
      npts = npts * int(line[0])

    for i in range(natms):
      headder.append(f.readline().strip('\n'))

    #Check if MO labels (may not be for singular cub files)
    line = f.readline().strip('\n')
    if 'E' in line:
      headder.append(' 4 1 1 1 1')
      val = np.asarray([float(v) for v in line.split()])
    else:
      headder.append(line)

    cube = np.asarray([float(v) for s in f for v in s.split()])
    if not val is None:
      #In case there was not a label
      cube = np.insert(cube, 0, val)
    
    if not len(cube) == npts:
      print("INITIAL CUBE: ", cube[0])
      print("END CUBE: ", cube[-1])
      print("FL: ", flNM)
      raise ValueError("NPTS: %i not equal to length of cube: %i?" \
                        %( npts, len(cube) ) )

  return (headder,cube)

def checkHeadder(headder1, headder):
  #Ensure both .cube files same:

  #NAtms:
  natm1 = int(headder1[2][0])
  natm2 = int(headder[2][0])
  if not natm1 == natm2:
      raise ValueError("Headders unequal number points (cur: %i, og: %i)" \
                        %( natm2, natm1 ) )
      return False

  #Atm Types:
  atms1 = np.asarray([ int(headder1[y].split()[0]) for y in range(6,natm1+7)])
  atms2 = np.asarray([ int(headder[y].split()[0]) for y in range(6,natm2+7)])
  test = atms1 == atms2
  if not test.all():
      raise ValueError("Headders unequal Atom Types")
      return False

  #Atm Positions:
  atms1 = []
  atms2 = []
  for y in range(6,natm1+7):
    atms1.append([ float(x) for x in headder1[y].split()[1:4]])
    atms2.append([ float(x) for x in headder[y].split()[1:4]])
  atms1 = np.asarray(atms1)
  atms2 = np.asarray(atms2)
  test = atms1 == atms2
  if not test.all():
      raise ValueError("Headders unequal Atom Positions")
      return False

  #Origin:
  orig1 = np.asarray([float(x) for x in headder1[2][1:4]])
  orig2 = np.asarray([float(x) for x in headder[2][1:4]])
  test = orig1 == orig2
  if not test.all():
      raise ValueError("Headders unequal origins")
      return False

  #Grid Points:
  npts1 = int(headder1[3][0]) * int(headder1[4][0]) * int(headder1[5][0])
  npts2 = int(headder[3][0]) * int(headder[4][0]) * int(headder[5][0])
  if not npts1 == npts2:
      raise ValueError("Headders unequal numper points (cur: %i, og: %i)" \
                        %( npts2, npts1 ) )
      return False
  
  #Steps X
  step1 = np.asarray([float(x) for x in headder1[3][1:4]])
  step2 = np.asarray([float(x) for x in headder[3][1:4]])
  test = step1 == step2
  if not test.all():
      raise ValueError("Headders unequal X-steps")
      return False
  
  #Steps Y
  step1 = np.asarray([float(x) for x in headder1[4][1:4]])
  step2 = np.asarray([float(x) for x in headder[4][1:4]])
  test = step1 == step2
  if not test.all():
      raise ValueError("Headders unequal Y-steps")
      return False
  
  #Steps Z
  step1 = np.asarray([float(x) for x in headder1[5][1:4]])
  step2 = np.asarray([float(x) for x in headder[5][1:4]])
  test = step1 == step2
  if not test.all():
      raise ValueError("Headders unequal Z-steps")
      return False

  return True


def write_to_file(headder, cube, fileNameOut="cub.cube"):
  '''
  Write single cube file to fileNameOut, need selected cube if you want the
  title to carry over from the combined cube file
  '''
  with open(fileNameOut, 'w') as f:
    for i in range(2):
      print(headder[i],file=f)
    
    natms = int(headder[2][0])
    if natms < 1:
      natms = natms * -1
    print(" ", *headder[2], file=f)
    
    for i in range(3,6):
      print(" ", *headder[i], file=f)
    
    for i in range(natms + 1):
      print(" ", *headder[6+i].split(), file=f)
   
    nz = int(headder[5][0])
    lineidx = 0
    for i in range(len(cube)):
      lineidx += 1
      print( "%.5e " % cube[i], file=f, end='')
      if (lineidx == 6):
        lineidx = 0
        print('', file=f)
      if (i % nz == nz-1):
        lineidx = 0
        print('', file=f)

if __name__ == '__main__':
  '''
  Recombine separate a/b r/i files into a file which can be read/split by 
  cubeSep.py (with the j25p=True flag).  Assumes that your files are arranged
  as fileRoot{a,b}{r,i}.cube.
  
**NOTE** at current time only .cube files with SINGLE CUBES contained within
  are supported.

  Usage: python cubeJoin.py fileRoot
  '''
  fileSuffix = '.cube'
  fileRoot = sys.argv[-1]
  cubOut = None
  headInit= None

  spins   = ['a', 'b']
  numbers = ['r', 'i']
  for i, spin in enumerate(spins):
    for j, number in enumerate(numbers):
      flNM = fileRoot + spin + number + fileSuffix
      if not os.path.exists(flNM):
        #File Exists?
        raise ValueError("File %s doesn't exist?" % flNM)
      
      headder, cub = readCub(flNM)
      if cubOut is None:
        #Need to initialize
        nX = int(headder[3][0])
        nY = int(headder[4][0])
        nZ = int(headder[5][0])
        cubOut = np.asarray([0. for x in range(nX * nY * nZ * 4)])
        headInit = headder
      
      elif not checkHeadder(headInit, headder):
        #Ensure same system (npts, natms, atm pos, step sizes)
        raise ValueError("File %s has different (grid/atoms) than others?" \
                          % flNM )
                          
      cubOut[(2 * i)+j::4] = cub

  #Write Output
  write_to_file(headder, cubOut, fileRoot+'_comb.cube')
