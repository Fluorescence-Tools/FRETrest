# -*- coding: utf-8 -*-
import argparse
import LabelLib as ll
import numpy as np 
import mdtraj as md
import argparse
import os
import json
import re
import sys


def main():
  parser = argparse.ArgumentParser(description='Create a FRET restraint file for AMBER and update pseudo atom positions.')
  parser.add_argument('-t','--top', required=True, type=str,
		      help='Input topology file path')
  parser.add_argument('-r','--restin', required=True, type=str,
		      help='Input restart file path')
  parser.add_argument('-j','--json', required=True, type=str,
		      help='FRET-restraint file path in .fps.json format')
  parser.add_argument('--restout', required=True, type=str,
		      help='Output restart file path with updated pseudoatom positions')
  parser.add_argument('--fout', required=True, type=str,
		      help='Output restraints file path (AMBER NMR restraints format, DISANG)')
  parser.add_argument('--force', required=True, type=float,
		      help='Maximum force applied between pseudoatoms in piconewtos')
  # Optional argument
  parser.add_argument('--chi2', type=str,
		      help='Name of the chi2[r] to consider. If this option is provided, Labelling Positions, not used by the given chi2[r] will be skipped.')
  parser.add_argument('--resoffset', type=int, default=0,
		      help='Residue numbering offset between the .fps.json and topology.')
  args = parser.parse_args()
  
  top=args.top
  inRestartPath=args.restin
  jsonPath=args.json
  outRestartPath=args.restout
  outDisang=args.fout
  maxForce=args.force
  chi2Name=args.chi2
  resSeqOffset=args.resoffset
  for path in [top,inRestartPath,jsonPath]:
    if not os.path.isfile(path):
      parser.error("file {} does not exist".format(path))
      
  with open(jsonPath) as jsonFile:    
    jdata = json.load(jsonFile)
    
  selDistList=selectedDistances(jdata,chi2Name)
  if selDistList is None:
    parser.error("evaluator "+chi2Name+" is not found in "+jsonPath)
  selLPs=selectedLPs(jdata,selDistList)
  
  print('#loading trajectory')
  frame=md.load(inRestartPath,top=top)[0]
  duIds=frame.topology.select('name DU')
  lpNames=sorted(selLPs.keys())
  if len(duIds) != len(lpNames):
    print('ERROR! Number of pseudoatoms in topology ({}) does not match to the number of labelling positions ({}).'.format(len(duIds),len(lpNames)))
    return
  
  #frame.image_molecules(inplace=True)
  
  #update pseudoatom positions
  print('#update pseudoatom positions')
  avs={}
  for ilp,lpName in enumerate(lpNames):
    av=getAV(frame,selLPs[lpName],resSeqOffset)
    if av is None:
      print('ERROR! Could not calculate av for labelling position '+lpName)
      return
    avs[lpName]=av
    duId=duIds[ilp]
    if np.max(av.grid)>0.0:
      mp=avMP(av)*0.1
      frame.xyz[0,duId,:]=mp
    else:
      print('ERROR! Calculation resulted in an empty AV for position {}.'.format(lpName))
      return
    print(lpName+' done!')
   
  restraints=[]
  #FRET restraints
  print('#FRET restraints')
  for dist in selDistList:
    sys.stdout.write(dist+' ... ')
    sys.stdout.flush()
    
    lp1name=jdata["Distances"][dist]["position1_name"]
    lp2name=jdata["Distances"][dist]["position2_name"]
    lp1Index=lpNames.index(lp1name)
    lp2Index=lpNames.index(lp2name)
    duId1=duIds[lp1Index]
    duId2=duIds[lp2Index]
    rda=float(jdata["Distances"][dist]["distance"])
    rtype=jdata["Distances"][dist]["distance_type"]
    R0=jdata["Distances"][dist]["Forster_radius"]
    rmp=RmpFromRda(avs[lp1name],avs[lp2name],rda,rtype,R0=R0)
    rest=AmberRestraint()
    rest.iat1=duId1+1
    rest.iat2=duId2+1
    rest.r2=rmp
    rest.r3=rest.r2
    rest.r1=rest.r2-float(jdata["Distances"][dist]["error_neg"])
    rest.r4=rest.r3+float(jdata["Distances"][dist]["error_pos"])
    rest.rk2=maxForce/(2.0*69.4786*(rest.r2-rest.r1)) #69.4786 pN = 1 kcal/mol Angstrom
    rest.rk3=maxForce/(2.0*69.4786*(rest.r4-rest.r3))
    rest.comment='{} ({}) <--> {} ({}) '.format(lp1name,rest.iat1,lp2name,rest.iat2)
    restraints.append(rest)
    
    print(' done!')
  
  #anchor restraints
  print('#anchor restraints')
  for ilp,lpName in enumerate(lpNames):
    sys.stdout.write('Position '+lpName+'. ')
    sys.stdout.flush()
    
    vmdSel=selLPs[lpName]['anchor_atoms']
    if len(vmdSel)==0:
      print("ERROR! Empty anchor atoms selection mask for position "+lpName)
      return
    anchorSel=selVmd2Mdtraj(vmdSel, resSeqOffset)
    
    selExprStr=frame.topology.select_expression(anchorSel)
    selExprStr=selExprStr.replace('[atom.index for atom in topology.atoms if ','')[:-1]
    print('Selection: '+selExprStr+'')
      
    ancIds=frame.topology.select(anchorSel)
    if len(ancIds)==0:
      print("ERROR! No anchor atoms selected for position "+lpName+":")
      print(anchorSel)
      return
    duId=duIds[ilp]
    for ancId in ancIds:
      ancXYZ=frame.xyz[0,ancId,:]
      mp=frame.xyz[0,duId,:]
      dist=np.sqrt(np.sum(np.square(ancXYZ-mp)))
      rest=AmberRestraint()
      rest.iat1=duId+1
      rest.iat2=ancId+1
      rest.r2=dist*10.0
      rest.r3=rest.r2
      rest.r1=rest.r2-1.0
      rest.r4=rest.r3+1.0
      rest.rk2=2.0*maxForce/(2.0*69.4786*(rest.r2-rest.r1)) #69.4786 pN = 1 kcal/mol Angstrom
      rest.rk3=rest.rk2
      resid2=frame.topology.atom(ancId).residue.resSeq
      atname2=frame.topology.atom(ancId).name
      rest.comment='{} ({}) <--> {}@{} ({})'.format(lpName,rest.iat1,resid2,atname2,rest.iat2)
      restraints.append(rest)

    print('Done! Atoms selected: {}\n'.format(len(ancIds)))

  with open(outDisang, "w") as text_file:
    for rest in restraints:
      text_file.write(rest.formString())
  
  saveRestart(outRestartPath,frame,inRestartPath)
 
class AmberRestraint:
  iat1=-100
  iat2=-100
  r1=0.0
  r2=0.0
  r3=0.0
  r4=0.0
  rk2 = 0.0
  rk3 = 0.0
  comment= ''
  def formString(self):
    return '#{}\n&rst iat = {}, {}, r1 = {:.3f}, r2 = {:.3f}, r3 = {:.3f}, r4 = {:.3f}, rk2 = {:.5f}, rk3 = {:.5f},\n/\n'.format(
      self.comment, self.iat1, self.iat2, self.r1, self.r2, self.r3, self.r4, self.rk2, self.rk3)

def readVelocities(inRestartPath):
  finRest = open(inRestartPath,'r')
  inlines=finRest.readlines()
  n_at=int(inlines[1].split(' ')[0])
  
  if len(inlines) < (3+n_at):
    #No velocities in the restart file
    return None, None, None
  
  vel=np.zeros([n_at,3])
  iLineVel=2+int(n_at/2.0+0.5)
  vel=np.empty([n_at,3])
  for i in range(0,n_at):
    if i%2==0:
      vel[i] = [float(x) for x in inlines[iLineVel].split()[:3]]
    else:
      vel[i] = [float(x) for x in inlines[iLineVel].split()[3:]]
      iLineVel+=1
  cell_length=np.array([[float(x) for x in inlines[-1].split()[:3]]])
  cell_angles=np.array([[float(x) for x in inlines[-1].split()[3:]]])
  return vel, cell_length, cell_angles

def saveRestart(outPath,frame,inRestartPath): 
  xyz=frame.xyz[0,:,:]*10.0
  cell_lengths=frame.unitcell_lengths*10.0
  cell_angles=frame.unitcell_angles
  time=frame.time[0]
  
  #read velocities
  vel, cell_length_rest, cell_angles_rest=readVelocities(inRestartPath)
  if vel is None:
    print('Information: No velocities in the restart file.')
    frame.save_amberrst7(outPath)
    return
  if len(vel) != frame.n_atoms:
    print('ERROR! Number of atoms in the frame ({}) and restart file ({}) do not match.'.format(frame.n_atoms,len(vel)))
    return

  if np.absolute(cell_length_rest-cell_lengths).max()>0.0001:
    print('ERROR! Cell length in the frame ({}) and restart file ({}) do not match.'.format(cell_lengths,cell_length_rest))
    return
  if np.absolute(cell_angles-cell_angles_rest).max()>0.0001:
    print('ERROR! Cell angles in the frame ({}) and restart file ({}) do not match.'.format(cell_angles,cell_angles_rest))
    return
  
  out=open(outPath, 'w')
  out.write('Amber restart file written by FRETrest\n')
  out.write('%5d%15.7e\n' % (frame.n_atoms, time))
  fmt = '%12.7f%12.7f%12.7f'

  #coordinates
  for i in range(frame.n_atoms):
      acor = xyz[i, :]
      out.write(fmt % (acor[0], acor[1], acor[2]))
      if i % 2 == 1: out.write('\n')
  if frame.n_atoms % 2 == 1: out.write('\n')
  #velocities
  for i in range(frame.n_atoms):
    avel = vel[i, :]
    out.write(fmt % (avel[0], avel[1], avel[2]))
    if i % 2 == 1: out.write('\n')
  if frame.n_atoms % 2 == 1: out.write('\n')
  if cell_lengths is not None:
      out.write(fmt % (cell_length_rest[0,0], cell_length_rest[0,1],
				cell_length_rest[0,2]))
      out.write(fmt % (cell_angles_rest[0,0], cell_angles_rest[0,1],
				cell_angles_rest[0,2]) + '\n')
  out.flush()
  
def av2points(grid):
  area = grid.shape[0] * grid.shape[1]
  nx, ny, nz = grid.shape
  ox, oy, oz = grid.originXYZ
  dx = grid.discStep
  g = np.array(grid.grid).reshape((nx, ny, nz),order='F')
  
  n_points=(g>0.0).sum()
  points=np.empty([n_points,4])
  i=0
  for iz in range(nz):
    for iy in range(ny):
      for ix in range(nx):
        val = g[ix, iy, iz]
        if val<=0.0:
          continue       
        x = ix * dx + ox
        y = iy * dx + oy
        z = iz * dx + oz
        points[i]=np.array([x,y,z,val])
        i+=1
  return points

def Rmp(av1,av2, transVec=None):
  mp1=avMP(av1)
  if transVec is None:
    transVec=np.zeros(3)
  mp2=avMP(av2)+transVec
  return np.sqrt(np.sum(np.square(mp2-mp1)))
  
def Rda2Efficiency(rda,R0):
  return 1.0/(np.power(rda/R0,6)+1.0)

def Efficiency2Rda(E,R0):
  return R0*np.power((1.0/E)-1.0,1.0/6.0)
  
def Rda(av1,av2, rtype, transVec=None, R0=None):
  if rtype=='Rmp':
    return Rmp(av1,av2, transVec)
  
  p1=av2points(av1)
  p2=av2points(av2)
  
  if transVec is not None:
    p2+=np.append(transVec,np.zeros(1))

  rand1=np.random.randint(0,len(p1),20000)
  rand2=np.random.randint(0,len(p2),20000)
  
  rdas=np.sqrt(np.sum(np.square(p2[rand2,:3]-p1[rand1,:3]),axis=1))
  weights=p2[rand2,3]*p1[rand1,3]
  
  if rtype=='RDAMean':
    return np.average(rdas,weights=weights)
  elif rtype=='RDAMeanE':
    aveE=np.average(Rda2Efficiency(rdas,R0),weights=weights)
    return Efficiency2Rda(aveE,R0)
  else:
    print('ERROR! Unknown distance type: "{}"'.format(rtype))
    return None
  
def RmpFromRda(av1, av2, rda, rtype, R0=None, tolerance=0.2, maxIter=10):
  mp1=avMP(av1)
  mp2=avMP(av2)
  drNorm=(mp2-mp1)/np.sqrt(np.sum(np.square(mp2-mp1)))
  
  shift=np.zeros(3)
  curRda=Rda(av1,av2,rtype,R0=R0)
  dev=rda-curRda
  for it in range(maxIter):
    shift+=drNorm*dev
    curRda=Rda(av1,av2,rtype,shift,R0=R0)
    dev=rda-curRda
    if abs(dev)<tolerance:
      break
  if abs(dev)>tolerance:
    print('ERROR! RmpFromRda() could not converge! Achieved deviation: {}, tolerance: {}'.format(abs(dev),tolerance))
    return None
  return np.sqrt(np.sum(np.square(mp2+shift-mp1)))
  
def selectedLPs(jdata,selDistList):
  selLPs={}
  for dist in selDistList:
    lp1name=jdata["Distances"][dist]["position1_name"]
    lp2name=jdata["Distances"][dist]["position2_name"]
    selLPs[lp1name]=jdata["Positions"][lp1name]
    selLPs[lp2name]=jdata["Positions"][lp2name]
  return selLPs

def selectedDistances(jdata,chi2Name):
  selDistList=[]
  if chi2Name is not None:
    if chi2Name in jdata['χ²']:
      selDistList=list(jdata['χ²'][chi2Name]['distances'])
    elif chi2Name in jdata['χᵣ²']:
      selDistList=list(jdata['χᵣ²'][chi2Name]['distances'])
    else:
      return None #error
  else:
    selDistList=list(jdata['Distances'].keys())
  return sorted(selDistList)
    
def chain2index(chain):
  return str(ord(chain[0])-ord('A'))

def attachmentString(lp, resSeqOffset=0):
  resi=int(lp['residue_seq_number'])+resSeqOffset
  resn=lp['residue_name']
  at=lp['atom_name']
  
  s=''
  chain=lp["chain_identifier"]
  if len(chain)>0:
    chain=chain2index(chain)
    s+='chainid {} and '.format(chain)
  s+='resSeq {} and resname {} and name {}'.format(resi,resn,at)
  return s

def selVmd2Mdtraj(sel, resSeqOffset=0):
  sel=sel.replace('resid', 'resSeq')
  p=re.compile('chain ([A-Z])')
  sel=p.sub(lambda m: 'chain '+chain2index(m.group(1)),sel)
  
  p=re.compile('resSeq ([0-9]+)')
  sel=p.sub(lambda m: 'resSeq {}'.format(int(m.group(1))+resSeqOffset),sel)
  
  p=re.compile('(resSeq [0-9]+ to )([0-9]+)')
  sel=p.sub(lambda m: '{}{}'.format(m.group(1),int(m.group(2))+resSeqOffset),sel)
  return sel

def keepString(lp, resSeqOffset=0):
  strip=selVmd2Mdtraj(lp["strip_mask"], resSeqOffset)  
  keep='not ('+strip+')'
  return keep

def avMP(grid):
  area = grid.shape[0] * grid.shape[1]
  nx, ny, nz = grid.shape
  ox, oy, oz = grid.originXYZ
  dx = grid.discStep
  g = np.array(grid.grid).reshape((nx, ny, nz),order='F')
  
  sumG = 0.0
  mp=np.zeros(3)
  for iz in range(nz):
    for iy in range(ny):
      for ix in range(nx):
        val = g[ix, iy, iz]
        if val<=0.0:
          continue
        sumG+=val
        
        x = ix * dx + ox
        y = iy * dx + oy
        z = iz * dx + oz
        mp+=np.array([x,y,z])*val
  mp=mp/sumG
  return mp

def getAV(fr,lp, resSeqOffset=0):
  avType=lp['simulation_type']
  if avType != 'AV1':
    print('ERROR! Simulation type is not supported: '+avType)
    return None
  linker_length=float(lp['linker_length'])
  linker_width=float(lp['linker_width'])
  dye_radius=float(lp['radius1'])
  disc_step=float(lp['simulation_grid_resolution'])
  attachSel=attachmentString(lp,resSeqOffset)
  try:
    iAttach=fr.topology.select(attachSel)[0]
  except IndexError:
    print('ERROR! Could not find the specified attachment atom in the topology:\n'+attachSel)
    return None
  xyzAttach=fr.xyz[0][iAttach]*10.0
  keepSel=keepString(lp, resSeqOffset)
  frSliced=fr.atom_slice(fr.topology.select(keepSel))
  radii=np.empty(frSliced.n_atoms)
  for iat in range(frSliced.n_atoms):
    radii[iat]=frSliced.topology.atom(iat).element.radius*10.0
  xyzr=np.vstack([frSliced.xyz[0].T*10.0,radii])
  
  limDist=float(lp['allowed_sphere_radius'])
  for i in range(frSliced.n_atoms):
    dist=np.sqrt(np.sum(np.square(xyzr[:3,i]-xyzAttach)))
    if dist < limDist:
      xyzr[3][i]=0.0

  av1 = ll.dyeDensityAV1(xyzr, xyzAttach, linker_length,linker_width,dye_radius, disc_step)
  return av1

def savePqr(fileName, grid):
  with open(fileName, "w") as out:
      area = grid.shape[0] * grid.shape[1]
      nx, ny, nz = grid.shape
      ox, oy, oz = grid.originXYZ
      dx = grid.discStep
      g = np.array(grid.grid).reshape((nx, ny, nz),order='F')
      
      iat = 0
      for iz in range(nz):
        for iy in range(ny):
          for ix in range(nx):
            val = g[ix, iy, iz]
            if val <= 0.0:
              continue
            
            iat += 1
            resi = int(iat / 10)
            
            x = ix * dx + ox
            y = iy * dx + oy
            z = iz * dx + oz
            
            sz = 'ATOM{0: 7}   AV  AV{1: 6}{2:12.1f}{3:8.1f}{4:8.1f}{5:8.2f}{6:7.3f}\n'
            sz = sz.format(iat, resi, x, y, z, val, dx * 0.5)
            out.write(sz)

if __name__ == "__main__":
    main()