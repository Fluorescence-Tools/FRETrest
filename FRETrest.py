import argparse
import LabelLib as ll
import numpy as np 
import mdtraj as md
import argparse
import os
import json
import re

def main():
  parser = argparse.ArgumentParser(description='Create a FRET restraint file for AMBER and update pseudo atom positions.')
  parser.add_argument('-t','--top', required=True, type=str,
		      help='Input topology file path')
  parser.add_argument('-r','--restin', required=True, type=str,
		      help='Input restart file path')
  parser.add_argument('-j','--json', required=True, type=str,
		      help='FRET-restraint file path in .fps.json format')
  #parser.add_argument('--restout', required=True, type=str,
		      #help='Output restart file path with updated pseudoatom positions')
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
  #outRestartPath=args.restout
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
  
  frame=md.load(inRestartPath,top=top)[0]
  duIds=frame.topology.select('name DU')
  lpNames=sorted(selLPs.keys())
  if len(duIds) != len(lpNames):
    print('ERROR! Number of pseudoatoms in topology ({}) does not match to the number of labelling positions ({}).'.format(len(duIds),len(lpNames)))
    return
  
  frame.image_molecules(inplace=True)
  
  #update pseudoatom positions
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
  
  restraints=[]
  #FRET restraints
  for dist in selDistList:
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
    rest.rk3=maxForce/(2.0*69.4786*(rest.r3-rest.r4))
    restraints.append(rest)
  
  #anchor restraints
  for ilp,lpName in enumerate(lpNames):
    vmdSel=selLPs[lpName]['anchor_atoms']
    anchorSel=selVmd2Mdtraj(vmdSel, resSeqOffset)
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
      rest.iat1=ancId+1
      rest.iat2=duId+1
      rest.r2=dist*10.0
      rest.r3=rest.r2
      rest.r1=rest.r2-1.0
      rest.r4=rest.r3+1.0
      rest.rk2=2.0*maxForce/(2.0*69.4786*(rest.r2-rest.r1)) #69.4786 pN = 1 kcal/mol Angstrom
      rest.rk3=rest.rk2
      restraints.append(rest)

  with open(outDisang, "w") as text_file:
    for rest in restraints:
      text_file.write(rest.formString())
 
class AmberRestraint:
  iat1=-100
  iat2=-100
  r1=0.0
  r2=0.0
  r3=0.0
  r4=0.0
  rk2 = 0.0
  rk3 = 0.0  
  def formString(self):
    return '&rst iat = {}, {}, r1 = {:.3f}, r2 = {:.3f}, r3 = {:.3f}, r4 = {:.3f}, rk2 = {:.5f}, rk3 = {:.5f},\n/\n'.format(
      self.iat1, self.iat2, self.r1, self.r2, self.r3, self.r4, self.rk2, self.rk3)

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
  elif rtype=='RDAMean':
    aveE=np.average(Rda2Efficiency(rdas,R0),weights=weights)
    return Efficiency2Rda(aveE,R0)
  else:
    print('ERROR! Unknown distance type: "{}"'.format(rtype))
    return None
  
def RmpFromRda(av1, av2, rda, rtype, R0=None, tolerance=0.2):
  mp1=avMP(av1)
  mp2=avMP(av2)
  drNorm=(mp2-mp1)/np.sqrt(np.sum(np.square(mp2-mp1)))
  
  shift=np.zeros(3)
  curRda=Rda(av1,av2,rtype)
  dev=rda-curRda
  while abs(dev)>tolerance:
    shift+=drNorm*dev
    curRda=Rda(av1,av2,rtype,shift,R0=R0)
    dev=rda-curRda
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
  return selDistList
    
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