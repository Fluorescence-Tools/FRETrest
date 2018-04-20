import LabelLib as ll
import numpy as np 
import mdtraj as md
import argparse
import os
import json
import re

def main():
  #input pdb [+ chi2_name] -> pdb with PAs
  parser = argparse.ArgumentParser(description='Append Pseudoatoms, representing Accessible Volume mean positions')
  parser.add_argument('-p','--pdbin', required=True, type=str,
		      help='Input PDB path')
  parser.add_argument('-o','--pdbout', required=True, type=str,
		      help='Output PDB path')
  parser.add_argument('-j','--json', required=True, type=str,
		      help='FRET-restraint file path in .fps.json format')
  # Optional argument
  parser.add_argument('--chi2', type=str,
		      help='Name of the chi2[r] to consider. If this option is provided, Labelling Positions, not used by the given chi2[r] will be skipped.')
  parser.add_argument('--avprefix', type=str,
		      help='If this option is provided, a file named {avprefix}{lp_name}.pqr will be saved for each generated AV.')
  args = parser.parse_args()
  
  inPath=args.pdbin
  outPath=args.pdbout
  jsonPath=args.json
  chi2Name=args.chi2
  avPrefix=args.avprefix
  
  if not os.path.isfile(inPath):
    parser.error("pdbin must be an existing PDB file")
  if not os.path.isfile(jsonPath):
    parser.error("json must be an existing .fps.json file")
  if avPrefix is None:
    avPrefix=''
  
  with open(jsonPath) as jsonFile:    
    jdata = json.load(jsonFile)
  
  selDistList=[]
  if chi2Name is not None:
    if chi2Name in jdata['χ²']:
      selDistList=list(jdata['χ²'][chi2Name]['distances'])
    elif chi2Name in jdata['χᵣ²']:
      selDistList=list(jdata['χᵣ²'][chi2Name]['distances'])
    else:
      parser.error("evaluator "+chi2Name+" is not found in "+jsonPath)
  else:
    selDistList=list(jdata['Distances'].keys())
  
  selLPs={}
  for dist in selDistList:
    lp1name=jdata["Distances"][dist]["position1_name"]
    lp2name=jdata["Distances"][dist]["position2_name"]
    selLPs[lp1name]=jdata["Positions"][lp1name]
    selLPs[lp2name]=jdata["Positions"][lp2name]
  
  fr=md.load_frame(inPath,0)
  
  resName='DU'
  atName='DU'
  el=md.element.Element(2,atName,'', 0, 0)
  topMP=md.Topology()
  
  xyzMP=np.zeros([1, len(selLPs), 3])
  resSeqLast=fr.topology.residue(fr.topology.n_residues-1).resSeq
  for ilp,lpName in enumerate(selLPs):
    chain=topMP.add_chain()
    res=topMP.add_residue(resName,chain,resSeqLast+ilp+1)
    topMP.add_atom(atName,el,res)
    av=getAV(fr,selLPs[lpName])
    if np.max(av.grid)<=0.0:
      print('ERROR! Could not build an AV for position {}. Is it buried?'.format(lpName))
      return
    xyzMP[0,ilp,:]=avMP(av)*0.1
    if len(avPrefix)>0:
      savePqr(avPrefix+lpName+'.pqr',av)
  trajMP=md.Trajectory(xyz=xyzMP, topology=topMP)
  out=fr.stack(trajMP)
  out.save_pdb(outPath)
    
def chain2index(chain):
  return str(ord(chain[0])-ord('A'))

def attachmentString(lp):
  resi=int(lp['residue_seq_number'])
  resn=lp['residue_name']
  at=lp['atom_name']
  
  s=''
  chain=lp["chain_identifier"]
  if len(chain)>0:
    chain=chain2index(chain)
    s+='chainid {} and '.format(chain)
  s+='resSeq {} and resname {} and name {}'.format(resi,resn,at)
  return s

def keepString(lp):
  strip=lp["strip_mask"].replace('resid', 'resSeq')
  p=re.compile('chain ([A-Z])')
  strip=p.sub(lambda m: 'chain '+chain2index(m.group(1)),strip)
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

def getAV(fr,lp):
  linker_length=float(lp['linker_length'])
  linker_width=float(lp['linker_width'])
  dye_radius=float(lp['radius1'])
  disc_step=float(lp['simulation_grid_resolution'])
  attachSel=attachmentString(lp)
  iAttach=fr.topology.select(attachSel)[0]
  xyzAttach=fr.xyz[0][iAttach]*10.0
  
  keepSel=keepString(lp)
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