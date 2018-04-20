import FRETrest as fr
import numpy as np 
import mdtraj as md
import argparse
import os
import json

def main():
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
  selDistList=fr.selectedDistances(jdata,chi2Name)
  if selDistList is None:
    parser.error("evaluator "+chi2Name+" is not found in "+jsonPath)
  
  selLPs=fr.selectedLPs(jdata,selDistList)
  frame=md.load_frame(inPath,0)
  
  resName='DU'
  atName='DU'
  el=md.element.Element(2,atName,'', 0, 0)
  topMP=md.Topology()
  
  xyzMP=np.zeros([1, len(selLPs), 3])
  resSeqLast=frame.topology.residue(frame.topology.n_residues-1).resSeq
  lpNames=sorted(selLPs.keys())
  for ilp,lpName in enumerate(lpNames):
    chain=topMP.add_chain()
    res=topMP.add_residue(resName,chain,resSeqLast+ilp+1)
    topMP.add_atom(atName,el,res)
    av=fr.getAV(frame,selLPs[lpName])
    if np.max(av.grid)<=0.0:
      print('ERROR! Could not build an AV for position {}. Is it buried?'.format(lpName))
      return
    xyzMP[0,ilp,:]=fr.avMP(av)*0.1
    if len(avPrefix)>0:
      fr.savePqr(avPrefix+lpName+'.pqr',av)
  trajMP=md.Trajectory(xyz=xyzMP, topology=topMP)
  out=frame.stack(trajMP)
  out.save_pdb(outPath)
  
if __name__ == "__main__":
    main()