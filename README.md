# FRETrest
Helper scripts for FRET-restrained MD simulations. Generates AMBER restraint files (DISANG).
The tool consists of two scripts: [placeAVmp.py](placeAVmp.py) and [FRETrest.py](FRETrest.py). 

**placeAVmp.py** adds dummy atoms to a PDB. Each dummy atom corresponds to a mean position of a FRET label. For the script to work, an input PDB and a FRET configuration file (.fps.json, can be generated by Olga software) must be provided.
```
python3 placeAVmp.py -p 148l_noH.pdb -o 148l_PA.pdb -j t4l.fps.json --chi2 'C3 χ²'
```
where `148l_noH.pdb` is the input PDB, `148l_PA.pdb` is the path to resulting PDB with dummy atoms added, `t4l.fps.json` is the FRET settings file and `C3 χ²` is the name of the relevant χ² evaluator from .fps.json. Labeling positions, that are present in the .fps.json, but are not relevant to the specified χ² will be ommited from the resulting PDB.

**FRETrest.py** generates AMBER restraint (DISANG) files, just like the NMR distance restraints. It also generates an updated restart file, so that if there are any inconsistencies between the conformation of the macromolecule and dummy atom positions, positions of the dummy atoms are adjusted accordingly.
```
python3 FRETrest.py -t 148l_watio_hmr.prmtop -r Equil/26_md_nvt_red_pme_11.restrt -j t4l.fps.json --chi2 'C3 χ²' --fout prod_0001.f --restout prod_0000.restrt --force 50 --resoffset -1
```
Here `148l_watio_hmr.prmtop` is an AMBER topology file, `26_md_nvt_red_pme_11.restrt` is an input restart file, `prod_0001.f` is the path to the resulting restraints file, `prod_0000.restrt` is the path to the resulting adjusted restart file. `--force 50` stands for the maximum inter-dummy force of 50 piconewtons. ` --resoffset -1` option should be provided if the residue numbering in .fps.json is shifted by one residue compared to the .prmtop file. This can happen if e.g. residue numbering starts from "1" in the .fps.json, but in the .prmtop it starts from "0", as is usual. The offset can be any integer number, positive or negative.

Complete step by step usage examples are available at [examples/T4L/t4l.sh](examples/T4L/t4l.sh) and [examples/hGBP1/hGBP1.sh](examples/hGBP1/hGBP1.sh). To use the scripts you would need a working installation of [AmberTools], and python libraries [LabelLib] and [mdtraj].

[AmberTools]: http://ambermd.org/GetAmber.php#ambertools
[LabelLib]: https://github.com/Fluorescence-Tools/labellib
[mdtraj]: https://github.com/Fluorescence-Tools/labellib
