#!/bin/bash

PDB=m1_107-29r.pdb
FORCE=50
JSON=screen_m1_chi2.fps.json
CHI2="chi2_all"
DIR="."

basename=$(basename "${PDB}")
basename="${basename%.*}"

mkdir -p "${DIR}/Equil"
#append DUmmy atoms to a PDB. DU atoms represent AV mean positions.
python3 /home/dimura/workspace/FRETrest/placeAVmp.py -p "${PDB}" -o "${DIR}/${basename}_PA.pdb" -j "${JSON}" --chi2 "${CHI2}" #--avprefix 'AV_'

#prepare AMBER input files
tmpfl=`mktemp /tmp/prep-leap-XXXX`
cat <<EOF > $tmpfl
source leaprc.protein.ff14SB
set default PBradii mbondi2
source leaprc.water.tip3p

loadamberparams DUM.frcmod
loadOff DUM.lib

PDB = loadpdb "${DIR}/${basename}_PA.pdb"
addions PDB NA 0
addions PDB CL 0
solvateoct PDB TIP3PBOX 11.0
saveamberparm PDB ${DIR}/${basename}_watio.prmtop ${DIR}/${basename}_watio.inpcrd
quit
EOF

tleap -f $tmpfl
rm $tmpfl
mv leap.log "${DIR}/"

#create a topology with HMR
tmpfl=`mktemp /tmp/prep-parmed-XXXX`
cat <<EOF > $tmpfl
parm ${DIR}/${basename}_watio.prmtop
HMassRepartition
parmout ${DIR}/${basename}_watio_hmr.prmtop
go
EOF
parmed -i $tmpfl
rm $tmpfl

#Run equilibration
cp Equil_template/* "${DIR}/Equil/"
RND_SEED=`od -A n -t d -N 2 /dev/urandom | tr -d ' '`
LASTRESI=`grep ' DU   DU ' "${DIR}/${basename}_PA.pdb" | tail -n 1 | awk '{print $6}'`
sed -i -e "s/ig     = 71277,/ig     = ${RND_SEED},/g" -e "s/RES 1 TEMPLATE_LASTRESI/RES 1 ${LASTRESI}/g" ${DIR}/Equil/*.in
(cd ${DIR}/Equil && ./equil.sh "../${basename}_watio_hmr.prmtop" "../${basename}_watio.inpcrd")

#Create an AMBER (DISANG) restraits file and update DU positions in the restart file.
python3 /home/dimura/workspace/FRETrest/FRETrest.py -t "${DIR}/${basename}_watio_hmr.prmtop" -r ${DIR}/Equil/26_md_nvt_red_pme_11.restrt -j "${JSON}" --chi2 "${CHI2}" --fout ${DIR}/prod.f --restout ${DIR}/prod_0000.restrt --force ${FORCE} --resoffset -1
ambpdb -p "${DIR}/${basename}_watio_hmr.prmtop" -c ${DIR}/prod_0000.restrt > ${DIR}/prod_0000.pdb

#Run production
cp prod.in ${JSON} "${DIR}/"
echo "Running production: "$( date --rfc-3339=s )
(cd ${DIR} && pmemd.cuda -O -i "prod.in" -o "prod_0001.out" -p "${basename}_watio_hmr.prmtop" -c "prod_0000.restrt" -r "prod_0001.restrt" -x "prod_0001.mdcrd")
echo "Done: "$( date --rfc-3339=s )
ambpdb -p "${DIR}/${basename}_watio_hmr.prmtop" -c ${DIR}/prod_0001.restrt > ${DIR}/prod_0001.pdb
