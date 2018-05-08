PRMTOP=$1
INPCRD=$2

#12_md_npt_ntr_pme is trippled to allow for reference updates
INLIST=(01_min_ntr_h_pme 02_min_ntr_l_pme 11_md_nvt_ntr_pme 12_md_npt_ntr_pme 12_md_npt_ntr_pme 12_md_npt_ntr_pme 21_md_nvt_red_pme 22_md_nvt_red_pme 23_md_nvt_red_pme 24_md_nvt_red_pme 25_md_nvt_red_pme 26_md_nvt_red_pme)

RESINIT="${INPCRD}"
listLen=${#INLIST[@]}
for (( i=0; i<${listLen}; i++ ));
do
  STEP=${INLIST[$i]}
  PREFIX="${STEP}_${i}"
  RESTOUT="${PREFIX}.restrt"
  EXE='pmemd.cuda'
  if (( $i == 0 )); then
    EXE='mpirun -n 6 pmemd.MPI'
  fi
  echo "Running equilibration step "$i"/"${listLen}": "$( date --rfc-3339=s )
  $EXE -O -i "${STEP}.in" -o "${PREFIX}.out" -p "$PRMTOP" -c "$RESINIT" -r "${RESTOUT}" -ref "${RESINIT}" -x "${PREFIX}.mdcrd"
  RESINIT="${RESTOUT}"
  if [ -f "${PREFIX}.mdcrd" ];then
    gzip -9 "${PREFIX}.mdcrd" &
  fi
done

wait

# mpirun -n 6 pmemd.MPI
#