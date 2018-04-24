PRMTOP=../148l_watio.prmtop 
RESINIT=../148l_watio.inpcrd
mpirun -n 6 pmemd.MPI -O -i 01_min_ntr_h_pme.in -o 01_min_ntr_h_pme.out -p $PRMTOP -c $RESINIT -r 01_min_ntr_h_pme.restrt -ref $RESINIT -x 01_min_ntr_h_pme.mdcrd
mpirun -n 6 pmemd.MPI -O -i 02_min_ntr_l_pme.in -o 02_min_ntr_l_pme.out -p $PRMTOP -c 01_min_ntr_h_pme.restrt -r 02_min_ntr_l_pme.restrt -ref $RESINIT -x 02_min_ntr_l_pme.mdcrd
mpirun -n 6 pmemd.MPI -O -i 11_md_nvt_ntr_pme.in -o 11_md_nvt_ntr_pme.out -p $PRMTOP -c 02_min_ntr_l_pme.restrt -r 11_md_nvt_ntr_pme.restrt -ref $RESINIT -x 11_md_nvt_ntr_pme.mdcrd 
