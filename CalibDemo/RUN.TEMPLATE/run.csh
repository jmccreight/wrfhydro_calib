#!/bin/csh
#
# LSF batch script to run the test MPI code
#
#BSUB -P NRAL0017                       # Project 99999999
#BSUB -x                                # exclusive use of node (not_shared)
#BSUB -n 16                      # number of total (MPI) tasks
#BSUB -R "span[ptile=16]"               # run a max of  tasks per node
#BSUB -J demo_calib      # job name
#BSUB -o wrf%J.out                      # output filename
#BSUB -e wrf%J.err                      # error filename
#BSUB -W 0:50                           # wallclock time
#BSUB -q regular                        # queue

mpirun.lsf ./wrf_hydro.exe
echo $?

set writeout=1

if ($writeout == 1) then
## clean up since wrf hydro dosent allow spec of an output folder.
set outFolder=OUTPUT
mkdir $outFolder
mkdir $outFolder/00_RUN
cp namelist.hrldas $outFolder/00_RUN/.
cp hydro.namelist $outFolder/00_RUN/.
cp CHANPARM.TBL $outFolder/00_RUN/.
cp HYDRO.TBL $outFolder/00_RUN/.
cp URBPARM.TBL $outFolder/00_RUN/.
cp DISTR_HYDRO_CAL_PARMS.TBL $outFolder/00_RUN/.
cp LAKEPARM.TBL $outFolder/00_RUN/.
cp VEGPARM.TBL $outFolder/00_RUN/.
cp GENPARM.TBL $outFolder/00_RUN/.
cp MPTABLE.TBL $outFolder/00_RUN/.
cp GWBUCKPARM.TBL $outFolder/00_RUN/.
cp SOILPARM.TBL $outFolder/00_RUN/.
## LSM
mv RESTART* $outFolder/.
mv *LDASOUT_DOMAIN* $outFolder/.
mv *LSMOUT_DOMAIN* $outFolder/.
## HYDRO
mv HYDRO*DOMAIN* $outFolder/.
mv diag_hydro* $outFolder/.
mv frxst_pts_out.txt* $outFolder/.
mv qstrmvolrt_accum.txt $outFolder/.
mv *CHANOBS_DOMAIN* $outFolder/.
mv *CHRTOUT_DOMAIN* $outFolder/.
mv *CHRTOUT_GRID* $outFolder/.
mv *RTOUT_DOMAIN* $outFolder/.
mv GW_*.txt $outFolder/.
mv qlink*.txt $outFolder/.
mv *LAKEOUT* $outFolder/.
endif

exit 0

