# This is a sample PBS script
# 
#   Request 1 CPUS
#PBS -l ncpus=1

#   Request 1 hour of CPU time 
#PBS -l cput=900:00:00

#   Request that regular output and terminal output go to the same file
#PBS -j oe

# Send job to the "slow" queue (potoff3-12)
#PBS -q batch

#	Send email


echo Running on host `hostname`
echo Time is `date`

#   The following is the body of the script. By default, 
#   PBS scripts execute in your home directory, not the 
#   directory from which they were submitted. The following
#   line places you in the directory from which the job
#   was submitted.
cd RUN_DIR/
echo Directory is `pwd`

# Run job
./MC_prog_icc_O3 > out_T_`grep System in.dat | awk '{print $3}'`_K_GCMC_u_`grep ChemPot in.dat | awk '{print $3*-1}'`_r`grep OutHistSettings in.dat | awk '{print $4$5}'`.log
