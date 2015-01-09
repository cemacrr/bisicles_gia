#!/bin/bash
#
# job script  for bluecrystal phase 3
#
#PBS -l walltime=12:00:00,nodes=@NODES:ppn=16
#request  nodes ,ppn processor per node
#
#PBS -j oe
# 
#get the job id no
export JOBNO="`echo $PBS_JOBID | sed s/.master.ic.cluster//`"

RUNDIR=/panfs/panasas01/geog/ggslc/intel/BISICLES/BISICLES/examples/MISMIPNG/lazyspinupshortdeep
cd $RUNDIR



export INFILEBASE=@INFILE
export INFILE=$INFILEBASE.$JOBNO
echo "INFILE = $INFILE"
cp $INFILEBASE $INFILE

#work out what the latest checkpoint file is (if it exists)
if test -n "$(find . -maxdepth 1 -name 'chk*@SMOD.@LEVlev.@SOLVER.*.2d.hdf5' -print -quit)"
    then
    LCHK=`ls -th chk*@SMOD.@LEVlev.@SOLVER.*.2d.hdf5 | head -n 1`
    echo "" >> $INFILE #ensure line break
    echo "amr.restart_file=$LCHK" >> $INFILE
    echo "" >> $INFILE #ensure line break
fi

export MYEX="/panfs/panasas01/geog/ggslc/intel/BISICLES/BISICLES/code/exec2D/driver2d.Linux.64.mpiCC.ifort.DEBUG.OPT.MPI.ex"


export CMD="$MYEX $INFILE"
echo "CMD = $CMD"
#
echo "PBS_NODEFILE = $PBS_NODEFILE"
nodelist=`cat $PBS_NODEFILE`

#
echo "JOBNO = $JOBNO" 
echo "nodelist = $nodelist"
#
#   name the configuration file
export CONFILE="$RUNDIR/ib.$JOBNO.conf"
echo "CONFILE =  $CONFILE"
#
#   put nodes in configuration file
#
#number of cores on each node
export cpn=16


for iii in `cat $PBS_NODEFILE | uniq` 
do 
  j=0
  while [ $j -lt $cpn ]
  do
        echo $iii >> $CONFILE
	let "j=j+1"
  done
done

#
#   get the number of processors
NUMPROC=`cat $CONFILE | wc -l`
export NUMPROC
echo $NUMPROC processors
#
#    run job
#
cm-launcher mpirun -np $NUMPROC -x PYTHONPATH=$RUNDIR/../ -machinefile $CONFILE -x CH_TIMER=1  $CMD;
