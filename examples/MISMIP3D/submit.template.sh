#!/bin/bash
#
#PBS -l walltime=@WALLTIME,nodes=@NODES:ppn=8
#request  nodes ,ppn processor per node
#
#PBS -j oe
# 
#get the job id no
export JOBNO="`echo $PBS_JOBID | sed s/.master.ic.cluster//`"
#create job directory
export BASEDIR=/gpfs/cluster/geog/ggslc/MISMIP3D/
export MYDIR=$BASEDIR/lev@LEVEL/@SMOD/$JOBNO
mkdir -p $MYDIR
export MYEX="$BASEDIR/so_bisicles  inputs.mismip3D.stnd.@SMOD.l@LEVEL"
cd $MYDIR
ln -s $BASEDIR/@SMOD.128.500mm.hdf5
ln -s $BASEDIR/@SMOD.256.500mm.hdf5
ln -s $BASEDIR/@SMOD.512.500mm.hdf5
ln -s $BASEDIR/@SMOD.1024.500mm.hdf5
ln -s $BASEDIR/@SMOD.2048.500mm.hdf5

cp $BASEDIR/inputs.mismip3D.stnd.@SMOD.l@LEVEL $MYDIR
#
echo $PBS_NODEFILE
nodelist=`cat $PBS_NODEFILE`
#
echo $JOBNO 
echo nodelist $nodelist
#
#   name the configuration file
export CONFILE="$MYDIR/ib.$JOBNO.conf"
echo config file $CONFILE
#
#   put nodes in configuration file
#
#number of cores on each node
export cpn=@CPN
#
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
mpirun -np $NUMPROC -machinefile $CONFILE -x CH_TIMER=1  $MYEX



