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
export BASEDIR=/gpfs/cluster/geog/ggslc/PineIsland/
export MYDIR=$BASEDIR/lev@LEVEL/@SMOD/$JOBNO
mkdir -p $MYDIR
export MYEX="$BASEDIR/simple_bisicles pigv5.1km.@SMOD.l@LEVEL.config"
cd $MYDIR
ln -s $BASEDIR/pigv5.1km.nc
cp $BASEDIR/pigv5.1km.@SMOD.l@LEVEL.config $MYDIR
cp $BASEDIR/inputs.pigv5.1km.@SMOD.l@LEVEL $MYDIR
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



