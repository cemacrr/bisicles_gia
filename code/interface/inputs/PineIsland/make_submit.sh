function mksubmit() {
walltime=$1
nodes=$2
cpn=$3
level=$4
smod=$5
ext=$smod.l$level
sed -e s/@WALLTIME/$walltime/ -e s/@NODES/$nodes/ -e s/@CPN/$cpn/ -e s/@LEVEL/$level/ -e s/@SMOD/$smod/ submit.template.sh > pigv5.1km.$ext.sh
}


for smod in l1l2 ssa; 
do
  mksubmit 6:00:00 1 4 0 $smod
  mksubmit 12:00:00 1 4 1 $smod
  mksubmit 24:00:00 1 4 2 $smod
  mksubmit 24:00:00 1 8 3 $smod
  mksubmit 24:00:00 2 8 4 $smod
  mksubmit 24:00:00 4 8 5 $smod
done
