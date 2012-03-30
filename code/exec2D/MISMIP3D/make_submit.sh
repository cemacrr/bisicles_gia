function mksubmit() {
walltime=$1
nodes=$2
cpn=$3
level=$4
smod=$5
ext=$smod.l$level
sed -e s/@WALLTIME/$walltime/ -e s/@NODES/$nodes/ -e s/@CPN/$cpn/ -e s/@LEVEL/$level/ -e s/@SMOD/$smod/ submit.template.sh > mismip3D.stnd.$ext.sh
}

function mksubmitp075() {
walltime=$1
nodes=$2
cpn=$3
level=$4
smod=$5
ext=$smod.l$level
sed -e s/@WALLTIME/$walltime/ -e s/@NODES/$nodes/ -e s/@CPN/$cpn/ -e s/@LEVEL/$level/ -e s/@SMOD/$smod/ mismip3D.p075.template.sh > mismip3D.p075.$ext.sh
}




for smod in l1l2 ssa; 
do
  mksubmit 6:00:00 1 4 0 $smod
  mksubmit 12:00:00 1 4 1 $smod
  mksubmit 24:00:00 1 4 2 $smod
  mksubmit 24:00:00 1 4 3 $smod
  mksubmit 24:00:00 1 8 4 $smod
  mksubmit 24:00:00 1 8 5 $smod

  mksubmitp075 6:00:00 1 4 0 $smod
  mksubmitp075 12:00:00 1 4 1 $smod
  mksubmitp075 24:00:00 1 4 2 $smod
  mksubmitp075 24:00:00 1 4 3 $smod
  mksubmitp075 24:00:00 2 4 4 $smod
  mksubmitp075 50:00:00 2 4 5 $smod
  mksubmitp075 100:00:00 4 4 6 $smod
  mksubmitp075 200:00:00 4 4 7 $smod
  mksubmitp075 400:00:00 4 4 8 $smod


done
