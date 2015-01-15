getcre()
{
  case $smod in
    l1l2) cre=L1L2;;
    ssa)  cre=GlensLaw;;
    *) echo "unknown stress model"
  esac
}

getnodes()
{
  case $lev in
    5)  NODES=1;;
    *) NODES=1
  esac
}

getsolveinterval()
{
  case $lev in
    7) SOLVEINTERVAL=16;;
    6) SOLVEINTERVAL=16;;
    5) SOLVEINTERVAL=8;;
    4) SOLVEINTERVAL=8;;
    3) SOLVEINTERVAL=4;;
    2) SOLVEINTERVAL=2;;
    1) SOLVEINTERVAL=1;;
    *) SOLVEINTERVAL=1
  esac  
}

for solver in Chombo PETSc fullPETSc
do
    for smod in l1l2 ssa;
    do
	for lev in 0 1 2 3 4 5 6; 
	do
	    tagcap=$(( lev - 1 )) 
	    getcre
	    getnodes
	    getsolveinterval
	    NAME=shortdeep.spin.$smod.l$lev.$solver
	    INFILE=inputs.$NAME
	    sed -e s/@SOLVEINTERVAL/$SOLVEINTERVAL/  -e s/#$solver// -e s/@SOLVER/$solver/ -e s/@NAME/$NAME/ -e s/@MAXLEVEL/$lev/ -e s/@TAGCAP/$tagcap/ -e s/@SMOD/$smod/ -e s/@CRE/$cre/ inputs.shortdeep.template > $INFILE
	    sed -e s/@SOLVEINTERVAL/$SOLVEINTERVAL/ -e s/@SOLVER/$solver/  -e s/@SMOD/$smod/ -e s/@NODES/$NODES/ -e s/@INFILE/$INFILE/ -e s/@LEV/$lev/  job.newblue.template.sh > job.newblue.$NAME.sh
	done
    done
done

