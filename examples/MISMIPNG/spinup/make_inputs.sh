getcre()
{
  case $smod in
    l1l2) cre=L1L2;;
    ssa)  cre=GlensLaw;;
    *) echo "unknown stress model"
  esac
}




for solver in Chombo PETSc
do
    for smod in l1l2 ssa;
    do
	for lev in 0 1 2 3 4 5; 
	do
	    tagcap=$(( lev - 1 )) 
	    getcre
	    
	    INFILE=inputs.ghgChannel.spin.$smod.l$lev.$solver
	    sed -e s/#$solver// -e s/@MAXLEVEL/$lev/ -e s/@TAGCAP/$tagcap/ -e s/@SMOD/$smod/ -e s/@CRE/$cre/ inputs.ghgChannel.template > $INFILE
	done
    done
done

