getcre()
{
  case $smod in
    l1l2) cre=L1L2;;
    ssa)  cre=GlensLaw;;
    *) echo "unknown stress model"
  esac
}



for method in implicit neglect explicit; do
    for smod in l1l2; do
	for lev in 0 1 2 3 4 5; 
	do
	    tagcap=$(( lev - 1 )) 
	    getcre
	    sed -e s/@METHOD/$method/ -e s/@MAXLEVEL/$lev/ -e s/@TAGCAP/$tagcap/ -e s/@SMOD/$smod/ -e s/@CRE/$cre/ inputs.mismip3D.spinup.template > inputs.mismip3D.spinup.$method.$smod.$lev"lev"
	done
    done
done

