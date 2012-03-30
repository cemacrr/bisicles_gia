for lev in 0 1 2 3 4 5; 
do
    tagcap=$(( lev - 1 )) 
    sed -e s/@MAXLEVEL/$lev/ -e s/@TAGCAP/$tagcap/ inputs.mismip3D.stnd.template > inputs.mismip3D.stnd.l$lev
    sed -e s/@MAXLEVEL/$lev/ mismip3D.stnd.template.config > mismip3D.stnd.l$lev.config
done
