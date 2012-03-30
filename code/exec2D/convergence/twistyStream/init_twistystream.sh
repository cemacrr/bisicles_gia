#!/bin/sh
EXECFILE=../../driver2d.Linux.64.g++.gfortran.OPT.ex
COMPAREEXEC=compare2d
INFILE_TEMPLATE=inputs.iceStream.template
INFILE_BASE=inputs.iceStream
COMPARE_TEMPLATE=inputs.compare.template
SCRIPTDIR=../../scripts
TEMPFILE=temp.out

CRE=L1L2
NLAYER=16
NSMOOTH=32
NPLUS=16
TAGFACTOR=0.25

#function to set tagging values
gettagval()
{
   case $RES in
     0032) TAGVAL="1000";;
     0064) TAGVAL="250";; 
     0128) TAGVAL="62.5";; 
     0256) TAGVAL="15.625";; 
     0512) TAGVAL="3.90625";; 
     1024) TAGVAL="0.9765625";; 
     2048) TAGVAL="0.244140625";; 
     *) echo "Unanticipated RES val for Tagging values";;
   esac
}

#function to set tagging values
gettagval3lev()
{
   case $RES in
     0032) TAGVAL="250";; 
     0064) TAGVAL="62.5";; 
     0128) TAGVAL="15.625";; 
     0256) TAGVAL="3.90625";; 
     0512) TAGVAL="0.9765625";; 
     1024) TAGVAL="0.244140625";; 
     *) echo "Unanticipated RES val for Tagging values";;
   esac
}

NREF=2
MAXLEVEL=0
CRSERES="0032"
FINESTRES=2048
COMPAREFILE=doCompare.single
RICHCOMPAREFILE=doRichardsonCompare.single
RUNFILE=doRuns.single
if [ -f $RUNFILE ]; then
  echo "deleting $RUNFILE"
  rm $RUNFILE
fi
if [ -f $RICHCOMPAREFILE ]; then
  echo "deleting  $RICHCOMPAREFILE"
  rm  $RICHCOMPAREFILE
fi
if [ -f $COMPAREFILE ]; then
  echo "deleting  $COMPAREFILE"
  rm  $COMPAREFILE
fi
echo "generating single-level inputs"
for RES in 0032 0064 0128 0256 0512 1024 2048  
do
    YRES=$RES
    gettagval 
    echo $RES "-- tagval = "  $TAGVAL;
    of=$INFILE_BASE.$CRE.$RES.l$MAXLEVEL
    sed -e s/@NSMOOTH/$NSMOOTH/ -e s/@CRE/$CRE/ -e s/@NLAYER/$NLAYER/ -e s/@RES/$RES/ -e s/@YRES/$RES/ -e s/@NREF/$NREF/ -e s/@TAGVAL/$TAGVAL/ -e s/@NREF1/$NREF/ -e s/@MAXLEVEL/$MAXLEVEL/  $INFILE_TEMPLATE > $of
    rcof=inputs.compare.Richardson.$CRSERES.l$MAXLEVEL
    sed -e s/@CRE/$CRE/ -e s/@FINERES/$RES/ -e s/@CRSERES/$CRSERES/ -e s/@NREF/$NREF/ -e s/@MAXLEVEL/$MAXLEVEL/  $COMPARE_TEMPLATE > $rcof
    cof=inputs.compare.$CRSERES.l$MAXLEVEL
    sed -e s/@CRE/$CRE/ -e s/@FINERES/$FINESTRES/ -e s/@CRSERES/$CRSERES/ -e s/@NREF/$NREF/ -e s/@MAXLEVEL/$MAXLEVEL/  $COMPARE_TEMPLATE > $cof
 
    outfile="run.$CRE.$RES.l$MAXLEVEL"
    innerConvergename="solverConverge/stream.$CRE.$RES.l$MAXLEVEL"
    outerConvergename="solverConverge/stream.$CRE.$RES.l$MAXLEVEL.outer"
    runcommand="$EXECFILE $of > $outfile"
    echo "echo \"doing $RES run\" " >> $RUNFILE
    echo $runcommand >> $RUNFILE
    echo "$SCRIPTDIR/innerJFNK.awk < $outfile > $TEMPFILE " >> $RUNFILE
    echo "$SCRIPTDIR/a.out $TEMPFILE  $innerConvergename" >> $RUNFILE
    echo "$SCRIPTDIR/jfnk.awk < $outfile > $outerConvergename " >> $RUNFILE
    if [ $RES!="0032" ]; then
      richcomparecommand="$COMPAREEXEC $rcof"
      echo $richcomparecommand >> $RICHCOMPAREFILE
      comparecommand="$COMPAREEXEC $cof"
      echo $comparecommand >> $COMPAREFILE
    fi
chmod +x $COMPAREFILE
chmod +x $RICHCOMPAREFILE
chmod +x $RUNFILE


CRSERES=$RES
NSMOOTH=$((NSMOOTH + NPLUS))
done 

COMPAREFILE=doCompare.2Ref
RUNFILE=doRuns.2Ref
if [ -f $RUNFILE ]; then
  echo "deleting $RUNFILE"
  rm $RUNFILE
fi
if [ -f $COMPAREFILE ]; then
  echo "deleting  $COMPAREFILE"
  rm  $COMPAREFILE
fi
FINERES=2048
NREF=2
NSMOOTH=32
MAXLEVEL=1
echo "generating nRef = 2 inputs"
for RES in 0032 0064 0128 0256 0512 
do
    YRES=$RES
    gettagval 
    echo $RES "-- tagval = "  $TAGVAL;
    of=$INFILE_BASE.$CRE.$RES.r$NREF.l$MAXLEVEL
    sed -e s/@NSMOOTH/$NSMOOTH/ -e s/@CRE/$CRE/ -e s/@NLAYER/$NLAYER/ -e s/@RES/$RES/ -e s/@YRES/$RES/ -e s/@NREF/$NREF/ -e s/@NREF1/$NREF/ -e s/@TAGVAL/$TAGVAL/ -e s/@MAXLEVEL/$MAXLEVEL/  $INFILE_TEMPLATE > $of
    cof=inputs.compare.$RES.r$NREF.l$MAXLEVEL
    sed -e s/@CRE/$CRE/ -e s/@FINERES/$FINERES/ -e s/@CRSERES/$RES/ -e s/@NREF/$NREF/ -e s/@MAXLEVEL/$MAXLEVEL/  $COMPARE_TEMPLATE > $cof
 
    outfile="run.$CRE.$RES.r$NREF.l$MAXLEVEL"
    innerConvergename="solverConverge/stream.$CRE.$RES.r$NREF.l$MAXLEVEL"
    outerConvergename="solverConverge/stream.$CRE.$RES.r$NREF.l$MAXLEVEL.outer"
    runcommand="$EXECFILE $of > $outfile"
    echo "echo \"doing $RES run\" " >> $RUNFILE
    echo $runcommand >> $RUNFILE
    echo "$SCRIPTDIR/innerJFNK.awk < $outfile > $TEMPFILE " >> $RUNFILE
    echo "$SCRIPTDIR/a.out $TEMPFILE  $innerConvergename" >> $RUNFILE
    echo "$SCRIPTDIR/jfnk.awk < $outfile > $outerConvergename " >> $RUNFILE

    comparecommand="$COMPAREEXEC $cof"
    echo $comparecommand >> $COMPAREFILE
CRSERES=$RES
NSMOOTH=$((NSMOOTH + NPLUS))
done 
chmod +x $COMPAREFILE
chmod +x $RUNFILE


COMPAREFILE=doCompare.3Lev
RUNFILE=doRuns.3Lev
if [ -f $RUNFILE ]; then
  echo "deleting $RUNFILE"
  rm $RUNFILE
fi
if [ -f $COMPAREFILE ]; then
  echo "deleting  $COMPAREFILE"
  rm  $COMPAREFILE
fi
FINERES=2048
NREF=2
NSMOOTH=32
MAXLEVEL=2
echo "generating three-level inputs" 
for RES in 0032 0064 0128 0256 
do
    YRES=$RES
    gettagval3lev 
    echo $RES "-- tagval = "  $TAGVAL;
    of=$INFILE_BASE.$CRE.$RES.r$NREF.l$MAXLEVEL
    sed -e s/@NSMOOTH/$NSMOOTH/ -e s/@CRE/$CRE/ -e s/@NLAYER/$NLAYER/ -e s/@RES/$RES/ -e s/@YRES/$RES/ -e s/@NREF/$NREF/ -e s/@NREF1/$NREF/ -e s/@TAGVAL/$TAGVAL/ -e s/@MAXLEVEL/$MAXLEVEL/  $INFILE_TEMPLATE > $of
    cof=inputs.compare.$RES.r$NREF.l$MAXLEVEL
    sed -e s/@CRE/$CRE/ -e s/@FINERES/$FINERES/ -e s/@CRSERES/$RES/ -e s/@NREF/$NREF/ -e s/@MAXLEVEL/$MAXLEVEL/  $COMPARE_TEMPLATE > $cof
 
    outfile="run.$CRE.$RES.r$NREF.l$MAXLEVEL"
    innerConvergename="solverConverge/stream.$CRE.$RES.r$NREF.l$MAXLEVEL"
    outerConvergename="solverConverge/stream.$CRE.$RES.r$NREF.l$MAXLEVEL.outer"
    runcommand="$EXECFILE $of > $outfile"
    echo "echo \"doing $RES run\" " >> $RUNFILE
    echo $runcommand >> $RUNFILE
    echo "$SCRIPTDIR/innerJFNK.awk < $outfile > $TEMPFILE " >> $RUNFILE
    echo "$SCRIPTDIR/a.out $TEMPFILE  $innerConvergename" >> $RUNFILE
    echo "$SCRIPTDIR/jfnk.awk < $outfile > $outerConvergename " >> $RUNFILE

    comparecommand="$COMPAREEXEC $cof"
    echo $comparecommand >> $COMPAREFILE

CRSERES=$RES
NSMOOTH=$((NSMOOTH + NPLUS))
done 
chmod +x $COMPAREFILE
chmod +x $RUNFILE

COMPAREFILE=doCompare.4Ref
RUNFILE=doRuns.4Ref
if [ -f $RUNFILE ]; then
  echo "deleting $RUNFILE"
  rm $RUNFILE
fi
if [ -f $COMPAREFILE ]; then
  echo "deleting  $COMPAREFILE"
  rm  $COMPAREFILE
fi
FINERES=2048
NREF=4
NSMOOTH=32
MAXLEVEL=1
echo "generating nRef = 4 inputs" 
for RES in 0032 0064 0128 0256 
do
    YRES=$RES
    gettagval 
    echo $RES "-- tagval = "  $TAGVAL;
    of=$INFILE_BASE.$CRE.$RES.$NREF.l$MAXLEVEL
    sed -e s/@NSMOOTH/$NSMOOTH/ -e s/@CRE/$CRE/ -e s/@NLAYER/$NLAYER/ -e s/@RES/$RES/ -e s/@YRES/$RES/ -e s/@NREF/$NREF/ -e s/@NREF1/$NREF/ -e s/@TAGVAL/$TAGVAL/ -e s/@MAXLEVEL/$MAXLEVEL/  $INFILE_TEMPLATE > $of
    cof=inputs.compare.$RES.r$NREF.l$MAXLEVEL
    sed -e s/@CRE/$CRE/ -e s/@FINERES/$FINERES/ -e s/@CRSERES/$RES/ -e s/@NREF/$NREF/ -e s/@MAXLEVEL/$MAXLEVEL/  $COMPARE_TEMPLATE > $cof
 
    outfile="run.$CRE.$RES.r$NREF.l$MAXLEVEL"
    innerConvergename="solverConverge/stream.$CRE.$RES.r$NREF.l$MAXLEVEL"
    outerConvergename="solverConverge/stream.$CRE.$RES.r$NREF.l$MAXLEVEL.outer"
    runcommand="$EXECFILE $of > $outfile"
    echo "echo \"doing $RES run\" " >> $RUNFILE
    echo $runcommand >> $RUNFILE
    echo "$SCRIPTDIR/innerJFNK.awk < $outfile > $TEMPFILE " >> $RUNFILE
    echo "$SCRIPTDIR/a.out $TEMPFILE  $innerConvergename" >> $RUNFILE
    echo "$SCRIPTDIR/jfnk.awk < $outfile > $outerConvergename " >> $RUNFILE

    comparecommand="$COMPAREEXEC $cof"
    echo $comparecommand >> $COMPAREFILE

CRSERES=$RES
NSMOOTH=$((NSMOOTH + NPLUS))
done 
chmod +x $COMPAREFILE
chmod +x $RUNFILE

exit 0


