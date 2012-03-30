#!/bin/csh -f
touch Antarctica/.tempDummy
foreach i (Antarctica/Antarctica_5km.config Antarctica/Antarctica_5km_withshelves.config Antarctica/Antarctica-noShelf.config Greenland/gis_1km.config  Greenland/gis_5km.config Greenland/gis-1km.zeroH.txt)
if (-e $i) cp -f $i  ../../../../gc1/parallel/bin/;
end;
rm -f ../../lib/src/BoxTools/.tempDummy


