#!/bin/csh -f
touch Antarctica/.tempDummy
foreach i (Antarctica_5km.config Antarctica_5km_withshelves.config Antarctica-noShelf.config)
cp -f  ../../../../gc1/parallel/bin/$i Antarctica/;
end;
foreach i (gis_1km.config  gis_5km.config gis-1km.zeroH.txt)
cp -f ../../../../gc1/parallel/bin/$i Greenland/;
end;


