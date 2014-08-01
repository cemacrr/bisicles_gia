function [ lo,hi ] = amr_query_domain_corners( amrID, level )
%extract the domain corners for a given level grid
load_libamrfile();
status = libpointer('int32Ptr',-1);
levelp=libpointer('int32Ptr',level);
lop=libpointer('int32Ptr',[0,0]);
hip=libpointer('int32Ptr',[1,1]);
calllib('libamrfile','amr_query_domain_corners',status,lop,hip,amrID,levelp);
amr_error(status);
lo = lop.value;
hi = hip.value;
end

