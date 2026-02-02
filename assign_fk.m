si=3;
%Outter layer
ibo=ncpml+1+si-1; %index=i, b=begin, O=outter -> ibO
ieo=nx-ncpml-si+1; %the but e=end -> ieO, and so on...
jbo=3;
jeo=nz-ncpml-si+1;
%Middle layer
ibm=ncpml+1+si;
iem=nx-ncpml-si;
jbm=3;
jem=nz-ncpml-si;
%Inner layer
ibi=ncpml+1+si+1;
iei=nx-ncpml-si-1;
jbi=3;
jei=nz-ncpml-si-1;

if (INJECTIONS_ORDER==2)
    %Assign fields to arrays
    vx_rm=fk_rm.fk_rm.vx.fk;
    vz_rm=fk_rm.fk_rm.vz.fk;
    txx_rm=fk_rm.fk_rm.sxx.fk;
    tzz_rm=fk_rm.fk_rm.szz.fk;
    txz_rm=fk_rm.fk_rm.sxz.fk;
    
    vx_lm=fk_lm.fk_lm.vx.fk;
    vz_lm=fk_lm.fk_lm.vz.fk;
    txx_lm=fk_lm.fk_lm.sxx.fk;
    tzz_lm=fk_lm.fk_lm.szz.fk;
    txz_lm=fk_lm.fk_lm.sxz.fk;
    
    vx_bm=fk_bm.fk_bm.vx.fk;
    vz_bm=fk_bm.fk_bm.vz.fk;
    txx_bm=fk_bm.fk_bm.sxx.fk;
    tzz_bm=fk_bm.fk_bm.szz.fk;
    txz_bm=fk_bm.fk_bm.sxz.fk;
end

if (INJECTIONS_ORDER==4)
    %Assign fields to arrays
    vx_rm=fk_rm.fk_rm.vx.fk;
    vz_rm=fk_rm.fk_rm.vz.fk;
    txx_rm=fk_rm.fk_rm.sxx.fk;
    tzz_rm=fk_rm.fk_rm.szz.fk;
    txz_rm=fk_rm.fk_rm.sxz.fk;
    
    vx_lm=fk_lm.fk_lm.vx.fk;
    vz_lm=fk_lm.fk_lm.vz.fk;
    txx_lm=fk_lm.fk_lm.sxx.fk;
    tzz_lm=fk_lm.fk_lm.szz.fk;
    txz_lm=fk_lm.fk_lm.sxz.fk;
    
    vx_bm=fk_bm.fk_bm.vx.fk;
    vz_bm=fk_bm.fk_bm.vz.fk;
    txx_bm=fk_bm.fk_bm.sxx.fk;
    tzz_bm=fk_bm.fk_bm.szz.fk;
    txz_bm=fk_bm.fk_bm.sxz.fk;
    
    vx_ri=fk_ri.fk_ri.vx.fk;
    vz_ri=fk_ri.fk_ri.vz.fk;
    txx_ri=fk_ri.fk_ri.sxx.fk;
    tzz_ri=fk_ri.fk_ri.szz.fk;
    txz_ri=fk_ri.fk_ri.sxz.fk;
    
    vx_li=fk_li.fk_li.vx.fk;
    vz_li=fk_li.fk_li.vz.fk;
    txx_li=fk_li.fk_li.sxx.fk;
    tzz_li=fk_li.fk_li.szz.fk;
    txz_li=fk_li.fk_li.sxz.fk;
    
    vx_bi=fk_bi.fk_bi.vx.fk;
    vz_bi=fk_bi.fk_bi.vz.fk;
    txx_bi=fk_bi.fk_bi.sxx.fk;
    tzz_bi=fk_bi.fk_bi.szz.fk;
    txz_bi=fk_bi.fk_bi.sxz.fk;
    
    vx_ro=fk_ro.fk_ro.vx.fk;
    vz_ro=fk_ro.fk_ro.vz.fk;
    txx_ro=fk_ro.fk_ro.sxx.fk;
    tzz_ro=fk_ro.fk_ro.szz.fk;
    txz_ro=fk_ro.fk_ro.sxz.fk;
    
    vx_lo=fk_lo.fk_lo.vx.fk;
    vz_lo=fk_lo.fk_lo.vz.fk;
    txx_lo=fk_lo.fk_lo.sxx.fk;
    tzz_lo=fk_lo.fk_lo.szz.fk;
    txz_lo=fk_lo.fk_lo.sxz.fk;
    
    vx_bo=fk_bo.fk_bo.vx.fk;
    vz_bo=fk_bo.fk_bo.vz.fk;
    txx_bo=fk_bo.fk_bo.sxx.fk;
    tzz_bo=fk_bo.fk_bo.szz.fk;
    txz_bo=fk_bo.fk_bo.sxz.fk;
end
