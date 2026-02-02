%Find max and min velocity 
vs_min=min(min(vs)); 
vp_min=min(min(vp)); 
vp_max=max(max(vp)); 
vs_max=max(max(vs)); 
cmax=max(vp_max,vs_max); 
cmin=min(vs_min,vp_min);
cmin=vs_min;
fmax=2*max(f0);
%Theoretical stability criteria for FD-SG 4 order
if (op==4)
dhstab=cmin/(8*fmax); 
dtstab=dhstab/((7/6)*sqrt(2)*cmax);
end
if (op==2)
dhstab=cmin/(12*fmax); 
dtstab=dhstab/(1*sqrt(2)*cmax);
end
display(sprintf('* The recommended dh value is %.3f[m] or smaller, you specified %.3f[m]',dhstab,dh))
display(sprintf('* The recommended dt value is %.3f[s] or smaller, you specified %.3f[s]',dtstab,dt))
display(sprintf('* Points per minimun wavelength (ppw)= %.2f',cmin/(fmax*dh)))

if(dt>dtstab)
    display('Choose smaller dt')
    return
else
    display('--------------------------------')
    display('Setup is ok, starting simulation')
    display('--------------------------------')
end