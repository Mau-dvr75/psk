%Free surface 
op=4;
fdoh=op/2; 
dthalbe=dt/2; 
dh24=1/dh; 

if op==2
   j=3; 
   hc=[1,1];
   in=-1; 
   fn=-1; 
end

if op==4
    j=3; 
    hc=[9/8,-1/24];
    in=0;
    fn=0;
end 

if op==8
   j=5
   hc=[1225.0/1024.0,-245.0/3072.0,49.0/5120.0,-5.0/7168.0];
   in=2;
   fn=2;
end 
%jbm:jem;
for i=3:nx-3%ncpml+1:nx-ncpml %3-fn
    szz(j,i)=0; 
    vxx=0; 
    vzz=0; 
    
    for m=1:fdoh
       % szz(j-m,i)=-szz(j+m,i); 
        sxz(j-m,i)=-sxz(j+m-1,i); 
      %  vxx=vxx+hc(m)*(vx(j,i+m-1)-vx(j,i-m));
      %  vzz=vzz+hc(m)*(vz(j+m-1,i)-vz(j-m,i)); 
    end 
%     vxx=vxx/dh; 
%     vzz=vzz/dh;
% 
%     %Apply CPML 
%     if(i<=ncpml && USE_PML_XMIN==1) 
%             mem_vxx(j,i)=b_x(i)*mem_vxx(j,i)+a_x(i)*vxx; 
%             vxx=vxx/K_x(i)+mem_vxx(j,i); 
%     end 
%     if(i>=nx-ncpml+1 && USE_PML_XMAX==1)
%             
%             h1=(i-nx+2*ncpml); 
%             h=i; 
%             mem_vxx(j,h1)=b_x(h1)*mem_vxx(j,h1)+a_x(h1)*vxx; 
%             vxx=vxx/K_x(h1)+mem_vxx(j,h1); 
%     end 
%     
   % vz(j-1,i)=vz(j,i)+(lamb(j,i)/(lamb(j,i)+2*mu(j,i)))*(vxx); 
  %vx(j-1,i)=vx(j+1,i)+vz(j-1,i+1)-vz(j-1,i)+vz(j,i+1)-vz(j,i);
  %trial
 % vz(j-2,i)=vz(j-1,i)+(vz(j+1,i)-vz(j,i))+(lamb(j,i)/(lamb(j,i)+2*mu(j,i)))*(vx(j+1,i)-vx(j+1,i-1)+vx(j-1,i)-vx(j-1,i-1));
  %  sxx(j,i)=sxx(j,i)-dt*((lamb(j,i)*lamb(j,i))/(lamb(j,i)+2*mu(j,i))*vxx+lamb(j,i)*vzz);
end 