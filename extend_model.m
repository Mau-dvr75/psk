function model=extend_model(model,nx,nz,ncpml)

modelext=zeros(nz+ncpml,nx+2*ncpml);


modelext(1:nz,ncpml+1:nx+ncpml)=model;
%vp(1:nz,ncpml+1:nx+ncpml)=vp_ne;
%rho(1:nz,ncpml+1:nx+ncpml)=rho_ne;

for i=1:nx+ncpml
modelext(nz+1:end,i)=modelext(nz,i);
%vp(nz+1:end,i)=vp(nz,i);
%rho(nz+1:end,i)=rho(nz,i);
end

for j=1:nz+ncpml
    modelext(j,1:ncpml)=modelext(j,ncpml+1);
    modelext(j,nx+ncpml+1:end)=modelext(j,ncpml+nx);
    
%     vp(j,1:ncpml)=vp(j,ncpml+1);
%     vp(j,nx+ncpml+1:end)=vp(j,ncpml+nx);
%     
%     rho(j,1:ncpml)=rho(j,ncpml+1);
%     rho(j,nx+ncpml+1:end)=rho(j,ncpml+nx);
end 

model=modelext; 


end 