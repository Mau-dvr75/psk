%Initializing dynamic fields 

%%stresses 
sxx=zeros(Nz,Nx);
sxx0=sxx; 
szz=zeros(Nz,Nx);
szz0=szz;
sxz=zeros(Nz,Nx);
sxz0=sxz; 
szx=zeros(Nz,Nx);
szx0=szx;

Fsxx=zeros(Nz,Nx);
Fszz=zeros(Nz,Nx);
Fsxz=zeros(Nz,Nx);
%Fszx=zeros(Nz,Nx);

%velocities 
vx=zeros(Nz,Nx);
vx0=vx; 
vz=zeros(Nz,Nx);
vz0=vz; 
Fvx=zeros(Nz,Nx);
Fvz=zeros(Nz,Nx);

vx_obs=zeros(nrec,nt);
ux_obs=zeros(nrec,nt);
vz_obs=zeros(nrec,nt);
uz_obs=zeros(nrec,nt);

vx_syn=zeros(nrec,nt);
ux_syn=zeros(nrec,nt);
vz_syn=zeros(nrec,nt);
uz_syn=zeros(nrec,nt);

%displacements
ux=zeros(Nz,Nx); 
uz=zeros(Nz,Nx); 
Fux=zeros(Nz,Nx); 
Fuz=zeros(Nz,Nx); 

%spatial derivatives
sxx_x=zeros(Nz,Nx); 
sxz_z=zeros(Nz,Nx); 
sxz_x=zeros(Nz,Nx);
szz_z=zeros(Nz,Nx);
vx_x=zeros(Nz,Nx);
vz_z=zeros(Nz,Nx);
vx_z=zeros(Nz,Nx);
vz_x=zeros(Nz,Nx);

%Memory variables for CPML (for forward simulation)
mem_vxx=zeros(Nz,Nx);
mem_vxz=zeros(Nz,Nx);
mem_vzx=zeros(Nz,Nx);
mem_vzz=zeros(Nz,Nx);
mem_sxx_x=zeros(Nz,Nx);
mem_szz_z=zeros(Nz,Nx);
mem_sxz_x=zeros(Nz,Nx);
mem_sxz_z=zeros(Nz,Nx);

%for adjoint simulation 
Amem_vxx=zeros(Nz,Nx);
Amem_vxz=zeros(Nz,Nx);
Amem_vzx=zeros(Nz,Nx);
Amem_vzz=zeros(Nz,Nx);
Amem_sxx_x=zeros(Nz,Nx);
Amem_szz_z=zeros(Nz,Nx);
Amem_sxz_x=zeros(Nz,Nx);
Amem_sxz_z=zeros(Nz,Nx);

%%stresses 
Asxx=zeros(Nz,Nx);
Aszz=zeros(Nz,Nx);
Asxz=zeros(Nz,Nx);
Aszx=zeros(Nz,Nx);

%velocities 
Avx=zeros(Nz,Nx);
Avz=zeros(Nz,Nx);

%displacements
Aux=zeros(Nz,Nx); 
Auz=zeros(Nz,Nx);

%Kernels
%rho
interaction.rhox=0; 
interaction.rhoz=0; 
K.rhox=zeros(Nz,Nx);
K.rhoz=zeros(Nz,Nx); 
K_rho=zeros(Nz,Nx); 
%Mu
interaction.mu=0;
K_mu=zeros(Nz,Nx); 
%Lambda
interaction.lamda=0;
K_lambda=zeros(Nz,Nx); 

if (store==1)
    vx_forward=zeros(nt/sfe,Nz,Nx);
    vz_forward=zeros(nt/sfe,Nz,Nx);
    % % displacement
    ux_forward=zeros(nt/sfe,Nz,Nx);
    uz_forward=zeros(nt/sfe,Nz,Nx);
else
    vx_forward=0;
    vz_forward=0;
    % % displacement
    ux_forward=0;
    uz_forward=0;
end

%Auxiliary fields to compute acceleration 
vxt=zeros(Nz,Nx); 
vzt=zeros(Nz,Nx);
sxxt=zeros(Nz,Nx); %time derivative of stresses
szzt=zeros(Nz,Nx);
sxzt=zeros(Nz,Nx);

%pseudo hessian 
Hrho=zeros(Nz,Nx);
Hlam=zeros(Nz,Nx); 
Hmu=zeros(Nz,Nx); 



