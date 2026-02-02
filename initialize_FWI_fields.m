misfit_total=0; 
misfit_sub=0;
%Kernels
%rho
interaction.rhox=0; 
interaction.rhoz=0; 
K.rhox=zeros(Nz,Nx);
K.rhoz=zeros(Nz,Nx); 
Ktotal.rho=zeros(Nz,Nx); 
Ktotal_rho=zeros(Nz,Nx);
%Mu
interaction.mu=0;
Ktotal.mu=zeros(Nz,Nx); 
Ktotal_mu=zeros(Nz,Nx);
%Lambda
interaction.lamda=0;
Ktotal.lambda=zeros(Nz,Nx); 

Ktotal_lambda=zeros(Nz,Nx);
prectotalkg=zeros(Nz,Nx);
prectotalrho=zeros(Nz,Nx);
prectotalmu=zeros(Nz,Nx);
prectotallam=zeros(Nz,Nx);


 beta_rho1=0;
 beta_rho2=0;
 beta_mu=0;
 beta_lambda=0;
 beta_vs2=0;
 beta_vp2=0;