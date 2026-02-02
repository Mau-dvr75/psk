muhxhz=zeros(Nz,Nx); 


%%%%%%%%%%%%%%%%%---Averaging parameters on staggered grid---%%%%%%%%%%%%%%%%%
%rho_half_x
for j=1:Nz
    for i=1:Nx-1
        rhohx(j,i)=0.5*(rho(j,i)+rho(j,i+1));  
    end 
end 
rhohx(:,Nx)=rho(:,Nx);

%rho_half_z
for j=1:Nz-1
    for i=1:Nx
        rhohz(j,i)=0.5*(rho(j,i)+rho(j+1,i)); 
    end 
end 
rhohz(Nz,:)=rho(Nz,:);

%mu_half_x_half_z
for j=1:Nz-1
    for i=1:Nx-1
        %muhxhz(j,i)=0.25*(mu(j,i)+mu(j+1,i)+mu(j,i+1)+mu(j+1,i+1));
        muhxhz(j,i)=4/(1/mu(j,i)+1/mu(j+1,i)+1/mu(j,i+1)+1/mu(j+1,i+1));
        
    end
end 
muhxhz(Nz,:)=mu(Nz,:);
muhxhz(:,Nx)=mu(:,Nx);

