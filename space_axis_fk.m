%Space axis (properly staggered) for FK-2D method
%Remember that origin is at (0,0)m which coreesponds to
%lower left corner. For this example, free surface is assumed to be at
%upper corners, i.e. 200,000 m = 200km

%Si=3 cause I am adding 2 extra nodes in Scat field zone
%This Si could be just 1
%Outter layer
ibo=ncpml+1+si-1; %index=i, b=begin, O=outter -> ibO
ieo=Nx-ncpml-si+1; %the but e=end -> ieO, and so on...
jbo=3;
jeo=Nz-ncpml-si+1;
%Middle layer
ibm=ncpml+1+si;
iem=Nx-ncpml-si;
jbm=3;
jem=Nz-ncpml-si;
%Inner layer
ibi=ncpml+1+si+1;
iei=Nx-ncpml-si-1;
jbi=3;
jei=Nz-ncpml-si-1;

total_thickness=sum(h_fk);

%Normal stresses
for i=1:Nx
    xns(i)=(i-ibm)*dh;
end

for i=1:Nz
    zns(i)=total_thickness-(i-jbm)*dh;
end

%Shear stress
for i=1:Nx
    xss(i)=(i-ibm)*dh+dh/2;
end

for i=1:Nz
    zss(i)=total_thickness-(i-jbm)*dh-dh/2;
end


%Horizontal velocity
xvx=xss;
zvx=zns;

%Vertical velocity
xvz=xns;
zvz=zss;

% %Let's take a look at the axis
% figure(1)
% [X,Y]=meshgrid(xns,zns);
% plot(X,Y,'bs','MarkerSize',18,'MarkerFaceColor','b')
% hold on
% [X,Y]=meshgrid(xss,zss);
% plot(X,Y,'ro','MarkerSize',18,'MarkerFaceColor','r')
% [X,Y]=meshgrid(xvx,zvx);
% plot(X,Y,'g>','MarkerSize',18,'MarkerFaceColor','g')
% [X,Y]=meshgrid(xvz,zvz);
% plot(X,Y,'kv','MarkerSize',18,'MarkerFaceColor','k')
% %set(gca, 'YDir','reverse')
% axis equal tight;
% xlim([-1000 1000])
% ylim([ total_thickness-1000 total_thickness+1000])
% grid on;






x_axis.sxx=xns;
x_axis.szz=xns;
x_axis.sxz=xss;
x_axis.vx=xvx;
x_axis.vz=xvz;

z_axis.sxx=zns;
z_axis.szz=zns;
z_axis.sxz=zss;
z_axis.vx=zvx;
z_axis.vz=zvz;
