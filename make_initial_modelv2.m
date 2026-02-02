for i=1:length(h_fk)
    h_fkm(i)=sum(h_fk(1:i));
end

nLayers = length(al_fk);

%Create 2D velocity model
vp_true= zeros(nz, nx);
vs_true= zeros(nz, nx);
rho_true= zeros(nz, nx);
% Depth axis (z = 0 at grid point 3)
z = ((1:nz) - 3) * dh;  % km, positive downward

% Assign velocities
for i = 1:nz
    % Find which layer the current depth belongs to
    layer_idx = find(z(i) < h_fkm, 1, 'first');

    if isempty(layer_idx)
        % Below last interface → last layer
        vp_true(i,:) = al_fk(end);
        vs_true(i,:) = be_fk(end);
        rho_true(i,:) = rho_fk(end);
    else
        vp_true(i,:) = al_fk(layer_idx);
        vs_true(i,:) = be_fk(layer_idx);
        rho_true(i,:) = rho_fk(layer_idx);
    end
end
vp=vp_true;
vs=vs_true;
rho=rho_true;

% ini=round(bz/2)-5;
% fin=nx-round(bz/2)+4;
% 
% intvp=vp(:,ini:fin);
% intvs=vs(:,ini:fin);
% intrho=rho(:,ini:fin);
% 
% sigma=8; pts=0;
% a=filter_2Dfield(intvp,sigma,pts);
% b=filter_2Dfield(intvs,sigma,pts);
% c=filter_2Dfield(intrho,sigma,pts);
% 
% vp(:,ini:fin)=a;
% vs(:,ini:fin)=b;
% rho(:,ini:fin)=c;

vp=extend_model(vp,nx,nz,ncpml);
vs=extend_model(vs,nx,nz,ncpml);
rho=extend_model(rho,nx,nz,ncpml);

sigma=10; pts=0;
vp_sm=filter_2Dfield(vp,sigma,pts);
vs_sm=filter_2Dfield(vs,sigma,pts);
rho_sm=filter_2Dfield(rho,sigma,pts);

x = 0:Nx-1;
z = 0:Nz-1;

[XX, ZZ] = meshgrid(x, z);
%Define inner heterogeneous region
x1 = ibi+3; x2 = iei-4;   % km
z1 = -floor(bz/2);  z2 = jei-3;   % km
 
inner_mask = (XX >= x1 & XX <= x2 & ZZ >= z1 & ZZ <= z2);
 
%Create heterogeneous perturbation (±3%)
% dvp_max = 0.03;
% perturb = dvp_max * sin(2*pi*XX/80) .* sin(2*pi*ZZ/60);
% perturb(~inner_mask) = 0;
 
%Smooth cosine taper
taper_width = floor(bz/2)-1; % km  %if no taper at top, make this equal to z1 
taper_widthz=10;%round(taper_width/2);
dx_dist = min(abs(XX - x1), abs(XX - x2));
dz_dist = min(abs(ZZ - z1), abs(ZZ - z2));
d = min(dx_dist, dz_dist);
 
taper = ones(size(XX));
idx = d < taper_width;
taper(idx) = 0.5 * (1 - cos(pi * d(idx)/taper_width));
 
taper(~inner_mask) = 0;

W = make_2d_taper(x, z, x1+taper_width, x2-taper_width, z2-taper_widthz, taper_width,taper_widthz);
taper=W;

%Final 2D velocity model
vp = (1 - taper) .* vp + taper .* vp_sm;
vs = (1 - taper) .* vs + taper .* vs_sm;
rho = (1 - taper) .* rho + taper .* rho_sm;

figure
imagesc(vp)
colormap jet
axis equal tight 
hold on 
plot(rec_x_id,rec_z_id,'wo')
xline(ncpml+1,'r')
xline(Nx-ncpml,'r')
yline(Nz-ncpml,'r')
xline(ibi,'w')
xline(iei,'w')
yline(jei,'w')

% figure
% imagesc(taper) 
% colormap jet
% hold on 
% plot(rec_x_id,rec_z_id,'wo')
% xline(ncpml+1,'r')
% xline(Nx-ncpml,'r')
% yline(Nz-ncpml,'r')
% xline(ibi,'w')
% xline(iei,'w')
% yline(jei,'w')
% axis equal tight 

vp=vp(nz0:nzf,nx0:nxf);
vs=vs(nz0:nzf,nx0:nxf);
rho=rho(nz0:nzf,nx0:nxf);
 
% figure
% imagesc(vp)
% colormap jet
% axis equal tight 
% hold on 
% plot(rec_x_id,rec_z_id,'wo')
% xline(ncpml+1,'r')
% xline(Nx-ncpml,'r')
% yline(Nz-ncpml,'r')
% xline(ibi,'w')
% xline(iei,'w')
% yline(jei,'w')
