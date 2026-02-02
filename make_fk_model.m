
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

% % Plot velocity model
% figure;
% subplot(3,1,1)
% %imagesc((0:nx-1)*dh/1000, z/1000, vp_true);
% imagesc(vp);
% hold on 
% plot(rec_x_id,rec_z_id,'-ow')
% axis equal tight
% colorbar;
% %set(gca, 'YDir', 'normal');  % depth increasing downward
% xlabel('x (km)');
% ylabel('Depth (km)');
% title('2D background Vp model');
% 
% subplot(3,1,2)
% %imagesc((0:nx-1)*dh/1000, z/1000, vs_true);
% imagesc(vs);
% axis equal tight
% colorbar;
% %set(gca, 'YDir', 'normal');  % depth increasing downward
% xlabel('x (km)');
% ylabel('Depth (km)');
% title('2D background Vs model');
% 
% subplot(3,1,3)
% %imagesc((0:nx-1)*dh/1000, z/1000, rho_true);
% imagesc(rho);
% axis equal tight
% colorbar;
% %set(gca, 'YDir', 'normal');  % depth increasing downward
% xlabel('x (km)');
% ylabel('Depth (km)');
% title('2D background rho model');
% 
