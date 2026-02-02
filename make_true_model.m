%True model
%Marmoussi 2 model (local)
% f=fopen('Elastic_models/marmousi_II_marine.vs');
% vs_true=fread(f,'float32');
% vs_true=reshape(vs_true,[nz nx]);
% f=fopen('Elastic_models/marmousi_II_marine.vp');
% vp_true=fread(f,'float32');
% vp_true=reshape(vp_true,[nz nx]);
% f=fopen('Elastic_models/marmousi_II_marine.rho');
% rho_true=fread(f,'float32');
% rho_true=reshape(rho_true,[nz nx]);
%TELESEISMIC CONFIGURATION (layered model)
%first layer
% x=1:1:nx;
% a=10;
% t=200;
% f=1/t;
% d=50;
% h=63;
% y=a*sin(2*pi*f*(x-d))+h;
% for i=1:nx
%     if(i<=50 || i>=nx-50)
%         y(i)=h;
%     end
% end
% for j=1:nz
%     for i=1:nx
%         if(j<round(y(i)))
%             vp_true(j,i)=5800.0883999999996;
%             vs_true(j,i)=3198.5574000000001;
%             rho_true(j,i)=2600;
%         else
%             %Second layer (30km)
%             vp_true(j,i)=8079.9750000000004;
%             vs_true(j,i)=4485.3476000000001;
%             rho_true(j,i)=3380;
%         end
%     end
% end
if (plane_wave)
vp_true=ones(nz,nx)*al_fk(1);
vs_true=ones(nz,nx)*be_fk(1);
rho_true=ones(nz,nx)*rho_fk(1);
%first layer (20km)
vp_true(round(h_fk(1)/dh)+3:end,:)=al_fk(2);
vs_true(round(h_fk(1)/dh)+3:end,:)=be_fk(2);
rho_true(round(h_fk(1)/dh)+3:end,:)=rho_fk(2);
%second layer (15km)
vp_true(round(h_fk(1)/dh)+round(h_fk(2)/dh)+3:end,:)=al_fk(3);
vs_true(round(h_fk(1)/dh)+round(h_fk(2)/dh)+3:end,:)=be_fk(3);
rho_true(round(h_fk(1)/dh)+round(h_fk(2)/dh)+3:end,:)=rho_fk(3);
end

%A SYNTHETIC KINDA REAL SUBDUCTION MODEL 
%ACAP stations is at node 48
%15° DIPPING SLAB from ACAP 75 km inland
%i=48:48+84; or 1:84+48*2? 
%topslab=m
%bottslab=
%then it becomes planar for 150 km 
%84+48+167
x=132:132+167; 
y=43; 
x=[299 366];
y=[43 228];
dipslab=(y(2)-y(1))/(x(2)-x(1));
x=1:nx; 
tds=dipslab*x-782.6; 
bds=dipslab*x-658;

x=[-18 132];
y=[43 88];
dipslab=(y(2)-y(1))/(x(2)-x(1));
x=1:nx; 
bds30=dipslab*x+48;
%%%%%%%%CUTE SLAB MODEL
% for j=1:nz
%     for i=1:nx
%         %15° dipping slab 
%         if i>=33 && i<=132
%             if j<round(bds30(i)) && j>=43
%                 vp_true(j,i)=al_fk(3)+al_fk(3)*0.1; 
%                 vs_true(j,i)=be_fk(3)+be_fk(3)*0.1; 
%                 rho_true(j,i)=rho_fk(3)+rho_fk(3)*0.1; 
%             end
% 
%         end
%         %Horizontal slab
%         if i>=132 && i<=132+167
%             vp_true(43:43+45,i)=al_fk(3)+al_fk(3)*0.1;
%             vs_true(43:43+45,i)=be_fk(3)+be_fk(3)*0.1;
%             rho_true(43:43+45,i)=rho_fk(3)+rho_fk(3)*0.1;
%         end
%         %70° dipping slab
%          if i>=254 && i<=366
%              if j>=round(tds(i)) && j<=round(bds(i)) && j>43
%                  vp_true(j,i)=al_fk(3)+al_fk(3)*0.1;
%                  vs_true(j,i)=be_fk(3)+be_fk(3)*0.1; 
%                 rho_true(j,i)=rho_fk(3)+rho_fk(3)*0.1; 
%              end
%          end
% %             if i<=299-45 && i<=321
% %                if j<=round(tds(i)) && j>=round(bds(i))
% %                    vp_true(j,i)=al_fk(3)+al_fk(3)*0.05;
% %                end
% %            end         
%     end
% end
%%%%%%%%END OF CUTE SLAB MODEL
% vp_true=ones(nz,nx)*al_fk(1);
% vs_true=ones(nz,nx)*be_fk(1);
% rho_true=ones(nz,nx)*rho_fk(1);
%Second layer (15km)
%vp_true(round(h_fk(1)/dh)+3:end,:)=al_fk(2);
%vs_true(round(h_fk(1)/dh)+3:end,:)=be_fk(2);

% figure
% imagesc(vp_true)
% hold on 
% %plot(x,tds,'ko')
% %plot(x,bds,'ko')
% %plot(x,bds30,'ko')
% plot(rec_x_id,rec_z_id,'wo')
% axis equal tight
% %Dipping 19 layer
% k=15; %distance away from upper  layer
% for i=15:3:nx-15 %10 - 10 moves laterally the slab
%     vp_true(60+k:140+k,i)=8079.975+8079.975*0.05;
%     vs_true(60+k:140+k,i)=4485.347+4485.437*0.05;
%     rho_true(60+k:140+k,i)=3380+3380*0.05;
%     
%     if(i>=1 && i<nx)
%         vp_true(60+k:140+k,i+1)=8079.975+8079.975*0.05;
%         vs_true(60+k:140+k,i+1)=4485.347+4485.437*0.05;
%         rho_true(60+k:140+k,i+1)=3380+3380*0.05;
%         vp_true(60+k:140+k,i-1)=8079.975+8079.975*0.05;
%         vs_true(60+k:140+k,i-1)=4485.347+4485.437*0.05;
%         rho_true(60+k:140+k,i-1)=3380+3380*0.05;
%     end
%     %         if(i>1 && i<nx)
%     %             vp_true(100+k:230+k,i-1)=8079.975+8079.975*0.05;
%     %            vs_true(100+k:230+k,i-1)=4485.347+4485.437*0.05;
%     %             rho_true(100+k:230+k,i-1)=3380+3380*0.05;
%     %         end
%     k=k+2; %this factor changes the angle of slab
% end

%Max and min elastic parameters values for hard constraints:
vp_min=min(min(vp_true));%Vp
vp_max=max(max(vp_true));
vp_mean=(vp_min+vp_max)/2;
vs_min=min(min(vs_true));%Vs
vs_max=max(max(vs_true));
vs_mean=(vs_min+vs_max)/2;
rho_min=min(min(rho_true));%Rho
rho_max=max(max(rho_true));
rho_mean=(rho_min+rho_max)/2;