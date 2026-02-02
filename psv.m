function [vx_syn,vz_syn,ux_syn,uz_syn,vx_forward,vz_forward,ux_forward,uz_forward,Fvx,Fux,Fvz,Fuz,Fsxx,Fszz,Fsxz,BD_vx,BD_vz,BD_sxx,BD_szz,BD_sxz,XS,ZS,NSP,prec]=psv(vp,vs,rho,nx,nz,dh,dt,nt,stf,xs,zs,nsrc,xr,zr,nrec,FREE_SURF,saveseis,IT_DISPLAY,sfe,makeobs,jf,mode,animation,kf,a_x,a_x_half,b_x,b_x_half,K_x,K_x_half,a_z,a_z_half,b_z,b_z_half,K_z,K_z_half,ncpml,store,op,field_order,fk_lo,fk_lm,fk_li,fk_ro,fk_rm,fk_ri,fk_bo,fk_bm,fk_bi,plane_wave,INJECTIONS_ORDER,makemovie,f0,STF_TYPE,iter)
%Solve the foward 2D PSV problem with FD(2,4) with staggered grid
si=0;
% if STF_TYPE==5 || plane_wave==0
%     f0=f0(kf);
% end
%makemovie=true;
%movie of simulation 
if (makemovie==1)
    mname=strcat('wavefieldvx_shot_',num2str(jf),'_',num2str(f0(kf)),'hz.mat');
    
    obj=VideoWriter(mname);
    obj.Quality=100;
    obj.FrameRate=10;
    open(obj)
end
%X axis
x0=0;
Nx=nx+2*ncpml;% Numero total de nodos en x, incluyendo nodos con absorbencia
nx=Nx;
x=(x0-ncpml*dh):dh:(x0+(Nx-ncpml-1)*dh); % malla con nodos CPML
nx0=ncpml+1; % posicion real de x0
nxf=Nx-ncpml; % posicion real de xf
%Z axis
z0=0;% coordenada inicial real sin CPML en z
%nz=nz;% Numero total de nodos en z sin absorbencia
Nz=nz+ncpml;% Numero total de nodos en z, incluyendo nodos con absorbencia,a_x,a_x_half,b_x,b_x_half,K_x,K_x_half,a_z,a_z_half,b_z,b_z_half,K_z,K_z_half
nz=Nz;
z=(z0-ncpml*dh):dh:(z0+(Nz-ncpml-1)*dh); % malla con nodos CPML
nz0=1; % posicion real de z0
nzf=Nz-ncpml;

%Setting some constants
dtdx=dt/dh;
a0=9/8;
a1=-1/24;

% if(jf==1 && mode==0)
%     courant=max(max(vp))*dt*sqrt(1/dh^2+1/dh^2);
%     display(sprintf('The courant number is %f',courant))
% end
if(mode==0)
    courant=max(max(vp))*dt*sqrt(1/dh^2+1/dh^2);
    display(sprintf('FWD_SIM %i, courant number is %f',jf,courant))
end

for i=1:nx
    for j=1:nz
        mu(j,i)=rho(j,i)*vs(j,i)*vs(j,i);
        lamb(j,i)=rho(j,i)*vp(j,i)*vp(j,i)-2.*mu(j,i);
    end
end
%Obtaining effective parameters
make_stagger;
%Assign fk variables to arrays
if (plane_wave==true)
    assign_fk;
end
if(mode==1)
    courant=max(max(vp))*dt*sqrt(1/dh^2+1/dh^2);
    display(sprintf('ADJ_SIM %i, courant number is %f',jf,courant))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Arrays geometry for storing all borders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We save the borders according to the next geometry
%  TTTTTTTTTTTTT
%  L           R
%  L           R
%  L           R
%  L           R
%  BBBBBBBBBBBBB
if (plane_wave==true)
    hst=1;%Saving factor
    fs=3; %free surface node
    %x and z node index for top traces
    xst_TOP=ibi:hst:iei;
    zst_TOP=(fs)*ones(1,length(xst_TOP));
    %x and z node index for bottom traces
    xst_BOTTOM=ibi:hst:iei;
    zst_BOTTOM=(jei)*ones(1,length(xst_BOTTOM));
    %x and z node index for left traces
    zst_LEFT=(fs):hst:(jei-1);
    xst_LEFT=(ibi)*ones(1,length(zst_LEFT));
    %x and z node index for left traces
    zst_RIGHT=(fs):hst:(jei-1);
    xst_RIGHT=(iei)*ones(1,length(zst_RIGHT)); 
else
    hst=1;%Saving factor
    fs=3; %free surface node
    %x and z node index for top traces
    xst_TOP=nx0:hst:nxf;
    zst_TOP=(fs)*ones(1,length(xst_TOP));
    %x and z node index for bottom traces
    xst_BOTTOM=nx0:hst:nxf;
    zst_BOTTOM=(nzf)*ones(1,length(xst_BOTTOM));
    %x and z node index for left traces
    zst_LEFT=(fs+1):hst:(nzf-1);
    xst_LEFT=(nx0)*ones(1,length(zst_LEFT));
    %x and z node index for right traces
    zst_RIGHT=(fs+1):hst:(nzf-1);
    xst_RIGHT=(nxf)*ones(1,length(zst_RIGHT));
end
nx_tr=length(xst_BOTTOM);%Length of traces in x direction
nz_tr=length(zst_LEFT);%Length of traces in z direction

%Array for storing traces
%vx
seis_vx_TOP=zeros(nt,nx_tr);
seis_vx_BOT=zeros(nt,nx_tr);
seis_vx_LEF=zeros(nt,nz_tr);
seis_vx_RIG=zeros(nt,nz_tr);
%vz
seis_vz_TOP=zeros(nt,nx_tr);
seis_vz_BOT=zeros(nt,nx_tr);
seis_vz_LEF=zeros(nt,nz_tr);
seis_vz_RIG=zeros(nt,nz_tr);
%sxx
seis_sxx_TOP=zeros(nt,nx_tr);
seis_sxx_BOT=zeros(nt,nx_tr);
seis_sxx_LEF=zeros(nt,nz_tr);
seis_sxx_RIG=zeros(nt,nz_tr);
%szz
seis_szz_TOP=zeros(nt,nx_tr);
seis_szz_BOT=zeros(nt,nx_tr);
seis_szz_LEF=zeros(nt,nz_tr);
seis_szz_RIG=zeros(nt,nz_tr);
%sxz
seis_sxz_TOP=zeros(nt,nx_tr);
seis_sxz_BOT=zeros(nt,nx_tr);
seis_sxz_LEF=zeros(nt,nz_tr);
seis_sxz_RIG=zeros(nt,nz_tr);

%Initial conditions
initialize_dynamic_fields;
preckg=zeros(nz,nx);
%Main loop time
for it=1:nt
    %it
    if(mode==0 && plane_wave==false)
        vz(zs,xs)=vz(zs,xs)+stf(it);
    end
    
    if(mode==1)
        for i=1:nrec
            vx(zr(i),xr(i))=vx(zr(i),xr(i))+stf.x(i,it);
            vz(zr(i),xr(i))=vz(zr(i),xr(i))+stf.z(i,it);
        end
    end
    
    if (store==0)
        %Traces for wavefiled reconstruction
        %x direction
        for ss=1:nx_tr
            seis_vx_TOP(it,ss)=vx(zst_TOP(ss),xst_TOP(ss));
            seis_vx_BOT(it,ss)=vx(zst_BOTTOM(ss),xst_BOTTOM(ss));
            seis_vz_TOP(it,ss)=vz(zst_TOP(ss),xst_TOP(ss));
            seis_vz_BOT(it,ss)=vz(zst_BOTTOM(ss),xst_BOTTOM(ss));
            seis_sxx_TOP(it,ss)=sxx(zst_TOP(ss),xst_TOP(ss));
            seis_sxx_BOT(it,ss)=sxx(zst_BOTTOM(ss),xst_BOTTOM(ss));
            seis_szz_TOP(it,ss)=szz(zst_TOP(ss),xst_TOP(ss));
            seis_szz_BOT(it,ss)=szz(zst_BOTTOM(ss),xst_BOTTOM(ss));
            seis_sxz_TOP(it,ss)=sxz(zst_TOP(ss),xst_TOP(ss));
            seis_sxz_BOT(it,ss)=sxz(zst_BOTTOM(ss),xst_BOTTOM(ss));
        end
        %z direction
        for ss=1:nz_tr %recorremos el # de receptores (= # de trazas)
            seis_vx_LEF(it,ss)=vx(zst_LEFT(ss),xst_LEFT(ss));
            seis_vx_RIG(it,ss)=vx(zst_RIGHT(ss),xst_RIGHT(ss));
            seis_vz_LEF(it,ss)=vz(zst_LEFT(ss),xst_LEFT(ss));
            seis_vz_RIG(it,ss)=vz(zst_RIGHT(ss),xst_RIGHT(ss));
            seis_sxx_LEF(it,ss)=sxx(zst_LEFT(ss),xst_LEFT(ss));
            seis_sxx_RIG(it,ss)=sxx(zst_RIGHT(ss),xst_RIGHT(ss));
            seis_szz_LEF(it,ss)=szz(zst_LEFT(ss),xst_LEFT(ss));
            seis_szz_RIG(it,ss)=szz(zst_RIGHT(ss),xst_RIGHT(ss));
            seis_sxz_LEF(it,ss)=sxz(zst_LEFT(ss),xst_LEFT(ss));
            seis_sxz_RIG(it,ss)=sxz(zst_RIGHT(ss),xst_RIGHT(ss));
        end
    end
    
    if(op==2)
        for j=3:nz-3
            for i=3:nx-3
                sxx_x=sxx(j,i+1)-sxx(j,i);
                szz_z=szz(j+1,i)-szz(j,i);
                sxz_x=sxz(j,i)-sxz(j,i-1);
                sxz_z=sxz(j,i)-sxz(j-1,i);
                
                mem_sxx_x(j,i)=b_x_half(i)*mem_sxx_x(j,i)+a_x_half(i)*sxx_x;
                sxx_x=sxx_x/K_x_half(i)+mem_sxx_x(j,i);
                
                mem_sxz_x(j,i)=b_x(i)*mem_sxz_x(j,i)+a_x(i)*sxz_x;
                sxz_x=sxz_x/K_x(i)+mem_sxz_x(j,i);
                
                mem_szz_z(j,i)=b_z_half(j)*mem_szz_z(j,i)+a_z_half(j)*szz_z;
                szz_z=szz_z/K_z_half(j)+mem_szz_z(j,i);
                
                mem_sxz_z(j,i)=b_z(j)*mem_sxz_z(j,i)+a_z(j)*sxz_z;
                sxz_z=sxz_z/K_z(j)+mem_sxz_z(j,i);
                vx(j,i)=vx(j,i)+dtdx/rhohx(j,i)*(sxx_x+sxz_z);
                vz(j,i)=vz(j,i)+dtdx/rhohz(j,i)*(sxz_x+szz_z);
                
            end
        end
        free_vel;
        %%%%%%---------Updating stresses------%%%%%
        
        for j=3:nz-3
            for i=3:nx-3
                vxx=vx(j,i)-vx(j,i-1);
                vzz=vz(j,i)-vz(j-1,i);
                vzx=vz(j,i+1)-vz(j,i);
                vxz=vx(j+1,i)-vx(j,i);
                
                mem_vxx(j,i)=b_x(i)*mem_vxx(j,i)+a_x(i)*vxx;
                vxx=vxx/K_x(i)+mem_vxx(j,i);
                
                mem_vzx(j,i)=b_x_half(i)*mem_vzx(j,i)+a_x_half(i)*vzx;
                vzx=vzx/K_x_half(i)+mem_vzx(j,i);
                
                mem_vzz(j,i)=b_z(j)*mem_vzz(j,i)+a_z(j)*vzz;
                vzz=vzz/K_z(j)+mem_vzz(j,i);
                
                mem_vxz(j,i)=b_z_half(j)*mem_vxz(j,i)+a_z_half(j)*vxz;
                vxz=vxz/K_z_half(j)+mem_vxz(j,i);
                
                sxx(j, i) = sxx(j, i) + dtdx * ( lamb(j,i) * (vxx + vzz) + 2.0 * mu(j,i) * vxx );
                szz(j, i) = szz(j, i) + dtdx * ( lamb(j,i) * (vxx + vzz) + 2.0 * mu(j,i) * vzz );
                sxz(j, i) = sxz(j, i) + dtdx * (  muhxhz(j,i) * (vzx + vxz) );
                
            end
        end
        free_tau;  
    end
    
    if(op==4)
        for j=3:nz-3
            for i=3:nx-3
                if (j==3)
                    sxx_x=a0*(sxx(j,i+1)-sxx(j,i))+a1*(sxx(j,i+2)-sxx(j,i-1));
                    szz_z=szz(j+1,i)-szz(j,i);
                    sxz_x=a0*(sxz(j,i)-sxz(j,i-1))+a1*(sxz(j,i+1)-sxz(j,i-2));
                    sxz_z=sxz(j,i)-sxz(j-1,i);
                else
                sxx_x=a0*(sxx(j,i+1)-sxx(j,i))+a1*(sxx(j,i+2)-sxx(j,i-1));
                szz_z=a0*(szz(j+1,i)-szz(j,i))+a1*(szz(j+2,i)-szz(j-1,i));
                sxz_x=a0*(sxz(j,i)-sxz(j,i-1))+a1*(sxz(j,i+1)-sxz(j,i-2));
                sxz_z=a0*(sxz(j,i)-sxz(j-1,i))+a1*(sxz(j+1,i)-sxz(j-2,i));
                end
                mem_sxx_x(j,i)=b_x_half(i)*mem_sxx_x(j,i)+a_x_half(i)*sxx_x;
                sxx_x=sxx_x/K_x_half(i)+mem_sxx_x(j,i);
                
                mem_sxz_x(j,i)=b_x(i)*mem_sxz_x(j,i)+a_x(i)*sxz_x;
                sxz_x=sxz_x/K_x(i)+mem_sxz_x(j,i);
                
                mem_szz_z(j,i)=b_z_half(j)*mem_szz_z(j,i)+a_z_half(j)*szz_z;
                szz_z=szz_z/K_z_half(j)+mem_szz_z(j,i);
                
                mem_sxz_z(j,i)=b_z(j)*mem_sxz_z(j,i)+a_z(j)*sxz_z;
                sxz_z=sxz_z/K_z(j)+mem_sxz_z(j,i);
                vx(j,i)=vx(j,i)+dtdx/rhohx(j,i)*(sxx_x+sxz_z);
                vz(j,i)=vz(j,i)+dtdx/rhohz(j,i)*(sxz_x+szz_z);
                
            end
        end
        free_vel;
        %%%%%%---------Updating stresses------%%%%%
        
        for j=3:nz-3
            for i=3:nx-3
                if(j==3)
                    vxx=a0*(vx(j,i)-vx(j,i-1))+a1*(vx(j,i+1)-vx(j,i-2));
                    vzz=vz(j,i)-vz(j-1,i);
                    vzx=a0*(vz(j,i+1)-vz(j,i))+a1*(vz(j,i+2)-vz(j,i-1));
                    vxz=vx(j+1,i)-vx(j,i);
                else
                vxx=a0*(vx(j,i)-vx(j,i-1))+a1*(vx(j,i+1)-vx(j,i-2));
                vzz=a0*(vz(j,i)-vz(j-1,i))+a1*(vz(j+1,i)-vz(j-2,i));
                vzx=a0*(vz(j,i+1)-vz(j,i))+a1*(vz(j,i+2)-vz(j,i-1));
                vxz=a0*(vx(j+1,i)-vx(j,i))+a1*(vx(j+2,i)-vx(j-1,i));
                end
                mem_vxx(j,i)=b_x(i)*mem_vxx(j,i)+a_x(i)*vxx;
                vxx=vxx/K_x(i)+mem_vxx(j,i);
                
                mem_vzx(j,i)=b_x_half(i)*mem_vzx(j,i)+a_x_half(i)*vzx;
                vzx=vzx/K_x_half(i)+mem_vzx(j,i);
                
                mem_vzz(j,i)=b_z(j)*mem_vzz(j,i)+a_z(j)*vzz;
                vzz=vzz/K_z(j)+mem_vzz(j,i);
                
                mem_vxz(j,i)=b_z_half(j)*mem_vxz(j,i)+a_z_half(j)*vxz;
                vxz=vxz/K_z_half(j)+mem_vxz(j,i);
                
                sxx(j, i) = sxx(j, i) + dtdx * ( lamb(j,i) * (vxx + vzz) + 2.0 * mu(j,i) * vxx );
                szz(j, i) = szz(j, i) + dtdx * ( lamb(j,i) * (vxx + vzz) + 2.0 * mu(j,i) * vzz );
                sxz(j, i) = sxz(j, i) + dtdx * (  muhxhz(j,i) * (vzx + vxz) );
                
            end
        end
        free_tau;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %PlANE WAVE SOURCE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(plane_wave==true && INJECTIONS_ORDER==2)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %LEFT SIDE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        i=ibm;
        for j=jbm:jem-1
            vx(j,i)=vx(j,i)+dtdx/rhohx(j,i)*(-txx_lm(it,j-jbm+1));
            vz(j,i)=vz(j,i)+dtdx/rhohz(j,i)*(-txz_lm(it,j-jbm+1));
            szz(j,i)=szz(j,i)+dtdx*lamb(j,i)*(-vx_lm(it,j-jbm+1));
            sxx(j,i)=sxx(j,i)+dtdx*(lamb(j,i)+2*mu(j,i))*(-vx_lm(it,j-jbm+1));
            sxz(j,i)=sxz(j,i)+dtdx*muhxhz(j,i)*(-vz_lm(it,j-jbm+1));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %BOTTOM SIDE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         j=jem;
         for i=ibm+1:iem-1
            vx(j,i)=vx(j,i)+dtdx/rhohx(j,i)*(+txz_bm(it,i-ibm));
             vz(j,i)=vz(j,i)+dtdx/rhohz(j,i)*(+tzz_bm(it,i-ibm));
             sxx(j,i)=sxx(j,i)+dtdx*lamb(j,i)*(+vz_bm(it,i-ibm));
             szz(j,i)=szz(j,i)+dtdx*(lamb(j,i)+2*mu(j,i))*(+vz_bm(it,i-ibm));
             sxz(j,i)=sxz(j,i)+dtdx*muhxhz(j,i)*(+vx_bm(it,i-ibm));
         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %RIGHT SIDE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        i=iem;
        for j=jbm:jem-1
            vx(j,i)=vx(j,i)+dtdx/rhohx(j,i)*(+txx_rm(it,j-jbm+1));
            vz(j,i)=vz(j,i)+dtdx/rhohz(j,i)*(+txz_rm(it,j-jbm+1));
            szz(j,i)=szz(j,i)+dtdx*lamb(j,i)*(+vx_rm(it,j-jbm+1));
            sxx(j,i)=sxx(j,i)+dtdx*(lamb(j,i)+2*mu(j,i))*(+vx_rm(it,j-jbm+1));
            sxz(j,i)=sxz(j,i)+dtdx*muhxhz(j,i)*(+vz_rm(it,j-jbm+1));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %CORNERS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %----BOTTOM LEFT CORNER----%
        j=jem;
        i=ibm;
        szz(jem,ibm)=szz(jem,ibm)+dtdx*lamb(jem,ibm)*(-vx_lm(it,end));
        sxx(jem,ibm)=sxx(jem,ibm)+dtdx*(lamb(jem,ibm)+2*mu(jem,ibm))*(-vx_lm(it,end));
        sxz(jem,ibm)=sxz(jem,ibm)+dtdx*muhxhz(jem,ibm)*(+vx_lm(it,end));
        vx(jem,ibm)=vx(jem,ibm)+dtdx/rhohx(jem,ibm)*(+txz_lm(it,end));
       vx(jem,ibm)=vx(jem,ibm)+dtdx/rhohx(jem,ibm)*(-txx_lm(it,end));
         %----BOTTOM RIGHT CORNER----%
        j=jem;
        i=iem;
        sxx(jem,iem)=sxx(jem,iem)+dtdx*lamb(jem,iem)*(+vz_rm(it,end));
        szz(jem,iem)=szz(jem,iem)+dtdx*(lamb(jem,iem)+2*mu(jem,iem))*(+vz_rm(it,end));
        szz(jem,iem)=szz(jem,iem)+dtdx*lamb(jem,iem)*(+vx_rm(it,end));
        sxx(jem,iem)=sxx(jem,iem)+dtdx*(lamb(jem,iem)+2*mu(jem,iem))*(+vx_rm(it,end));
        vx(jem,iem)=vx(jem,iem)+dtdx/rhohx(jem,iem)*(+txx_rm(it,end));
        vz(jem,iem)=vz(jem,iem)+dtdx/rhohz(jem,iem)*(+tzz_rm(it,end));
    end
    
    if (plane_wave==true && INJECTIONS_ORDER==4)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %LEFT SIDE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Reference j is at middle layer
         for j=jbm:jem-1
            %middle
%             vx_li=
%             vz_li=
%             txx_li=
%             tzz_li=
%             txz_li=
    
            vx(j,ibm)=vx(j,ibm)+dtdx/rhohx(j,ibm)*(-a0*txx_li(it,j-jbm+1)-a1*txx_lo(it,j-jbm+1));
            vz(j,ibm)=vz(j,ibm)+dtdx/rhohz(j,ibm)*(-a0*txz_li(it,j-jbm+1)-a1*txz_li(it,j-jbm+1));
            sxx(j,ibm)=sxx(j,ibm)+dtdx*(lamb(j,ibm)+2*mu(j,ibm))*(-a0*vx_li(it,j-jbm+1)-a1*vx_li(it,j-jbm+1));
            szz(j,ibm)=szz(j,ibm)+dtdx*lamb(j,ibm)*(-a0*vx_li(it,j-jbm+1)-a1*vx_li(it,j-jbm+1));
            sxz(j,ibm)=sxz(j,ibm)+dtdx*muhxhz(j,ibm)*(-a0*vz_li(it,j-jbm+1)-a1*vz_lo(it,j-jbm+1));
            %inner
            vx(j,ibi)=vx(j,ibi)+dtdx/rhohx(j,ibi)*(-a1*txx_li(it,j-jbm+1));
            vz(j,ibi)=vz(j,ibi)+dtdx/rhohz(j,ibi)*(-a1*txz_lo(it,j-jbm+1));
            sxx(j,ibi)=sxx(j,ibi)+dtdx*(lamb(j,ibi)+2*mu(j,ibi))*(-a1*vx_lo(it,j-jbm+1));
            szz(j,ibi)=szz(j,ibi)+dtdx*lamb(j,ibi)*(-a1*vx_lo(it,j-jbm+1));
            sxz(j,ibi)=sxz(j,ibi)+dtdx*muhxhz(j,ibi)*(-a1*vz_li(it,j-jbm+1));
            %outter
            vx(j,ibo)=vx(j,ibo)+dtdx/rhohx(j,ibo)*(-a1*txx_li(it,j-jbm+1));
            vz(j,ibo)=vz(j,ibo)+dtdx/rhohz(j,ibo)*(-a1*txz_li(it,j-jbm+1));
            sxx(j,ibo)=sxx(j,ibo)+dtdx*(lamb(j,ibo)+2*mu(j,ibo))*(-a1*vx_li(it,j-jbm+1));
            szz(j,ibo)=szz(j,ibo)+dtdx*lamb(j,ibo)*(-a1*vx_li(it,j-jbm+1));
            sxz(j,ibo)=sxz(j,ibo)+dtdx*muhxhz(j,ibo)*(-a1*vz_li(it,j-jbm+1));
         end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %RIGHT SIDE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Reference j corresponds to middle layer
        for j=jbm:jem-1
            %middle
            vx(j,iem)=vx(j,iem)+dtdx/rhohx(j,iem)*(+a0*txx_rm(it,j-jbm+1)+a1*txx_ri(it,j-jbm+1));
            vz(j,iem)=vz(j,iem)+dtdx/rhohz(j,iem)*(+a0*txz_rm(it,j-jbm+1)+a1*txz_ro(it,j-jbm+1));
            sxx(j,iem)=sxx(j,iem)+dtdx*(lamb(j,iem)+2*mu(j,iem))*(+a0*vx_rm(it,j-jbm+1)+a1*vx_ro(it,j-jbm+1));
            szz(j,iem)=szz(j,iem)+dtdx*lamb(j,iem)*(+a0*vx_rm(it,j-jbm+1)+a1*vx_ro(it,j-jbm+1));
            sxz(j,iem)=sxz(j,iem)+dtdx*muhxhz(j,iem)*(+a0*vz_rm(it,j-jbm+1)+a1*vz_ri(it,j-jbm+1));
            %inner
            vx(j,iei)=vx(j,iei)+dtdx/rhohx(j,iei)*(+a1*txx_rm(it,j-jbm+1));
            vz(j,iei)=vz(j,iei)+dtdx/rhohz(j,iei)*(+a1*txz_rm(it,j-jbm+1));
            sxx(j,iei)=sxx(j,iei)+dtdx*(lamb(j,iei)+2*mu(j,iei))*(+a1*vx_rm(it,j-jbm+1));
            szz(j,iei)=szz(j,iei)+dtdx*lamb(j,iei)*(+a1*vx_rm(it,j-jbm+1));
            sxz(j,iei)=sxz(j,iei)+dtdx*muhxhz(j,iei)*(+a1*vz_ro(it,j-jbm+1));
            %outter
            vx(j,ieo)=vx(j,ieo)+dtdx/rhohx(j,ieo)*(+a1*txx_rm(it,j-jbm+1));
            vz(j,ieo)=vz(j,ieo)+dtdx/rhohz(j,ieo)*(+a1*txz_ri(it,j-jbm+1));
            sxx(j,ieo)=sxx(j,ieo)+dtdx*(lamb(j,ieo)+2*mu(j,ieo))*(+a1*vx_ri(it,j-jbm+1));
            szz(j,ieo)=szz(j,ieo)+dtdx*lamb(j,ieo)*(+a1*vx_ri(it,j-jbm+1));
            sxz(j,ieo)=sxz(j,ieo)+dtdx*muhxhz(j,ieo)*(+a1*vz_rm(it,j-jbm+1));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %BOTTOM SIDE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         for i=ibi:iei
            %middle
            vx(jem,i)=vx(jem,i)+dtdx/rhohx(jem,i)*(+a0*txz_bm(it,i-ibi+3)+a1*txz_bo(it,i-ibi+3));
            vz(jem,i)=vz(jem,i)+dtdx/rhohz(jem,i)*(+a0*tzz_bm(it,i-ibi+3)+a1*tzz_bm(it,i-ibi+3));
            sxx(jem,i)=sxx(jem,i)+dtdx*lamb(jem,i)*(+a0*vz_bm(it,i-ibi+3)+a1*vz_bo(it,i-ibi+3));
            szz(jem,i)=szz(jem,i)+dtdx*(lamb(jem,i)+2*mu(j,iem))*(+a0*vz_bm(it,i-ibi+3)+a1*vz_bo(it,i-ibi+3));
            sxz(jem,i)=sxz(jem,i)+dtdx*muhxhz(jem,i)*(+a0*vx_bm(it,i-ibi+3)+a1*vx_bm(it,i-ibi+3));
            %inner
            vx(jei,i)=vx(jei,i)+dtdx/rhohx(jei,i)*(+a1*txz_bm(it,i-ibi+3));
            vz(jei,i)=vz(jei,i)+dtdx/rhohz(jei,i)*(+a1*tzz_bo(it,i-ibi+3));
            sxx(jei,i)=sxx(jei,i)+dtdx*lamb(jei,i)*(+a1*vz_bm(it,i-ibi+3));
            szz(jei,i)=szz(jei,i)+dtdx*(lamb(jei,i)+2*mu(j,iem))*(+a1*vz_bm(it,i-ibi+3));
            sxz(jei,i)=sxz(jei,i)+dtdx*muhxhz(jei,i)*(+a1*vx_bo(it,i-ibi+3));
            %outter
            vx(jeo,i)=vx(jeo,i)+dtdx/rhohx(jeo,i)*(+a1*txz_bm(it,i-ibi+3));
            vz(jeo,i)=vz(jeo,i)+dtdx/rhohz(jeo,i)*(+a1*tzz_bm(it,i-ibi+3));
            sxx(jeo,i)=sxx(jeo,i)+dtdx*lamb(jeo,i)*(+a1*vz_bm(it,i-ibi+3)); 
            szz(jeo,i)=szz(jeo,i)+dtdx*(lamb(jeo,i)+2*mu(j,iem))*(+a1*vz_bm(it,i-ibi+3)); 
            sxz(jeo,i)=sxz(jeo,i)+dtdx*muhxhz(jeo,i)*(+a1*vx_bm(it,i-ibi+3)); 
         end
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %CORNERS
         %There are 14 (+2 because "ns" are in the same node) corner points 
         %to correct :S 
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %----BOTTOM LEFT CORNERS----% %A lot of points :S 
         j=jem; 
         i=ibi; 
         vx(j,i)=vx(j,i)+dtdx/rhohx(j,i)*(-a1*txx_li(it,end-1));
         sxx(j,i)=sxx(j,i)+dtdx*(lamb(j,i)+2*mu(j,i))*(-a1*vx_li(it,end-1));
         szz(j,i)=szz(j,i)+dtdx*lamb(j,i)*(-a1*vx_li(it,end-1));
         
         j=jei; 
         i=ibm; 
         vx(j,i)=vx(j,i)+dtdx/rhohx(j,i)*(+a1*txz_bm(it,2));
         sxz(j,i)=sxz(j,i)+dtdx*muhxhz(j,i)*(+a1*vx_bo(it,2));
         
         j=jem; 
         i=ibi; 
         vx(j,i)=vx(j,i)+dtdx/rhohx(j,i)*(-a1*txx_li(it,end-1));
         sxx(j,i)=sxx(j,i)+dtdx*(lamb(j,i)+2*mu(j,i))*(-a1*vx_lo(it,end-1));
         szz(j,i)=szz(j,i)+dtdx*lamb(j,i)*(-a1*vx_lo(it,end-1));
       
         j=jeo; 
         i=ibm; 
         vx(j,i)=vx(j,i)+dtdx/rhohx(j,i)*(+a1*txz_bm(it,2));
         sxz(j,i)=sxz(j,i)+dtdx*muhxhz(j,i)*(+a1*vx_bm(it,2));
         
         j=jem;
         i=ibm; 
         
         sxx(j,i)=sxx(j,i)+dtdx*(lamb(j,i)+2*mu(j,i))*(-a0*vx_li(it,end-1)-a1*vx_li(it,end-1));
         szz(j,i)=szz(j,i)+dtdx*lamb(j,i)*(-a0*vx_li(it,end-1)-a1*vx_li(it,end-1));
         sxz(j,i)=sxz(j,i)+dtdx*muhxhz(j,i)*(+a0*vx_bm(it,2)+a1*vx_bm(it,2));
         vx(j,i)=vx(j,i)+dtdx/rhohx(j,i)*(-a0*txx_li(it,end-1)-a1*txx_lo(it,end-1));
         vx(j,i)=vx(j,i)+dtdx/rhohx(j,i)*(+a0*txz_bm(it,2)+a1*txz_bo(it,2));
       
         
         %----BOTTOM RIGHT CORNERS----% %A lot of points :S 
         j=jem; 
         i=iei; 
         sxx(j,i)=sxx(j,i)+dtdx*(lamb(j,i)+2*mu(j,i))*(+a1*vx_bm(it,end-1)); 
         szz(j,i)=szz(j,i)+dtdx*lamb(j,i)*(+a1*vx_bm(it,end-1));
         vx(j,i)=vx(j,i)+dtdx/rhohx(j,i)*(+a1*txx_bm(it,end));

         j=jei; 
         i=iem; 
         sxx(j,i)=sxx(j,i)+dtdx*lamb(j,i)*(+a1*vz_bm(it,end-1)); 
         szz(j,i)=szz(j,i)+dtdx*(lamb(j,i)+2*mu(j,i))*(+a1*vz_bm(it,end-1)); 
         vz(j,i)=vz(j,i)+dtdx/rhohz(j,i)*(+a1*tzz_bo(it,end-1)); 
         
         j=jem; 
         i=ieo; 
         vx(j,i)=vx(j,i)+dtdx/rhohx(j,i)*(+a1*txx_bm(it,end-1)); 
         sxx(j,i)=sxx(j,i)+dtdx*(lamb(j,i)+2*mu(j,i))*(+a1*vx_bm(it,end-2)); 
         szz(j,i)=szz(j,i)+dtdx*lamb(j,i)*(+a1*vx_bm(it,end-2)); 
         
         j=jeo; 
         i=iem; 
         vz(j,i)=vz(j,i)+dtdx/rhohz(j,i)*(+a1*tzz_bm(it,end-1)); 
         sxx(j,i)=sxx(j,i)+dtdx*lamb(j,i)*(+a1*vz_bm(it,end-1)); 
         szz(j,i)=szz(j,i)+dtdx*(lamb(j,i)+2*mu(j,iem))*(+a1*vz_bm(it,end-1)); 
         
         j=jem; 
         i=iem; 
         vx(j,i)=vx(j,i)+dtdx/rhohx(j,i)*(+a0*txx_bm(it,end-1)+a1*txx_bm(it,end-1));
         vz(j,i)=vz(j,i)+dtdx/rhohz(j,i)*(+a0*tzz_bm(it,end-1)+a1*tzz_bm(it,end-1));
         sxx(j,i)=sxx(j,i)+dtdx*lamb(j,i)*(+a0*vz_bm(it,end-1)+a1*vz_bo(it,end-1));
         szz(j,i)=szz(j,i)+dtdx*(lamb(j,i)+2*mu(j,i))*(+a0*vz_bm(it,end-1)+a1*vz_bo(it,end-1));
         sxx(j,i)=sxx(j,i)+dtdx*(lamb(j,i)+2*mu(j,i))*(+a0*vx_rm(it,end-1)+a1*vx_ro(it,end-1));
         szz(j,i)=szz(j,i)+dtdx*lamb(j,i)*(+a0*vx_rm(it,end-1)+a1*vx_ro(it,end-1));
          
    end
    
    %Compute displacement field via numerical integration
    ux=ux+vx.*dt;
    uz=uz+vz.*dt;
    
    %Compute Kaelin and guiton (200...) illumination factor
    preckg= (vz + vx).^2 + preckg;
    
    %Compute pseudo-hessian matrix (Shin et al. 2001)
    vxt=time_der(dt,vx0,vx);
    vx0=vx; 
    vzt=time_der(dt,vz0,vz);
    vz0=vz;
    sxxt=time_der(dt,sxx0,sxx);
    sxx0=sxx;
    szzt=time_der(dt,szz0,szz);
    szz0=szz;
    sxzt=time_der(dt,sxz0,sxz);
    sxz0=sxz;
    
    %pseudo hessian
    Hrho=Hrho+(vxt.^2+vzt.^2)./rho.^2;
    Hlam=Hlam+(sxxt+szzt).^2./(2*(lamb+mu).^2);
    Hmu=Hmu+(sxxt+szzt).^2./(2*(lamb+mu).^2)+(sxxt-szzt).^2./(2*mu.^2)+sxzt.^2./mu.^2;
    %Hlam=zeros(Nz,Nx);
    %Hmu=zeros(Nz,Nx);
    
    
    if (mod(it,sfe)==0 && store==1)
        %Storing time reversed forward velocity and displacement fields
        vx_forward(nt/sfe+1-it/sfe,:,:)=vx(:,:);
        vz_forward(nt/sfe+1-it/sfe,:,:)=vz(:,:);
        % displacement
        ux_forward(nt/sfe+1-it/sfe,:,:)=ux(:,:);
        uz_forward(nt/sfe+1-it/sfe,:,:)=uz(:,:);
    end
    
    
%     totv=sqrt(vx.^2+vz.^2);
%     totu=sqrt(ux.^2+uz.^2);
    if (mod(it,IT_DISPLAY) == 0 && animation==1)
        %fprintf('Time step: %d \t %.4f s\n',it, it*dt);
        %u=sqrt(ux3.^2 + uz3.^2);
        figure(2)
        %imagesc(vx(1:end,1:end))
        imagesc([0:nx+1]*dh/1000, [0:nz+1]*dh/1000, vz(1:end,:)); colorbar; colormap jet;
        axis equal tight;
        hold on, plot(xr*dh/1000,zr*dh/1000,'wo')
        %imagesc(totu(50:end,:));
        colorbar; colormap jet; axis equal tight
        title(['Stepvx = ',num2str(it),' Time: ',sprintf('%.4f',it*dt),' sec']);
        %xlabel('x[km]'); ylabel('z[km]'); axis equal tight;
        drawnow;
        if(makemovie==1)
            frame = getframe(gcf);
            writeVideo(obj,frame);
        end
    end
    
    %Make observed seismograms
    if (makeobs==1)
        for jr=1:nrec
            vx_obs(jr,it)=vx(zr(jr),xr(jr));
            ux_obs(jr,it)=ux(zr(jr),xr(jr));
            vz_obs(jr,it)=vz(zr(jr),xr(jr));
            uz_obs(jr,it)=uz(zr(jr),xr(jr));
        end
    end
    
    %Seismic records
    for jr=1:nrec
        vx_syn(jr,it)=vx(zr(jr),xr(jr));
        ux_syn(jr,it)=ux(zr(jr),xr(jr));
        vz_syn(jr,it)=vz(zr(jr),xr(jr));
        uz_syn(jr,it)=uz(zr(jr),xr(jr));
    end
    
    
end %end of loop time

%Hrho=Hrho*(1/dt^2);
%Hlam=Hlam*(1/dt^2);
%Hmu=Hmu*(1/dt^2);

prec.kg=preckg; 
prec.rho=Hrho; 
prec.lam=Hlam; 
prec.mu=Hmu; 

%Last forward field frames
Fvx=vx;
Fvz=vz;
Fux=ux;
Fuz=uz;
Fsxx=sxx;
Fszz=szz;
Fsxz=sxz;
%Rearranging arrays
BD_vx=[seis_vx_TOP,seis_vx_RIG,seis_vx_BOT,seis_vx_LEF];
BD_vz=[seis_vz_TOP,seis_vz_RIG,seis_vz_BOT,seis_vz_LEF];
BD_sxx=[seis_sxx_TOP,seis_sxx_RIG,seis_sxx_BOT,seis_sxx_LEF];
BD_szz=[seis_szz_TOP,seis_szz_RIG,seis_szz_BOT,seis_szz_LEF];
BD_sxz=[seis_sxz_TOP,seis_sxz_RIG,seis_sxz_BOT,seis_sxz_LEF];
XS=[ xst_TOP , xst_RIGHT , xst_BOTTOM , xst_LEFT ];
ZS=[ zst_TOP , zst_RIGHT , zst_BOTTOM , zst_LEFT];
NSP=size(BD_vx,2);



%Saving seismograms
if (saveseis==1 && makeobs==1 && plane_wave==0)
    nameux=strcat('Observed_data/ux_shot',num2str(jf),'_iter',num2str(iter),'_',num2str(f0(kf)),'hz.mat');
    nameuz=strcat('Observed_data/uz_shot',num2str(jf),'_iter',num2str(iter),'_',num2str(f0(kf)),'hz.mat');
    namevx=strcat('Observed_data/vx_shot',num2str(jf),'_iter',num2str(iter),'_',num2str(f0(kf)),'hz.mat');
    namevz=strcat('Observed_data/vz_shot',num2str(jf),'_iter',num2str(iter),'_',num2str(f0(kf)),'hz.mat');
    
    save(namevx,'vx_obs');
    save(namevz,'vz_obs');
    save(nameux,'ux_obs');
    save(nameuz,'uz_obs');
end
if (saveseis==1 && makeobs==1 && STF_TYPE(jf)==5)
    nameux=strcat('Observed_data/ux_green_shot',num2str(jf),'_iter',num2str(iter),'_',num2str(f0(kf)),'hz.mat');
    nameuz=strcat('Observed_data/uz_green_shot',num2str(jf),'_iter',num2str(iter),'_',num2str(f0(kf)),'hz.mat');
    namevx=strcat('Observed_data/vx_green_shot',num2str(jf),'_iter',num2str(iter),'_',num2str(f0(kf)),'hz.mat');
    namevz=strcat('Observed_data/vz_green_shot',num2str(jf),'_iter',num2str(iter),'_',num2str(f0(kf)),'hz.mat');
    
    save(namevx,'vx_obs');
    save(namevz,'vz_obs');
    save(nameux,'ux_obs');
    save(nameuz,'uz_obs');
end
if (saveseis==1 && makeobs==0 && plane_wave==0)
    nameux=strcat('Synthetic_data/ux_shot',num2str(jf),'_iter',num2str(iter),'_',num2str(f0(kf)),'hz.mat');
    nameuz=strcat('Synthetic_data/uz_shot',num2str(jf),'_iter',num2str(iter),'_',num2str(f0(kf)),'hz.mat');
    namevx=strcat('Synthetic_data/vx_shot',num2str(jf),'_iter',num2str(iter),'_',num2str(f0(kf)),'hz.mat');
    namevz=strcat('Synthetic_data/vz_shot',num2str(jf),'_iter',num2str(iter),'_',num2str(f0(kf)),'hz.mat');
    
    save(namevx,'vx_syn');
    save(namevz,'vz_syn');
    save(nameux,'ux_syn');
    save(nameuz,'uz_syn');
end

if (saveseis==1 && makeobs==0 && plane_wave==1)
    nameux=strcat('Inversion_results/synthetic_data/ux_shot',num2str(jf),'_iter',num2str(iter),'_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat');
    nameuz=strcat('Inversion_results/synthetic_data/uz_shot',num2str(jf),'_iter',num2str(iter),'_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat');
    namevx=strcat('Inversion_results/synthetic_data/vx_shot',num2str(jf),'_iter',num2str(iter),'_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat');
    namevz=strcat('Inversion_results/synthetic_data/vz_shot',num2str(jf),'_iter',num2str(iter),'_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat');
    
    save(namevx,'vx_syn');
    save(namevz,'vz_syn');
    save(nameux,'ux_syn');
    save(nameuz,'uz_syn');
end

if (saveseis==1 && makeobs==1 && plane_wave==1)
    nameux=strcat('Observed_data/ux_shot',num2str(jf),'_iter',num2str(iter),'_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat');
    nameuz=strcat('Observed_data/uz_shot',num2str(jf),'_iter',num2str(iter),'_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat');
    namevx=strcat('Observed_data/vx_shot',num2str(jf),'_iter',num2str(iter),'_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat');
    namevz=strcat('Observed_data/vz_shot',num2str(jf),'_iter',num2str(iter),'_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat');
    
    save(namevx,'vx_obs');
    save(namevz,'vz_obs');
    save(nameux,'ux_obs');
    save(nameuz,'uz_obs');
end

if(makemovie==1)
    obj.close()
end


end %end of function
