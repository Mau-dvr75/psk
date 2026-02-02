clear all; clc; 
input_parameters2;
tic
make_fk_model; %make true model
%parpool(2); 
if (makeobs==1)
    maxiter=1;
    ekf=length(f0);
    vp=vp_true;
    vs=vs_true;
    rho=rho_true;
end

check_stability;

if(plane_wave==true && compute_initial_field==true)
    %add cpml nodes to models
    if STF_TYPE(1)==4 || STF_TYPE(1)==1 || STF_TYPE(1)==3 
        ekf=length(f0)-1; %because bandpass filtering.
    end
    vp=extend_model(vp,nx,nz,ncpml);
    vs=extend_model(vs,nx,nz,ncpml);
    rho=extend_model(rho,nx,nz,ncpml);
    FK_fields;
    vp=vp(nz0:nzf,nx0:nxf);
    vs=vs(nz0:nzf,nx0:nxf);
    rho=rho(nz0:nzf,nx0:nxf);
end
if plane_wave
    if STF_TYPE(1)==4
        ekf=length(f0)-1; %because bandpass filtering
    end
end
%Loop over frequencies
for kf=1:ekf
    display((sprintf('Working with data filtered from %f to %f [Hz]',f0(kf),f0(kf+1))))
    %make corresponding stf
    stf=make_stf(nt,dt,f0(kf),1,3,0); %dummy
    %Compute cpml arrays once and for all for each frequency
    [a_x,a_x_half,b_x,b_x_half,K_x,K_x_half,a_z,a_z_half,b_z,b_z_half,K_z,K_z_half]=cpml(ncpml,dh,Nx,Nz,(f0(kf)+f0(kf+1))/2,dt,vp_max);
    %loop over iterations 
    for iter=1:maxiter
        if makeobs==0
        display((sprintf('FWI iteration %i of %i',iter,maxiter)))
        end
        initialize_FWI_fields;
        %add cpml nodes to models
        vp=extend_model(vp,nx,nz,ncpml);
        vs=extend_model(vs,nx,nz,ncpml);
        rho=extend_model(rho,nx,nz,ncpml);
        %loop over sources (can be paralelized)
        for jf=1:1
            if(plane_wave==true) %load fk injection fields
                if stf_method(jf)==1
                    method='average';
                elseif stf_method(jf)==2
                    method='pratt';
                end
                if STF_TYPE(jf)==4
                    namestf=['./data_korea/stations/arrays/',num2str(noa),'/STF/',num2str(f0(ekf)),'-',num2str(f0(ekf+1)),'Hz_event',num2str(jf),catalog,'_dt',num2str(0.01),method,'.mat'];%stf is at 0.01 samples
                    if exist(namestf,'file')==0
                        continue %if no source, skip
                    end
                end
                if (INJECTIONS_ORDER==2)
                    fk_bi=[];fk_bo=[];
                    fk_ri=[];fk_ro=[];
                    fk_lo=[];fk_li=[];
                    fk_lm=load(strcat('Plane_sources/',num2str(angleforced(jf)),'_','fk_lm','_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat'));
                    fk_rm=load(strcat('Plane_sources/',num2str(angleforced(jf)),'_','fk_rm','_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat'));
                    fk_bm=load(strcat('Plane_sources/',num2str(angleforced(jf)),'_','fk_bm','_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat'));
                elseif (INJECTIONS_ORDER==4)
                    fk_lo=load(strcat('Plane_sources/',num2str(angleforced(jf)),'_','fk_lo','_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat'));
                    fk_lm=load(strcat('Plane_sources/',num2str(angleforced(jf)),'_','fk_lm','_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat'));
                    fk_ro=load(strcat('Plane_sources/',num2str(angleforced(jf)),'_','fk_ro','_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat'));
                    fk_rm=load(strcat('Plane_sources/',num2str(angleforced(jf)),'_','fk_rm','_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat'));
                    fk_ri=load(strcat('Plane_sources/',num2str(angleforced(jf)),'_','fk_ri','_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat'));
                    fk_bo=load(strcat('Plane_sources/',num2str(angleforced(jf)),'_','fk_bo','_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat'));
                    fk_bm=load(strcat('Plane_sources/',num2str(angleforced(jf)),'_','fk_bm','_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat'));
                    fk_bi=load(strcat('Plane_sources/',num2str(angleforced(jf)),'_','fk_bi','_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat'));
                end
            else %set to 0 FK fields (dummy values)
                fk_bi=[];fk_bm=[];fk_bo=[];
                fk_ri=[];fk_rm=[];fk_ro=[];
                fk_lo=[];fk_lm=[];fk_li=[];
            end
            %Solve the forward problem
            [vx_syn,vz_syn,ux_syn,uz_syn,vx_fw,vz_fw,ux_fw,uz_fw,Fvx,Fux,Fvz,Fuz,Fsxx,Fszz,Fsxz,BD_vx,BD_vz,BD_sxx,BD_szz,BD_sxz,XS,ZS,NSP,prec]=psv(vp,vs,rho,nx,nz,dh,dt,nt,stf,source_x_id(jf),source_z_id(jf),nsrc,rec_x_id,rec_z_id,nrec,FREE_SURF,saveseis,IT_DISPLAY,sfe,makeobs,jf,0,animation,kf,a_x,a_x_half,b_x,b_x_half,K_x,K_x_half,a_z,a_z_half,b_z,b_z_half,K_z,K_z_half,ncpml,0,op,field_order,fk_lo,fk_lm,fk_li,fk_ro,fk_rm,fk_ri,fk_bo,fk_bm,fk_bi,plane_wave,INJECTIONS_ORDER,makemovie,f0,STF_TYPE,iter);
        end
    end %end of iterations
    
    if (makeobs==1)
        vp=vp(nz0:nzf,nx0:nxf);
        vs=vs(nz0:nzf,nx0:nxf);
        rho=rho(nz0:nzf,nx0:nxf);
    end
    
end  %end of stages

if (analytical==true)
    analyticalFK;
end
toc