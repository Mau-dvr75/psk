for kf=1:ekf
    for jf=1:1
        %jf
        %ne=length(SOURCE_TYPE);
        if stf_method(jf)==1
            method='average';
        elseif stf_method(jf)==2
            method='pratt';
        end
        if STF_TYPE(jf)==4
        %namestf=['./data_mase/STF/STF_',num2str(f0(ekf)),'-',num2str(f0(ekf+1)),'Hz_event',num2str(jf),catalog,'_dt',num2str(dt),method,'.mat'];
        %namestf=['./data_mase/STF/STF_nf_event',num2str(jf),catalog,'_dt',num2str(dt),'average_amp.mat'];
        namestf=['./data_korea/stations/arrays/',num2str(noa),'/STF/',num2str(f0(ekf)),'-',num2str(f0(ekf+1)),'Hz_event',num2str(jf),catalog,'_dt',num2str(0.01),method,'.mat'];%stf is at 0.01 samples

        if exist(namestf)==0
             continue;
        end
        end
%             return
%         else
%             continue
%         end
        for i=1:Nx
            for j=1:Nz
                mu(j,i)=rho(j,i)*vs(j,i)*vs(j,i);
                lamb(j,i)=rho(j,i)*vp(j,i)*vp(j,i)-2.*mu(j,i);
            end
        end
        
        %for jf=1:ne
        %Now in radians
        angleforce=angleforced(jf)*pi/180;
        %Ray parameter of lower layer (half-space)
        if (SOURCE_TYPE(jf)==1)
            p=sin(angleforce)/al_fk(nlayer_fk);
        end
        if (SOURCE_TYPE(jf)==2)
            p=sin(angleforce)/be_fk(nlayer_fk);
        end
        HDUR=1./f0; %Pulses' duration (if predifined STF is used)
        tstart=0;
        stag=true;
        %deltat=dt;
        %%%%%%%%%%%%%%
        if resample_order
        nstepfk=round(nt/resample_order);
        dt_fk=dt*resample_order; %for FK simulation 
        else
           nstepfk=nt; 
           dt_fk=dt;
        end
        %nstepfk = nstep;%max(1+round(50*f0/deltat),nstep);
        t=0:dt:dt*(nt-1);
        tfk=0:dt_fk:dt_fk*(nstepfk-1);
        %Quick check if wave comes from left
        if (angleforced(jf)>0 && angleforced(jf)<=90)
            %Just a quick test to see if analytical FK field is ok
            ib=ibm;
            ie=ibm;
            jb=jem;
            je=jem;
            flag1="lc";
            flag2="all";%ns=normal stresses, ss=shear streses; vx=vx, vz=vz; all=all
            [vxu,vzu,txxu,tzzu,txzu]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
            ib=rec_x_id(1);
            ie=ib; 
            jb=rec_z_id(1);
            je=jb; 
           [vxt,vzt,txxt,tzzt,txzt]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
            %itmod=finddelay(vxu,vxt);
            %d3=sqrt(165000^2+((rec_x_id(1)-ibm)*dh)^2);
            
            
            %temp=data(i,abs(round(teodelayid(i)))+1:end,ic);
%             vxu=[zeros(1,abs(itmod)) vxu'];
%             vxu=vxu(1:nstepfk); 
            %data_alligned(i,:,ic)=temp;
        end
        %Quick check if wave comes from right
        if (angleforced(jf)>270)
            %Just a quick test to see if analytical FK field is ok
            ib=iem;
            ie=iem;
            jb=jem;
            je=jem;
            flag1="rc";
            flag2="all";%ns=normal stresses, ss=shear streses; vx=vx, vz=vz; all=all
            [vxu,vzu,txxu,tzzu,txzu]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
            ib=rec_x_id(end);
            ie=ib;
            jb=rec_z_id(end);
            je=jb;
            [vxt,vzt,txxt,tzzt,txzt]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);

        end
        
        if STF_TYPE(jf)==1 || STF_TYPE(jf)==2 || STF_TYPE(jf)==3
        display((sprintf('Creating initial wavefield for plane wave %i of %i at freq %f',jf,nsrc,f0(kf))))
        elseif STF_TYPE(jf)==4 
            display((sprintf('Creating initial wavefield for plane wave %i of %i, filtering from %f %f [Hz]',jf,nsrc,f0(kf),f0(kf+1))))
        elseif STF_TYPE(jf)==5
            display((sprintf('Creating Green function wavefield for plane wave %i of %i',jf,nsrc)))
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %SECOND ORDER FIELD
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(field_order==2)
            parfor k=1:12
                
                %         fk_l=[];
                %         fk_b=[];
                %         fk_r=[];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %---LEFT BOUNDARY
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(k==1)
                    ib=ibm;
                    ie=ibm;
                    jb=jbm;
                    je=jem;
                    flag1="left";
                    flag2="ns";%ns=normal stresses, ss=shear streses; vx=vx, vz=vz; all=all
                    [~,~,txx_lm,tzz_lm,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txx_lm',txx_lm);
                    parsave_fk('Plane_sources/tzz_lm',tzz_lm);
                end
                if (k==2)
                    ib=ibm;
                    ie=ibm;
                    jb=jbm;
                    je=jem;
                    flag1="left";
                    flag2="ss";
                    [~,~,~,~,txz_lm]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txz_lm',txz_lm);
                end
                if (k==3)
                    ib=ibm;
                    ie=ibm;
                    jb=jbm;
                    je=jem;
                    flag1="left";
                    flag2="vx";
                    [vx_lm,~,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vx_lm',vx_lm);
                end
                if (k==4)
                    ib=ibm;
                    ie=ibm;
                    jb=jbm;
                    je=jem;
                    flag1="left";
                    flag2="vz";
                    [~,vz_lm,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vz_lm',vz_lm);
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %---RIGHT BOUNDARY
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(k==5)
                    ib=iem;
                    ie=iem;
                    jb=jbm;
                    je=jem;
                    flag1="right";
                    flag2="ns";%ns=normal stresses, ss=shear streses; vx=vx, vz=vz; all=all
                    [~,~,txx_rm,tzz_rm,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txx_rm',txx_rm);
                    parsave_fk('Plane_sources/tzz_rm',tzz_rm);
                end
                if(k==6)
                    ib=iem;
                    ie=iem;
                    jb=jbm;
                    je=jem;
                    flag1="right";
                    flag2="ss";
                    [~,~,~,~,txz_rm]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txz_rm',txz_rm);
                end
                if(k==7)
                    ib=iem;
                    ie=iem;
                    jb=jbm;
                    je=jem;
                    flag1="right";
                    flag2="vx";
                    [vx_rm,~,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vx_rm',vx_rm);
                end
                if(k==8)
                    ib=iem;
                    ie=iem;
                    jb=jbm;
                    je=jem;
                    flag1="right";
                    flag2="vz";
                    [~,vz_rm,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vz_rm',vz_rm);
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %---BOTTOM BOUNDARY
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(k==9)
                    ib=ibm+1;
                    ie=iem-1;
                    jb=jem;
                    je=jem;
                    flag1="bottom"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="ns";
                    [~,~,txx_bm,tzz_bm,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txx_bm',txx_bm);
                    parsave_fk('Plane_sources/tzz_bm',tzz_bm);
                end
                if(k==10)
                    ib=ibm+1;
                    ie=iem-1;
                    jb=jem;
                    je=jem;
                    flag1="bottom"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="ss";
                    [~,~,~,~,txz_bm]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txz_bm',txz_bm);
                end
                if(k==11)
                    ib=ibm+1;
                    ie=iem-1;
                    jb=jem;
                    je=jem;
                    flag1="bottom"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="vx";
                    [vx_bm,~,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vx_bm',vx_bm);
                end
                if(k==12)
                    ib=ibm+1;
                    ie=iem-1;
                    jb=jem;
                    je=jem;
                    flag1="bottom"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="vz";
                    [~,vz_bm,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vz_bm',vz_bm);
                end
            end
            txx_lm=load('Plane_sources/txx_lm');
            tzz_lm=load('Plane_sources/tzz_lm');
            txz_lm=load('Plane_sources/txz_lm');
            vx_lm=load('Plane_sources/vx_lm');
            vz_lm=load('Plane_sources/vz_lm');
            if(resample_order)
                txx_lm.fk(round(nt/resample_order)+1:end,:)=[];
                txx_lm.fk=resample(txx_lm.fk,resample_order,1);
                tzz_lm.fk(round(nt/resample_order)+1:end,:)=[];
                tzz_lm.fk=resample(tzz_lm.fk,resample_order,1);
                txz_lm.fk(round(nt/resample_order)+1:end,:)=[];
                txz_lm.fk=resample(txz_lm.fk,resample_order,1);
                vx_lm.fk(round(nt/resample_order)+1:end,:)=[];
                vx_lm.fk=resample(vx_lm.fk,resample_order,1);
                vz_lm.fk(round(nt/resample_order)+1:end,:)=[];
                vz_lm.fk=resample(vz_lm.fk,resample_order,1);
            end
            fk_lm=[];
            fk_lm.vx=vx_lm;
            fk_lm.vz=vz_lm;
            fk_lm.sxx=txx_lm;
            fk_lm.szz=tzz_lm;
            fk_lm.sxz=txz_lm;
            name=strcat('Plane_sources/',num2str(angleforced(jf)),'_','fk_lm','_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat');
            save(name,'fk_lm');
            
            txx_rm=load('Plane_sources/txx_rm');
            tzz_rm=load('Plane_sources/tzz_rm');
            txz_rm=load('Plane_sources/txz_rm');
            vx_rm=load('Plane_sources/vx_rm');
            vz_rm=load('Plane_sources/vz_rm');
            if(resample_order)
                txx_rm.fk(round(nt/resample_order)+1:end,:)=[];
                txx_rm.fk=resample(txx_rm.fk,resample_order,1);
                tzz_rm.fk(round(nt/resample_order)+1:end,:)=[];
                tzz_rm.fk=resample(tzz_rm.fk,resample_order,1);
                txz_rm.fk(round(nt/resample_order)+1:end,:)=[];
                txz_rm.fk=resample(txz_rm.fk,resample_order,1);
                vx_rm.fk(round(nt/resample_order)+1:end,:)=[];
                vx_rm.fk=resample(vx_rm.fk,resample_order,1);
                vz_rm.fk(round(nt/resample_order)+1:end,:)=[];
                vz_rm.fk=resample(vz_rm.fk,resample_order,1);
            end
            fk_rm=[];
            fk_rm.vx=vx_rm;
            fk_rm.vz=vz_rm;
            fk_rm.sxx=txx_rm;
            fk_rm.szz=tzz_rm;
            fk_rm.sxz=txz_rm;
            name=strcat('Plane_sources/',num2str(angleforced(jf)),'_','fk_rm','_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat');
            save(name,'fk_rm');
            
            txx_bm=load('Plane_sources/txx_bm');
            tzz_bm=load('Plane_sources/tzz_bm');
            txz_bm=load('Plane_sources/txz_bm');
            vx_bm=load('Plane_sources/vx_bm');
            vz_bm=load('Plane_sources/vz_bm');
            if(resample_order)
                txx_bm.fk(round(nt/resample_order)+1:end,:)=[];
                txx_bm.fk=resample(txx_bm.fk,resample_order,1);
                tzz_bm.fk(round(nt/resample_order)+1:end,:)=[];
                tzz_bm.fk=resample(tzz_bm.fk,resample_order,1);
                txz_bm.fk(round(nt/resample_order)+1:end,:)=[];
                txz_bm.fk=resample(txz_bm.fk,resample_order,1);
                vx_bm.fk(round(nt/resample_order)+1:end,:)=[];
                vx_bm.fk=resample(vx_bm.fk,resample_order,1);
                vz_bm.fk(round(nt/resample_order)+1:end,:)=[];
                vz_bm.fk=resample(vz_bm.fk,resample_order,1);
            end
            fk_bm=[];
            fk_bm.vx=vx_bm;
            fk_bm.vz=vz_bm;
            fk_bm.sxx=txx_bm;
            fk_bm.szz=tzz_bm;
            fk_bm.sxz=txz_bm;
            name=strcat('Plane_sources/',num2str(angleforced(jf)),'_','fk_bm','_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat');
            save(name,'fk_bm');
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %FOURTH ORDER FIELD
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(field_order==4)
            parfor k=1:36
                %fk_lo=[];
                %fk_lm=[];
                %fk_li=[];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %---LEFT BOUNDARY
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Outter layer
                if(k==1)
                    ib=ibo;
                    ie=ibo;
                    jb=jbo;
                    je=jeo;
                    flag1="left"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="ns";%ns=normal stresses, ss=shear streses; vx=vx, vz=vz; all=all
                    [~,~,txx_lo,tzz_lo,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txx_lo',txx_lo);
                    parsave_fk('Plane_sources/tzz_lo',tzz_lo);
                end
                if(k==2)
                    ib=ibo;
                    ie=ibo;
                    jb=jbo;
                    je=jeo;
                    flag1="left";
                    flag2="ss";
                    [~,~,~,~,txz_lo]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txz_lo',txz_lo);
                end
                if(k==3)
                    ib=ibo;
                    ie=ibo;
                    jb=jbo;
                    je=jeo;
                    flag1="left";
                    flag2="vx";
                    [vx_lo,~,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vx_lo',vx_lo);
                end
                if(k==4)
                    ib=ibo;
                    ie=ibo;
                    jb=jbo;
                    je=jeo;
                    flag1="left";
                    flag2="vz";
                    [~,vz_lo,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vz_lo',vz_lo);
                end
                %Middle layer
                if(k==5)
                    ib=ibm;
                    ie=ibm;
                    jb=jbm;
                    je=jeo;
                    flag1="left";
                    flag2="ns";%ns=normal stresses, ss=shear streses; vx=vx, vz=vz; all=all
                    [~,~,txx_lm,tzz_lm,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txx_lm',txx_lm);
                    parsave_fk('Plane_sources/tzz_lm',tzz_lm);
                end
                if(k==6)
                    ib=ibm;
                    ie=ibm;
                    jb=jbm;
                    je=jeo;
                    flag1="left";
                    flag2="ss";
                    [~,~,~,~,txz_lm]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txz_lm',txz_lm);
                end
                if(k==7)
                    ib=ibm;
                    ie=ibm;
                    jb=jbm;
                    je=jeo;
                    flag1="left";
                    flag2="vx";
                    [vx_lm,~,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vx_lm',vx_lm);
                end
                if(k==8)
                    ib=ibm;
                    ie=ibm;
                    jb=jbm;
                    je=jeo;
                    flag1="left";
                    flag2="vz";
                    [~,vz_lm,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vz_lm',vz_lm);
                end
                %Inner layer
                if(k==9)
                    ib=ibi;
                    ie=ibi;
                    jb=jbi;
                    je=jeo;
                    flag1="left";
                    flag2="ns";%ns=normal stresses, ss=shear streses; vx=vx, vz=vz; all=all
                    [~,~,txx_li,tzz_li,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txx_li',txx_li);
                    parsave_fk('Plane_sources/tzz_li',tzz_li);
                end
                if(k==10)
                    ib=ibi;
                    ie=ibi;
                    jb=jbi;
                    je=jeo;
                    flag1="left";
                    flag2="ss";
                    [~,~,~,~,txz_li]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txz_li',txz_li);
                end
                if(k==11)
                    ib=ibi;
                    ie=ibi;
                    jb=jbi;
                    je=jeo;
                    flag1="left";
                    flag2="vx";
                    [vx_li,~,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vx_li',vx_li);
                end
                if(k==12)
                    ib=ibi;
                    ie=ibi;
                    jb=jbi;
                    je=jeo;
                    flag1="left";
                    flag2="vz";
                    [~,vz_li,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vz_li',vz_li);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %---RIGHT BOUNDARY
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Outter layer
                if(k==13)
                    ib=ieo;
                    ie=ieo;
                    jb=jbo;
                    je=jeo;
                    flag1="right";
                    flag2="ns";%ns=normal stresses, ss=shear streses; vx=vx, vz=vz; all=all
                    [~,~,txx_ro,tzz_ro,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txx_ro',txx_ro);
                    parsave_fk('Plane_sources/tzz_ro',tzz_ro);
                end
                if(k==14)
                    ib=ieo;
                    ie=ieo;
                    jb=jbo;
                    je=jeo;
                    flag1="right";
                    flag2="ss";
                    [~,~,~,~,txz_ro]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txz_ro',txz_ro);
                end
                if(k==15)
                    ib=ieo;
                    ie=ieo;
                    jb=jbo;
                    je=jeo;
                    flag1="right";
                    flag2="vx";
                    [vx_ro,~,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vx_ro',vx_ro);
                end
                if(k==16)
                    ib=ieo;
                    ie=ieo;
                    jb=jbo;
                    je=jeo;
                    flag1="right";
                    flag2="vz";
                    [~,vz_ro,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vz_ro',vz_ro);
                end
                %Middle layer
                if(k==17)
                    ib=iem;
                    ie=iem;
                    jb=jbm;
                    je=jeo;
                    flag1="right"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="ns";%ns=normal stresses, ss=shear streses; vx=vx, vz=vz; all=all
                    [~,~,txx_rm,tzz_rm,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txx_rm',txx_rm);
                    parsave_fk('Plane_sources/tzz_rm',tzz_rm);
                end
                if(k==18)
                    ib=iem;
                    ie=iem;
                    jb=jbm;
                    je=jeo;
                    flag1="right"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="ss";
                    [~,~,~,~,txz_rm]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txz_rm',txz_rm);
                end
                if(k==19)
                    ib=iem;
                    ie=iem;
                    jb=jbm;
                    je=jeo;
                    flag1="right"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="vx";
                    [vx_rm,~,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vx_rm',vx_rm);
                end
                if(k==20)
                    ib=iem;
                    ie=iem;
                    jb=jbm;
                    je=jeo;
                    flag1="right"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="vz";
                    [~,vz_rm,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vz_rm',vz_rm);
                end
                %Inner layer
                if(k==21)
                    ib=iei;
                    ie=iei;
                    jb=jbi;
                    je=jeo;
                    flag1="right"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="ns";%ns=normal stresses, ss=shear streses; vx=vx, vz=vz; all=all
                    [~,~,txx_ri,tzz_ri,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txx_ri',txx_ri);
                    parsave_fk('Plane_sources/tzz_ri',tzz_ri);
                end
                if(k==22)
                    ib=iei;
                    ie=iei;
                    jb=jbi;
                    je=jeo;
                    flag1="right"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="ss";
                    [~,~,~,~,txz_ri]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txz_ri',txz_ri);
                end
                if(k==23)
                    ib=iei;
                    ie=iei;
                    jb=jbi;
                    je=jeo;
                    flag1="right"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="vx";
                    [vx_ri,~,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vx_ri',vx_ri);
                end
                if(k==24)
                    ib=iei;
                    ie=iei;
                    jb=jbi;
                    je=jeo;
                    flag1="right"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="vz";
                    [~,vz_ri,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vz_ri',vz_ri);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %---BOTTOM BOUNDARY
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Outter layer
                if(k==25)
                    ib=ibo;
                    ie=ieo;
                    jb=jeo;
                    je=jeo;
                    flag1="bottom"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="ns";%ns=normal stresses, ss=shear streses; vx=vx, vz=vz; all=all
                    [~,~,txx_bo,tzz_bo,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txx_bo',txx_bo);
                    parsave_fk('Plane_sources/tzz_bo',tzz_bo);
                end
                if(k==26)
                    ib=ibo;
                    ie=ieo;
                    jb=jeo;
                    je=jeo;
                    flag1="bottom"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="ss";
                    [~,~,~,~,txz_bo]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txz_bo',txz_bo);
                end
                if(k==27)
                    ib=ibo;
                    ie=ieo;
                    jb=jeo;
                    je=jeo;
                    flag1="bottom"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="vx";
                    [vx_bo,~,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vx_bo',vx_bo);
                end
                if(k==28)
                    ib=ibo;
                    ie=ieo;
                    jb=jeo;
                    je=jeo;
                    flag1="bottom"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="vz";
                    [~,vz_bo,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vz_bo',vz_bo);
                end
                %Middle layer
                if(k==29)
                    ib=ibo;
                    ie=ieo;
                    jb=jem;
                    je=jem;
                    flag1="bottom"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="ns";%ns=normal stresses, ss=shear streses; vx=vx, vz=vz; all=all
                    [~,~,txx_bm,tzz_bm,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txx_bm',txx_bm);
                    parsave_fk('Plane_sources/tzz_bm',tzz_bm);
                end
                if(k==30)
                    ib=ibo;
                    ie=ieo;
                    jb=jem;
                    je=jem;
                    flag1="bottom"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="ss";
                    [~,~,~,~,txz_bm]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txz_bm',txz_bm);
                end
                if(k==31)
                    ib=ibo;
                    ie=ieo;
                    jb=jem;
                    je=jem;
                    flag1="bottom"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="vx";
                    [vx_bm,~,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vx_bm',vx_bm);
                end
                if(k==32)
                    ib=ibo;
                    ie=ieo;
                    jb=jem;
                    je=jem;
                    flag1="bottom"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="vz";
                    [~,vz_bm,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vz_bm',vz_bm);
                end
                %Inner layer
                if(k==33)
                    ib=ibo;
                    ie=ieo;
                    jb=jei;
                    je=jei;
                    flag1="bottom"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="ns";%ns=normal stresses, ss=shear streses; vx=vx, vz=vz; all=all
                    [~,~,txx_bi,tzz_bi,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txx_bi',txx_bi);
                    parsave_fk('Plane_sources/tzz_bi',tzz_bi);
                end
                if(k==34)
                    ib=ibo;
                    ie=ieo;
                    jb=jei;
                    je=jei;
                    flag1="bottom"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="ss";
                    [~,~,~,~,txz_bi]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/txz_bi',txz_bi);
                end
                if(k==35)
                    ib=ibo;
                    ie=ieo;
                    jb=jei;
                    je=jei;
                    flag1="bottom"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="vx";
                    [vx_bi,~,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vx_bi',vx_bi);
                end
                if(k==36)
                    ib=ibo;
                    ie=ieo;
                    jb=jei;
                    je=jei;
                    flag1="bottom"; %left,right,bottom, lc=left corner, rc=right corner
                    flag2="vz";
                    [~,vz_bi,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    parsave_fk('Plane_sources/vz_bi',vz_bi);
                end
            end
            
            fk_lo=[];
            fk_lm=[];
            fk_li=[];
            
            txx_lo=load('Plane_sources/txx_lo');
            tzz_lo=load('Plane_sources/tzz_lo');
            txz_lo=load('Plane_sources/txz_lo');
            vx_lo=load('Plane_sources/vx_lo');
            vz_lo=load('Plane_sources/vz_lo');
            if(resample_order)
                txx_lo.fk(round(nt/resample_order)+1:end,:)=[];
                txx_lo.fk=resample(txx_lo.fk,resample_order,1);
                tzz_lo.fk(round(nt/resample_order)+1:end,:)=[];
                tzz_lo.fk=resample(tzz_lo.fk,resample_order,1);
                txz_lo.fk(round(nt/resample_order)+1:end,:)=[];
                txz_lo.fk=resample(txz_lo.fk,resample_order,1);
                vx_lo.fk(round(nt/resample_order)+1:end,:)=[];
                vx_lo.fk=resample(vx_lo.fk,resample_order,1);
                vz_lo.fk(round(nt/resample_order)+1:end,:)=[];
                vz_lo.fk=resample(vz_lo.fk,resample_order,1);
            end
            fk_lo.vx=vx_lo;
            fk_lo.vz=vz_lo;
            fk_lo.sxx=txx_lo;
            fk_lo.szz=tzz_lo;
            fk_lo.sxz=txz_lo;
            name=strcat('Plane_sources/',num2str(angleforced(jf)),'_','fk_lo','_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat');
            save(name,'fk_lo');
            
            txx_lm=load('Plane_sources/txx_lm');
            tzz_lm=load('Plane_sources/tzz_lm');
            txz_lm=load('Plane_sources/txz_lm');
            vx_lm=load('Plane_sources/vx_lm');
            vz_lm=load('Plane_sources/vz_lm');
            if(resample_order)
                txx_lm.fk(round(nt/resample_order)+1:end,:)=[];
                txx_lm.fk=resample(txx_lm.fk,resample_order,1);
                tzz_lm.fk(round(nt/resample_order)+1:end,:)=[];
                tzz_lm.fk=resample(tzz_lm.fk,resample_order,1);
                txz_lm.fk(round(nt/resample_order)+1:end,:)=[];
                txz_lm.fk=resample(txz_lm.fk,resample_order,1);
                vx_lm.fk(round(nt/resample_order)+1:end,:)=[];
                vx_lm.fk=resample(vx_lm.fk,resample_order,1);
                vz_lm.fk(round(nt/resample_order)+1:end,:)=[];
                vz_lm.fk=resample(vz_lm.fk,resample_order,1);
            end
            fk_lm.vx=vx_lm;
            fk_lm.vz=vz_lm;
            fk_lm.sxx=txx_lm;
            fk_lm.szz=tzz_lm;
            fk_lm.sxz=txz_lm;
            name=strcat('Plane_sources/',num2str(angleforced(jf)),'_','fk_lm','_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat');
            save(name,'fk_lm');
            
            txx_li=load('Plane_sources/txx_li');
            tzz_li=load('Plane_sources/tzz_li');
            txz_li=load('Plane_sources/txz_li');
            vx_li=load('Plane_sources/vx_li');
            vz_li=load('Plane_sources/vz_li');
            if(resample_order)
                txx_li.fk(round(nt/resample_order)+1:end,:)=[];
                txx_li.fk=resample(txx_li.fk,resample_order,1);
                tzz_li.fk(round(nt/resample_order)+1:end,:)=[];
                tzz_li.fk=resample(tzz_li.fk,resample_order,1);
                txz_li.fk(round(nt/resample_order)+1:end,:)=[];
                txz_li.fk=resample(txz_li.fk,resample_order,1);
                vx_li.fk(round(nt/resample_order)+1:end,:)=[];
                vx_li.fk=resample(vx_li.fk,resample_order,1);
                vz_li.fk(round(nt/resample_order)+1:end,:)=[];
                vz_li.fk=resample(vz_li.fk,resample_order,1);
            end
            fk_li.vx=vx_li;
            fk_li.vz=vz_li;
            fk_li.sxx=txx_li;
            fk_li.szz=tzz_li;
            fk_li.sxz=txz_li;
            name=strcat('Plane_sources/',num2str(angleforced(jf)),'_','fk_li','_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat');
            save(name,'fk_li');
            
            fk_ro=[];
            txx_ro=load('Plane_sources/txx_ro');
            tzz_ro=load('Plane_sources/tzz_ro');
            txz_ro=load('Plane_sources/txz_ro');
            vx_ro=load('Plane_sources/vx_ro');
            vz_ro=load('Plane_sources/vz_ro');
            if(resample_order)
                txx_ro.fk(round(nt/resample_order)+1:end,:)=[];
                txx_ro.fk=resample(txx_ro.fk,resample_order,1);
                tzz_ro.fk(round(nt/resample_order)+1:end,:)=[];
                tzz_ro.fk=resample(tzz_ro.fk,resample_order,1);
                txz_ro.fk(round(nt/resample_order)+1:end,:)=[];
                txz_ro.fk=resample(txz_ro.fk,resample_order,1);
                vx_ro.fk(round(nt/resample_order)+1:end,:)=[];
                vx_ro.fk=resample(vx_ro.fk,resample_order,1);
                vz_ro.fk(round(nt/resample_order)+1:end,:)=[];
                vz_ro.fk=resample(vz_ro.fk,resample_order,1);
            end
            fk_ro.vx=vx_ro;
            fk_ro.vz=vz_ro;
            fk_ro.sxx=txx_ro;
            fk_ro.szz=tzz_ro;
            fk_ro.sxz=txz_ro;
            name=strcat('Plane_sources/',num2str(angleforced(jf)),'_','fk_ro','_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat');
            save(name,'fk_ro');
            
            fk_rm=[];
            txx_rm=load('Plane_sources/txx_rm');
            tzz_rm=load('Plane_sources/tzz_rm');
            txz_rm=load('Plane_sources/txz_rm');
            vx_rm=load('Plane_sources/vx_rm');
            vz_rm=load('Plane_sources/vz_rm');
            if(resample_order)
                txx_rm.fk(round(nt/resample_order)+1:end,:)=[];
                txx_rm.fk=resample(txx_rm.fk,resample_order,1);
                tzz_rm.fk(round(nt/resample_order)+1:end,:)=[];
                tzz_rm.fk=resample(tzz_rm.fk,resample_order,1);
                txz_rm.fk(round(nt/resample_order)+1:end,:)=[];
                txz_rm.fk=resample(txz_rm.fk,resample_order,1);
                vx_rm.fk(round(nt/resample_order)+1:end,:)=[];
                vx_rm.fk=resample(vx_rm.fk,resample_order,1);
                vz_rm.fk(round(nt/resample_order)+1:end,:)=[];
                vz_rm.fk=resample(vz_rm.fk,resample_order,1);
            end
            fk_rm.vx=vx_rm;
            fk_rm.vz=vz_rm;
            fk_rm.sxx=txx_rm;
            fk_rm.szz=tzz_rm;
            fk_rm.sxz=txz_rm;
            name=strcat('Plane_sources/',num2str(angleforced(jf)),'_','fk_rm','_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat');
            save(name,'fk_rm');
            
            fk_ri=[];
            txx_ri=load('Plane_sources/txx_ri');
            tzz_ri=load('Plane_sources/tzz_ri');
            txz_ri=load('Plane_sources/txz_ri');
            vx_ri=load('Plane_sources/vx_ri');
            vz_ri=load('Plane_sources/vz_ri');
            if(resample_order)
                txx_ri.fk(round(nt/resample_order)+1:end,:)=[];
                txx_ri.fk=resample(txx_ri.fk,resample_order,1);
                tzz_ri.fk(round(nt/resample_order)+1:end,:)=[];
                tzz_ri.fk=resample(tzz_ri.fk,resample_order,1);
                txz_ri.fk(round(nt/resample_order)+1:end,:)=[];
                txz_ri.fk=resample(txz_ri.fk,resample_order,1);
                vx_ri.fk(round(nt/resample_order)+1:end,:)=[];
                vx_ri.fk=resample(vx_ri.fk,resample_order,1);
                vz_ri.fk(round(nt/resample_order)+1:end,:)=[];
                vz_ri.fk=resample(vz_ri.fk,resample_order,1);
            end
            fk_ri.vx=vx_ri;
            fk_ri.vz=vz_ri;
            fk_ri.sxx=txx_ri;
            fk_ri.szz=tzz_ri;
            fk_ri.sxz=txz_ri;
            name=strcat('Plane_sources/',num2str(angleforced(jf)),'_','fk_ri','_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat');
            save(name,'fk_ri');
            
            fk_bo=[];
            txx_bo=load('Plane_sources/txx_bo');
            tzz_bo=load('Plane_sources/tzz_bo');
            txz_bo=load('Plane_sources/txz_bo');
            vx_bo=load('Plane_sources/vx_bo');
            vz_bo=load('Plane_sources/vz_bo');
            fk_bo.vx=vx_bo;
            fk_bo.vz=vz_bo;
            fk_bo.sxx=txx_bo;
            fk_bo.szz=tzz_bo;
            fk_bo.sxz=txz_bo;
            name=strcat('Plane_sources/',num2str(angleforced(jf)),'_','fk_bo','_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat');
            save(name,'fk_bo');
            
            fk_bm=[];
            txx_bm=load('Plane_sources/txx_bm');
            tzz_bm=load('Plane_sources/tzz_bm');
            txz_bm=load('Plane_sources/txz_bm');
            vx_bm=load('Plane_sources/vx_bm');
            vz_bm=load('Plane_sources/vz_bm');
            if(resample_order)
                txx_bm.fk(round(nt/resample_order)+1:end,:)=[];
                txx_bm.fk=resample(txx_bm.fk,resample_order,1);
                tzz_bm.fk(round(nt/resample_order)+1:end,:)=[];
                tzz_bm.fk=resample(tzz_bm.fk,resample_order,1);
                txz_bm.fk(round(nt/resample_order)+1:end,:)=[];
                txz_bm.fk=resample(txz_bm.fk,resample_order,1);
                vx_bm.fk(round(nt/resample_order)+1:end,:)=[];
                vx_bm.fk=resample(vx_bm.fk,resample_order,1);
                vz_bm.fk(round(nt/resample_order)+1:end,:)=[];
                vz_bm.fk=resample(vz_bm.fk,resample_order,1);
            end
            fk_bm.vx=vx_bm;
            fk_bm.vz=vz_bm;
            fk_bm.sxx=txx_bm;
            fk_bm.szz=tzz_bm;
            fk_bm.sxz=txz_bm;
            name=strcat('Plane_sources/',num2str(angleforced(jf)),'_','fk_bm','_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat');
            save(name,'fk_bm');
            
            fk_bi=[];
            txx_bi=load('Plane_sources/txx_bi');
            tzz_bi=load('Plane_sources/tzz_bi');
            txz_bi=load('Plane_sources/txz_bi');
            vx_bi=load('Plane_sources/vx_bi');
            vz_bi=load('Plane_sources/vz_bi');
            if(resample_order)
                txx_bi.fk(round(nt/resample_order)+1:end,:)=[];
                txx_bi.fk=resample(txx_bi.fk,resample_order,1);
                tzz_bi.fk(round(nt/resample_order)+1:end,:)=[];
                tzz_bi.fk=resample(tzz_bi.fk,resample_order,1);
                txz_bi.fk(round(nt/resample_order)+1:end,:)=[];
                txz_bi.fk=resample(txz_bi.fk,resample_order,1);
                vx_bi.fk(round(nt/resample_order)+1:end,:)=[];
                vx_bi.fk=resample(vx_bi.fk,resample_order,1);
                vz_bi.fk(round(nt/resample_order)+1:end,:)=[];
                vz_bi.fk=resample(vz_bi.fk,resample_order,1);
            end
            fk_bi.vx=vx_bi;
            fk_bi.vz=vz_bi;
            fk_bi.sxx=txx_bi;
            fk_bi.szz=tzz_bi;
            fk_bi.sxz=txz_bi;
            name=strcat('Plane_sources/',num2str(angleforced(jf)),'_','fk_bi','_',num2str(f0(kf)),'-',num2str(f0(kf+1)),'hz.mat');
            save(name,'fk_bi');
            
        end
    end
end