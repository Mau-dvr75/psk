    seisA_vx=zeros(nt,length(rec_x_id));
    seisA_vz=zeros(nt,length(rec_x_id));
    if (makeobs==1)
        vp=extend_model(vp,nx,nz,ncpml);
        vs=extend_model(vs,nx,nz,ncpml);
        rho=extend_model(rho,nx,nz,ncpml);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ANALYTICAL VS NUMERICAL COMPARISON
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     nsrc=19;
%     angleforced=fwi_info.inc_angle;
    for kf=1:ekf
        for jf=1:1
            for i=1:Nx
                for j=1:Nz
                    mu(j,i)=rho(j,i)*vs(j,i)*vs(j,i);
                    lamb(j,i)=rho(j,i)*vp(j,i)*vp(j,i)-2.*mu(j,i);
                end
            end
            ne=length(SOURCE_TYPE);
            
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
            HDUR=1./f0;
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
            for ssz=1:1
                for ssx=1:length(rec_x_id)
                    ib=rec_x_id(ssx);
                    ie=ib;
                    jb=rec_z_id(ssz);
                    je=jb;
                    flag1="an";
                    flag2="vx";
                    [vx_an,~,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    flag2="vz";
                    [~,vz_an,~,~,~]=fk_field(al_fk,be_fk,mul_fk,h_fk,nlayer_fk,HDUR(kf),p,x0_source(kf,jf),z0_source(jf),tstart,dt_fk,nstepfk,stag,lamb,mu,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE(jf),STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);
                    
                     vx_an=resample(vx_an,resample_order,1);
                     vz_an=resample(vz_an,resample_order,1);
                    
                    seisA_vx(:,ssx,ssz)=vx_an(1:nt);
                    seisA_vz(:,ssx,ssz)=vz_an(1:nt);
                end
            end
        end
    end
    if(saveseis==1)
        namevxa=strcat('Analytical/vx_analytic_OP',num2str(op),'_FO',num2str(field_order),'_IO',num2str(INJECTIONS_ORDER),'shot_',num2str(jf));
        namevza=strcat('Analytical/vz_analytic_OP',num2str(op),'_FO',num2str(field_order),'_IO',num2str(INJECTIONS_ORDER),'shot_',num2str(jf));

        save(namevxa,'seisA_vx');
        save(namevza,'seisA_vz');
        %  save(nameux,'u_obsx');
        %  save(nameuz,'u_obsz');
    end
    
    