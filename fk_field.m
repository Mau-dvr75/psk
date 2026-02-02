function [vx,vz,txx,tzz,txz]=fk_field(al,be,mu,H,nlayer,Tg,p,x0,z0,t0,dt,npts,staggered,lamb,mul,ib,ie,jb,je,flag1,flag2,x_axis,z_axis,SOURCE_TYPE,STF_TYPE,f0,kf,catalog,jf,resample_order,stf_method,ft,delay,noa);%,NF_FOR_STORING,NPOW_FOR_FFT,NP_RESAMP,DF_FK);
%I've changed the sign in VZ and TXZ componenents to make it consistent with
%my coordinate axis system. Otherwise, it is the same as SPECFEM2D-FK

if(flag1=="left")
    np=length(jb:je);
end

if(flag1=="right")
    np=length(jb:je);
end

if(flag1=="bottom")
    np=length(ib:ie);
end

if(flag1=="lc")
    np=1;
end

if(flag1=="rc")
    np=1;
end

if(flag1=="an")
    np=1;
end
%Allocate wavefields
vx=zeros(npts,np);
vz=zeros(npts,np);
txx=zeros(npts,np);
tzz=zeros(npts,np);
txz=zeros(npts,np);

%Common freq information for all fields

% df    = DF_FK;
% nf2   = NF_FOR_STORING+1;   %! number of positive frequency sample points
% nf    = 2*NF_FOR_STORING;   %! number of total frequencies after symetrisation
% npts2 = nf;
%
% %!! VM VM recompute new values for new way to do
% npow = ceil(log(npts2*1.0)/log(2.0));
% npts2 = 2^npow;
% NPOW_FOR_FFT = npow;
% dt_fk = 1/(df*(npts2-1));
% % !! number of points for resmpled vector
% npoints2 = NP_RESAMP*(npts2-1)+1;
%
% lnpts2=ceil(log(npts2*1)/log(2));
% npts2=2^lnpts2;
% %nf=npts2;
%nf2=nf/2+1; %Positive frequency points
% df=1/(dt*nf);
% fmax=1/(2*dt);
lnpts2=ceil(log(npts*1)/log(2));
npts2=2^lnpts2;
nf=npts2;
nf2=nf/2+1; %Positive frequency points
df=1/(dt*nf);
fmax=1/(2*dt);
fvec=zeros(1,nf2);
STF_TYPE=STF_TYPE(jf);
stf_method=stf_method(jf);
%SOURCE_TYPE=SOURCE_TYPE(jf);
dtm=dt/resample_order;

for i=1:nf2
    fvec(i)=(i-1)*df;
end

if STF_TYPE==4
    if stf_method==1
        method='average';
    elseif stf_method==2
        method='pratt';
    end
    %     if ft==1 || ft==2
    %         ftname='Hz_';
    %     else
    %         ftname='nf_';
    %     end
    %
    stf=zeros(1,nf2);
    cof1=f0(kf);
    cof2=f0(kf+1);
    if ft==0;
        %namestf=['./data_mase/STF/STF_',num2str(cof1),'-',num2str(cof2),'Hz_event',num2str(jf),catalog,'_dt',num2str(dtm),method,'.mat'];
        namestf=['./data_korea/stations/arrays/',num2str(noa),'/STF/',num2str(cof1),'-',num2str(cof2),'Hz_event',num2str(jf),catalog,'_dt',num2str(0.01),method,'.mat'];%stf is at 0.01 samples
    else
        %namestf=['./data_mase/STF/STF_nf_event',num2str(jf),catalog,'_dt',num2str(dtm),method,'.mat'];
    end
    stf_temp=load(namestf);
    stf_temp=stf_temp.stf_final;
    % if stf_method==1
    % stf_temp=stf_temp.z;
    % end
    stf_temp=downsample(stf_temp,4);%stf at 0.01, FD dt is dt (0.04)

    fs=1/dtm;
    if ft==1
        fob=6;
        [B,A] = butter(fob,[cof1,cof2]/(fs/2),'bandpass');
        [zhi,phi,khi] = butter(fob,[cof1,cof2]/(fs/2),'bandpass');
        soshi = zp2sos(zhi,phi,khi);
        % stf_temp=filter(B,A,stf_temp);
        stf_temp=sosfilt(soshi,stf_temp);
    elseif ft==2
        bpfilt = designfilt( 'bandpassfir' , ...
            'FilterOrder' ,512, 'CutoffFrequency1' ,cof1, ...
            'CutoffFrequency2' ,cof2, 'SampleRate' ,fs);
        stf_temp=filter(bpfilt,stf_temp);
        % stf_temp=stf_temp(mean(grpdelay(bpfilt)):end);
    end

    corrfac=[290 0 250 0 290 290 250 290 0 0 0 0 290 250 0 290 290 0 0];
    %delay=corrfac(jf);
    %stf_temp=cumtrapz(stf_temp); %hereeeeeeeeeee
    %stf_temp=-diff(stf_temp);
    % stf_temp=stf_temp(corrfac:end);
    delay=1;
    stf_temp=stf_temp(delay:end);
    win=tukeywin(length(stf_temp),0.01);
    %stf_t=stf_av.*win;

    % stf_temp=stf_temp-mean(stf_temp);
    % stf_temp=detrend(stf_temp);
    % stf_temp=stf_temp.*win;
    %
    %     figure
    %     plot(stf_temp)
    %     legend actstfusedinFK

    ls=length(stf_temp);
    cl=round(npts*resample_order);
    stf_temp=[stf_temp;zeros(cl-ls,1)];
    stf_temp=downsample(stf_temp,resample_order);
    if length(stf_temp)>length(stf)
        stf(1:length(stf_temp))=stf_temp;
    else
        stf=stf_temp;
    end
    STF=fft(stf);
end
if STF_TYPE==5
    %GOOD FOR GF
    fo=5;
    % d = designfilt( 'lowpassfir','FilterOrder' ,250,'PassbandFrequency' ,1.8, ...      % Frequency constraints
    %           'StopbandFrequency' ,2,'DesignMethod' , 'ls' ,'PassbandWeight' ,1, ...           % Design method options
    %           'StopbandWeight',2,'SampleRate' ,(1/(0.02)/2))                % Sample rate
    %      [B,A] = butter(8,[5]/(1/(0.02)/2),'low'); %0.02 dt simulation
    %      %[B,A] = butter(fo,[1]/(1/(dt)/2),'low');
    %      stf=zeros(1,10002);
    %      stf(25)=1;
    %      stf=filtfilt(B,A,stf);
    %      stf=downsample(stf,3);
    %     STF=fft(stf);
    %
    %      stf=zeros(1,nf2);
    %%%--for dh=500, dt=0.02
%     stf=make_stf(10002,0.02,2,1,1,0);
%     stf=stf(1:10002);
%     stf=downsample(stf,3);
%     STF=fft(stf);
    %%%---for fh=900, dy=0.04
    %5001 with dt=.04;
    f0gf=1;
    stf=make_stf(10001,0.04,f0gf,1,1,0);
    stf=stf(1:10001);
    stf=downsample(stf,resample_order);
    STF=fft(stf);
    % % %
    
    %pprooo
    %      stf=zeros(1,nf2);
    %      stf_temp=make_stf(npts,dt,0.8,1,3,0);
    %      stf(1:npts)=stf_temp;
    %      %stf=downsample(stf,2);
    %      STF=fft(stf);
end

%   if(flag1=="lc")
%
%       display(sprintf('Entering the FK synthetics program'))
%       display(sprintf('Number of points used for FFT %f= ',npts2))
%       display(sprintf('Number of samples stored for FK solution = %f', NF_FOR_STORING))
%       display(sprintf('Total time length used for FK= %f ', t0+(npts2-1)*dt_fk))
%       display(sprintf('  FK time step       = %f', dt_fk))
%       display(sprintf('  FK frequency step  = %f ', df))
%       display(sprintf('  power of 2 for FFT = %d ', npow))
%   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALL FIELDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (flag2=="all")
    %We assume by default that nodes coordinates are in normal stress
    %positions
    xcoord=x_axis.sxx;
    zcoord=z_axis.sxx;
    nvar=5;
    % display(sprintf('FK_2D program.Computing all fields in %s boundary portion',flag1))
    %fprintf('Enter fk synthetics program, calculating all fields')
    %     lnpts2=ceil(log(npts*1)/log(2));
    %     npts2=2^lnpts2;
    %     nf=npts2;
    %     nf2=nf/2+1; %Positive frequency points
    %     df=1/(dt*nf);
    %     fmax=1/(2*dt);
    %     fvec=zeros(1,nf2);
    %     for i=1:nf2
    %         fvec(i)=(i-1)*df;
    %     end
    
    if (SOURCE_TYPE==1)
        i=complex(0,1);
        C_3=i*p*al(nlayer); %amplitude of incoming P in bott layer
        eta_p=sqrt(1/al(nlayer)^2-p^2); %vertical slowness for lowest layer
        
    end
    if (SOURCE_TYPE==2)
        % p=sin(angleforce)/be_fk(nlayer_fk);
        C_3= p*be(nlayer);
        eta_p=sqrt(1/be(nlayer)^2-p^2);
    end
    
    
    coeff=zeros(2,nf2);
    field_f=zeros(nf,nvar);
    field=zeros(npts2,nvar);
    %dtmp=zeros(1,npts);
    
    %wave coefficients in bottom layer for al freqs
    for i=1:nf2
        om=2*pi*fvec(i);
        [N_mat]=compute_N_rayleigh(al,be,mu,H,nlayer,om,p,sum(H(1:nlayer-1)));
        a=N_mat(3,2);
        b=N_mat(3,4);
        c=N_mat(4,2);
        d=N_mat(4,4);
        delta_mat=a*d-b*c;
        coeff(1,i)=-(d*N_mat(3,3)-b*N_mat(4,3))/delta_mat*C_3;
        coeff(2,i)=-(-c*N_mat(3,3)+a*N_mat(4,3))/delta_mat*C_3;
    end
    
    ip=1;
    for i=ib:ie
        for j=jb:je
            
            %X and Z coordinates of point
            x=xcoord(i);
            z=zcoord(j);
            
            if(flag1=="lc")
                x;
                z;
            end
            
            %Defining parameters for each point
            xil=mul(j,i)/(lamb(j,i)+2*mul(j,i));
            xi1=1-2.0*xil;
            xim=(1-xil)*mul(j,i);
            field_f=0;
            
            %Delay time
            if (flag1=="lc")
                tdelay=p*(x-x0)+eta_p*(0-z0);
                display(['Delay in left corner ',num2str(tdelay),' [s]'])
            end
            if (flag1=="rc")
                tdelay=p*(x-x0)+eta_p*(0-z0);
                display(['Delay in right corner ',num2str(tdelay),' [s]'])
            end
            
            
            for ii=1:nf2
                %Angular freq. vector
                om=2*pi*fvec(ii);
                if STF_TYPE==1
                    %Gaussian pulse
                    stf_coeff=exp(-(om*Tg/2)^2);
                end
                if STF_TYPE==2
                    %1st derivative of gaussian
                    stf_coeff=exp(-(om*Tg/2)^2)*om*complex(0,1);
                end
                if STF_TYPE==3
                    %1st derivative of gaussian
                    stf_coeff=exp(-(om*Tg/2)^2)*om*complex(0,1)*om*complex(0,1);
                end
                if STF_TYPE==4 || STF_TYPE==5
                    stf_coeff=STF(ii);
                end
                stf_coeff=stf_coeff*exp(complex(0,-1)*om*tdelay);
                [N_mat]=compute_N_rayleigh(al,be,mu,H,nlayer,om,p,z);
                
                %Velocities
                dx_f=N_mat(1,2)*coeff(1,ii)+N_mat(1,4)*coeff(2,ii)+N_mat(1,3)*C_3; %y1
                dz_f=N_mat(2,2)*coeff(1,ii)+N_mat(2,4)*coeff(2,ii)+N_mat(2,3)*C_3; %y3
                field_f(ii,1)=stf_coeff*dx_f*complex(0,-1)*complex(0,om);
                field_f(ii,2)=stf_coeff*dz_f*complex(0,om);
                
                %Stresses
                txz_f=N_mat(3,2)*coeff(1,ii)+N_mat(3,4)*coeff(2,ii)+N_mat(3,3)*C_3;%y4
                tzz_f=N_mat(4,2)*coeff(1,ii)+N_mat(4,4)*coeff(2,ii)+N_mat(4,3)*C_3;%y6
                field_f(ii,3)=stf_coeff*om*p*(xi1*tzz_f-4*xim*dx_f);% T_xx
                field_f(ii,4)=stf_coeff*om*p*txz_f*complex(0,-1);% T_xz
                field_f(ii,5)=stf_coeff*om*p*tzz_f;% T_zz
            end
            
            %Construct the negative part of FFT
            
            for ii=2:nf2-1
                field_f(nf+2-ii,:) = conj(field_f(ii,:));
            end
            
            
            %Inverse Fourier Transform
            for j=1:nvar
                field(:,j)=ifft(field_f(:,j))/dt;
            end
            
            %Stack velocity by deltat/2 time
            %             if (staggered==1)
            %                 for ii = 1:2
            %                     field(1:npts-1,ii)=(field(1:npts-1,ii)+field(2:npts,ii))/2;
            %                 end
            %             end
            
            vx(1:npts,ip)=real(field(1:npts,1));
            vz(1:npts,ip)=-real(field(1:npts,2));
            
            %Stack stress by deltat time
            %             if (staggered==1)
            %                 for ii = 3:5
            %                     field(1:npts-1,ii)=field(2:npts,ii);
            %                 end
            %             end
            if (staggered==1)
                for ii = 3:5
                    field(1:npts-1,ii)=(field(1:npts-1,ii)+field(2:npts,ii))/2;
                end
            end
            
            
            txx(1:npts,ip)=real(field(1:npts,3));
            txz(1:npts,ip)=-real(field(1:npts,4));
            tzz(1:npts,ip)=real(field(1:npts,5));
            
            %             if(mod(ip,100)==0)
            %                 display(sprintf('It has been computed %d of %d',ip,np))
            %             end
            
            ip=ip+1;
        end
    end
    
    %name=strcat(flag1,'_',flag2);
    
    
    
    %display(sprintf('Done'))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NORMAL STRESSES FIELDS ONLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (flag2=="ns")
    xcoord=x_axis.sxx;
    zcoord=z_axis.sxx;
    nvar=2;
    %display(sprintf('FK_2D program.Computing sxx and szz in %s boundary portion',flag1))
    %fprintf('Enter fk synthetics program, calculating normal stresses')
    %     lnpts2=ceil(log(npts*1)/log(2));
    %     npts2=2^lnpts2;
    %     nf=npts2;
    %     nf2=nf/2+1; %Positive frequency points
    %     df=1/(dt*nf);
    %     fmax=1/(2*dt);
    %     fvec=zeros(1,nf2);
    %     for i=1:nf2
    %         fvec(i)=(i-1)*df;
    %     end
    
    %     i=complex(0,1);
    %     C_3=i*p*al(nlayer); %amplitude of incoming P in bott layer
    %     eta_p=sqrt(1/al(nlayer)^2-p^2); %vertical slowness for lowest layer
    if (SOURCE_TYPE==1)
        i=complex(0,1);
        C_3=i*p*al(nlayer); %amplitude of incoming P in bott layer
        eta_p=sqrt(1/al(nlayer)^2-p^2); %vertical slowness for lowest layer
        
    end
    if (SOURCE_TYPE==2)
        % p=sin(angleforce)/be_fk(nlayer_fk);
        C_3= p*be(nlayer);
        eta_p=sqrt(1/be(nlayer)^2-p^2);
    end
    coeff=zeros(2,nf2);
    field_f=zeros(nf,nvar);
    field=zeros(npts2,nvar);
    
    %wave coefficients in bottom layer for al freqs
    for i=1:nf2
        om=2*pi*fvec(i);
        [N_mat]=compute_N_rayleigh(al,be,mu,H,nlayer,om,p,sum(H(1:nlayer-1)));
        a=N_mat(3,2);
        b=N_mat(3,4);
        c=N_mat(4,2);
        d=N_mat(4,4);
        delta_mat=a*d-b*c;
        coeff(1,i)=-(d*N_mat(3,3)-b*N_mat(4,3))/delta_mat*C_3;
        coeff(2,i)=-(-c*N_mat(3,3)+a*N_mat(4,3))/delta_mat*C_3;
    end
    
    ip=1;
    for i=ib:ie
        for j=jb:je
            
            %X and Z coordinates of point
            x=xcoord(i);
            z=zcoord(j);
            
            %Defining parameters for each point
            xil=mul(j,i)/(lamb(j,i)+2*mul(j,i));
            xi1=1-2.0*xil;
            xim=(1-xil)*mul(j,i);
            field_f=0;
            
            %Delay time
            tdelay=p*(x-x0)+eta_p*(0-z0);
            
            for ii=1:nf2
                %Angular freq. vector
                om=2*pi*fvec(ii);
                if STF_TYPE==1
                    %Gaussian pulse
                    stf_coeff=exp(-(om*Tg/2)^2);
                end
                if STF_TYPE==2
                    %1st derivative of gaussian
                    stf_coeff=exp(-(om*Tg/2)^2)*om*complex(0,1);
                end
                if STF_TYPE==3
                    %1st derivative of gaussian
                    stf_coeff=exp(-(om*Tg/2)^2)*om*complex(0,1)*om*complex(0,1);
                end
                if STF_TYPE==4 || STF_TYPE==5
                    stf_coeff=STF(ii);
                end
                stf_coeff=stf_coeff*exp(complex(0,-1)*om*tdelay);
                
                [N_mat]=compute_N_rayleigh(al,be,mu,H,nlayer,om,p,z);
                
                %Stresses
                dx_f=N_mat(1,2)*coeff(1,ii)+N_mat(1,4)*coeff(2,ii)+N_mat(1,3)*C_3; %y1
                tzz_f=N_mat(4,2)*coeff(1,ii)+N_mat(4,4)*coeff(2,ii)+N_mat(4,3)*C_3;%y6
                field_f(ii,1)=stf_coeff*om*p*(xi1*tzz_f-4*xim*dx_f);% T_xx
                field_f(ii,2)=stf_coeff*om*p*tzz_f;% T_zz
            end
            
            %Construct the negative part of FFT
            for ii=2:nf2-1
                field_f(nf+2-ii,:) = conj(field_f(ii,:));
            end
            
            %Inverse Fourier Transform
            for j=1:nvar
                field(:,j)=ifft(field_f(:,j))/dt;
                %field(:,5)=ifft(field_f(:,5))/dt;
            end
            
            %Stack stress by deltat time
            if (staggered==1)
                for ii = 1:nvar
                    %  field(1:npts-1,ii)=field(2:npts,ii);
                    field(1:npts-1,ii)=(field(1:npts-1,ii)+field(2:npts,ii))/2;
                end
            end
            
            txx(1:npts,ip)=real(field(1:npts,1));
            tzz(1:npts,ip)=real(field(1:npts,2));
            
            %             if(mod(ip,100)==0)
            %                 display(sprintf('It has been computed %d of %d',ip,np))
            %             end
            
            ip=ip+1;
        end
    end
    %display(sprintf('Done'))
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SHEAR STRESS FIELD ONLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag2=="ss")
    xcoord=x_axis.sxz;
    zcoord=z_axis.sxz;
    nvar=1;
    %display(sprintf('FK_2D program.Computing sxz in %s boundary portion',flag1))
    % fprintf('Enter fk synthetics program, calculating shear stresses only')
    %     lnpts2=ceil(log(npts*1)/log(2));
    %     npts2=2^lnpts2;
    %     nf=npts2;
    %     nf2=nf/2+1; %Positive frequency points
    %     df=1/(dt*nf);
    %     fmax=1/(2*dt);
    %     fvec=zeros(1,nf2);
    %     for i=1:nf2
    %         fvec(i)=(i-1)*df;
    %     end
    
    %     i=complex(0,1);
    %     C_3=i*p*al(nlayer); %amplitude of incoming P in bott layer
    %     eta_p=sqrt(1/al(nlayer)^2-p^2); %vertical slowness for lowest layer
    if (SOURCE_TYPE==1)
        i=complex(0,1);
        C_3=i*p*al(nlayer); %amplitude of incoming P in bott layer
        eta_p=sqrt(1/al(nlayer)^2-p^2); %vertical slowness for lowest layer
        
    end
    if (SOURCE_TYPE==2)
        %p=sin(angleforce)/be_fk(nlayer_fk);
        C_3= p*be(nlayer);
        eta_p=sqrt(1/be(nlayer)^2-p^2);
    end
    coeff=zeros(2,nf2);
    field_f=zeros(nf,nvar);
    field=zeros(npts2,nvar);
    %dtmp=zeros(1,npts);
    
    %wave coefficients in bottom layer for al freqs
    for i=1:nf2
        om=2*pi*fvec(i);
        [N_mat]=compute_N_rayleigh(al,be,mu,H,nlayer,om,p,sum(H(1:nlayer-1)));
        a=N_mat(3,2);
        b=N_mat(3,4);
        c=N_mat(4,2);
        d=N_mat(4,4);
        delta_mat=a*d-b*c;
        coeff(1,i)=-(d*N_mat(3,3)-b*N_mat(4,3))/delta_mat*C_3;
        coeff(2,i)=-(-c*N_mat(3,3)+a*N_mat(4,3))/delta_mat*C_3;
    end
    
    ip=1;
    for i=ib:ie
        for j=jb:je
            
            %X and Z coordinates of point
            x=xcoord(i);
            z=zcoord(j);
            
            %Defining parameters for each point
            xil=mul(j,i)/(lamb(j,i)+2*mul(j,i));
            xi1=1-2.0*xil;
            xim=(1-xil)*mul(j,i);
            field_f=0;
            
            %Delay time
            tdelay=p*(x-x0)+eta_p*(0-z0);
            
            for ii=1:nf2
                %Angular freq. vector
                om=2*pi*fvec(ii);
                if STF_TYPE==1
                    %Gaussian pulse
                    stf_coeff=exp(-(om*Tg/2)^2);
                end
                if STF_TYPE==2
                    %1st derivative of gaussian
                    stf_coeff=exp(-(om*Tg/2)^2)*om*complex(0,1);
                end
                if STF_TYPE==3
                    %1st derivative of gaussian
                    stf_coeff=exp(-(om*Tg/2)^2)*om*complex(0,1)*om*complex(0,1);
                end
                if STF_TYPE==4 || STF_TYPE==5
                    stf_coeff=STF(ii);
                end
                stf_coeff=stf_coeff*exp(complex(0,-1)*om*tdelay);
                
                [N_mat]=compute_N_rayleigh(al,be,mu,H,nlayer,om,p,z);
                
                
                %Stresses
                txz_f=N_mat(3,2)*coeff(1,ii)+N_mat(3,4)*coeff(2,ii)+N_mat(3,3)*C_3;%y4
                field_f(ii)=stf_coeff*om*p*txz_f*complex(0,-1);% T_xz
            end
            
            %Construct the negative part of FFT
            for ii=2:nf2-1
                field_f(nf+2-ii) = conj(field_f(ii));
            end
            
            
            %Inverse Fourier Transform
            field(:)=ifft(field_f(:))/dt;
            
            %Stack stress by deltat time
            if (staggered==1)
                %field(1:npts-1)=field(2:npts);
                field(1:npts-1,1)=(field(1:npts-1,1)+field(2:npts,1))/2;
            end
            
            %txx(1:npts,ip)=field(1:npts,3);
            txz(1:npts,ip)=-real(field(1:npts,1));
            %tzz(1:npts,ip)=field(1:npts,5);
            
            %             if(mod(ip,100)==0)
            %                 display(sprintf('It has been computed %d of %d',ip,np))
            %             end
            
            ip=ip+1;
        end
    end
    %display(sprintf('Done'))
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VX FIELD ONLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (flag2=="vx")
    xcoord=x_axis.vx;
    zcoord=z_axis.vx;
    nvar=1;
    
    %display(sprintf('FK_2D program.Computing vx in %s boundary portion',flag1))
    %fprintf('Enter fk synthetics program, calculating all fields')
    %     lnpts2=ceil(log(npts*1)/log(2));
    %     npts2=2^lnpts2;
    %     nf=npts2;
    %     nf2=nf/2+1; %Positive frequency points
    %     df=1/(dt*nf);
    %     fmax=1/(2*dt);
    %     fvec=zeros(1,nf2);
    %     for i=1:nf2
    %         fvec(i)=(i-1)*df;
    %     end
    
    %     i=complex(0,1);
    %     C_3=i*p*al(nlayer); %amplitude of incoming P in bott layer
    %     eta_p=sqrt(1/al(nlayer)^2-p^2); %vertical slowness for lowest layer
    if (SOURCE_TYPE==1)
        i=complex(0,1);
        C_3=i*p*al(nlayer); %amplitude of incoming P in bott layer
        eta_p=sqrt(1/al(nlayer)^2-p^2); %vertical slowness for lowest layer
        
    end
    if (SOURCE_TYPE==2)
        % p=sin(angleforce)/be_fk(nlayer_fk);
        C_3= p*be(nlayer);
        eta_p=sqrt(1/be(nlayer)^2-p^2);
    end
    
    coeff=zeros(2,nf2);
    field_f=zeros(nf,nvar);
    field=zeros(npts2,nvar);
    %dtmp=zeros(1,npts);
    
    %wave coefficients in bottom layer for al freqs
    for i=1:nf2
        om=2*pi*fvec(i);
        [N_mat]=compute_N_rayleigh(al,be,mu,H,nlayer,om,p,sum(H(1:nlayer-1)));
        a=N_mat(3,2);
        b=N_mat(3,4);
        c=N_mat(4,2);
        d=N_mat(4,4);
        delta_mat=a*d-b*c;
        coeff(1,i)=-(d*N_mat(3,3)-b*N_mat(4,3))/delta_mat*C_3;
        coeff(2,i)=-(-c*N_mat(3,3)+a*N_mat(4,3))/delta_mat*C_3;
    end
    
    ip=1;
    for i=ib:ie
        for j=jb:je
            
            %X and Z coordinates of point
            x=xcoord(i);
            z=zcoord(j);
            
            %Defining parameters for each point
            xil=mul(j,i)/(lamb(j,i)+2*mul(j,i));
            xi1=1-2.0*xil;
            xim=(1-xil)*mul(j,i);
            field_f=0;
            
            %Delay time
            tdelay=p*(x-x0)+eta_p*(0-z0);
            
            for ii=1:nf2
                %Angular freq. vector
                om=2*pi*fvec(ii);
                if STF_TYPE==1
                    %Gaussian pulse
                    stf_coeff=exp(-(om*Tg/2)^2);
                end
                if STF_TYPE==2
                    %1st derivative of gaussian
                    stf_coeff=exp(-(om*Tg/2)^2)*om*complex(0,1);
                end
                if STF_TYPE==3
                    %1st derivative of gaussian
                    stf_coeff=exp(-(om*Tg/2)^2)*om*complex(0,1)*om*complex(0,1);
                end
                if STF_TYPE==4 || STF_TYPE==5
                    stf_coeff=STF(ii);
                end
                stf_coeff=stf_coeff*exp(complex(0,-1)*om*tdelay);
                
                [N_mat]=compute_N_rayleigh(al,be,mu,H,nlayer,om,p,z);
                
                %Velocities
                dx_f=N_mat(1,2)*coeff(1,ii)+N_mat(1,4)*coeff(2,ii)+N_mat(1,3)*C_3; %y1
                field_f(ii)=stf_coeff*dx_f*complex(0,-1)*complex(0,om);
            end
            
            %Construct the negative part of FFT
            
            for ii=2:nf2-1
                field_f(nf+2-ii) = conj(field_f(ii));
            end
            
            
            %Inverse Fourier Transform
            
            field(:)=ifft(field_f(:))/dt;
            
            
            %Stack velocity by deltat/2 time
            %             if (staggered==1)
            %
            %                     field(1:npts-1)=(field(1:npts-1)+field(2:npts))/2;
            %
            %             end
            
            vx(1:npts,ip)=real(field(1:npts));
            
            %             if(mod(ip,100)==0)
            %                 display(sprintf('It has been computed %d of %d',ip,np))
            %             end
            
            ip=ip+1;
        end
    end
    %display(sprintf('Done'))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VZ FIELD ONLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (flag2=="vz")
    xcoord=x_axis.vz;
    zcoord=z_axis.vz;
    nvar=1;
    
    %display(sprintf('FK_2D program.Computing vz in %s boundary portion',flag1))
    %fprintf('Enter fk synthetics program, calculating all fields')
    %     lnpts2=ceil(log(npts*1)/log(2));
    %     npts2=2^lnpts2;
    %     nf=npts2;
    %     nf2=nf/2+1; %Positive frequency points
    %     df=1/(dt*nf);
    %     fmax=1/(2*dt);
    %     fvec=zeros(1,nf2);
    %     for i=1:nf2
    %         fvec(i)=(i-1)*df;
    %     end
    
    %     i=complex(0,1);
    %     C_3=i*p*al(nlayer); %amplitude of incoming P in bott layer
    %     eta_p=sqrt(1/al(nlayer)^2-p^2); %vertical slowness for lowest layer
    if (SOURCE_TYPE==1)
        i=complex(0,1);
        C_3=i*p*al(nlayer); %amplitude of incoming P in bott layer
        eta_p=sqrt(1/al(nlayer)^2-p^2); %vertical slowness for lowest layer
        
    end
    if (SOURCE_TYPE==2)
        % p=sin(angleforce)/be_fk(nlayer_fk);
        C_3= p*be(nlayer);
        eta_p=sqrt(1/be(nlayer)^2-p^2);
    end
    
    coeff=zeros(2,nf2);
    field_f=zeros(nf,nvar);
    field=zeros(npts2,nvar);
    %dtmp=zeros(1,npts);
    
    %wave coefficients in bottom layer for al freqs
    for i=1:nf2
        om=2*pi*fvec(i);
        [N_mat]=compute_N_rayleigh(al,be,mu,H,nlayer,om,p,sum(H(1:nlayer-1)));
        a=N_mat(3,2);
        b=N_mat(3,4);
        c=N_mat(4,2);
        d=N_mat(4,4);
        delta_mat=a*d-b*c;
        coeff(1,i)=-(d*N_mat(3,3)-b*N_mat(4,3))/delta_mat*C_3;
        coeff(2,i)=-(-c*N_mat(3,3)+a*N_mat(4,3))/delta_mat*C_3;
    end
    
    ip=1;
    for i=ib:ie
        for j=jb:je
            
            %X and Z coordinates of point
            x=xcoord(i);
            z=zcoord(j);
            
            %Defining parameters for each point
            xil=mul(j,i)/(lamb(j,i)+2*mul(j,i));
            xi1=1-2.0*xil;
            xim=(1-xil)*mul(j,i);
            field_f=0;
            
            %Delay time
            tdelay=p*(x-x0)+eta_p*(0-z0);
            
            for ii=1:nf2
                %Angular freq. vector
                om=2*pi*fvec(ii);
                if STF_TYPE==1
                    %Gaussian pulse
                    stf_coeff=exp(-(om*Tg/2)^2);
                end
                if STF_TYPE==2
                    %1st derivative of gaussian
                    stf_coeff=exp(-(om*Tg/2)^2)*om*complex(0,1);
                end
                if STF_TYPE==3
                    %1st derivative of gaussian
                    stf_coeff=exp(-(om*Tg/2)^2)*om*complex(0,1)*om*complex(0,1);
                end
                if STF_TYPE==4 || STF_TYPE==5
                    stf_coeff=STF(ii);
                end
                stf_coeff=stf_coeff*exp(complex(0,-1)*om*tdelay);
                
                [N_mat]=compute_N_rayleigh(al,be,mu,H,nlayer,om,p,z);
                
                %Velocities
                dz_f=N_mat(2,2)*coeff(1,ii)+N_mat(2,4)*coeff(2,ii)+N_mat(2,3)*C_3; %y3
                field_f(ii,1)=stf_coeff*dz_f*complex(0,om);
                
            end
            
            %Construct the negative part of FFT
            
            for ii=2:nf2-1
                field_f(nf+2-ii) = conj(field_f(ii));
            end
            
            
            %Inverse Fourier Transform
            field(:)=ifft(field_f(:))/dt;
            
            %Stack velocity by deltat/2 time
            %             if (staggered==1)
            %
            %                     field(1:npts-1)=(field(1:npts-1)+field(2:npts))/2;
            %
            %             end
            
            vz(1:npts,ip)=-real(field(1:npts));
            
            %             if(mod(ip,100)==0)
            %                 display(sprintf('It has been computed %d of %d',ip,np))
            %             end
            
            ip=ip+1;
        end
    end
    %display(sprintf('Done'))
end


end

%ANOTHER SOURCE (IS WRONG)
%Ricker pulse
%tp=0.1;
%ts=2;
%wp=2*pi/tp;
%b=om./wp;
%stf_coeff=((-tp/(pi)^(1/2))*b.^2.*exp(-b.^2).*exp(-1i*om*ts));





