vxa=load('./Analytical/vx_analytic_OP4_FO2_IO2shot_1.mat');
vxa=vxa.seisA_vx;
vza=load('./Analytical/vz_analytic_OP4_FO2_IO2shot_1.mat');
vza=vza.seisA_vz;

vxs=load('./Observed_data/vx_shot1_iter1_0.03-0.25hz.mat');
vxs=vxs.vx_obs;
vzs=load('./Observed_data/vz_shot1_iter1_0.03-0.25hz.mat');
vzs=vzs.vz_obs;

nr=20;
figure
subplot(2,1,1)
plot(t,vxa(:,nr))
hold on;
plot(t,vxs(nr,:))
legend VXanalytical VXnumerical
title(['Receiver # ',num2str(nr)])
subplot(2,1,2)
plot(t,vza(:,nr))
hold on;
plot(t,vzs(nr,:))
legend VZanalytical VZnumerical

[nt,nrecp]=size(vza);
t=0:dt:(nt-1)*dt;
afac=7000;fac=1;ww=1;
figure
for ic=1:2
    for i=1:nrecp
        subplot(1,2,ic)
        if ic==1
            plot(t(1:nt),vza(:,i)*afac+dgc(i)/fac,'r','LineWidth',ww)
            hold on;
            tr=vzs(i,:);
            plot(t(1:length(tr)),tr*afac+dgc(i)/fac,'k','LineWidth',ww)
            xlim([0 220])
            title (['Vertical velocity'])
        end
        if ic==2
            plot(t(1:nt),vxa(:,i)*afac+dgc(i)/fac,'r','LineWidth',ww)
            hold on;
            tr=vxs(i,:);
            plot(t(1:length(tr)),tr*afac+dgc(i)/fac,'k','LineWidth',ww)
            xlim([0 220])
            % ylim([-5 550])
            title (['Horizontal velocity'])
            legend analytical modeled
        end
    end
end