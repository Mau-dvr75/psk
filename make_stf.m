function stf=make_stf(nt,dt,f0,amp,type,plot_source)
    t=0:dt:dt*(nt-1);
    switch type 
%         case 'delta'
%          stf=zeros(1,length(t)); 
%          stf(1)=1/dt; 
        case 1
            t0=1.2/f0; %source time duration (and shift) in seconds
            a = pi*pi*f0*f0;
            stf=exp(-a.*(t-t0).^2);
         
        case 2
            t0=1.2/f0; %source time shift in seconds 
            a = pi*pi*f0*f0;
            stf = -2*a*(t-t0).*exp(-a.*(t-t0).^2);
        
        case 3
            %third derivative of a gaussian
%             da=pi*f0;
%             t0=1.5*sqrt(6.)/(pi*f0); %shifts in seconds 
%             a=pi*f0*(t-t0);
%             a2=(pi*f0*(t-t0)).^2;
%             ppx=-(0.5/(da.*da)).*exp(-a2);
%             px=(a/da).*exp(-a2);
%             stf=(1-2.*a2).*exp(-a2);
              tshift=20;   
             ts=1/f0;    
            tau = pi .* (t - 1.5 .* ts - tshift) ./ (1.5 * ts);
            stf = (1.0 - 4.0 .* tau .* tau) .*exp(-2.0 .* tau .* tau);
            case 4
                stf=zeros(1,length(t));
                %namestf=['./data_mase/STF/STF_',num2str(cof1),'-',num2str(cof2),'Hz_event',num2str(jf),catalog,'_dt0.01pratt.mat'];

        otherwise 
            error 'STF type not recognised' 
    end 
    stf=stf/max(stf); 
    stf=amp*stf; 
    
    [f,am]=check_frequency_spectrum(stf,dt); 
    if(plot_source==1);
    figure 
    subplot(211)
    plot(t,stf)
    hold on;
    title (['STF ', type])
    xlabel 't[s]'
    ylabel 'Source amp'
    grid on; 
    subplot(212)
    plot(f,am)
    hold on; 
    grid on; 
    title 'STF spectrum' 
    xlabel 'f[hz]'
    ylabel 'amp'
    end
end 