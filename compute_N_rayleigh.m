function [Nmat]=compute_N_rayleigh(alpha,beta,mu,H,nlayer,om,p,z)
hola=zeros(4,4);
Nmat=complex(hola,0);
Emat=complex(hola,0);  
Player=complex(hola,0);
hola=zeros(1,nlayer);
eta_al=complex(hola,0);
eta_be=complex(hola,0); 
nu_al=complex(hola,0); 
nu_be=complex(hola,0);

for i=1:nlayer
   eta_al(i)=-complex(0,1)*sqrt(1/(alpha(i)^2)-p^2); %i vertical slowness 
   eta_be(i)=-complex(0,1)*sqrt(1/(beta(i)^2)-p^2); %i vertical slowness 
   nu_al(i)=om*eta_al(i);
   nu_be(i)=om*eta_be(i); %i*vertical wavenumber
   gamma(i)=2*beta(i)^2*p^2;
   gamma1(i)=1-1/gamma(i);
    
end

  %! note Emat is not omega dependent
  Emat(1,1) =  eta_be(nlayer)/p;
  Emat(1,2) = -Emat(1,1);
  Emat(1,3) = 1;
  Emat(1,4) = 1;
  Emat(2,1) = 1;
  Emat(2,2) = 1;
  Emat(2,3) =  eta_al(nlayer)/p;
  Emat(2,4) = -Emat(2,3);

  Emat(3,1) = 2*mu(nlayer)*gamma1(nlayer);
  Emat(3,2) = Emat(3,1);
  Emat(3,3) = 2*mu(nlayer)*eta_al(nlayer)/p;
  Emat(3,4) = -Emat(3,3);
  Emat(4,1) = 2*mu(nlayer)*eta_be(nlayer)/p;
  Emat(4,2) = -Emat(4,1);
  Emat(4,3) = Emat(3,1);
  Emat(4,4) = Emat(3,1);

  if(z>sum(H(1:nlayer-1)))
      error('z point is above free surface')
  end

  if (z<=0)
      Gmat=zeros(4,4); 
      Gmat(1,1)=exp(nu_be(nlayer)*z); 
      Gmat(2,2)=exp(-nu_be(nlayer)*z); 
      Gmat(3,3)=exp(nu_al(nlayer)*z); 
      Gmat(4,4)=exp(-nu_al(nlayer)*z); 
      Nmat=Emat*Gmat; 
  else
      hh=H; 
      ilayer=nlayer; 
      for j=nlayer-1:-1:1
          if(z<=sum(H(j:nlayer-1)))
              ilayer=j; 
              break
          end
      end 
      hh(ilayer+1:nlayer-1)=H(ilayer+1:nlayer-1); 
      hh(ilayer)=z-sum(H(ilayer+1:nlayer-1));
      if(hh(ilayer)<0)
          error('Thickness must be greater than 0')
      end 
      Nmat=Emat; 
      
      for j=nlayer-1:-1:ilayer
          c1=nu_al(j)*hh(j);
          ca=(exp(c1)+exp(-c1))/2; 
          sa=(exp(c1)-exp(-c1))/2;
          xa=eta_al(j)*sa/p; 
          ya=p*sa/eta_al(j);
          c2=nu_be(j)*hh(j);
          cb=(exp(c2)+exp(-c2))/2; 
          sb=(exp(c2)-exp(-c2))/2;
          xb=eta_be(j)*sb/p; 
          yb=p*sb/eta_be(j);
          g1=gamma1(j); 
          mul=mu(j) ;

          Player(1,1) = ca-g1*cb;
          Player(1,2) = xb-g1*ya;
          Player(1,3) = (ya-xb)/(2*mul);
          Player(1,4) = (cb-ca)/(2*mul);
          Player(2,1) = xa-g1*yb;
          Player(2,2) = cb-g1*ca;
          Player(2,3) = (ca-cb)/(2*mul);
          Player(2,4) = (yb-xa)/(2*mul);
          Player(3,1) = 2*mul*(xa-g1^2*yb);
          Player(3,2) = 2*mul*g1*(cb-ca);
          Player(3,3) = ca-g1*cb;
          Player(3,4) = g1*yb-xa;
          Player(4,1) = 2*mul*g1*(ca-cb);
          Player(4,2) = 2*mul*(xb-g1^2*ya);
          Player(4,3) = g1*ya-xb;
          Player(4,4) = cb-g1*ca;
          
          Nmat=gamma(j)*Player*Nmat; 
         
      end 
  end

end 