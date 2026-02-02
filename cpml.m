function [a_x,a_x_half,b_x,b_x_half,K_x,K_x_half,a_z,a_z_half,b_z,b_z_half,K_z,K_z_half]=cpml(ncpml,dh,nx,nz,f0,dt,vp_max);


USE_PML_XMIN=true;
USE_PML_XMAX=true; 
USE_PML_ZMIN=false; 
USE_PML_ZMAX=true; 


cp=7500;


K_MAX_PML=1;
NPOWER = 2.0;
ALPHA_MAX_PML=2.0*pi*(f0/2.0);
%ncpml=20; 
Rcoef = 0.001;

d_x=zeros(nx,1); 
d_x_half=zeros(nx,1);
K_x=ones(nx,1); 
K_x_half=ones(nx,1); 
alpha_x=zeros(nx,1); 
alpha_x_half=zeros(nx,1); 
a_x=zeros(nx,1); 
a_x_half=zeros(nx,1); 

d_z=zeros(nz,1); 
d_z_half=zeros(nz,1);
K_z=ones(nz,1); 
K_z_half=ones(nz,1); 
alpha_z=zeros(nz,1); 
alpha_z_half=zeros(nz,1); 
a_z=zeros(nz,1); 
a_z_half=zeros(nz,1); 


thickness_PML_x=ncpml*dh;
thickness_PML_z=ncpml*dh;

%! reflection coefficient (INRIA report section 6.1) http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  Rcoef = 0.001d0;

%! check that NPOWER is okay
 

%! compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  d0_x = - (NPOWER + 1) * cp * log(Rcoef) / (2.d0 * thickness_PML_x);
  d0_z = - (NPOWER + 1) * cp * log(Rcoef) / (2.d0 * thickness_PML_z);


%! damping in the X direction

%! origin of the PML layer (position of right edge minus thickness, in meters)
  xoriginleft = thickness_PML_x;
  xoriginright = (nx-1)*dh - thickness_PML_x;

  for i = 1:nx

%! abscissa of current grid point along the damping profile
    xval=dh*(i-1);

%!---------- left edge
    if (USE_PML_XMIN==1) 

%! define damping profile at the grid points
      abscissa_in_PML = xoriginleft - xval;
      if (abscissa_in_PML >= 0) 
        abscissa_normalized = abscissa_in_PML / thickness_PML_x;
        d_x(i) = d0_x * abscissa_normalized^NPOWER;
%! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized^NPOWER;
        alpha_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized);
      end

%! define damping profile at half the grid points
      abscissa_in_PML = xoriginleft - (xval + dh/2.d0);
      if (abscissa_in_PML >= 0) 
        abscissa_normalized = abscissa_in_PML / thickness_PML_x;
        d_x_half(i) = d0_x * abscissa_normalized^NPOWER;
%! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_x_half(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized^NPOWER;
        alpha_x_half(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized);
      end

      end

%!---------- right edge
    if (USE_PML_XMAX==1) 

%! define damping profile at the grid points
      abscissa_in_PML = xval - xoriginright;
      if (abscissa_in_PML >= 0) 
        abscissa_normalized = abscissa_in_PML / thickness_PML_x;
        d_x(i) = d0_x * abscissa_normalized^NPOWER;
%! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized^NPOWER;
        alpha_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized);
      end

%! define damping profile at half the grid points
      abscissa_in_PML = xval + dh/2.d0 - xoriginright;
      if (abscissa_in_PML >= 0) 
        abscissa_normalized = abscissa_in_PML / thickness_PML_x;
        d_x_half(i) = d0_x * abscissa_normalized^NPOWER;
%! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_x_half(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized^NPOWER;
        alpha_x_half(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized);
      end

    end

%! just in case, for -5 at the end
    if (alpha_x(i) < 0) 
        alpha_x(i) = 0;
    end 
    if (alpha_x_half(i) < 0) 
        alpha_x_half(i) = 0;
    end

    b_x(i) = exp(- (d_x(i) / K_x(i) + alpha_x(i)) * dt);
    b_x_half(i) = exp(- (d_x_half(i) / K_x_half(i) + alpha_x_half(i)) * dt);

%! this to avoid division by zero outside the PML
    if (abs(d_x(i)) > 1.d-6) 
        a_x(i) = d_x(i) * (b_x(i) - 1.d0) / (K_x(i) * (d_x(i) + K_x(i) * alpha_x(i)));
    end
    if (abs(d_x_half(i)) > 1.d-6) 
        a_x_half(i) = d_x_half(i)*(b_x_half(i) - 1.d0) / (K_x_half(i) * (d_x_half(i) + K_x_half(i) * alpha_x_half(i)));
    end
    end

%! damping in the Y direction

%! origin of the PML layer (position of right edge minus thickness, in meters)
   zoriginbottom = thickness_PML_z;
   zorigintop = nz*dh - thickness_PML_z;

  for j = 1:nz

%! abscissa of current grid point along the damping profile
    zval = dh*(j-1);

%!---------- top edge
    if (USE_PML_ZMIN==1) 

%! define damping profile at the grid points
      abscissa_in_PML = zoriginbottom - zval;
      if (abscissa_in_PML >= 0) 
        abscissa_normalized = abscissa_in_PML / thickness_PML_z;
        d_z(j) = d0_z * abscissa_normalized^NPOWER;
%! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_z(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized^NPOWER;
        alpha_z(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized);
      end

%! define damping profile at half the grid points
      abscissa_in_PML = zoriginbottom - (zval + dh/2.d0);
      if (abscissa_in_PML >= 0) 
        abscissa_normalized = abscissa_in_PML / thickness_PML_z;
        d_z_half(j) = d0_z * abscissa_normalized^NPOWER;
%! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_z_half(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized^NPOWER;
        alpha_z_half(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized);
      end

    end

%!---------- bottom edge
    if (USE_PML_ZMAX) 

%! define damping profile at the grid points
      abscissa_in_PML = zval - zorigintop;
      if (abscissa_in_PML >= 0) 
        abscissa_normalized = abscissa_in_PML / thickness_PML_z;
        d_z(j) = d0_z * abscissa_normalized^NPOWER;
%! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_z(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized^NPOWER;
        alpha_z(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized);
      end

%! define damping profile at half the grid points
      abscissa_in_PML = zval + dh/2.d0 - zorigintop;
      if (abscissa_in_PML >= 0)
        abscissa_normalized = abscissa_in_PML / thickness_PML_z;
        d_z_half(j) = d0_z * abscissa_normalized^NPOWER;
%! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_z_half(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized^NPOWER;
        alpha_z_half(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized);
      end

    end

    b_z(j) = exp(- (d_z(j) / K_z(j) + alpha_z(j)) * dt);
    b_z_half(j) = exp(- (d_z_half(j) / K_z_half(j) + alpha_z_half(j)) * dt);

%! this to avoid division by zero outside the PML
    if (abs(d_z(j)) > 1.d-6) 
        a_z(j) = d_z(j) * (b_z(j) - 1.d0) / (K_z(j) * (d_z(j) + K_z(j) * alpha_z(j)));
    end
    if (abs(d_z_half(j)) > 1.d-6) 
        a_z_half(j) = d_z_half(j)*(b_z_half(j) - 1.d0) / (K_z_half(j) * (d_z_half(j) + K_z_half(j) * alpha_z_half(j)));
    end

  end
  
  
end