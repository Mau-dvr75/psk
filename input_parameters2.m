%--------------------------------------------------------------------------
%Model dimensions
%--------------------------------------------------------------------------
lx = 423*1000;  %length in x [m]
lz = 210*1000;  %depht in [m]
lb=30000; %buffer zone length (just one side) in [m].
dh=800; %Grid spacing [m]
%Number of grid points
nx=round(lx/dh)+1;
nz=round(lz/dh)+3+2;%+3 because of free surface,+2 is because we reduce the number of nodes in z dir
bz=round(lb/dh)*2+1;% %width of buffer zone *2, add this to psv->assing_fk :(
nx=nx+bz;
nz=nz;%+round(bz/2);
op=4; %FD spatial operator order
noa=2004; %number of array we are working with
%--------------------------------------------------------------------------
%Space axis configuration considering CPML
%--------------------------------------------------------------------------
ncpml=13;
%X axis
x0=0;
Nx=nx+2*ncpml;% Numero total de nodos en x, incluyendo nodos con absorbencia
x=(x0-ncpml*dh):dh:(x0+(Nx-ncpml-1)*dh); % malla con nodos CPML
nx0=ncpml+1; % posicion real de x0
nxf=Nx-ncpml; % posicion real de xf
%Z axis
z0=0;% coordenada inicial real sin CPML en z
%nz=nz;% Numero total de nodos en z sin absorbencia
Nz=nz+ncpml;% Numero total de nodos en z, incluyendo nodos con absorbencia
z=(z0-ncpml*dh):dh:(z0+(Nz-ncpml-1)*dh); % malla con nodos CPML
nz0=1; % posicion real de z0
nzf=Nz-ncpml;
%--------------------------------------------------------------------------
%Time discretization
%--------------------------------------------------------------------------
%*T/dt must be integer to store the fw field correctly
T=220.0;%Duration of simulation. Must be long enough to compute the FK fourier integral properly. 
dt=0.04; %dt for FD simulation
nt=round(T/dt);
sfe=10; %save the fw field every sfe time-steps fro gradients calculation
nt=sfe*round(nt/sfe);
dso=dt/0.01; %order of data downsampling. Original data is at 0.01;
FREE_SURF=true; %Currently dummy. Change this manually in the cpml  function
%--------------------------------------------------------------------------
%PLANE WAVE SOURCE SETUP
%--------------------------------------------------------------------------
plane_wave=true; %Plane wave sources?
INJECTIONS=1;
INJECTIONS_ORDER=2;
compute_initial_field=true; %when computing initial field, check the order (width of si)
field_order=2;
%resample_fk_fields=true;
resample_order=11; %the FK method do not require a small dt -> we can resample it
if(plane_wave==0)
    si=0;
else
    si=3;%si is for the injection surface. It has to be a little far away from cpml nodes
    %----------------------------------------------------------------------
    %FK background model definition
    %----------------------------------------------------------------------
    al_fk= [5.80 6.50 8.04 8.04 8.11 8.23 8.23]*1000; %in m/s
    be_fk= [3.46 3.85 4.48 4.49 4.50 4.51 4.51]*1000; %in m/s
    rho_fk=[2.72 2.92 3.33 3.35 3.38 3.41 3.41]*1000; %in km/m^3
    %Mu lame parameter for each layer
    mul_fk=[rho_fk.*be_fk.*be_fk];
    %Thickness of each layer
    % h_fk=[20000 15000 165000 0];% in m
    h_fk=[20 15 42.5 42.5 45 45 0]*1000;% in [m]
    nlayer_fk=length(h_fk);%Number of layers

    %Create boundary coordinates to compute analytical solution
    space_axis_fk; %zvx(3) must be total thickness ;)
    %--EARTHQUAKE DATA INFO--%
    catalog='best'; %earthquake catalog to read. Dummy 
    f0=[0.03 0.25];%corner frequencies of the BP filter. The last value is a dummy
    ft=0;%Filter type-> 0=nofilter 1=butt 2=FIR1
    fob=6; %butt order
    delay=1;% 356 para 0.03-0.1
    %P wace incidence angle in degrees, 0:vertical, 90:horizontal from
    %left, 360-IA=incidence from right
    angleforced=load(['./data_korea/stations/arrays/',num2str(noa),'/inc_angle.mat']);
    angleforced=angleforced.angleforced;
    SOURCE_TYPE=[1,1]; %1=P / 2=SV wavefront
    STF_TYPE=[4,1]; %1=Gaussian/ 2=1st derivative of gaussian/ 3=ricker/ 4=user's time-domainSTF/
    stf_method=[1,1]; %1=averaging, 2=deconvolution (Pratt)
    nsrc=length(SOURCE_TYPE);
    source_x_id=zeros(1,nsrc); %Dummy array
    source_z_id=zeros(1,nsrc); %Dummy array
end
%--------------------------------------------------------------------------
%RECEIVER CONFIGURATION
%--------------------------------------------------------------------------
all_stations_info=readtable(['./data_korea/stations/arrays/',num2str(noa),'/stations_proj.txt']);
stlo=all_stations_info.Var2;
stla=all_stations_info.Var3;
dgc=zeros(1,length(stla));
for i=1:length(stla)
    dgc(i)=distance('gc',[stla(1) stlo(1)],[stla(i) stlo(i)]);
    dgc(i)=deg2km(dgc(i)); %km
end
xrec=(nx0+si+round(bz/2))*dh+dgc*1000;
rdepth=0; %rec depth in meters
zrec=ones(1,length(xrec))*rdepth+3*dh-3*dh;
[rec_x_id, rec_z_id, nthings] = compute_indices(xrec,zrec,dh,Nx,Nz);
nrec=nthings;
%dltl=xns(rec_x_id(1)); dltr=xns(rec_x_id(end));%earthquakes from left %from right
dltl=xns(ibm); dltr=xns(iem);
%Source origin coordinates
x0_source=[dltl+210000 dltl]; %This controls the delay of incoming plane wavefield
z0_source=[0,0];%z0 is always 0, -x0 if wave from left,+x0 from right.
%--------------------------------------------------------------------------
%GRAPHIC SETUP
%--------------------------------------------------------------------------
IT_DISPLAY=100;
animation=true; %can't display animation if parfor is used
makemovie=false;
%--------------------------------------------------------------------------
%OUTPUT SEISMOGRAMS
%--------------------------------------------------------------------------
saveseis=true;
makeobs=true;%turn this on when computing Green Functions
analytical=true;
