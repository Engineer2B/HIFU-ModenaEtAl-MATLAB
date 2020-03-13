function [Pressure, p_bone_longT,  Vel_bone_s_x, Vel_bone_s_y,...
	Vel_bone_s_z,vl1,vl2,vl3,eps12,eps13, eps23]= Process_rays(Raylist, object)
% Process_rays Process all the rays from on trd
% called 256 times
% Each group of rays has the same phase!
%

[ xmin, xmax, ymin,ymax, zmin, zmax, Nx, Ny, Nz, dx, dy, dz,...
	xx, yy, zz, xxb,yyb, zzb ] = Define_table();
material=Define_material();
mem={Raylist.path}';%save in cell the path of the Rays
pat=cell2mat(mem);%convert in a matrix


groups=unique(pat, 'rows');%find the groups
[sizegroups,~]=size(groups); %how many groups
Vel_bone_s_x=zeros(Nx,Ny,Nz);%x comp of the velocity due to shear waves
Vel_bone_s_y=zeros(Nx,Ny,Nz);%y comp of the velocity due to shear waves
Vel_bone_s_z=zeros(Nx,Ny,Nz);%z comp of the velocity due to shear waves

Pressure=zeros(Nx,Ny,Nz);%Pressure in soft
p_bone_longT=zeros(Nx,Ny,Nz); %give back the intesity matrix for the bone

%full int
vl1=zeros(Nx,Ny,Nz);%x comp of the velocity due to shear waves
vl2=zeros(Nx,Ny,Nz);%y comp of the velocity due to shear waves
vl3=zeros(Nx,Ny,Nz);%z comp of the velocity due to shear waves
eps12=zeros(Nx,Ny,Nz);%x comp of the velocity due to shear waves
eps13=zeros(Nx,Ny,Nz);%y comp of the velocity due to shear waves
eps23=zeros(Nx,Ny,Nz);%z comp of the velocity due to shear waves
% Definition of the parameters here 
% TODO: change and connect with define_material
fr=1.2e6;
omega=2*pi*fr; 
cS= 1995; % m/s, speed of shear waves in bone
cL= 3736; % m/s, speed of longit waves in bone
rho_s=2025; % kg/m3, density of bone
alphaS=280; % /m, attenuation of shear waves in bone, in Nepers
alphaL=190; % /m, attenuation of longit waves in bone, in Nepers
kS=omega/cS;
kL=omega/cL;
% v1=zeros(Nx,Ny,Nz); %Intensity in soft
% v2=zeros(Nx,Ny,Nz); %Phase in soft
% v3=zeros(Nx,Ny,Nz);%Number of rays in soft in each cube


[mu,eta]=FSolvePars(rho_s,omega,cS,alphaS);
% next: compute p1=lambda+2 mu and p2=xi+4 eta/3:
[p1,p2]=FSolvePars(rho_s,omega,cL,alphaL);
lambda=p1-2*mu;
xi=p2-4*eta/3;
% Define other constants

% for solid (bone), longitudinal:
Const2=omega*((lambda+2*mu)*kL+(xi+4*eta/3)*omega*alphaL)/2;
Const3=1/Const2;
Const4=(-1i*omega)*(1i*kL-alphaL); % complex constant
% for solid (bone), shear:
Const5=(mu*omega*kS+eta*omega^2*alphaS)/2;
Const6=1/Const5;
Const7=(-1i*omega)*(1i*kS-alphaS);% complex constant


for ct=1:sizegroups %process the rays in the same group
    
    tf=(ismember(pat,groups(ct,:),'rows')); %find the indices of the rows which are equal to the groups
    group_now=Raylist(tf); %tf is a logical vector!
    %i call intensity with group_now
    [Intensity,PhaseOneEl,NraysOneEl,int_bone_long, int_bone_shear,...
			PhaseLong,PhaseShear,Nlong, Nshear,polx,poly,polz, kvs1,kvs2,kvs3, kvl1,kvl2,kvl3] = ...
        Calc_intensity(group_now, object);

    NraysOneEl_tot=max(NraysOneEl,1); %if no rays pass through the cube. To avoid dicision by 0
    PhaseOneEl_t=PhaseOneEl./NraysOneEl_tot; %average phase
    
    NraysOneEl_l=max(Nlong,1);%if no rays pass through the cube. To avoid dicision by 0
    PhaseOneEl_l=PhaseLong./ NraysOneEl_l;%average phase
    
    NraysOneEl_s=max(Nshear,1);%if no rays pass through the cube. To avoid dicision by 0
    PhaseOneEl_s=PhaseShear./NraysOneEl_s;%average phase
 
    
    p_bone_long=sqrt(int_bone_long).*exp(1i*PhaseOneEl_l);%pressure bone due to long
    Vel_bone_shear_x=sqrt(int_bone_shear).*exp(1i*PhaseOneEl_s).*polx;%calculation of the velocities due to shear waves
    Vel_bone_shear_y=sqrt(int_bone_shear).*exp(1i*PhaseOneEl_s).*poly;%calculation of the velocities due to shear waves
    Vel_bone_shear_z=sqrt(int_bone_shear).*exp(1i*PhaseOneEl_s).*polz;%calculation of the velocities due to shear waves
    Pressure_one_group=sqrt(Intensity*2*material.muscle.z).*exp(1i*PhaseOneEl_t); %pressure muscle
    
    Pressure=Pressure+Pressure_one_group; %Sum up pressure in soft of different groups (that creates interference)
    
    p_bone_longT=p_bone_longT+p_bone_long;%Sum up pressure in bone of different groups (that creates interference)
    
    Vel_bone_s_x=Vel_bone_s_x+Vel_bone_shear_x;
    Vel_bone_s_y=Vel_bone_s_y+Vel_bone_shear_y;%Sum up velocity in bone of different groups (that creates interference)
    Vel_bone_s_z=Vel_bone_s_z+Vel_bone_shear_z;
    
    % new method
    B1=sqrt(Const6*(int_bone_shear.*1e4)); % amplitude shear ray in bone
    A1=sqrt(Const3*(int_bone_long.*1e4)); % amplitude longit ray in bone
    %A1=0;
    Is=int_bone_shear;
    Il=int_bone_long;
%     vh1=vh1-1i*omega*B1.*exp(1i*PhaseOneEl_s).*polx-1i*omega*A1.*exp(1i*PhaseOneEl_l).*kvl1;
%     vh2=vh2-1i*omega*B1.*exp(1i*PhaseOneEl_s).*poly-1i*omega*A1.*exp(1i*PhaseOneEl_l).*kvl2;
%     vh3=vh3-1i*omega*B1.*exp(1i*PhaseOneEl_s).*polz-1i*omega*A1.*exp(1i*PhaseOneEl_l).*kvl3;
    
    vl1=vl1+A1*Const4.*exp(1i*PhaseOneEl_l).*kvl1.^2+B1*Const7.*exp(1i*PhaseOneEl_s).*kvs1.*polx;
    vl2=vl2+A1*Const4.*exp(1i*PhaseOneEl_l).*kvl2.^2+B1*Const7.*exp(1i*PhaseOneEl_s).*kvs2.*poly;
    vl3=vl3+A1*Const4.*exp(1i*PhaseOneEl_l).*kvl3.^2+B1*Const7.*exp(1i*PhaseOneEl_s).*kvs3.*polz;
    eps12=eps12+A1*Const4.*exp(1i*PhaseOneEl_l).*2.*kvl1.*kvl2+B1.*Const7.*exp(1i*PhaseOneEl_s).*(kvs1.*poly+kvs2.*polx);
    eps13=eps13+A1*Const4.*exp(1i*PhaseOneEl_l).*2.*kvl1.*kvl3+B1.*Const7.*exp(1i*PhaseOneEl_s).*(kvs1.*polz+kvs3.*polx);
    eps23=eps23+A1*Const4.*exp(1i*PhaseOneEl_l).*2.*kvl2.*kvl3+B1.*Const7.*exp(1i*PhaseOneEl_s).*(kvs2.*polz+kvs3.*poly);

end

