function [PowerLoss,PowerLossL,PowerLossS,Ptot,Q]= RayTracer(material,object, xTrd,yTrd,zTrd,Ntrd, Nrays, phase_matrix)
%==========================================================================
% MAIN PROGRAM -> calls 256 times Calculation_Pressure
%==========================================================================

[ xmin, xmax, ymin,ymax, zmin, zmax, Nx, Ny, Nz, dx, dy, dz, xx, yy, zz, xxb,yyb, zzb ] = Define_table();


Ptot=0;
Pressure=zeros(Nx,Ny,Nz); %Pressure in SOFT tissue
Pressure_bone_long=zeros(Nx,Ny,Nz);%Pressure in bone due to longitudinal waves

vx_shear=zeros(Nx,Ny,Nz); % x component velocity due to shear waves
vy_shear=zeros(Nx,Ny,Nz); % y component velocity due to shear waves
vz_shear=zeros(Nx,Ny,Nz); % z component velocity due to shear waves
 
%full interference

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


%
% If you want to work in parallel uncomment and change the for to a parfor
%
% myCluster= parcluster('local') %UNCOMMENT FOR PARALLEL POOL
% myCluster.NumWorkers=12        %UNCOMMENT FOR PARALLEL POOL
% parpool('local', 12)           %UNCOMMENT FOR PARALLEL POOL

Ntrd=1;
for trd = 1:Ntrd % loop over all transducer elements
    [Pressure1, Ptot1,Int_bone_long,  velx,vely,velz,vl1_1,vl2_1,vl3_1,eps12_1,eps13_1, eps23_1]=CalculationPressure(Nrays,trd, phase_matrix, object,xTrd,yTrd,zTrd,material);
    Pressure=Pressure1+Pressure; 
    Ptot=Ptot1+Ptot; %sum the power produced total
    Pressure_bone_long=Pressure_bone_long+Int_bone_long; %pressure in bone due to long waves
    vx_shear=vx_shear+velx; %VELOCITY due to shear waves propagation
    vy_shear=vy_shear+vely;
    vz_shear=vz_shear+velz;
    
    %new full int
    vl1=vl1+vl1_1;%x comp of the velocity due to shear waves
    vl2=vl2+vl2_1;%y comp of the velocity due to shear waves
    vl3=vl3+vl3_1;%z comp of the velocity due to shear waves
    eps12=eps12+eps12_1;%x comp of the velocity due to shear waves
    eps13=eps13+eps13_1;%y comp of the velocity due to shear waves
    eps23=eps23+eps23_1;%z comp of the velocity due to shear waves
    
    if mod(trd,20)==0
        fprintf('trans el %3d of %3d\n', trd,Ntrd);
    end;
end
% PowerLoss:  PowerLoss in soft tissue
% PowerLossL: PowerLoss in bone due to long waves
% PowerLossS: PowerLoss in bone due to shear waves

PowerLoss= 2* material.muscle.attenuation*((abs(Pressure.^2))/(2*material.muscle.z)); %due to soft tissue
PowerLossL=2* material.bone.attenuationl*((abs(Pressure_bone_long.^2))); %due to long

Ix_shear=abs(vx_shear).^2; %From velocity to Intensity 
Iy_shear=abs(vy_shear).^2;
Iz_shear=abs(vz_shear).^2;

I_tot_shear=Ix_shear+Iy_shear+Iz_shear;
% V_tot=sqrt(V_tot);
% V_tot=Vtot.^2;
PowerLossS=2* material.bone.attenuations*I_tot_shear;

phiv1=angle(vl1);
phiv2=angle(vl2);
phiv3=angle(vl3);

Q=(xi+4*eta/3).*(abs(vl1).^2+abs(vl2).^2+abs(vl3).^2)/2 + ...
eta* (abs(eps12).^2 + abs(eps13).^2 +abs(eps23).^2)/2+...
(xi-2*eta/3) *(abs(vl1).*abs(vl2).*cos(phiv1-phiv2)+ ...
abs(vl1).*abs(vl3).*cos(phiv1-phiv3)+abs(vl2).*abs(vl3).*cos(phiv2-phiv3));

























