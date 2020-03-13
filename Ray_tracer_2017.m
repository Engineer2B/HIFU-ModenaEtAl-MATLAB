
%==========================================================================
% Ray_tracer_2017
% Version including interference in soft tissue and bone
%==========================================================================

%==========================================================================
% Definition general parameters
%==========================================================================
tic()
material= Define_material(); % definition of the materials
object=Object_Definition(); % Definition of the object
[ xmin, xmax, ymin,ymax, zmin, zmax, Nx, Ny, Nz, dx, dy, dz, xx, yy, zz, xxb,yyb, zzb ] = Define_table(); %definition table
RequiredTotalPower = 25; % Watt
[xTrd,yTrd,zTrd] = Transd_position( 0, 0, 0, 0,0,0 ); %definition of trd position-> Transd_position( alpha, beta, gamma, xd,yd,zd )
Ntrd = length(xTrd); % nr of transducers
Nrays = 20000; % nr of rays PER transducer element

% first 3 zeros of besselj(1,x)=J1(x):
r1 = 3.8317;
% SELECT LOBES TO INCLUDE:
% sin_theta_max=r1/ka;  % then:  k*a*sin(theta_max)=r1
sin_theta_max=r1/(material.oil.ka*4); % smaller, only for shape of focal point
theta_max=asin(sin_theta_max);
% compute the fraction of the power that an individual transducer element
% radiates between 0 and theta_max
f=@(theta) sin(theta).*(2*besselj(1,material.oil.ka*sin(theta))./(material.oil.ka*sin(theta))).^2;
% QUESTION: Why start at 1e-8 and not 0?
% QUESTION: Why use quad and not integral?
IntTot=quad(f,0.00000001,pi/2); % integral over theta from 0 to pi/2
IntRestr=quad(f,0.00000001,theta_max);
powerfr=IntRestr/IntTot;

% Deefinition of the parameters for the case in which sonication cell is 4
% QUESTION: Why is this needed?
Power_total_son_cell_4=0;
PowerLoss_son_cell4=zeros(Nx,Ny,Nz); %power loss in soft tissue for sonication cell of 4 mm
PowerLossL_son_cell4=zeros(Nx,Ny,Nz);%power loss in bone(longitudinal) for sonication cell of 4 mm
PowerLossS_son_cell4=zeros(Nx,Ny,Nz);%power loss in bone(shear) for sonication cell of 4 mm

son_cell_size=4; %sonication cell (2 or 4)

%==========================================================================
% Computation of the power produced
%==========================================================================
%
% For a 2mm sonication cell if you want to use the parallel pool, please go
% to RayTracer.m 
%
switch son_cell_size
    case 2 %sonication cell of 2 mm -> no electronic steering applied
    xd=0; %change xd,yd, zd if you want to apply electronic steering is applied
    yd=0;
    zd=0;
    
    phase_matrix=Define_el_steering(xTrd,yTrd,zTrd, xd, yd, zd); %matrix in which the starting phases for each trd are collected
    % RayTracer is the main program which calls ->  CalculationPressure
    % (256 times)
    % which calls -> UpdateRays ( each trd: Nrays Times)-> ProcessRays (256 times)
    [PowerLoss,PowerLossL,PowerLossS,Ptot,Q]= RayTracer(material,object, xTrd,yTrd,zTrd, Ntrd,Nrays, phase_matrix);
    
    %================================================================================
    % Corrections!
    %
    % RequiredTotalPower is the total acoustic power for the transducer
    % RequiredTotalPower*powerfraction is the part of the total power, emitted
    % within theta_max, Ptot is te actual emitted power by the rays within
    % theta_max. Therefore we scale by PowerCorrection, and divide by dx*dy*dz to
    % obtain densities
    %================================================================================
    PowerCorrectionFactor=powerfr*RequiredTotalPower/Ptot;
    PowerLoss1=PowerLoss*(PowerCorrectionFactor);
    Intensity_bone=PowerLossL+PowerLossS;
    Intensity_bone=Intensity_bone*(PowerCorrectionFactor);%/dx*dy*dz;
    PowerLoss= Intensity_bone+ PowerLoss1;
    Q=Q*(PowerCorrectionFactor*1e-6);
    case 4 %application of electronic steering. the system focus in eight different points defined by pt.
        tic()
        pt=[0.2,0; sqrt(0.2),sqrt(0.2);0,0.2;-sqrt(0.2),sqrt(0.2);-0.2,0;-sqrt(0.2),-sqrt(0.2);0,-0.2;sqrt(0.2),-sqrt(0.2)];
        xd=0; %change xd if the electronic steering is applied in the x-direction
       
%         myCluster = parcluster('local') %parallel pool -> change the for in a parfor if you want to use it
%         myCluster.NumWorkers = 12
%         parpool('local',12)

        tic()
         %summing up the heat production for each point pt 
        for j=1:length(pt)
            phase_matrix=Define_el_steering(xTrd,yTrd,zTrd, xd,pt(j,1), pt(j,2));
            [PowerLoss_m,PowerLossL,PowerLossS,Ptot]= RayTracer(material,object, xTrd,yTrd,zTrd, Ntrd,Nrays, phase_matrix);
            Power_total_son_cell_4=Power_total_son_cell_4+Ptot;
            PowerLoss_son_cell4=PowerLoss_son_cell4+PowerLoss_m;
            PowerLossL_son_cell4=PowerLossL+PowerLossL_son_cell4;
            PowerLossS_son_cell4=PowerLossS+PowerLossS_son_cell4;
        end
        t1=toc();
        disp('time to compute the power production for the 4 mm son cell')
        disp(toc())
        PowerCorrectionFactor=powerfr*RequiredTotalPower/Power_total_son_cell_4;
        PowerLoss1=PowerLoss_son_cell4*PowerCorrectionFactor;
        Intensity_bone= PowerLossL_son_cell4 + PowerLossS_son_cell4;
        Intensity_bone= Intensity_bone * PowerCorrectionFactor;%/dx*dy*dz;
        PowerLoss= Intensity_bone+PowerLoss1;
    
end


%================================================================================
% Visualization
%================================================================================

% plot for fixed z=z0
z0=0;
[zmin,iz]=min(abs(zz-z0));
z1=zz(iz);
fprintf('nearest occuring height value z=%7.3f\n',z1);
Powerz=squeeze(PowerLoss(:,:,iz));
[x3,y3]=meshgrid(xx,yy);
figure
surfc(x3,y3,Powerz');
xlabel('x (cm)');
ylabel('y (cm)');
zlabel('Q (W/cm^3)');
st2=['heat production old at z=',num2str(z1),' in watt/cm^3'];
title(st2);

% % plot for fixed x=x0
% x0=-0.75;
% [xmin,ix]=min(abs(xx-x0));
% x1=xx(ix);
% fprintf('nearest occuring x value x=%7.3f\n',x1);
% Powerx=squeeze(PowerLoss(ix,:,:));
% [x4,y4]=meshgrid(yy,zz);
% figure
% surfc(x4,y4,Powerx');
% xlabel('y (cm)');
% ylabel('z (cm)');
% zlabel('Q (W/cm^3)');
% st3=['heat production at x=',num2str(x1),' in watt/cm^3'];
% title(st3);

% plot for fixed z=z0
z0=0;
[zmin,iz]=min(abs(zz-z0));
z1=zz(iz);
fprintf('nearest occuring height value z=%7.3f\n',z1);
Powerz=squeeze(Q(:,:,iz));
[x3,y3]=meshgrid(xx,yy);
figure
surfc(x3,y3,Powerz');
xlabel('x (cm)');
ylabel('y (cm)');
zlabel('Q (W/cm^3)');
st2=['heat production new at z=',num2str(z1),' in watt/cm^3'];
title(st2);



