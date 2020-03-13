function [Refl_long,Refl_shear,Transm_long]= B2MReflTransmCoef_ShearInp(alpha_in,rho_liquid,rho_solid,c_long_liquid,c_long_solid,c_shear_solid);
% BEREKENING REFLECTIE TRANSMISSIE COEFF VAN SOLID (SHEAR INPUT VERTICAL POLARIZATION, LONGIT en SHEAR REFLECTIE)
% NAAR LIQUID (LONGIT OUTPUT) 
% MET MATRIX VERGEL OPLOSSEN VOLGENS AULD 
% LEVERT  AMPLITUDE REFL/TRANSMISSIE COEFF (PLOT 1)
% EN POWER REFL EM TRANSMISSIE COEFF (PLOT 2)
% versie volgens Auld
% Voorbeeld: solid naar liquid/Marrow of  Aluminium naar Water
% alpha_in in radialen !!

% Lamé pars solid:
mu_solid=rho_solid*c_shear_solid^2;
lambda_solid=rho_solid*c_long_solid^2-2*mu_solid;

% Lamé pars liquid:
lambda_liquid=rho_liquid*c_long_liquid^2;

% wave numbers: k c= omega=2 pi f, dus k=2 pi f/c
%%f=10^4; %% f= 1.4 MHz
% value of f arbitrary, select f value such that elements in N matrix of approx same size:
f=1e-4; 
k_long_liquid=2*pi*f/c_long_liquid;
k_long_solid=2*pi*f/c_long_solid;
k_shear_solid=2*pi*f/c_shear_solid;


% compute critical incoming angles:
%alpha_long1_crit= asin(c_shear_solid/c_long_solid); % kritieke hoek voor longit reflectie naar solid terug
%% alpha_long1_crit_degr=alpha_long1_crit*180/pi;

%  alpha_long2_crit= asin(c_shear_solid/c_long_liquid);
%  c_shear_solid/c_long_liquid>1 dus beta_l (long transmissie) < alpha_in (shear incoming)
%  dus geen kritieke hoek voor longit transmissie naar liquid
%  alpha_long2_crit_degr=alpha_long_crit*180/pi;

%fprintf('critical angle shear in -> long reflect = %7.3f \n',alpha_long1_crit_degr);



% compute alpha_long (reflectie longit): Snellius: sin(alpha_in)/sin(alpha_long)=c_shear_solid/c_long_solid;
sinalpha_long=sin(alpha_in)*c_long_solid/c_shear_solid;
cosalpha_long=sqrt(1-sinalpha_long^2);  % pure imag if alpha_in>alpha_long1_crit

% compute beta_long  (refractie longit): Snellius: sin(alpha_in)/sin(beta_long)=c_shear_solid/c_long_liquid;
sinbeta_long=sin(alpha_in)*c_long_liquid/c_shear_solid;
cosbeta_long=sqrt(1-sinbeta_long^2);  % always real

% refraction leads to longit wave, reflection to longit and shear(vertically polarized) wave
% matrix:   
N=[ sin(alpha_in), cosalpha_long,   cosbeta_long;
    2*mu_solid*k_shear_solid*sin(alpha_in)*cos(alpha_in) ,2*mu_solid*k_long_solid*cosalpha_long^2+lambda_solid*k_long_solid,-lambda_liquid*k_long_liquid;
     k_shear_solid*((sin(alpha_in))^2-(cos(alpha_in))^2), 2*k_long_solid*sinalpha_long*cosalpha_long,0
 ];
% right hand side:
b=[-sin(alpha_in);2*mu_solid*k_shear_solid*sin(alpha_in)*cos(alpha_in)           ; -k_shear_solid*((sin(alpha_in))^2-(cos(alpha_in))^2)];

% oplossen:
xx=N^(-1)*b;
% power coeff:
Refl_shear=real(abs(xx(1))^2); 
Refl_long=real(abs(xx(2))^2  * c_long_solid/c_shear_solid* cosalpha_long/cos(alpha_in));
Transm_long=real(abs(xx(3))^2  * c_long_liquid*rho_liquid/(c_shear_solid*rho_solid) * cosbeta_long/cos(alpha_in));
% klaar

    
      
     



