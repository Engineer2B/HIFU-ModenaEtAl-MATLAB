function [Refl_long,Refl_shear,Transm_long]= B2MReflTransmCoef_LongInp(alpha_in,rho_liquid,rho_solid,c_long_liquid,c_long_solid,c_shear_solid);
% BEREKENING REFLECTIE TRANSMISSIE COEFF VAN SOLID (LONGIT INPUT, LONGIT en SHEAR REFLECTIE)
% NAAR LIQUID (LONGIT OUTPUT) 
% MET MATRIX VERGEL OPLOSSEN VOLGENS AULD 
% Voorbeeld: Bone naar Muscle/Marrow of  Aluminium naar Water
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
%% alpha_long_crit= asin(c_long_solid/c_long_liquid); c_long_liquid/c_long_solid<0.5, dus geen kritieke hoek voor
%% beta_long, de gebroken long wave (in liquid) is er dus altijd
%% alpha_long_crit_degr=alpha_long_crit*180/pi;

%% alpha_shear_crit=asin(c_long_solid/c_shear_solid);  c_shear_solid/c_long_solid<1, dus geen kritieke hoek voor
%% beta_shear, de gereflecteerde shear wave (in solid) bestaat dus altijd
%% alpha_shear_crit_degr=alpha_shear_crit*180/pi;
%%fprintf('critical angle long = %7.3f critical angle shear= %7.3f \n',alpha_long_crit_degr,alpha_shear_crit_degr);

% compute beta_long: Snellius: sin(alpha_in)/sin(beta_long)=c_long_solid/c_long_liquid;
sinbeta_long=sin(alpha_in)*c_long_liquid/c_long_solid;
cosbeta_long=sqrt(1-sinbeta_long^2);  % pure imag if alpha_in>alpha_long_crit

% compute beta_shear: Snellius: sin(alpha_in)/sin(beta_shear)=c_long_solid/c_shear_solid;
sinbeta_shear=sin(alpha_in)*c_shear_solid/c_long_solid;
cosbeta_shear=sqrt(1-sinbeta_shear^2);  % pure imag if alpha_in>alpha_shear_crit


%% fprintf('alpha_in_degrees= %7.3f sinbeta_long=%7.3f  sinbeta_shear=%7.3f\n',alpha_in*180/pi,sinbeta_long,sinbeta_shear);
% refraction leads to longit wave, reflection to longit and shear wave
       
% matrix:   
N=[ cos(alpha_in), sinbeta_shear,   cosbeta_long;
   (lambda_solid+2*mu_solid*(cos(alpha_in))^2)*k_long_solid , 2*mu_solid*k_shear_solid*cosbeta_shear*sinbeta_shear,  -lambda_liquid*k_long_liquid;
   2*mu_solid*k_long_solid*cos(alpha_in)*sin(alpha_in) , mu_solid*k_shear_solid *(sinbeta_shear^2-cosbeta_shear^2) , 0         
  ];
% right hand side:
b=[cos(alpha_in);-k_long_solid*(lambda_solid+2*mu_solid*(cos(alpha_in))^2); 2*mu_solid*k_long_solid*cos(alpha_in)*sin(alpha_in)];

% oplossen amplitude  coeff:
xx=N^(-1)*b;
 
% power coeff:
% coeff zijn altijd real, er is geen critical angle      
Refl_long=abs(xx(1))^2; 
Refl_shear=abs(xx(2))^2  * c_shear_solid/c_long_solid* cosbeta_shear/cos(alpha_in);
Transm_long=abs(xx(3))^2  * c_long_liquid*rho_liquid/(c_long_solid*rho_solid) * cosbeta_long/cos(alpha_in);

 
