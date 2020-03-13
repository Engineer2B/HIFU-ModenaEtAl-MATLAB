function [Refl, Transm_long,Transm_shear]= M2BReflTransmCoef(alpha_in,rho_liquid,rho_solid,c_long_liquid,c_long_solid,c_shear_solid);
% BEREKENING REFLECTIE TRANSMISSIE COEF VAN LIQUID (LONGIT INPUT)  NAAR SOLID (LONG EN SHEAR OUTPUT) 
% alpha_in in radialen !!
% Voorbeeld: Muscle  naar Bone of Water naar Aluminium


% Lamé params solid:
mu_solid=rho_solid*c_shear_solid^2;
lambda_solid=rho_solid*c_long_solid^2-2*mu_solid;

% Lamé params liquid:
lambda_liquid=rho_liquid*c_long_liquid^2;

% compute critical angles:
alpha_long_crit =asin(c_long_liquid/c_long_solid);
alpha_shear_crit=asin(c_long_liquid/c_shear_solid);
alpha_long_crit_degr=alpha_long_crit*180/pi;
alpha_shear_crit_degr=alpha_shear_crit*180/pi;
% fprintf('critical angle long = %7.3f critical angle shear= %7.3f \n',alpha_long_crit_degr,alpha_shear_crit_degr);



% compute beta_long: Snellius: sin(alpha_in)/sin(beta_long)=c_long_liquid/c_long_solid;
sinbeta_long=sin(alpha_in)*c_long_solid/c_long_liquid;
cosbeta_long=sqrt(1-sinbeta_long^2);  % pure imag if alpha_in>alpha_long_crit

% compute beta_shear: Snellius: sin(alpha_in)/sin(beta_shear)=c_long_liquid/c_shear_solid;
sinbeta_shear=sin(alpha_in)*c_shear_solid/c_long_liquid;
cosbeta_shear=sqrt(1-sinbeta_shear^2);  % pure imag if alpha_in>alpha_shear_crit

 

%FORMULE VERSIE;
Z1=rho_liquid*c_long_liquid/cos(alpha_in);
ZL=rho_solid*c_long_solid/cosbeta_long;
Zs=rho_solid*c_shear_solid/cosbeta_shear;
Zeff=ZL*(cosbeta_shear^2-sinbeta_shear^2)^2 + Zs*(2*sinbeta_shear*cosbeta_shear)^2;
 

% REFLECTIE:     
Rdf=(Zeff-Z1)/(Zeff+Z1);
% dit is velocity potential refl coef, maar die is gelijk aan de uitwijkings of snelheids transm coefficient.

% nu power refl coeff berekenen ( zelfde impedanties vooe en na, en zelfde hoeken):
Refl=abs(Rdf)^2;
     

% TRANSMISSIE, LONGITUDINAL: 
Tdf=(rho_liquid/rho_solid)*2*ZL*(cosbeta_shear^2-sinbeta_shear^2)/(Zeff+Z1);
% correctie, formule hierboven geeft de velocity potential longit transm  coefficients,
% de matrix berekening volgwens Auld geeft de uitwijking of (gelijk!) snelheids refl en transm coefficients
% Tdf,velocity(Auld) = Tdf, vel potential * || k_shear_solid|| /||k_long_liquid|| 
% = Tdf, vel potential * c_long_liquid/c_long_solid

Tdf=Tdf* c_long_liquid/c_long_solid;  
%% was Tdf=Tdf* k_long_solid/k_long_liquid;  
% nu is Tdf de uitwijkings of snelheids transm coefficient, zoals in Auld

% nu power transmissie longit coeff berekenen ( incl correctie impedanties en hoeken):  
Tdpow=abs(Tdf)^2*(rho_solid*c_long_solid)/(rho_liquid*c_long_liquid)*cosbeta_long/cos(alpha_in);
% Merk op: Tdpow is reëel onder de critical angle, daaarboven is cosbeta_long zuiver imaginair en dus
% is ook Tdpow zuiver imaginair. Er is dan ook geen power transmissie in de long mode, dus we kunnen het reëele
% deel van Tdpow pakken:
Transm_long = real(Tdpow);


% TRANSMISSIE, SHEAR: 
Tsf=-(rho_liquid/rho_solid)*2*Zs*(2*sinbeta_shear*cosbeta_shear)/(Zeff+Z1);
% correctie, formule hierboven geeft de velocity potential shear transm  coefficients,
Tsf=Tsf* c_long_liquid/c_shear_solid;  
%% was: Tsf=Tsf*k_shear_solid/k_long_liquid;
% nu is Tdf de uitwijkings of snelheids transm coefficient, zoals in Auld

% nu power transmissie longit coeff berekenen ( incl correctie impedanties en hoeken):  
Tspow=abs(Tsf)^2 *(rho_solid*c_shear_solid)/(rho_liquid*c_long_liquid)*cosbeta_shear/cos(alpha_in);
% Merk op: Tspow is reëel onder de critical angle, daaarboven is cosbeta_shear zuiver imaginair en dus
% is ook Tspow zuiver imaginair. Er is dan ook geen power transmissie in de long mode, dus we kunnen het reëele
% deel van Tspow pakken:
Transm_shear = real(Tspow); 


% correctie, formules geven de velocity potential reflection and transmission coefficients,
% de matrix berekening volgwens Auld geeft de uitwijking of (gelijk!) snelheids refl en transm coefficients
% Ts,velocity(Auld) = Ts, vel potential * || k_shear_solid|| /||k_long_liquid||


     
% Merk op: Tspow is reëel onder de critical angle, daaarboven is cosbeta_shear zuiver imaginair en dus
% is ook Tspow zuiver imaginair. Er is dan ook geen power transmissie in de shear mode, dus we kunnen het reëele
% deel van Tspow pakken:
Transm_shear=real(Tspow);



