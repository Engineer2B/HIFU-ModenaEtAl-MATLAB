function [T, R]=CoefficientLL(v_in,n,c1,c2, rho1,rho2)
% refraction with incoming direction v_in, normal n , velocities c_in and c_out
% v_in, n and v_out must be normalised to length
%I checked and it seems correct!
v_in=v_in/norm(v_in);
n=n/norm(n);

if v_in'*n<0
   n=-n; % n must be directed into out
end;
% nbr=c2/c1;
% nb1=c1/c2;
% a=v_in-(v_in'*n)*n; 
% na=norm(a);
% sintheta2= nbr*na; 
% sintheta1=nb1*sintheta2;
% costheta_t=sqrt(1-sintheta2^2);
% costheta_i=sqrt(1-sintheta1^2);

costheta_i= v_in'*n;
%theta_i= acos(costheta_i);
sintheta_i=sqrt(1-costheta_i^2);
sintheta_t= (c2/c1)*sintheta_i;
costheta_t=sqrt(1-sintheta_t^2);



T= (4*rho1*c1*rho2*c2*costheta_t*costheta_i)/(rho2* c2*costheta_i + rho1* c1*costheta_t)^2;
R= ((rho2*c2*costheta_i-rho1*c1*costheta_t)^2)/(rho2* c2*costheta_i + rho1* c1*costheta_t)^2;
end

