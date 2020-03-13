function [refr, v_out]=refraction0(v_in,n,c_in,c_out);
% refraction with incoming direction v_in, normal n , velocities c_in and c_out
% v_in, n and v_out must be normalised to length
v_in=v_in/norm(v_in);
n=n/norm(n);
if v_in'*n<0
   n=-n; % n must be directed into out
end;
nbr=c_out/c_in;
a=v_in-(v_in'*n)*n; 
na=norm(a);
if na==0
   % incoming ray parallel to normal, i.e. perpendicular to surface
   refr=true; v_out=v_in;
else % incoming ray not perpendicular to surface
   % na=sin(theta1)
   sintheta2= nbr*na;  % Snellius: sintheta1/sintheta2= c_in/c_out
   if sintheta2<=1
       % refraction 
       costheta2=sqrt(1-sintheta2^2);
       refr=true; 
       v_out= costheta2*n+nbr*a;
   else % no refraction
       refr=false;
       v_out=[0;0;0];
   end;
end;

