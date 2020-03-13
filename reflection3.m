function  [v_out, vpol_dir_in,vpol_dir_refl, hor_pol_dir]=reflection3(v_in,n);
% reflection of shear wave with incoming direction v_in, polarization direction pol, normal n, 
% computes direction and pol directions of the incoming and reflected outgoing shear wave 
% v_in, n and v_out must be normalised to length
v_in=v_in/norm(v_in);
n=n/norm(n);
if v_in'*n>0
   n=-n; % n must be directed towards incoming and reflected wave
end;

a=v_in-(v_in'*n)*n; % hor component of n
na=norm(a);
if na==0
   % incoming ray parallel to normal, i.e. perpendicular to surface
   v_out=v_in;
   % in this case the polarization direction of the reflected shear wave is undefined,
   % but the amplitude of this reflected wave is zero
   vert_pol_dir=[0;0;0]; % arbitrary value
   hor_pol_dir=[0;0;0];  % arbitrary value
else % incoming ray not perpendicular to surface
   %%na=sin(theta1)
   sinalpha_in= na ; 
   cosalpha_in=sqrt(1-sinalpha_in^2);
   v_out= -v_in + 2*a;
   vpol_dir_in  =  sinalpha_in* n  +cosalpha_in/na * a;
   vpol_dir_refl= -sinalpha_in* n  +cosalpha_in/na * a;
   hor_pol_dir= cross(a,n)/sinalpha_in;
end;




