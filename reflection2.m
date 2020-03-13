function  [reflect_possible,v_out, pol_direction]=reflection2(v_in,n,c1,c2);
% reflection of wave with incoming direction v_in, normal n, speed c1
% outgoing wave with speed c2, pol_direction gives vertical polarization direction in case reflected wave is shear
% v_in, n and v_out must be normalised to length
v_in=v_in/norm(v_in);
n=n/norm(n);
if v_in'*n>0
   n=-n; % n must be directed towards incoming and reflected wave
end;
nbr=c2/c1;
a=v_in-(v_in'*n)*n; 
sinalpha_in=norm(a);
sinalpha_out=nbr*sinalpha_in;
if sinalpha_in==0
   % incoming ray parallel to normal, i.e. perpendicular to surface
   reflect_possible=true;
   v_out=v_in;
   % in this case the polarization direction of the reflected (shear) wave is undefined,
   % but the amplitude of this reflected wave is zero
    pol_direction=[0;0;0]; % arbitrary value
else % incoming ray not perpendicular to surface
   %%na=sin(theta1)
   %%sintheta2= nbr*na;  % Snellius: sintheta1/sintheta2= c1/c2
   if sinalpha_out>1
       reflect_possible=false;
       v_out=[0;0;0];
       pol_direction=[0;0;0]; 
   else
       reflect_possible=true;
       cosalpha_out=sqrt(1-sinalpha_out^2);
       v_out= cosalpha_out*n+nbr*a;
       pol_direction=-sinalpha_out * n + cosalpha_out/sinalpha_in * a;
   end;
end;





