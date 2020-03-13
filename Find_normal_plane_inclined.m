function normal = Find_normal_plane_inclined( beta, gamma)
%UNTITLED Gives back the normal of a plane giving the rotation angles beta and gamma 
%   Detailed explanation goes here
alpha=0; %always zero


R_alpha=[1,0,0;...
    0, cos(alpha), -sin(alpha);...
    0,sin(alpha), cos(alpha)];
     
R_beta=[cos(beta), 0,sin(beta);... %ydirection
         0, 1,0;...
         -sin(beta),0,cos(beta)];
 
R_gamma=[cos(gamma),-sin(gamma),0;... %zdirection
    sin(gamma),cos(gamma),0;...
    0,0,1];

V=[1,0,0];
normal=R_gamma*R_beta*R_alpha*[V(1); V(2);V(3)];
normal=normal/norm(normal);


end

