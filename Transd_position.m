function [ xt, yt, zt] = Transd_position( alpha, beta, gamma, xd,yd,zd )
%Transd_position : transd position with translation and rotation
% alpha not used! beta and gamma are the rotation keeping y and z fixed 
%respectevely. xd, yd, zd are the translation in x, y and z 
% first translation then rotation!

TransPos = TransPosV2(); Focusd=14;
%  Transpos = TransP(); Focusd=12;  % for V1 transducer with 12 cm Focal Length

xTrd = TransPos.x -Focusd; % unit cm
yTrd = TransPos.y;
zTrd = TransPos.z;

%rotation matrices

R_alpha=[1,0,0;...
    0, cos(alpha), -sin(alpha);...
    0,sin(alpha), cos(alpha)];

R_beta=[cos(beta), 0,sin(beta);...
    0, 1,0;...
    -sin(beta),0,cos(beta)];

R_gamma=[cos(gamma),-sin(gamma),0;...
    sin(gamma),cos(gamma),0;...
    0,0,1];

New_val= zeros(256,3);
for i=1:256
    newX= R_gamma*R_beta*R_alpha*[xTrd(i); yTrd(i);zTrd(i)];
    New_val(i,:)=newX;
end

xt= New_val(:,1)+xd;
yt= New_val(:,2)+yd;
zt= New_val(:,3)+zd;

end

