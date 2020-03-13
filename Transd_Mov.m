TransPos = TransPosV2(); Focusd=14;
%  Transpos = TransP(); Focusd=12;  % for V1 transducer with 12 cm Focal Length
close all
xTrd = TransPos.x -Focusd; % unit cm
yTrd = TransPos.y;
zTrd = TransPos.z;
scatter3(xTrd,yTrd,zTrd, 40, 'red', 'filled');
axis([-20 20 -20 20 -20 20])
 hold on


% euler's approach

% R_alpha=[cos(alpha), sin(alpha),0;...
%          -sin(alpha), cos(alpha),0;...
%          0,0,1];
%      
% R_beta=[1, 0,0;...
%          0, cos(beta),sin(beta);...
%          0,-sin(beta),cos(beta)];
%  
% R_gamma=[cos(gamma),sin(gamma),0;...
%     -sin(gamma),cos(gamma),0;...
%     0,0,1];
% 
% New_val= zeros(256,3);
% for i=1:256
%     newX= R_gamma*R_beta*R_alpha*[xTrd(i); yTrd(i);zTrd(i)];
%     New_val(i,:)=newX;
% end
% scatter3(New_val(:,1),New_val(:,2),New_val(:,3), 40, 'green', 'filled');
%  

%rotation matrices 

alpha=0;
beta=pi/6; % yaxis rotation
gamma=0; %z axis rotation

R_alpha=[1,0,0;...
    0, cos(alpha), -sin(alpha);...
    0,sin(alpha), cos(alpha)];
     
R_beta=[cos(beta), 0,sin(beta);...
         0, 1,0;...
         -sin(beta),0,cos(beta)];
 
R_gamma=[cos(gamma),-sin(gamma),0;...
    sin(gamma),cos(gamma),0;...
    0,0,1];
xTrd=xTrd+2; %translation

New_val= zeros(256,3);
for i=1:256
    newX= R_gamma*R_beta*R_alpha*[xTrd(i); yTrd(i);zTrd(i)];
    New_val(i,:)=newX;
end
New_val(:,1)=New_val(:,1)+1;
scatter3(New_val(:,1),New_val(:,2),New_val(:,3), 40, 'green', 'filled');



focal_1=[1,0,0]; %translation
focal_2=R_gamma*R_beta*R_alpha*[focal_1(1); focal_1(2);focal_1(3)];%rotation
hold on
plot3(focal_2(1), focal_2(2), focal_2(3), 'o')
hold on
plot3(0,0,0, '*')
axis equal

figure()
V=[1,0,0];
Z=R_gamma*R_beta*R_alpha*[V(1); V(2);V(3)];
Z=Z/norm(Z);
quiver3(0,0,0,V(1),V(2),V(3), 'r');
hold on
quiver3(0,0,0,Z(1),Z(2),Z(3), 'r');
xlabel('x')
ylabel('y')
zlabel('z')


