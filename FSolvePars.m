function [p_1,p_2]=FSolvePars(rho, omega, speed ,alpha)
%
k=omega/speed; % wave number
C=omega^2*rho/(alpha^2+ k^2);
D=sqrt(2)*C/(speed*sqrt(rho));
p_1=D^2-C;
p_2=sqrt(C^2-p_1^2)/omega;
%
end