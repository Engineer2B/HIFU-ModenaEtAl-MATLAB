function [nn, lambdaplane]= check_sign(dir, plane_normal)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

denominator=dir'*plane_normal;
% if abs(denominator)<1e-12
%     % ray is parallel to plane, no intersection
%     lambdaplane=1e12;
% else
%     lambdaplane= (plane_constant-start'*plane_normal)/denominator;
% end;
if denominator>0
    nn=-plane_normal;
else
    nn=plane_normal;
end;

end

