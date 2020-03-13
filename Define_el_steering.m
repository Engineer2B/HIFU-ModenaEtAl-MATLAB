function phase_matrix = Define_el_steering(xTrd, yTrd, zTrd, xd, yd, zd)
%DEFINE_EL_STEERING calculation of the E.steering for each trd element
%   Electronic steering is needed to:
% 1) create a sonication cell of 4 mm
% 2) Move the focal point
% 
% xd, yd, zd is the new location of the NEW focal point
% To fill phase matrix I need the new distances between the NEW focal point and
% and the trd elements (when the focal point is in 0,0,0 the distances are
% the same!) 
mat=Define_material();
distances= zeros(1,256);
phase_matrix=zeros(1,256);
% THIS IS NOT CORRECT, FOR A VERY SMALL SHIFT FROM THE FOCUS POINT, 
% THE PHASES SHOULD BE NEAR ZERO, WHICH THEY AREN'T.
if xd ~=0 || yd ~=0 || zd~=0 % if one of them is different from zero
    for i=1:256
        X = [xTrd(i),yTrd(i),zTrd(i);xd,yd,zd];
        distances(1,i)=pdist(X,'euclidean');%calculation of the new distances between the trd elements and the new focal point
        if i==1
            phase_matrix(1,i)=0;
        else
            phase_matrix(1,i)= mat.oil.k* distances(1,1) -  mat.oil.k*distances(1,i);
        end
    end
end
end