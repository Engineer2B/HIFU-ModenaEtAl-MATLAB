function [son, ptrd1, ray_lossless] = Update_rays(theta, phi, startpoint,I_init , normal,v2norm,v3,material,phase_initial,object) 
%Update_rays gives back the path of the rays
%called Nrays times each transducer

% Each Ray has 
% 1. start point
% 2. Vray: the vector of the direction
% 3. material: 1=lossless , 2=fat, 3=muscle, 4=bone, 5=marrow
% 4. initial phase
% 5. direction=> 1=towards the focal point , -1 towards the transducer
% 6. IO the initial intensity
% 7. end point
% 8. IF the final intensity
% material 1= lossless 2=fat 3=muscle 4= bone 5=marrow 


%%% they should be moveddddd!!!! remembeeer!


pat=zeros(1,20);% a ray can go/bounch through a material up to 20 times
pat(1,1)=1;
Rtheta=tan(theta);
Vray=normal+Rtheta*(cos(phi)*v2norm+sin(phi)*v3);
ray_lossless.start=startpoint;
ray_lossless.Vray=Vray/norm(Vray);
B=material.oil.ka*sin(theta);
if B==0
    Scaling=1;
else    
    Scaling= (2*besselj(1,B)/B)^2;
end;    

ray_lossless.I0=I_init*Scaling;
ray_lossless.material=1;
ray_lossless.direction=1;
ptrd1=I_init*Scaling;

[inters_point]=Intersection(ray_lossless,object, 1); % call intersection
ray_lossless.end= inters_point.pt;
ray_lossless.phase_initial= phase_initial;
ray_lossless.phase_final=phase_initial+material.oil.k*inters_point.d;
ray_lossless.IF=ray_lossless.I0;
ray_lossless.path=pat;
ray_lossless.next_object=inters_point.index_obj;
ray_lossless.next_kind=object(inters_point.index_obj).kind;
ray_lossless.actual_object=1;
ray_lossless.nn=inters_point.normal;
ray_lossless.actual_kind=0;


%======================================================================
% Generated Rays from ray_lossless
%======================================================================

son = Reflection_Refraction(ray_lossless,object); %generation of the sons




end

