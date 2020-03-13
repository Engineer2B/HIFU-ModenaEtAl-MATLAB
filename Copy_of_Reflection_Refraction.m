function son  = Reflection_Refraction( ray, configuration, object)
%Reflection_Reflaction Calculates the rays generated from a ray
%  Different Cases will be trated
son(3).xi=0;
mat=Define_material();
cnt=1;

%find the first zero of the path
idx = find(ray.path(1,:)==0, 1, 'first');



%
%Take it away
%==========================================================================
% Location of interfaces
%==========================================================================
% configuration 1: only muscle
type_end_lossless_conf1=1;
xloc_end_lossless_conf1=-6;

type_end_muscle_conf1=1; %flat interface
xloc_end_muscle_conf1=2; %initial point of interface
% configuration 2: only muscle
type_end_lossless_conf2=1;
xloc_end_lossless_conf2=-6;

type_end_fat_conf2=1;
xloc_end_fat_conf2=-4;

type_end_muscle_conf2=1;
xloc_end_muscle_conf2=1;

%configuration 5:

type_end_lossless_conf5=1;
xloc_end_lossless_conf5=-6;

type_end_muscle_conf5=1;
xloc_end_muscle_conf5=-3;

type_end_bone_conf5=1;
xloc_end_bone_conf5=-2;

type_end_marrow_conf5=1;
xloc_end_marrow_conf5=-1;

switch(configuration.type)
    case 1 % L-M
        %int0 not used
        int0=Interface(type_end_lossless_conf1, xloc_end_lossless_conf1);% interface L-F
        int1= Interface(type_end_muscle_conf1, xloc_end_muscle_conf1);% creation of the interface
    case 2 % L-F-M
        int0=Interface(type_end_lossless_conf2, xloc_end_lossless_conf2);% interface L-F
        int1=Interface(type_end_fat_conf2,xloc_end_fat_conf2);% interface F-M
        int2=Interface(type_end_muscle_conf2, xloc_end_muscle_conf2);% interface end muscle
    case 5 % L-M-B-M
        int0=Interface(type_end_lossless_conf5, xloc_end_lossless_conf5);% interface L-F
        int1=Interface(type_end_muscle_conf5,xloc_end_muscle_conf5);% interface F-M
        int2=Interface(type_end_bone_conf5, xloc_end_bone_conf5);% interface L-F
        int3=Interface(type_end_marrow_conf5,xloc_end_marrow_conf5);% interface F-M
        
end



switch (ray.material)
    case 1 % ray is in lossless
        %xend=2; %end point of interface
        if (configuration.type==1) %next mat is muscle and the focal is in muscle
            mat2_c=mat.muscle.c* 1e-2;
            mat2_alpha=mat.muscle.attenuation;
            son(cnt).material=3;
            path_r= ray.path;
            path_r(1, idx)=son(cnt).material;
            son(cnt).path=path_r;

            
        elseif (configuration.type==5)%next mat is muscle
            mat2_c=mat.muscle.c* 1e-2;
            mat2_alpha=mat.muscle.attenuation;
            son(cnt).material=2;
            path_r= ray.path;
            path_r(1, idx)=son(cnt).material;
            son(cnt).path=path_r;
            
        elseif (configuration.type==2| configuration.type==3|configuration.type==4)%next material is fat
            mat2_c=mat.fat.c* 1e-2;
            mat2_alpha=mat.fat.attenuation;
            son(cnt).material=2;
            path_r= ray.path;
            path_r(1, idx)=son(cnt).material;
            son(cnt).path=path_r;
            
        end
        
        nn= check_sign(ray.Vray/norm(ray.Vray), int0.normal); %check the normal of the plane
        [refr, v_out]=refraction0(ray.Vray,nn,mat.oil.c* 1e-2,mat2_c);
        if refr==1 %refraction possible
            son(cnt).start=ray.end;
            son(cnt).I0= ray.I0;
            son(cnt).Vray=v_out/norm(v_out);
            son(cnt).direction=1;
            inters_point= Intersection(son(cnt),object);  %ok for all configurations
            son(cnt).end= inters_point.pt;
            son(cnt).IF=ray.I0*exp(-2*mat2_alpha*inters_point.d);
            son(cnt).phase= ray.phase;
           
        end
        
        
    case 2 %the ray is in fat %ok the name of the interfaces!
        % generation of the refracted
        
        if ray.direction==1 % ray is going towards the focal point
            nn= check_sign(ray.Vray/norm(ray.Vray), int1.normal); %check the normal of the end plane of the current ray
            [T,R]=CoefficientLL(ray.Vray,nn,mat.fat.c* 1e-2,mat.muscle.c* 1e-2, mat.fat.attenuation,mat.muscle.attenuation);
            
            % the refracted ray
            [refr, v_out]=refraction0(ray.Vray,nn,mat.fat.c* 1e-2,mat.muscle.c* 1e-2);
            % the reflected ray
           
            if refr==1
                son(cnt).start=ray.end;
                son(cnt).I0= ray.IF * T;
                son(cnt).Vray=v_out/norm(v_out);
                inters_point= Intersection(son(cnt),int2);  %ok for all configurations
                son(cnt).end= inters_point.pt;
                son(cnt).IF=son(cnt).I0*exp(-2*mat.muscle.attenuation*inters_point.d);
                son(cnt).material=3;
                son(cnt).direction=+1;
                son(cnt).phase= ray.phase +inters_point.d* mat.fat.k;
                path_r= ray.path;
                path_r(1, idx)=son(cnt).material;
                son(cnt).path=path_r;
                cnt=cnt+1;
            end
            % the reflected ray
          
            v_out=reflection(ray.Vray,nn);
            son(cnt).start=ray.end;
            son(cnt).I0= ray.IF * R;
            son(cnt).Vray=v_out/norm(v_out);
            inters_point= Intersection(son(cnt),int0);  %ok for all configurations
            son(cnt).end= inters_point.pt;
            son(cnt).IF=son(cnt).I0*exp(-2*mat.fat.attenuation*inters_point.d);
            son(cnt).material=2;
            son(cnt).phase= ray.phase +inters_point.d* mat.fat.k;
            son(cnt).direction=-1;
            path_r= ray.path;
            path_r(1, idx)=son(cnt).material;
            son(cnt).path=path_r;
            cnt=cnt+1;
            
        else %ray is coming towards the transducer Generation Only the reflected!
            
            nn= check_sign(ray.Vray/norm(ray.Vray), int0.normal); %check the normal of the end plane of the current ray
            [T,R]=CoefficientLL(ray.Vray,nn,mat.fat.c* 1e-2,mat.oil.c* 1e-2, mat.fat.attenuation,mat.oil.attenuation);
            v_out=reflection(ray.Vray,nn);
            son(cnt).start=ray.end;
            son(cnt).I0= ray.IF * R;
            son(cnt).Vray=v_out/norm(v_out);
            inters_point= Intersection(son(cnt),int1);  %ok for all configurations
            son(cnt).end= inters_point.pt;
            son(cnt).IF=son(cnt).I0*exp(-2*mat.fat.attenuation*inters_point.d);
            son(cnt).material=2;
            son(cnt).phase= ray.phase +inters_point.d* mat.fat.k;
            son(cnt).direction=+1;
            path_r= ray.path;
            path_r(1, idx)=son(cnt).material;
            son(cnt).path=path_r;
            cnt=cnt+1;
        end
        
        
        
    case 3 %ray is in muscle
        
        
        
        
        
        
    case 4 %ray is in bone
        
end



%clean the son vector 
vec_mat={son.material}';
svid=cellfun(@isempty, vec_mat);
son(svid)=[];
end

