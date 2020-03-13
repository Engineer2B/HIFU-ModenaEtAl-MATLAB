function int = Intersection( ray, object, obj_ray)
%UNTITLED5 Calculates the intersection point between a ray and a interface
%   Detailed explanation goes here

[a,b]=size(object);
int(b).xi=0;
cnt=1;


for i=1:b
    if(object(i).activate==1)
        switch(object(i).interface.type)
            case 1 %flat or inclined
                int_object=find(object(i).int_index(2,:)==i);
                
                d= ((object(i).interface.point - ray.start)'*object(i).interface.normal)/(ray.Vray' * object(i).interface.normal);
                if d>0.1
                    int(cnt).pt= d*ray.Vray+ray.start; %what point
                    int(cnt).d= d;
                    
                    int(cnt).x_int=int(cnt).pt(1);%x of the intersection
                    plane_normal = object(i).interface.normal;
                    dir = ray.Vray/norm(ray.Vray);
                    denominator=dir'*plane_normal;
                    if denominator>0
                        nn=-plane_normal;
                    else
                        nn=plane_normal;
                    end;
                    int(cnt).normal=nn;
                    if ray.direction==1 || int_object==1
                        int(cnt).obj=object(i);
                        int(cnt).index_obj=i;
                        
                    else
                        int(cnt).obj=object(object(i).int_index(2,int_object-1));
                        int(cnt).index_obj=object(i).int_index(2,int_object-1);
                    end
                    cnt=cnt+1;
                end
                
            case 2 % the object is a cylinder
                
                int_object=find(object(i).int_index(2,:)==i);
                cylcentpoint=object(i).interface.center;
                cyldir= object(i).interface.dir/norm(object(i).interface.dir);
                cylradius= object(i).interface.radius;
                [d,nn]=  Intersect_cylinder(ray.start,ray.Vray,cylcentpoint,cyldir,cylradius);
                
%                 plane_normal=nn;
%                 dir = ray.Vray/norm(ray.Vray);
%                 denominator=dir'*plane_normal;
%                 if denominator>0
%                     nn=-plane_normal;
%                 else
%                     nn=plane_normal;
%                 end;
                
                
                
                if d>0.1
                    int(cnt).pt= d*ray.Vray+ray.start; %what point
                    int(cnt).d= d;
                    int(cnt).x_int=int(cnt).pt(1);%x of the intersection
                    int(cnt).normal=nn;
                    
                    if ray.direction==1 
                        int(cnt).obj=object(i);
                        int(cnt).index_obj=i;
                        if int(cnt).index_obj==obj_ray %if the ray passes through the cylinder without reaching marrow
                            int(cnt).obj=object(object(i).int_index(2,int_object-1));
                            int(cnt).index_obj=object(i).int_index(2,int_object-1);
                        end
                    else % ray is coming back to the trds
                        int(cnt).obj=object(object(i).int_index(2,int_object-1));
                        int(cnt).index_obj=object(i).int_index(2,int_object-1);
                    end
                    
                     if  int_object==1
                        int(cnt).obj=object(i);
                        int(cnt).index_obj=i;
                     end
                    
                    cnt=cnt+1;
                end
                
                
        end
        
    end
end


vec_x={int.x_int}';
%clean d

svid=cellfun(@isempty, vec_x); %clean the interfaces empty
int(svid)=[];



Myx=cell2mat(vec_x);

if ray.direction==1 %only the x bigger than the start
    st= Myx> ray.start(1); %take only positive d
    int=int(st);
    Myx=Myx(st);
    [m, index]=min(Myx);
else %only negative
    st= Myx< ray.start(1); %change the number here eh!
    int=int(st);
    Myx=Myx(st);
    [m, index]=max(Myx);
end

int= int(index);


