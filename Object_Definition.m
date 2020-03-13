function object=Object_Definition()
%UNTITLED definition of the objects
%   Detailed explanation goes here
%
material= Define_material(); % definition of the materials

object(1).name_material='oil';
object(1).material=1;
object(1).obj=1;
object(1).interface=Interface(1, -15, Find_normal_plane_inclined(0,0));
object(1).activate=1; %1= object is present 0= object is not
object(1).attenuation=material.(object(1).name_material).attenuation;
object(1).density=material.(object(1).name_material).density;
object(1).c=material.(object(1).name_material).c;
object(1).z=material.(object(1).name_material).z;
object(1).k=material.(object(1).name_material).k;
object(1).ka=material.(object(1).name_material).ka;
object(1).kind=0; %soft
object(1).xint=object(1).interface.xi; %a point of the interface


object(2).name_material='fat';
object(2).material=2;
object(2).obj=2; % this is the index in the struct!
object(2).interface=Interface(1, -6, Find_normal_plane_inclined(0,0));
object(2).activate=0; %1= object is present 0= object is not
object(2).attenuation=material.(object(2).name_material).attenuation;
object(2).density=material.(object(2).name_material).density;
object(2).c=material.(object(2).name_material).c;
object(2).z=material.(object(2).name_material).z;
object(2).k=material.(object(2).name_material).k;
object(2).ka=material.(object(2).name_material).ka;
object(2).kind=0; %soft
object(2).xint=object(2).interface.xi; %a point of the interface


object(3).name_material='muscle';
object(3).material=3;
object(3).obj=3; % this is the index in the struct!
object(3).interface=Interface(1, -6,Find_normal_plane_inclined(0,0)); %Interface(type, xloc, normal) xloc is where the interface is
object(3).activate=1; %1= object is present 0= object is not
object(3).attenuation=material.(object(3).name_material).attenuation;
object(3).density=material.(object(3).name_material).density;
object(3).c=material.(object(3).name_material).c;
object(3).z=material.(object(3).name_material).z;
object(3).k=material.(object(3).name_material).k;
object(3).ka=material.(object(3).name_material).ka;
object(3).kind=0; %soft
object(3).xint=object(3).interface.xi; %a point of the interface


object(4).name_material='bone';
object(4).material=4;
object(4).obj=4; % this is the index in the struct!
center.x=0; % center of the bone
center.y=0;
inclination_bone=0;
object(4).interface=Interface(1, -1,Find_normal_plane_inclined(0,0));
%object(4).interface=Interface(2, center,inclination_bone,1); %Interface(type, xloc, normal) xloc is where the interface is
object(4).activate=1; %1= object is present 0= object is not
object(4).attenuationl=material.(object(4).name_material).attenuationl;
object(4).attenuations=material.(object(4).name_material).attenuations;
object(4).density=material.(object(4).name_material).density;
object(4).clong=material.(object(4).name_material).clong;
object(4).cshear=material.(object(4).name_material).cshear;
object(4).zl=material.(object(4).name_material).zl;
object(4).zs=material.(object(4).name_material).zs;
object(4).ks=material.(object(4).name_material).ks;
object(4).kl=material.(object(4).name_material).kl;
object(4).kal=material.(object(4).name_material).kal;
object(4).kas=material.(object(4).name_material).kas;
object(4).kind=1; %bone
object(4).xint=object(4).interface.xi; %a point of the interface

object(5).name_material='muscle'; %marrow
object(5).material=5;
object(5).obj=5; % this is the index in the struct!
center.x=0; % center of the bone
center.y=0;
inclination_marrow=0;

object(5).interface=Interface(2, center,inclination_marrow,0.5); %Interface(type, xloc, normal) xloc is where the interface is
object(5).activate=0; %1= object is present 0= object is not
object(5).attenuation=material.(object(3).name_material).attenuation;
object(5).density=material.(object(3).name_material).density;
object(5).c=material.(object(3).name_material).c;
object(5).z=material.(object(3).name_material).z;
object(5).k=material.(object(3).name_material).k;
object(5).ka=material.(object(3).name_material).ka;
object(5).kind=0; %soft
object(5).xint=object(5).interface.xi; %a point of the interface


% here definition of the last material
% depending on what is the last one
% needs to arrange it automatically
vec_mat={object.activate}'; %check what objects are activated
Mymat=cell2mat(vec_mat);
st= Mymat == 1   ;
obj2=object(st);

vec_mat={obj2.xint}'; %check what is the bigger interface starting point
Mymat=cell2mat(vec_mat);
[~,b]=max(Mymat); % to find the last object
ind=obj2(b).obj; % I have to save what object is


[~,size_obj]=size(object);

object(size_obj+1).name_material='final';
object(size_obj+1).material=ind;
object(size_obj+1).obj=[]; % this is the index in the struct!
object(size_obj+1).interface=Interface(1, 2,[1;0;0]); %Interface(type, xloc, normal) xloc is where the interface is
object(size_obj+1).activate=1; %1= object is present 0= object is not
object(size_obj+1).density=material.(object(ind).name_material).density;
object(size_obj+1).xint=object(size_obj+1).interface.xi; %a point of the interface
if object(ind).kind==1 % it's a bone
    object(size_obj+1).attenuationl=material.(object(ind).name_material).attenuationl;
    object(size_obj+1).attenuations=material.(object(ind).name_material).attenuations;
    object(size_obj+1).clong=material.(object(ind).name_material).clong;
    object(size_obj+1).cshear=material.(object(ind).name_material).cshear;
    object(size_obj+1).zl=material.(object(ind).name_material).zl;
    object(size_obj+1).zs=material.(object(ind).name_material).zs;
    object(size_obj+1).ks=material.(object(ind).name_material).ks;
    object(size_obj+1).kl=material.(object(ind).name_material).kl;
    object(size_obj+1).kal=material.(object(ind).name_material).kal;
    object(size_obj+1).kas=material.(object(ind).name_material).kas;
    object(size_obj+1).kind=1; %bone
else
    object(size_obj+1).attenuation=material.(object(ind).name_material).attenuation;
    object(size_obj+1).c=material.(object(ind).name_material).c;
    object(size_obj+1).z=material.(object(ind).name_material).z;
    object(size_obj+1).k=material.(object(ind).name_material).k;
    object(size_obj+1).ka=material.(object(ind).name_material).ka;
    object(size_obj+1).kind=0; %soft
end

% vector of xint
[~,size_obj]=size(object);
z= zeros(2,size_obj);
ct=1;
for i=1:size_obj
    if object(i).activate==1
        z(1,ct)=object(i).xint;
        z(2,ct)=i;
        ct=ct+1;
    end
end
A=find(z(1,:)==0);
z(:,A)=[];
[Ordered_int,I] = sort(z(1,:)); %order the interfaces
%[~,sizez]=size(z);

sizez=ct-1;
int_index= zeros(2,sizez);
for i=1:sizez
    int_index(1,i)=Ordered_int(i);
    int_index(2,i)=z(2,I(i));
end
    
for i=1:size_obj
    if object(i).activate==1
       object(i).int_index=int_index;
    end
end

end

