function [Pressure,  Ptrd, Int_bone_long, velx,vely,velz,vl1,vl2,vl3,eps12,eps13, eps23]=CalculationPressure(Nrays,trd, phase_matrix, object,xTrd,yTrd,zTrd,material)
%========================================================================
% Loop called in the parfor (256 times).
%this function calls for every trd UpdateRays
%========================================================================
I_init     = 1; % DO NOT CHANGE !!!
r1=3.8317; %first zero bessels
sin_theta_max=r1/(material.oil.ka*4); % smaller, only for shape of focal point
cos_theta_max=sqrt(1-sin_theta_max^2);
AA=1-cos_theta_max;
Ptrd=0;
startpoint=[xTrd(trd);yTrd(trd);zTrd(trd)];
normal =-startpoint/norm(startpoint);
v2=[-normal(3);0;normal(1)]; % a vector perpendicular to the normal, so in the plane of the transd.
if norm(v2)<1e-12
    input('ERROR, vector v2 vanishes'); 
end;
v2norm=v2/norm(v2); % normalised v2
v3=cross(v2norm,normal); % another vector in the plane of the transd
%=========================================================================
% In this loop I calculate the Rays Refracted and reflected
%=========================================================================

numberMaxSons=15; %each ray in lossless can generate up to 15 rays sons
Raylist(Nrays*numberMaxSons).xi= 0; %Initialize my big struct with ray
cnt=1;

for nray=1: Nrays
    r=rand;
    theta=acos(1-AA*r);
    phi=2*pi*rand(1);%  other angle (azimuth)
    [son, ptr,yy]= Update_rays(theta, phi, startpoint,I_init,normal, v2norm,v3,material, phase_matrix(trd),object);% For each ray in lossless I generate the sons!
    [~,sizeS]=size(son);
    fn=fieldnames(son);% fields of son
    [s,~]=size(fn);
    for j=1:s  %copy the rays in the big Raylist
        [Raylist(cnt:cnt+sizeS-1).(fn{j})]=son.(fn{j}); %copying the fields
    end
    cnt=cnt+sizeS;
    Ptrd= ptr+ Ptrd;
end


%clean Raylist Time compsuming! But necessary!

vec_obj={Raylist.actual_object}';
svid=cellfun(@isempty, vec_obj);
Raylist(svid)=[];
%==========================================================================
% %use this part of the code to print rays -> take care it prints ALL rays
% % suitable only for Nrays=1
% % print Rays
% 
% scatter3(startpoint(1),startpoint(2),startpoint(3), 40, 'red', 'filled');
% hold on;
% Pt=[yy.start(1) ,yy.start(2),yy.start(3);...
%     yy.end(1) ,yy.end(2),yy.end(3)];
% line(Pt(:,1), Pt(:,2), Pt(:,3),'color', 'green')
% hold on;
% % theta= 0:0.01:2*pi;
% % r=1;
% % x= r*cos(theta);
% % y= r*sin(theta);
% % plot(0+x,0+y)
% [~,b]=size(Raylist);
% for i=1:b
%      hold on;
%      if abs(Raylist(i).start(1)- Raylist(i).end(1)) < 10
%     if i==1 || i==3
%         Pt=[Raylist(i).start(1) ,Raylist(i).start(2),Raylist(i).start(3);...
%             Raylist(i).end(1) ,Raylist(i).end(2),Raylist(i).end(3)];
%         line(Pt(:,1), Pt(:,2), Pt(:,3))
%     else
%         Pt=[Raylist(i).start(1) ,Raylist(i).start(2),Raylist(i).start(3);...
%             Raylist(i).end(1) ,Raylist(i).end(2),Raylist(i).end(3)];
%         line(Pt(:,1), Pt(:,2), Pt(:,3),'color', 'red')
%     end
%      end
%     hold on;
% end
% % hold on
%==========================================================================

%=========================================================================
% TAKE CARE -> please pay attention here. To speed up Change the Material
%in which object you want to calculate the power produced 
% Mat=3 muscle
% Mat=4 bone
% Mat= 2 fat
% Mat= 5 marrow
%=========================================================================

%check in what material/s I have to calculate the Intensity and Pressure
[~,b]=size(Raylist);
vec_obj={Raylist.actual_object}';
Mymat=cell2mat(vec_obj);
 %find the indices of the rows which are equal to the groups
st= Mymat == 4 |Mymat == 3    ; %change the number here eh! To apply this row correctly I need to use a clean Raylist
Raylist=Raylist(st);

[Pressure, Int_bone_long,  velx,vely,velz,vl1,vl2,vl3,eps12,eps13, eps23]=Process_rays(Raylist, object); %function which processes all the interesting rays of one transducer




