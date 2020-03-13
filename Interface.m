function Interfac = Interface(type, xloc, fn,radius)
%INTERFACE define the different interfaces
%   % TYPE
%   type= 1 flat or Inclined depending on fn
%     xloc= location of the interface, fn= the normal at the interface,
%     radius=0 (NO RADIUS!)
%     fn=normal to the interface
%
%   type= 2 cylinder 
%     xloc= struct cointaing x and y of the center of the cylinder, fn= the inclination of the cylinder,
%     radius= radius of the cylinder
%     fn= bone inclination



switch(type)
    case 1 %flat or inclined plane
        Interfac.type=1;
        Interfac.xi= xloc; 
        %Interfac.normal=[1;0;0]; %normal to the plane
      %  Interfac.normal=Find_normal_plane_inclined( bet, gamm);
        Interfac.normal= fn;
        Interfac.point=[xloc;0;0]; %point belonging to the plane
      %  Interfac.xf=[xend;0;0];
    
        
    % ... other cases
    case 2 % cylinder
        BoneIncl=fn;
        sinI=sin(BoneIncl/180*pi);
        cosI=cos(BoneIncl/180*pi);
        xcenter=xloc.x;
        ycenter=xloc.y;
        Interfac.type=2;
        Interfac.center=[xcenter; ycenter; 0];
        Interfac.dir=[sinI; 0;cosI];% Direction of the cylinder
        Interfac.radius=radius;%radius
        Interfac.xi= xcenter-radius; % I need that to order the interfaces!


        
        
       % ..
end



end

