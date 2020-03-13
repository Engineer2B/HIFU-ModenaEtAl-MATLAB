%[xmin, xmax, ymin,ymax, zmin, zmax, Nx, Ny, Nz, dx, dy, dz, xx, yy, zz, xxb,yyb, zzb ] = Define_table();
clear all

xmin=0;
xmax=1;
ymin=0;
ymax=1;
zmin=0;
zmax=1;
dx=1;
dy=1;
dz=1;
Nx=ceil((xmax-xmin)/dx ); % nr of steps in x direction
Ny=ceil(ymax-ymin)/dy ;
Nz=ceil(zmax-zmin)/dz;

xx=xmin-dx/2+[1:Nx]*dx; % x of cube centers
yy=ymin-dy/2+[1:Ny]*dy; % y of cube centers
zz=zmin-dz/2+[1:Nz]*dz; % z of cube centers

xxb=xmin+[0:Nx]*dx; % x of cube boundaries
yyb=ymin+[0:Ny]*dy; % x of cube boundaries
zzb=zmin+[0:Nz]*dz; % x of cube boundaries

sizeR=1;
% Ray_tot(1).end(1)=-0.5;
% Ray_tot(1).end(2)=-0.5+2*(1/6);
% Ray_tot(1).end(3)=-1+2.5*(1/6);

Ray_tot(1).start(1)=-0.5;
Ray_tot(1).start(2)=-0.5+2*(1/6);
Ray_tot(1).start(3)=-1+2.5*(1/6);
%
% Ray_tot(1).start(1)=0.5; %%check when a ray starts inside a cube
% Ray_tot(1).start(2)=0.5;
% Ray_tot(1).start(3)=0.5;


% %
% Ray_tot(1).end(1)=1.5;
% Ray_tot(1).end(2)=-0.5+2*(2.5/3);
% Ray_tot(1).end(3)=-1+2.5*(2.5/3);

% Ray_tot(1).start(1)=0.5;
% Ray_tot(1).start(2)=0.5;
% Ray_tot(1).start(3)=0.5;

Ray_tot(1).end(1)=0.5;%%check when a ray terminates inside a cube
Ray_tot(1).end(2)=0.5;
Ray_tot(1).end(3)=0.5;

Ray_tot(1).start=[Ray_tot(1).start(1);Ray_tot(1).start(2);Ray_tot(1).start(3)];
Ray_tot(1).Vray(1)=Ray_tot(1).end(1)-Ray_tot(1).start(1);
Ray_tot(1).Vray(2)=Ray_tot(1).end(2)-Ray_tot(1).start(2);
Ray_tot(1).Vray(3)=Ray_tot(1).end(3)-Ray_tot(1).start(3);

Ray_tot(1).Vray=[Ray_tot(1).Vray(1);Ray_tot(1).Vray(2);Ray_tot(1).Vray(3)];

for ri=1: sizeR
    if  Ray_tot(ri).start(1)<xmin 
        if Ray_tot(ri).Vray(1)>0
            lambda_x=(xxb-Ray_tot(ri).start(1))/Ray_tot(ri).Vray(1); % all positive elements
        else % Vray(1)<=0, ray cannot reach table region
            lambda_x=inf;
        end
    elseif Ray_tot(ri).start(1)>xmax
        if Ray_tot(ri).Vray(1)<0
            lambda_x=(xxb-Ray_tot(ri).start(1))/Ray_tot(ri).Vray(1);
            % all positive elements, BUT: decreasing sequence
            lambda_x=lambda_x(end:-1:1); % reverse
        else % Vray(1)>=0, ray cannot reach table region
            lambda_x=inf;
        end
    else % xmin<= startpoint(1)<= xmax, start inside table region
        if Ray_tot(ri).Vray(1)>0
            lambda_x=(xxb-Ray_tot(ri).start(1))/Ray_tot(ri).Vray(1); % may contain pos and neg values
            lambda_x=lambda_x(lambda_x>=0); % no neg values
            lambda_x=[0,lambda_x]; % add lambda=0, doubles will be removed later on
        elseif Ray_tot(ri).Vray(1)<0
            lambda_x=(xxb-Ray_tot(ri).start(1))/Ray_tot(ri).Vray(1); % may contain pos and neg values
            lambda_x=lambda_x(lambda_x>=0); % no neg values, BUT decreasing sequence
            lambda_x=lambda_x(end:-1:1); % reverse
            lambda_x=[0,lambda_x]; % add lambda=0, doubles will be removed later
        else % Vray(1) ==0
            lambda_x=0;
        end;
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if lambda_x(1)<inf % ray may pass through table region
        % Generate lambda_y, the series of lambdas such that
        % startpoint+lambda_y*Vray crosses a y-boundary between two cubes
        if  Ray_tot(ri).start(2)<ymin
            if Ray_tot(ri).Vray(2)>0
                lambda_y=(yyb-Ray_tot(ri).start(2))/Ray_tot(ri).Vray(2); % all positive elements
            else % Vray(2)<=0, ray cannot reach table region
                lambda_y=inf;
            end
        elseif Ray_tot(ri).start(2)>ymax
            if Ray_tot(ri).Vray(2)<0
                lambda_y=(yyb-Ray_tot(ri).start(2))/Ray_tot(ri).Vray(2);
                % all positive elements, BUT decreasing sequence
                lambda_y=lambda_y(end:-1:1); %reverse
            else % Vray(2)>=0, ray cannot reach table region
                lambda_y=inf;
            end
        else % ymin<= startpoint(2)<= ymax, start inside table region
            if Ray_tot(ri).Vray(2)>0
                lambda_y=(yyb-Ray_tot(ri).start(2))/Ray_tot(ri).Vray(2); % may contain pos and neg values
                lambda_y=lambda_y(lambda_y>=0); % no neg values
                lambda_y=[0,lambda_y]; % add lambda=0, doubles will be removed later on
            elseif Ray_tot(ri).Vray(2)<0
                lambda_y=(yyb-Ray_tot(ri).start(2))/Ray_tot(ri).Vray(2); % may contain pos and neg values
                lambda_y=lambda_y(lambda_y>=0); % no neg values, BUT decreasing sequence
                lambda_y=lambda_y(end:-1:1); % reverse
                lambda_y=[0,lambda_y];% add lambda=0, doubles will be removed later
            else % Vray(2) ==0
                lambda_y=0;
            end;
        end;
    else
        lambda_y=inf; % lambda_x(1)=inf, so ray does not pass through table region anyway
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if lambda_x(1)<inf && lambda_y(1)<inf % ray may pass through table region
        % Generate lambda_z, the seies of lambdas such that
        % startpoint+lambda_z*Vray crosses a z-boundary between two cubes
        if  Ray_tot(ri).start(3)<zmin
            if Ray_tot(ri).Vray(3)>0
                lambda_z=(zzb-Ray_tot(ri).start(3))/Ray_tot(ri).Vray(3); % all positive elements
            else % Vray(3)<=0, ray cannot reach table region
                lambda_z=inf;
            end
        elseif Ray_tot(ri).start(3)>zmax
            if Ray_tot(ri).Vray(3)<0
                lambda_z=(zzb-Ray_tot(ri).start(3))/Ray_tot(ri).Vray(3);
                % all positive elements, BUT decreasing sequence
                lambda_z=lambda_z(end:-1:1); % reverse
            else % Vray(2)>=0, ray cannot reach table region
                lambda_z=inf;
            end
        else % zmin<= startpoint(3)<= zmax, start inside table region
            if Ray_tot(ri).Vray(3)>0
                lambda_z=(zzb-Ray_tot(ri).start(3))/Ray_tot(ri).Vray(3); % may contain pos and neg values
                lambda_z=lambda_z(lambda_z>=0); % no neg values
                lambda_z=[0,lambda_z]; % add lambda=0, doubles will be removed later on
            elseif Ray_tot(ri).Vray(3)<0
                lambda_z=(zzb-Ray_tot(ri).start(3))/Ray_tot(ri).Vray(3); % may contain pos and neg values
                lambda_z=lambda_z(lambda_z>=0); % no neg values, BUT decreasing sequence
                lambda_z=lambda_z(end:-1:1); % reverse
                lambda_z=[0,lambda_z]; %add lambda=0, doubles will be removed later
            else % Vray(3) ==0
                lambda_z=0;
            end;
        end;
    else
        lambda_z=inf; % lambda_x(1)=inf or lambda_y(1)=inf,
        % so ray does not pass through table region anyway
    end;
    
    % now process the three lambda sequences
    Min_lambda = max([lambda_x(1),lambda_y(1),lambda_z(1)]);
    Max_lambda = min([lambda_x(end),lambda_y(end),lambda_z(end)]);
    % part of the ray between Min_lambda and Max_lambda is
    % in the table region;
    if Min_lambda<Max_lambda %ray passes through table region
        P_in=Ray_tot(ri).start+Min_lambda*Ray_tot(ri).Vray;
        P_out=Ray_tot(ri).start+Max_lambda*Ray_tot(ri).Vray;
        lambda_x=lambda_x(Min_lambda<=lambda_x);
        lambda_x_restr=lambda_x(lambda_x<=Max_lambda);
        lambda_y=lambda_y(Min_lambda<=lambda_y);
        lambda_y_restr=lambda_y(lambda_y<=Max_lambda);
        lambda_z=lambda_z(Min_lambda<=lambda_z);
        lambda_z_restr=lambda_z(lambda_z<=Max_lambda);
        lambda_interesting=unique(sort([lambda_x_restr,lambda_y_restr,lambda_z_restr]));
        
        for n=1:length(lambda_interesting)-1
            % Add a condition here about Pin and Pout
            lambda_1= lambda_interesting(n);
            lambda_2= lambda_interesting(n+1);
            lambda_12=(lambda_1+lambda_2)/2;
            P_in_cube=Ray_tot(ri).start+lambda_1*Ray_tot(ri).Vray;
            P_out_cube=Ray_tot(ri).start+lambda_2*Ray_tot(ri).Vray;
            ind=floor((Ray_tot(ri).start+lambda_12*Ray_tot(ri).Vray-[xmin;ymin;zmin])./[dx;dy;dz])+1;
           
            %
            % Possible condition: I can reduce it, to explain
            %
            if P_out_cube(1)> P_in_cube(1) %ray is going from left to right
                if P_out_cube(1) >  Ray_tot(ri).end(1)% ray ends in the cube
                    nlambda_2=(Ray_tot(ri).end- Ray_tot(ri).start')/Ray_tot(ri).Vray';
                end
            else %ray is going from right to left
                if P_out_cube(1) < Ray_tot(ri).end(1)% ray ends in the cube
                    nlambda_2=(Ray_tot(ri).end- Ray_tot(ri).start')/Ray_tot(ri).Vray';
                end
            end;
           % lambda_12=(lambda_1+lambda_2)/2;
            newPout=Ray_tot(ri).start+nlambda_2*Ray_tot(ri).Vray;
            %nlambda_2=(Ray_tot(ri).end - Ray_tot(ri).start)/Ray_tot(ri).Vray;
        end
    end
end

x=[
    0 1 1 0 0 % bottom
    0 1 1 0 0 % top
    0 0 0 0 0 % left
    1 1 1 1 1 % right
    ];
y=[
    0 0 1 1 0
    0 0 1 1 0
    0 1 1 0 0
    0 1 1 0 0
    ];
z=[
    0 0 0 0 0
    1 1 1 1 1
    0 0 1 1 0
    0 0 1 1 0
    ];



figure()
Pt=[Ray_tot(1).start(1) , Ray_tot(1).start(2),Ray_tot(1).start(3);Ray_tot(1).end(1), Ray_tot(1).end(2),Ray_tot(1).end(3)];
line(Pt(:,1), Pt(:,2), Pt(:,3),'LineWidth',2)
hold on
plot3(xxb, yyb,zzb, 'o')
hold on
line(x',y',z','color',[0 0 0],'LineWidth',1);
view(3);
hold on
plot3(P_out(1), P_out(2),P_out(3), 'o')

hold on
plot3(newPout(1), newPout(2),newPout(3), '*')
hold on
plot3(P_in_cube(1), P_in_cube(2),P_in_cube(3), '*')
axis equal