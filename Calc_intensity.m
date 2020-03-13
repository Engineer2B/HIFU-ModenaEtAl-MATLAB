function [Intensity,PhaseOneEl,NraysOneEl,Int_bone_long, Int_bone_shear,...
	Pl, PS, NL, NS, polx,poly,polz, ks1,ks2,ks3, kl1,kl2,kl3] = Calc_intensity(ray, object)
% Calc_intensity given ray (struct containing a number of rays) this
% function will calculate the intensiy, phase, nrays in each point of the
% grid


% velx,vely,velz are the component (x,y,z) of the average of the
% polarization direction
% The velocity will be calculate later

[ xmin, xmax, ymin,ymax, zmin, zmax, Nx, Ny, Nz, dx, dy,...
	dz, xx, yy, zz, xxb,yyb, zzb ] = Define_table();

Intensity=zeros(Nx,Ny,Nz); %Intensity in soft
PhaseOneEl=zeros(Nx,Ny,Nz); %Phase in soft
NraysOneEl=zeros(Nx,Ny,Nz);%Number of rays in soft in each cube

poly=zeros(Nx,Ny,Nz); %POLArization direction x 
polx=zeros(Nx,Ny,Nz); %POLArization direction y 
polz=zeros(Nx,Ny,Nz); %POLArization direction z 

Int_bone_long=zeros(Nx,Ny,Nz); %Intensity in bone due to long
Pl=zeros(Nx,Ny,Nz); %Phase in bone due to long
NL=zeros(Nx,Ny,Nz); %Number of rays in each cube long

Int_bone_shear=zeros(Nx,Ny,Nz);%Intensity in bone due to shear
PS=zeros(Nx,Ny,Nz);%Phase in bone due to shear
NS=zeros(Nx,Ny,Nz);%Number of rays in each cube shear

% full approach

ks1=zeros(Nx,Ny,Nz); %direction ray 
ks2=zeros(Nx,Ny,Nz); %direction ray 
ks3=zeros(Nx,Ny,Nz); %direction ray 

kl1=zeros(Nx,Ny,Nz); %direction ray 
kl2=zeros(Nx,Ny,Nz); %direction ray 
kl3=zeros(Nx,Ny,Nz); %direction ray 

% vec_mat={ray.actual_object}';
% svid=~cellfun(@isempty, vec_mat);
% ray(svid)=[];


[~,sizeR]=size(ray);
if (sizeR>0)
		for ri=1: sizeR
				
				if  ray(ri).start(1)<xmin
					%  if ray(ri).end(1)>xmin
						if ray(ri).Vray(1)>0
								lambda_x=(xxb-ray(ri).start(1))/ray(ri).Vray(1); % all positive elements
						else % Vray(1)<=0, ray cannot reach table region
								lambda_x=inf;
						end
						
						if ray(ri).end(1)<xmin % ray cannot reach the table
								lambda_x=inf;
						end
						
				elseif ray(ri).start(1)>xmax
						if ray(ri).Vray(1)<0
								lambda_x=(xxb-ray(ri).start(1))/ray(ri).Vray(1);
								% all positive elements, BUT: decreasing sequence
								lambda_x=lambda_x(end:-1:1); % reverse
						else % Vray(1)>=0, ray cannot reach table region
								lambda_x=inf;
						end
						
						 if ray(ri).end(1)>xmax % ray cannot reach the table
								lambda_x=inf;
						end
						
						
						
				else % xmin<= startpoint(1)<= xmax, start inside table region
						if ray(ri).Vray(1)>0
								lambda_x=(xxb-ray(ri).start(1))/ray(ri).Vray(1); % may contain pos and neg values
								lambda_x=lambda_x(lambda_x>=0); % no neg values
								lambda_x=[0,lambda_x]; % add lambda=0, doubles will be removed later on
						elseif ray(ri).Vray(1)<0
								lambda_x=(xxb-ray(ri).start(1))/ray(ri).Vray(1); % may contain pos and neg values
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
						if  ray(ri).start(2)<ymin
								if ray(ri).Vray(2)>0
										lambda_y=(yyb-ray(ri).start(2))/ray(ri).Vray(2); % all positive elements
								else % Vray(2)<=0, ray cannot reach table region
										lambda_y=inf;
								end
						elseif ray(ri).start(2)>ymax
								if ray(ri).Vray(2)<0
										lambda_y=(yyb-ray(ri).start(2))/ray(ri).Vray(2);
										% all positive elements, BUT decreasing sequence
										lambda_y=lambda_y(end:-1:1); %reverse
								else % Vray(2)>=0, ray cannot reach table region
										lambda_y=inf;
								end
						else % ymin<= startpoint(2)<= ymax, start inside table region
								if ray(ri).Vray(2)>0
										lambda_y=(yyb-ray(ri).start(2))/ray(ri).Vray(2); % may contain pos and neg values
										lambda_y=lambda_y(lambda_y>=0); % no neg values
										lambda_y=[0,lambda_y]; % add lambda=0, doubles will be removed later on
								elseif ray(ri).Vray(2)<0
										lambda_y=(yyb-ray(ri).start(2))/ray(ri).Vray(2); % may contain pos and neg values
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
						if  ray(ri).start(3)<zmin
								if ray(ri).Vray(3)>0
										lambda_z=(zzb-ray(ri).start(3))/ray(ri).Vray(3); % all positive elements
								else % Vray(3)<=0, ray cannot reach table region
										lambda_z=inf;
								end
						elseif ray(ri).start(3)>zmax
								if ray(ri).Vray(3)<0
										lambda_z=(zzb-ray(ri).start(3))/ray(ri).Vray(3);
										% all positive elements, BUT decreasing sequence
										lambda_z=lambda_z(end:-1:1); % reverse
								else % Vray(2)>=0, ray cannot reach table region
										lambda_z=inf;
								end
						else % zmin<= startpoint(3)<= zmax, start inside table region
								if ray(ri).Vray(3)>0
										lambda_z=(zzb-ray(ri).start(3))/ray(ri).Vray(3); % may contain pos and neg values
										lambda_z=lambda_z(lambda_z>=0); % no neg values
										lambda_z=[0,lambda_z]; % add lambda=0, doubles will be removed later on
								elseif ray(ri).Vray(3)<0
										lambda_z=(zzb-ray(ri).start(3))/ray(ri).Vray(3); % may contain pos and neg values
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
						%P_in=ray(ri).start+Min_lambda*ray(ri).Vray;
						%P_out=ray(ri).start+Max_lambda*ray(ri).Vray;
						lambda_x=lambda_x(Min_lambda<=lambda_x);
						lambda_x_restr=lambda_x(lambda_x<=Max_lambda);
						lambda_y=lambda_y(Min_lambda<=lambda_y);
						lambda_y_restr=lambda_y(lambda_y<=Max_lambda);
						lambda_z=lambda_z(Min_lambda<=lambda_z);
						lambda_z_restr=lambda_z(lambda_z<=Max_lambda);
						lambda_interesting=unique(sort([lambda_x_restr,lambda_y_restr,lambda_z_restr]));
						
						for n=1:length(lambda_interesting)-1
								lambda_1= lambda_interesting(n);
								lambda_2= lambda_interesting(n+1);
								lambda_12=(lambda_1+lambda_2)/2;
								P_in_cube=ray(ri).start+lambda_1*ray(ri).Vray;
								P_out_cube=ray(ri).start+lambda_2*ray(ri).Vray;
								%ray(ri).start+lambda_12*ray(ri).Vray;
								ind=floor((ray(ri).start+lambda_12*ray(ri).Vray-[xmin;ymin;zmin])./[dx;dy;dz])+1;
								%
								% Possible condition: I can reduce it, to explain
								%
%                 if ray(ri).actual_object==2 
%                     ffjdsjfkewf
%                 end
						if ray(ri).direction == 1 %ray is going from left to right
								if P_out_cube(1) >  ray(ri).end(1)% ray ends in the cube
										lambda_2=(ray(ri).end- ray(ri).start)'/ray(ri).Vray';
										lambda_12=(lambda_1+lambda_2)/2;
								end
						else %ray is going from right to left
								if P_out_cube(1) < ray(ri).end(1)% ray ends in the cube
										lambda_2=(ray(ri).end- ray(ri).start)'/ray(ri).Vray';
										lambda_12=(lambda_1+lambda_2)/2;
								end
						end;
						%CONDITION NECESSARY to avoid that the rays ending in the
						% middle of the table region will be considered to go through
						% and intersect the cubes 
						if ((ray(ri).direction==1 && P_in_cube(1)<= ray(ri).end(1) ) || ... 
										(ray(ri).direction==-1 && P_in_cube(1)>= ray(ri).end(1)) )  
								% find the correct k and alpha
							if object(ray(ri).actual_object).kind==0 %soft tissue
								alpha= object(ray(ri).actual_object).attenuation;
								kw=object(ray(ri).actual_object).k;
								Intensity(ind(1),ind(2),ind(3))=Intensity(ind(1),ind(2),ind(3))+ ...
										(ray(ri).I0*(exp(-2*alpha*lambda_12)))* ((abs(lambda_1-lambda_2))/(dx*dy*dz));
								PhaseOneEl(ind(1),ind(2),ind(3))=PhaseOneEl(ind(1),ind(2),ind(3))+kw*lambda_12;
								PhaseOneEl(ind(1),ind(2),ind(3))=PhaseOneEl(ind(1),ind(2),ind(3))+ray(ri).phase_initial;
								NraysOneEl(ind(1),ind(2),ind(3))=NraysOneEl(ind(1),ind(2),ind(3))+1;
										
										
							else if ray(ri).shear==0 %long
								alpha= object(ray(ri).actual_object).attenuationl;
								kw=object(ray(ri).actual_object).kl;
								% treat bone here
								Int_bone_long(ind(1),ind(2),ind(3))=Int_bone_long(ind(1),ind(2),ind(3))+...
										(ray(ri).I0*(exp(-2*alpha*lambda_12)))* ((abs(lambda_1-lambda_2))/(dx*dy*dz));
								Pl(ind(1),ind(2),ind(3))=Pl(ind(1),ind(2),ind(3))+kw*lambda_12;
								Pl(ind(1),ind(2),ind(3))=Pl(ind(1),ind(2),ind(3))+ray(ri).phase_initial;
								NL(ind(1),ind(2),ind(3))=NL(ind(1),ind(2),ind(3))+1;
								
								else
									alpha= object(ray(ri).actual_object).attenuations;
									kw=object(ray(ri).actual_object).ks;
									% treat bone here
									Int_bone_shear(ind(1),ind(2),ind(3))=Int_bone_shear(ind(1),ind(2),ind(3))+...
											(ray(ri).I0*(exp(-2*alpha*lambda_12)))* ((abs(lambda_1-lambda_2))/(dx*dy*dz));
									PS(ind(1),ind(2),ind(3))=PS(ind(1),ind(2),ind(3))+kw*lambda_12;
									PS(ind(1),ind(2),ind(3))=PS(ind(1),ind(2),ind(3))+ray(ri).phase_initial;
									NS(ind(1),ind(2),ind(3))=NS(ind(1),ind(2),ind(3))+1;

									polx(ind(1),ind(2),ind(3))=polx(ind(1),ind(2),ind(3))+ray(ri).polarization(1);% sum of the pol dir x
									poly(ind(1),ind(2),ind(3))=poly(ind(1),ind(2),ind(3))+ray(ri).polarization(2);% sum of the pol dir y
									polz(ind(1),ind(2),ind(3))=polz(ind(1),ind(2),ind(3))+ray(ri).polarization(3);% sum of the pol dir z

									ks1(ind(1),ind(2),ind(3))= ks1(ind(1),ind(2),ind(3))+ray(ri).Vray(1);% sum of the dir x
									ks2(ind(1),ind(2),ind(3))=ks2(ind(1),ind(2),ind(3))+ray(ri).Vray(2);% sum of the dir y
									ks3(ind(1),ind(2),ind(3))=ks3(ind(1),ind(2),ind(3))+ray(ri).Vray(3);% sum of the  dir z
								end
							end
						end
					end;
				end
		end
end



polx=polx./max(NS,1); %calculation of the average
poly=poly./max(NS,1);
polz=polz./max(NS,1);


velx2=polx.^2;
vely2=poly.^2;
velz2=polz.^2;

vnorm=velx2+vely2+velz2;
vnorm(vnorm==0)=1;
vnorm=sqrt(vnorm);

polx=polx./(vnorm); %normalization
poly=poly./(vnorm);
polz=polz./(vnorm);
%Vray long
kl1=kl1./max(NL,1); %calculation of the average
kl2=kl2./max(NL,1);
kl3=kl3./max(NL,1);

klx2=kl1.^2;
kly2=kl2.^2;
klz2=kl3.^2;

vnorm=klx2+kly2+klz2;
vnorm(vnorm==0)=1;
vnorm=sqrt(vnorm);

kl1=kl1./(vnorm); %normalization
kl2=kl2./(vnorm);
kl3=kl3./(vnorm);

%Vray shear
ks1=ks1./max(NS,1); %calculation of the average
ks2=ks2./max(NS,1);
ks3=ks3./max(NS,1);

ksx2=ks1.^2;
ksy2=ks2.^2;
ksz2=ks3.^2;

vnorm=ksx2+ksy2+ksz2;
vnorm(vnorm==0)=1;
vnorm=sqrt(vnorm);

ks1=ks1./(vnorm); %normalization
ks2=ks2./(vnorm);
ks3=ks3./(vnorm);
