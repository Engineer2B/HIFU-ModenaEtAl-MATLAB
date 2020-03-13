function [ xmin, xmax, ymin,ymax, zmin, zmax, Nx, Ny, Nz, dx, dy, dz, xx, yy, zz, xxb,yyb, zzb ] = Define_table()
%Definition of the table region

xmin =  -1-0.6;  % startpoint(1)+ FocalLength-2;  % cm
xmax =  -1 +0.6; % startpoint(1)+ FocalLength+2;  % cm
ymin=-0.4; %-2; % cm
ymax=0.4;  % cm
zmin=-0.4; % cm
zmax=0.4;  % cm

Nx=floor((xmax-xmin)/0.02) +1; % nr of steps in x direction
Ny=floor((ymax-ymin)/0.02) +1;
Nz=floor((zmax-zmin)/0.02) +1;

dx=(xmax-xmin)/Nx;
dy=(ymax-ymin)/Ny;
dz=(zmax-zmin)/Nz;

xx=xmin-dx/2+[1:Nx]*dx; % x of cube centers
yy=ymin-dy/2+[1:Ny]*dy; % y of cube centers
zz=zmin-dz/2+[1:Nz]*dz; % z of cube centers
xxb=xmin+[0:Nx]*dx; % x of cube boundaries 
yyb=ymin+[0:Ny]*dy; % x of cube boundaries 
zzb=zmin+[0:Nz]*dz; % x of cube boundaries 
end

