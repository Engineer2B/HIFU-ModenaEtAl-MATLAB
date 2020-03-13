function Mat =Define_material( )
% Define_material Definition of properties and boundaries
%   Detailed explanation goes here

fr=1.2e6;
omega=2*pi*fr; 
radius=0.35; % transducer diameter =7 mm, radius=3.5 mm=0.35 cm
% FocalLength=14; % cm
% Area=pi*radius^2; % area of transducer
% RequiredTotalPower=100; % watt

%% Lossless

% definition of the material parameters in Lossless
material.oil.c=1380*100; % speed of sound in gel, in: cm/s*100 
material.oil.density=1030e-6; % gel, density in kg/cm^3
material.oil.z=material.oil.c*material.oil.density;
material.oil.attenuation=0;% attenuation for pressure, per cm
material.oil.k= omega/material.oil.c; % unit: cm^-1
material.oil.ka=material.oil.k* radius;

% definition of the boundaries in LossLess


%% Muscle
material.muscle.c=1537*100; % speed of sound in gel, in: cm/s*100 
material.muscle.density=1010e-6; % gel, density in kg/cm^3
material.muscle.z=material.muscle.c*material.muscle.density;
material.muscle.attenuation=0.0576;% attenuation for pressure, per cm
material.muscle.k= omega/material.muscle.c; % unit: cm^-1
material.muscle.ka=material.muscle.k* radius;

%% definition Fat
material.fat.c=1450*100; % speed of sound in gel, in: cm/s*100 
material.fat.density=0.9094e-3; % gel, density in kg/cm^3
material.fat.z=material.fat.c*material.fat.density;
material.fat.attenuation=0.066;% attenuation for pressure, per cm
material.fat.k= omega/material.fat.c; % unit: cm^-1
material.fat.ka=material.fat.k* radius;

%% definition Bone
material.bone.clong=3736*100; % speed of sound in gel, in: cm/s*100 
material.bone.cshear=1995*100; % speed of sound in gel, in: cm/s*100 
material.bone.density=2025e-6;% gel, density in kg/cm^3
material.bone.attenuationl=1.9;% attenuation for pressure, per cm
material.bone.attenuations=2.8;% attenuation for pressure, per cm
material.bone.zl=material.bone.clong*material.bone.density;
material.bone.zs=material.bone.cshear*material.fat.density;
material.bone.kl= omega/material.bone.clong; % unit: cm^-1
material.bone.ks= omega/material.bone.cshear; % unit: cm^-1
material.bone.kas=material.bone.kl* radius;
material.bone.kal=material.bone.ks* radius;

Mat= material;
end

