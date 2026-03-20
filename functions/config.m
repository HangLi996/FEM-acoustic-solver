air.c = 340; % Speed of sound [m/s]
air.rho = 1.3; % Density of air [kg/m^3]

freq = 10:10:90;
omega = 2*pi*freq;
k = omega / air.c;

file.meshfolder = 'mesh/';
file.path = './';

conf.mesher = 'GMSH';