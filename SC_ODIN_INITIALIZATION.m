%% SPOOQY-1                
% 1 44332U 98067QH  19350.77721468  .00006274  00000-0  10054-3 0  9998
% 2 44332  51.6386 172.4798 0005656  48.9347 311.2132 15.54514188 2835547
%TLE epoch 16 Dec 2019 18:39:11.348
%SPOOQY-1 is a 3U cubesat lauched from ISS in 2019

%*************TO RUN IN SIMULINK ONLY*********************
%****************ADCS_Outline.slx*************************
clear
close all
longstr1 = '1 44332U 98067QH  19350.77721468  .00006274  00000-0  10054-3 0  9998';
longstr2 = '2 44332  51.6386 172.4798 0005656  48.9347 311.2132 15.54514188 2835547';
opsmode = 'i';
whichconst = 84;
%% simulation setup
duration = 5*60*90;%seconds
step = 20; %seconds for SGP4
%% targets
ground = [1266.06;-4295.81;4526.14]; %ITRF [X;Y;Z] %concordia
argentina = [1612.36;-4182.16;-4522.83]; %ITRF [X;Y;Z]
namibia = [5659.26;1238.85;-2659.2]; %ITRF [X;Y;Z;]
%% start time

    %TLE epoch start case
    year = 2019; %yyyy
    month = 12; %mm
    day = 16; %dd
    hour = 18; %HH
    minute = 39; %MM
    second = 11.348; %SS.FFF

    %interesting case when the sun and the magnetic field are in the same
    %direction, *also the eclipse case*
%     year = 2020; %yyyy
%     month = 1; %mm
%     day = 9; %dd
%     hour = 19; %HH
%     minute = 25; %MM
%     second = 0; %SS.FFF
    
[start_day, start_fraction] = jday(year,month,day,hour,minute,second);
%% for Lat and Lon in Visual
a=6378137.0;
finv=298.257222101;
f=1/finv;
b=a*(1-f);
e2=1-(1-f)^2;
%% SGP4 setup
[satrec] = SCtwoline2rv(longstr1, longstr2, opsmode, whichconst); %initialise for SGP4


sat_epoch = satrec.jdsatepoch + satrec.jdsatepochf;
JD_0 = start_day+start_fraction;% - 2400000.5% + (32.184+37)/60/60/24; %epoch time used in the simulink model
difference = (JD_0-sat_epoch)*24*60; %minutes

data = zeros(duration/step,7); %initialize the array
Orbit = zeros(duration/step,3);

for time=difference:step/60:difference+duration/60 %minutes
JD = time/60/24+sat_epoch;
[satrec, r, v] = sgp4(satrec,time); %translates tle into positions and velocities over the simulation period
timestep = round((time-difference)/step*60); %position in array, has to be an integer
data(1+timestep,:) = [JD r v];

%for visual
[secef, recef, vecef, aecef, DCM_ET] = teme2ecef([0;0;0],[r(1)*1000;r(2)*1000;r(3)*1000],[0;0;0],[0;0;0],(JD-2451545.0)/36525,JD,0,0.127631*pi/3600/180,0.247746*pi/3600/180,2);
X = recef(1);%vector in ecef
Y = recef(2);
Z = recef(3);
lon=atan2(Y,X);
e=e2*(a/b)^2;
p=sqrt(X.*X+Y.*Y);
r2=sqrt(p.*p+Z.*Z);
u=atan(b.*Z.*(1+e.*b./r2)./(a.*p));
lat=atan((Z+e.*b.*sin(u).^3)./(p-e2.*a.*cos(u).^3));

Orbit(1+timestep,:)=[JD lat*180/pi lon*180/pi];

end

breakpoints = data(:,1);
     datapx = data(:,2);
     datapy = data(:,3);
     datapz = data(:,4);
     datavx = data(:,5);
     datavy = data(:,6);
     datavz = data(:,7);
%% initial states
W0 =[1;0;1-1/240].*pi()/180; %initial angular rate rad/s
DCM_BT0 = eye(3);
%           [  0.6666667,  0.3333333, -0.6666667;...
%              0.3333333,  0.6666667,  0.6666667;...
%              0.6666667, -0.6666667,  0.3333333 ]; %initial satellite orientation
%% Sun almanac
load nut80.dat; %for sun
%% Moment of Inertia
Ixx = 0.04582275;
Iyy = 0.04561065;
Izz = 0.0085932349;
%product of inertia %kg*m^2
Ixy = 0.000947555;
Ixz = -0.000072545546;
Iyz = -0.000031658756;

I = [Ixx -Ixy -Ixz;...
     -Ixy Iyy -Iyz;...
     -Ixz -Iyz Izz] ;

%controller coded inertia
C_I = [0.045 0.0009 -0.00007;...
       0.0009 0.045 -0.00003;...
       -0.00007 -0.00003 0.009];
%% Disturbances
%Atmospheric Drag

A = 0.2*0.3; %m^2
Cd = 2.0; %dimensionless
L = [0.0000812 -0.0022136 -0.0225868]; %m
V_Vec = [7600 1 0]; %m/s / Not the real reference frame

%Gravity Gradient

Ra = 6721000;
G = 6.67408e-11; %Gravitational Constant (m^3*kg^-1*s^-2)
M = 5.9722e24;  %Mass of Earth (kg)

%Solar Radiation Pressure

As = 0.2*0.3 ; %Sunlit Area
c = 2.99792458e8; %Speed of Light (m/s)
q = 0.6; %Unitless Reflectance Factor (Varies from 0 to 1)
Ks = 1367; %Solar Constant (W/m^2)
Ls =0.01*[-3; 3; 5]; %Vector from center of mass to center of pressure

D=0.01*[1;1;1];%overall residual dipole
%% Magnetorquers
magnetorquer_x_saturation = 0.24;
magnetorquer_y_saturation = 0.24;
magnetorquer_z_saturation = 0.13;
magnetorquer_tau = 0.02;
%% Magnetometer

radii = [1;0.9;1.1];%Ellipsoid radiuses (Shape of the ellipsoid)
Ellrot = [   0.1022443, -0.0217327, -0.9945219;
            -0.7026129,  0.7061515, -0.0876649;
             0.7041884,  0.7077271,  0.0569303 ]; %Ellipsoid rotation
scale=[radii(1),0,0;
       0,radii(2),0;
       0,0,radii(3)]/mean(radii);
invA_mag = inv(Ellrot/scale*Ellrot'); %Magnetometer stretch (soft iron)
b_mag=[10;-15.5;-12.2]*1e-6; %Magnetometer Bias (hard iron)

% Calculated values from last calibration
load('MagnetometerBias.mat')%get the latest values for calibration parameters
% Amag=inv(invA_mag); %uncomment to get the perfect values
% bmag=b_mag; %uncomment to get the perfect values
%% Reaction Wheels
wheel_saturation = 0.23e-3; %Nm
wheel_tau = 2; %seconds
%% Controller
controller_sample_rate = 0.5; %seconds
Clockerror0=0.1;
Clockerror=Clockerror0;
ClockerrorRamp=3.5e-1;
%% Reset SImulink blocks
try
set_param('ADCS_Outline/Guidance&Determination/Mode Selection/Constant3','Value','0') %reset modeOverdrive
set_param('ADCS_Outline/VISUALISATION/Minimap','Value','0') %reset minimap
set_param('ADCS_Outline/Pulse Generator','Period','1') %reset sample time
set_param('ADCS_Outline/Real-Time Pacer','MaskValueString','inf')
catch
end
%% sun sensors in body frame

% sun_sensors = [ 0       1       0     ;...%1
%                -0.8944  0.4472  0     ;...%2
%                -0.2764  0.4472  0.8506;...%3
%                 0.7236  0.4472  0.5257;...%4
%                -0.2764  0.4472 -0.8506;...%5
%                 0.7236  0.4472 -0.5257;...%6
%                -0.7236 -0.4472 -0.5257;...%7
%                -0.7236 -0.4472  0.5257;...%8
%                 0.2764 -0.4472  0.8506;...%9
%                 0.8944 -0.4472  0     ;...%10
%                 0.2764 -0.4472 -0.8506;...%11
%                 0      -1       0     ];  %12