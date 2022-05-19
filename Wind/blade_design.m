%%
clc, clear all, close all
%% Loading Wind Data - Tehachapi 2011

%%%Loading Data
[dt] = xlsread('H:\Classes\MAE_219P\Tehachapi_2011.csv');

%%%Eliminating NaN rows
dt = dt(~isnan(dt(:,1)),:);

%%%Splitting Data
%%%Time as serial vector
time = datenum([dt(:,1:5) zeros(length(dt),1)]);

%%%Time step (seconds)
delta = (time(2) - time(1))*(24*60*60);

%%%Wind Speed (m/s)
speed_t = dt(:,8);

%% Constants

%%Betz Limit
betz = 16/27;

%%%Assumed density of air
dens = 1.225; %kg/m^3

%%%Design Power
P = 3*1000^2; %Watts

%%%Design Tip Speed
ts_design = 8;

%%%Peak efficiency
eff_design = 0.8; %

%%%Blade number
bn = 3;

%%%Cut in fraction
cutin = 0.05;

%%%Cut out fraction
cutout = 1.5;


%% Design parameters
v_design = 10; %m/s

a = 13; %Angle of attack

%% General turbine design parameters
%%%Design power densigy (w/m^2 swept area)
p_design = 0.5*betz*dens*eff_design*v_design^3;

%%%Swept area (m^2)
sa = P/p_design;

%%%Blade length
bl = (sa/pi)^(1/2);

%%%Rotational Speed (Radians per second)
rot_speed = ts_design*v_design/bl;

%% Coefficients of Lift/Drag
CL = 0.6 + 0.066*a - 0.001*a^2 - 3.8e-5*a^3;
CD = .12 - .11*cosd(a);
%% Loop across different parts of the blade

%%%Blade discretizaiton
disc = 50; %%%Number of blade chunks considerd in the analysis

%%%Declaring blade area variable
blade_area = 0;

%%%Discretized blade results
for ii = 1:disc
    %%%Average radius of current blade zone
    r_i(ii) = (ii-1)*(bl/disc) + (bl/disc)/2; %%%m
    
    %%%Relative wind speed
    w_i = ((v_design*2/3)^2 + (r_i(ii)*rot_speed)^2)^(1/2); %m/s
    
    %%%Local Betz Power Limit
    betz_i(ii) = betz*.5*dens*v_design^3*2*pi*r_i(ii)*(bl/disc);
    
    %%%Constant
    multi(ii) = bn*rot_speed*0.5*dens*w_i*r_i(ii)*(bl/disc)*(CL*2/3*v_design - CD*r_i(ii)*rot_speed);
    
    %%%Chord Length
    K(ii) = betz_i(ii)/multi(ii);
    
    %%%Approximate blade area
    blade_area = blade_area + bn*K(ii)*(bl/disc);
    
    setup_angle(ii) = a - 180/pi*atan(2*bl/3/ts_design/r_i(ii));
    
   
end

%% Plotting some results
close all
figure
hold on
plot(r_i,K,'LineWidth',2)
box on
grid on
ylabel('Chord Size (m)','FontSize',18)
xlabel('Blade Position from Hub (m)','FontSize',18)
hold off


figure
hold on
plot(r_i,setup_angle,'LineWidth',2)
box on
grid on
ylabel('Setup Angle (Degree)','FontSize',18)
xlabel('Blade Position from Hub (m)','FontSize',18)
hold off

%% trash bin
%  %%%lift force
%     lift(ii) = (1/2)*dens*w_i^2*CL*(bn*K(ii)*(bl/disc))*sind(setup_angle(ii));
%     
%     %%%Drag Force
%     drag(ii) = (1/2)*dens*w_i^2*CD*(bn*K(ii)*(bl/disc))*cosd(setup_angle(ii));
%     
%     %%%Torque
%     torque(ii) = (lift(ii) - drag(ii))*r_i(ii);
%     
%     %%%Power
%     power(ii) = torque(ii)*bn*rot_speed;