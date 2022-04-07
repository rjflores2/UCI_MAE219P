%%%Predicting radiation falling on a surface. This script introduces
%%%equation 12.15. This script assumes that each day is 12 hours long, sun
%%%rise is at 6 am, and sun set is 6pm. These are obvious sources of error
close all, clear all, close all

%% Defining constants

%%%Radiative power
P = 1000; % W/m^2

%%%Days of the year in decimal form: Each integer value corresponds to a single day
%%% one hour = 1/24
%%% one minute = 1/(24*60)
%%% 15 seconds = 1/(24*60*4)
d = 0:1/(24):365; %Date integer

%%% Length of each time step - must keep time step constant
delta_t = (d(2) - d(1))*(24*60*60); %S

%%%Start of each day in decimal form
day_start = 6/24;
%%%End of each day in decimal form
day_end = 18/24;

%% Surface orientation and location
%%%Panel Tilt - insert in terms of degrees, output is in radianas
eta = 35; %Degrees 

%%%Azimuthal orientation Zeta- insert in terms of degrees, output is in radianas
zeta = 180; %0 deg = north, 90 deg = east, 180 = south

%%%Latitude of location
lat = 33.7; %Degrees
%% Incident power throughout the year
close all
%%%Allocating space for power vector
P_incident = zeros(size(d));
angles = zeros(length(d),5);
for ii = 1:length(d)
    
    %%%Hour angle
    alpha = (360/24)*(rem(d(ii),1)*24 - 12); %Degrees
    %%%Solar declinaiton
    sd = 23.44*sind(360/365.25*(d(ii)-80)); %Degrees
    %%%Zenith Angle
    zen = acosd(sind(sd)*sind(lat) + cosd(sd)*cosd(lat)*cosd(alpha));
    %%%Azimuth
    az = atand(sind(alpha)./ ...
        (sind(lat)*cosd(alpha) - cosd(lat)*tand(sd)));
    
    %%%Adjusting az value using logic provided in Section 12.2.1
    if alpha>0 && tand(az) > 0
        az = az + 180;
        tr = 1;
    elseif alpha>0 && tand(az) < 0
        az = az + 360;
        tr = 2;
    elseif alpha < 0 && tand(az) > 0
        az = az + 0;
        tr = 3;
    else%if alpha < 0 && tand(az) < 0
        az = az + 180;
        tr = 4;
    end
    
    angles(ii,:) = [alpha sd  zen az tr];
    %%% if the sun is shining
    if  (cosd(eta)*cosd(zen) + sind(eta)*sind(zen)*cosd(az - zeta)) > 0
        %%%Incident power
        P_incident(ii) = P*(cosd(eta)*cosd(zen) + sind(eta)*sind(zen)*cosd(az - zeta));
    end
end

%%%Plotting annual incident radiation
figure
hold on
plot(d,P_incident)
xlim([0 365])
ylim([0 1050])
box on
% grid on
datetick('x','mmm')
set(gca,'FontSize',14)
ylabel('Incident Radiative Power (W/m^2)','FontSize',16)
set(gcf, 'Position',  [50, 50, 1000, 450])
hold off

%%%Reiman sum of incident energy
energy = sum(P_incident.*delta_t)... %%% J/m^2/year
    /1000 %Convert to kJ/m^2/year