%% Cantera Demo
clc, clear all, close all
%% Isentropic Compression in an open system from 1 atm to different pressure levels
%%%This example shows the isentropic compression of air from 1 atm, 300 K to
%%%different pressures ranging from 2 atm to 20 atm. State point 1 is the
%%%initial state (1 atm, 300 K), State point 2 is the state after compression.

%%%Declare gas model - GRI-Mech 3.0
%%%GRI-Mech 3.0 is an old set of chemical equilibrium and kinetics models
%%%designed to capture the combustion of methane. More info at: http://combustion.berkeley.edu/gri-mech/version30/text30.html
% gas = GRI30;
gas = Solution('GRI30.yaml');

%%%Find air species
io2 = speciesIndex(gas,'O2');
in2  = speciesIndex(gas,'N2');

%%%Number of species in gas model
nsp = nSpecies(gas);
%% Setting state 1 properties

%%%Temperature
T1 = 300; %Kelvin

%%%Pressure
P1 = 1; %atm - NOTE THAT CANTERA USES Pa AS THE BASE UNIT: https://cantera.org/documentation/dev/sphinx/html/yaml/general.html

%%%Air Species Concentrations
x = zeros(nsp,1);
x(io2,1) = 1;
x(in2,1) = 3.76;

%%%Setting gas at T and P
set(gas,'Temperature',T1,'Pressure',101325.0*P1,'MoleFractions',x);

%%%Initial enthalpy
h1 = enthalpy_mass(gas); % J/kg
%%%Initial entropy
s1 = entropy_mass(gas); % J/g*K

%% Evaluating compression

%%% Creating a for loop to cycle through different pressures
%%% - idx is a counter index that increases with each loop
%%% - P2 is the pressure range we are exmaining
idx = 1;
P2 = 2:20 %atm

for P = P2;  %atm
    %%%Resetting pressure at state 2 through an isentropic process
    set(gas,'S',s1,'Pressure',101325.0*P)
    
    %%%Extracting relevant thermodynamic info
    h2(idx,1) = enthalpy_mass(gas); %J/kg
    T2(idx,1) = temperature(gas); %K
    
    %%%Increasing the counter
    idx = idx + 1;
end

%% Plotting Results

%%%Plotting temperature
figure
hold on
plot(P2,T2,'LineWidth',2)
box on
grid on
set(gca,'FontSize',14)
ylabel('Temperature (Kelvin)','FontSize',16)
xlabel('Pressure After Compression (atm)','FontSize',16)
xlim([min(P2) max(P2)])
hold off

%%%Plotting specific work
figure
hold on
plot(P2,(h2 - h1)./1000,'LineWidth',2)
box on
grid on
set(gca,'FontSize',14)
ylabel('Specific Work (kJ/kg/s)','FontSize',16)
xlabel('Pressure After Compression (atm)','FontSize',16)
xlim([min(P2) max(P2)])
hold off