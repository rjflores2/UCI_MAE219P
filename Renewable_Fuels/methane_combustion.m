%% Cantera Demo
clc, clear all, close all
%% Combustion of Methane
%%% This example shows the chemical equilibrium associated with methane
%%% combustion at different air/fuel ratios
%%% State 1 is an air fuel mixture at 1 atm, 300K
%%% State 2 is after combustion in an open system
%%% Stoichiometric combustion of methane: (CH4 + 2O2 --> CO2 + 2H2O) --> Stoichiometric coefficient  a = 2
a_ng = 2;

%%%Declare gas model - GRI-Mech 3.0
%%%GRI-Mech 3.0 is an old set of chemical equilibrium and kinetics models
%%%designed to capture the combustion of methane. More info at: http://combustion.berkeley.edu/gri-mech/version30/text30.html
% gas = GRI30;
gas = Solution('GRI30.yaml');

%%%Find air species
io2 = speciesIndex(gas,'O2');
in2  = speciesIndex(gas,'N2');

%%%Species of complete combusiton
ico2 = speciesIndex(gas,'CO2');
ih2o = speciesIndex(gas,'H2O');

%%%Find fuel species 
ich4 = speciesIndex(gas,'CH4');

%%%Pollutant species - NOx
ino = speciesIndex(gas,'NO');
ino2 = speciesIndex(gas,'NO2');

%%%Species that will kill you - CO
ico = speciesIndex(gas,'CO');

%%%Number of species in gas model
nsp = nSpecies(gas);

%% Looping through states
phi = [0.5:0.025:2]; %%%Equivalence ratios to be examined

for ii = 1:length(phi)
    %%%Defining Mole Fractions
    x = zeros(nsp,1);
    %%%Setting species concentrations
    x(ich4) = 1;
    x(io2) = 2*(1/phi(ii));
    x(in2) = 2*3.76*(1/phi(ii));
    
    %%%Setting gas species
    set(gas,'Temperature',300,'Pressure',101325,'MoleFractions',x)
    
    %%%React the system at constant pressure and enthalpy
    equilibrate(gas,'HP');
    
    %%%Record items of interest
    T2(ii,:) = temperature(gas); %%%Temperature (K)
    x2(ii,:) = moleFractions(gas); %%%Species (%volume or %mol)
    
end
%% Plotting Results
close all

figure
hold on
plot(phi,T2,'LineWidth',2)
box on
grid on
set(gca,'FontSize',14)
ylabel('Temperature (Kelvin)','FontSize',16)
xlabel('Equivalent Ratio \phi','FontSize',16)
xlim([min(phi) max(phi)])
hold off

figure
hold on
plot(phi,x2(:,[io2 in2 ico2 ih2o]),'LineWidth',2)
box on
grid on
set(gca,'FontSize',14)
xlabel('Equivalent Ratio \phi','FontSize',16)
ylabel('Species Concentrations (% Volume)','FontSize',16)
title('Complete Combusiton Species','FontSize',16)
legend('O_2','N_2','CO_2','H_2O')
xlim([min(phi) max(phi)])
hold off

figure
hold on
plot(phi,x2(:,[ich4 ico]),'LineWidth',2)
box on
grid on
set(gca,'FontSize',14)
xlabel('Equivalent Ratio \phi','FontSize',16)
ylabel('Species Concentrations (% Volume)','FontSize',16)
title('Unused Fuel','FontSize',16)
legend('CH_4','CO')
xlim([min(phi) max(phi)])
hold off


figure
hold on
plot(phi,x2(:,[ino ino2]),'LineWidth',2)
box on
grid on
set(gca,'FontSize',14)
xlabel('Equivalent Ratio \phi','FontSize',16)
ylabel('Species Concentrations (% Volume)','FontSize',16)
title('O_3 Precursors (NO_x)','FontSize',16)
legend('NO','NO_2')
xlim([min(phi) max(phi)])
hold off