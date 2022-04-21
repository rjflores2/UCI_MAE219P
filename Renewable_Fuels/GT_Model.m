clear all, clc
%% Engine informaiton

%%%Compressor isentropic efficiency
comp_eff = 0.85;

%%%Turbine isentropic efficiency
turb_eff = 0.9;

%%%Pressure ratio
p_r = 16;

%%%Inlet airflow rate
mflow = 40; %kg/s
%% Combustion properties

%%%Methane combustion
%%%Stoichiometric coefficient a (CH4 + 2O2 --> CO2 + 2H2O)
a_ng = 2;
%%%Excess air
ea_ng = 2;

%%%H2 combustion
%%%Stoichiometric coefficient a  (H2 + 1/2 O2 --> H2O)
a_h2 = 0.5;

%%%Biogas combustion - Biogas concentrations normalized to 1 kmol fuel = CH4 + 2/3 CO2
%%%Stoichiometric coefficient a  (CH4 + 2/3 CO2 + 2O2 --> 5/3 CO2 + 2H2O)
a_bg = 2;

%%%Syngas combustion - Syngas concentrations normalized to 1 kmol fuel = 0.67 CO + 0.33 H2 + 0.11 CO2
%%%Stoichiometric coefficient a (0.67 CO + 0.33 H2 + 0.11 CO2 + 1 O2 --> 0.77 CO2 + 0.33 H2O
a_sg = 1;

%%%Excess air range
ea = [1:0.01:3];
%% Setting species info
%%%Gas properties at state 1
gas1 = IdealGasMix('gri30.yaml');

MW = molecularWeights(gas1);

nsp = nSpecies(gas1);

%%%Find air species
io2 = speciesIndex(gas1,'O2');
in2  = speciesIndex(gas1,'N2');

%%%Find fuel species
ih2 = speciesIndex(gas1,'H2');
ico = speciesIndex(gas1,'CO');
ico2 = speciesIndex(gas1,'CO2');
ich4 = speciesIndex(gas1,'CH4');

%%%Pollutant species - NOx
ino = speciesIndex(gas1,'NO');
ino2 = speciesIndex(gas1,'NO2');

%% Defining properties at State Point 1

%%%Inlet temperature
T1 = 300;

%%%Inlet Species
x = zeros(nsp,1);
x(io2,1) = 1;
x(in2,1) = 3.76;

%%%Setting gas properties at state 1
set(gas1,'Temperature',T1,'Pressure',101325.0,'MoleFractions',x);

%%%Inlet enthalpy
h1 = enthalpy_mass(gas1);
%%%Inlet entropy
s1 = entropy_mass(gas1);

%%%Molecular weight at state 1
MW_air = sum(MW.*(x./sum(x)));
%% Compression & Defining properties at State Point 2
%%%Increasing gas pressure at constant entropy
set(gas1,'Pressure',101325.0*p_r,'S',s1)
temperature(gas1)

%%%Compressor outlet enthalpy - isentropic compression
h2s = enthalpy_mass(gas1);

%%%Actual enthalpy  comp_eff = (h2s - h1) / (h2 - h1)
h2 = h1 + (h2s - h1)/comp_eff;

%%%Resetting gas at actual enthalpy
set(gas1,'Pressure',101325.0*p_r,'H',h2)

%%%Actual outlet temperature
T2 = temperature(gas1);

%% Combustion of Methane

%%% Combsution gases
gas2 = IdealGasMix('gri30.yaml');

%%Combustion reactants
x = zeros(nsp,1);
x(ich4,1) = 1;
x(io2,1) = (1+ea_ng)*a_ng*1;
x(in2,1) = (1+ea_ng)*a_ng*3.76;

rat_ch4 = x(ich4,1)/(x(io2,1) + x(in2,1));

%%%Setting combusiton gases
set(gas2,'Temperature',T2,'Pressure',101325.0*p_r,'MoleFractions',x);

%%%React mixture to equilibrium
equilibrate(gas2,'HP');
cp_ch4 = cp_mass(gas2)/1000;

%%%Combustion temperature
T3 = temperature(gas2);
%%%Species concentration
x3 = moleFractions(gas2);
%%%Enthalpy
h3 = enthalpy_mass(gas2);
%%%entropy
s3 = entropy_mass(gas2);

%%%Species
xch4 = moleFractions(gas2);

%%%Molecular weight at state 1
MW_air = sum(molecularWeights(gas1).*(xch4./sum(xch4)));

%%%Massflow rate of CH4
ch4_mflow = mflow/MW_air*rat_ch4*MW(ich4);

%%%Heat in -  LHV = 50 MJ/ kg CH4
heat_in_ch4 =ch4_mflow*50000;
%% Expansion & Definint State Point 4

%%%Resetting gas at actual enthalpy
set(gas2,'Pressure',101325.0,'S',s3)

%%%Turbine outlet enthalpy - isentropic compression
h4s = enthalpy_mass(gas2);

%%%Actual Enthalpy
h4 = h3 - (h3 - h4s)*turb_eff;


%%%Resetting gas at actual enthalpy
set(gas2,'Pressure',101325.0,'H',h4)

%%%Actual outlet temperature
T4 = temperature(gas2);

%% Methane turbine
%%%Work output
w_net = ((h1 - h2)*mflow + (h3 - h4)*(mflow + ch4_mflow))/1000
%%%Efficiency
efficiency = w_net/heat_in_ch4
(h3 - h4)*50/1000 + (h1 - h2)*50/1000

%%%Pollutant emissions
pollutant_concentrations(1,:) = [xch4(ico) xch4(ino)+xch4(ino2)];

%% Examining the combustion of H2 - Hw3-3
%%% Combsution gases
gas_h2 = IdealGasMix('gri30.yaml');

for ii = 1:length(ea)
    
    %%Combustion reactants
    x = zeros(nsp,1);
    x(ih2,1) = 1;
    x(io2,1) = (1 + ea(ii))*a_h2*1;
    x(in2,1) = (1 + ea(ii))*a_h2*3.76;
    
    %%%Setting combusiton gases
    set(gas_h2,'Temperature',T2,'Pressure',101325.0*p_r,'MoleFractions',x);
    %%%React mixture to equilibrium
    equilibrate(gas_h2,'HP');
    
    T3_h2(ii,1) = temperature(gas_h2);
    xh2(:,ii) = moleFractions(gas_h2);
end

%%%Which EA gets to the closest temperature?
[~,idx_h2] = min((T3_h2 - T3).^2);
% [~,idx_h2] = min(abs(T3_h2 - T3));
ea_h2 = ea(idx_h2);
xh2 = xh2(:,idx_h2);


%%%Pollutant emissions
pollutant_concentrations(2,:) = [xh2(ico) xh2(ino)+xh2(ino2)];



%% change in pollutant species
for ii = 1:size(pollutant_concentrations,1)
    pollutant_concentrations_delta(ii,:) = 100.*(pollutant_concentrations(ii,:) - pollutant_concentrations(1,:))./pollutant_concentrations(1,:);
end
%% Examining the operation of turbine using alternative fuels
%% Hydrogen
%%% Combsution gases
gas_h2 = IdealGasMix('gri30.yaml'); %%Combustion reactants
x = zeros(nsp,1);
x(ih2,1) = 1;
x(io2,1) = (1 + ea_h2)*a_h2*1;
x(in2,1) = (1 + ea_h2)*a_h2*3.76;

rat_h2 = x(ih2,1)/(x(io2,1) + x(in2,1));

%%%Setting combusiton gases
set(gas_h2,'Temperature',T2,'Pressure',101325.0*p_r,'MoleFractions',x);
%%%React mixture to equilibrium
equilibrate(gas_h2,'HP');
cp_h2 = cp_mass(gas_h2)/1000;


%%%Combustion temperature
T3_h2 = temperature(gas_h2);
%%%Species concentration
x3_h2 = moleFractions(gas_h2);
%%%Enthalpoy
h3_h2 = enthalpy_mass(gas_h2);
%%%entropy
s3_h2 = entropy_mass(gas_h2);

%%%Massflow rate of CH4
h2_mflow = mflow/MW_air*rat_h2*MW(ih2);

%%%Heat in -  LHV = 120 MJ/ kg H2
heat_in_h2 = h2_mflow*120000;

%%%Expansion in the turbine
%%%Resetting gas at actual enthalpy
set(gas_h2,'Pressure',101325.0,'S',s3_h2)

%%%Turbine outlet enthalpy - isentropic compression
h4s_h2 = enthalpy_mass(gas_h2);

%%%Actual Enthalpy
h4_h2 = h3_h2 - (h3_h2 - h4s_h2)*turb_eff;


%%%Resetting gas at actual enthalpy
set(gas_h2,'Pressure',101325.0,'H',h4_h2)

%%%Actual outlet temperature
T4_h2 = temperature(gas_h2);

%%%H2 turbine
w_net_h2 = ((h1 - h2)*mflow + (h3_h2 - h4_h2)*(mflow + h2_mflow))/1000
efficiency_h2 = w_net_h2/heat_in_h2