clc, clear all, close all

 %% Adjust concentration here to tackle Problem 1_10
 C = 1;
 
%% Constants

%%%Plank's Constant
h = 6.6260755e-34; %J*s

%%%Boltzmann Constant
k = 1.380658e-23; %J/K

%%%Speed of light
c = 2.99792456e8; %m/s

%%%Charge of electron
q = 1.60217733e-19; %C

%%%Radiaiton Source Temperature
T = 6000; %K

%%%Radiation source intensity
P = 1000*C; %W/m^2

%% Cell properties

%%%Band gaps %All in volts
Vg_1 = 1.1; %Volts - bandgap of lower energy semiconductor

%%%Cell Area
area = 0.0001; %m^2 - converted from cm^2 to m^2
% area = 1;

%%%Leakage current %%%Normalizing leakage current by area
Io = 0.5e-12/0.0001; %Amps, converted form pA to A
Io = 4e-10;


%%%Recombination current %%%Normalizing recombinaiton current by area
Ir = 4e-8/0.0001; %Camps, converted from microA to A
% Ir = 0; %Camps, converted from microA to A %Used to find efficiency difference due to recombination

%%%Cell Temperature
T_cell = 300;

%% Counting Electrons Incident on our material

%%%Fraction of photons at or above specific energy level
%%%Integral portion of Equaiton 14.25
fun1 = @(X) X.^2./(exp(X)-1);

%%%Integral limits
lim_1 = q*Vg_1/(k*T);
lim_2 = inf;

%%%Equation 14.25 - (%)
sig = 0.416*integral(fun1,lim_1,lim_2);

%%%Total number of photons using Equation 14.16
phi = P./(37.28e-24*T); %Photons/m^2/s

%%%Photons at or above specific energy level
phi_1 = sig*phi;%Photons/m^2/s

%%%Photons incident on panel
PHI = phi_1; %Photons/s

%%%Induced current -Equation 14.41
I_v = PHI*q; %C/s 

%% Current Sweep

%%%First step in voltage sweep - Voltage = 0
V = 0;
%%%Delta in voltage level
V_step = 0.001;
%%%Counter
idx = 1;
%%%Current at first time step
I(idx) = I_v ... %%%Induced current 
    - Ir*(exp(q*V/(2*k*T_cell)) - 1) ... %%%Recombination current
    - Io*(exp(q*V/(k*T_cell))-1); %%%Diode leakage current
 I_comp(idx,:) = [Ir*(exp(q*V(idx)/(2*k*T_cell)) - 1)  Io*(exp(q*V(idx)/(k*T_cell))-1)];

%%%Increasing counter for next step
idx = idx + 1;


%%%Increase and evaluate cell performance as long as current is positive
while  I(idx-1) > 0 %%%Keep going until negative current is induced
    
    
    %%%Increasing voltage for next step (Volts)
    V(idx) = V(idx-1) + V_step;
    
    %%%Check current at each step (Amps/m^2)
    I(idx) = I_v ... %%%Induced current
        - Ir*(exp(q*V(idx)/(2*k*T_cell)) - 1) ... %%%Recombination current
        - Io*(exp(q*V(idx)/(k*T_cell))-1); %%%Diode leakage current
    
    I_comp(idx,:) = [Ir*(exp(q*V(idx)/(2*k*T_cell)) - 1)  Io*(exp(q*V(idx)/(k*T_cell))-1)];
     
    %%%Increasing counter for next step
    idx = idx + 1;
    
end

%%%Adjusting current fron per unit area to cell area
I = I*area;

%%%Cell Power (watts)
p = V.*I;

%%%Cell efficiency (eff)
eff = p/(P*area);

%%%Plot V-i curve
figure
yyaxis left
hold on
plot(I,V)
plot(I,eff)
ylim([0 max(V)*1.05])
xlim([0 max(I)*1.05])
ylabel('Cell Voltage (V) & Efficiency','FontSize',18)
xlabel('Current (amps)','FontSize',18)


yyaxis right
plot(I,p)
ylabel('Power (W)','FontSize',18)
legend('Voltage','Efficiency','Power','Location','Best')
box on
grid on
hold off
