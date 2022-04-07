%%%The following script examines and integrates Plank's law to determine
%%%radiation intensity vs. wavelength and total solar intensity

%%%Cleaning current workspace
clc, clear all, close all
%% Defining constants 
%%%Planck's constant h
h = 6.67259e-34; %%J*s

%%%Boltzmann Constant k
k = 1.380658e-23; %J/K

%%%Speed of Light c
c = 2.99792458e8; %m/s

%%%Wavelength - wavelength vector stops arbitrarily at 20 thousand nm.
%%%Radiation is negligible at these wavelengths for temperatures of
%%%interest
wvlgth = [0.1 1:1:20000]/10^9; %m

%% calculating and plotting solar intensity as a function of wavelength
close all
%%%Temperature of interest going from blackbody radiation 
T = [3000 4000 5000 6000 ]; %K - All temperatures will show up as absolute 
%%%Running loop to predict radiative power intensity  for each temperature
%%%Equation the wavelength form taken from https://en.wikipedia.org/wiki/Planck%27s_law
%%%Multiplied by 4pi, or the solid angle of a sphere : https://en.wikipedia.org/wiki/Solid_angle
for i = 1:length(T)
    P(i,:) = 8.*pi.*c^2./(wvlgth.^5).*(1./(exp(h.*c./(k.*T(i).*wvlgth))-1));
    %%%String that goes into legend
    leg_str{i} = num2str(T(i));
end

%%%Plotting normalized results - normalization occurs in the plotting
%%%function. To plot absolute values, plot the commented line below
figure
hold on
% semilogy(wvlgth.*10^9,P./max(P,[],2),'LineWidth',2)
semilogy(wvlgth.*10^9,P,'LineWidth',2)
box on
set(gca,'FontSize',14)
xlabel('Wavelength (nm)','FontSize',16)
ylabel('Radiative Power','FontSize',16) %%%Units are not provided because this could be a normalized value (no units) or radiative power (W/m)
grid on
xlim([0 0.5.*10^4])
legend(leg_str)
hold off

%% Integrating Planck's law

%%%Specifying temperature of body
T = 3000; %K

%%%Defining material emissivity
em = 1; %%% 1 = black body
%%%Defining planck's law as a function where w is the input variable in
%%%meters
fun = @(w) em*8.*pi.*c^2./(w.^5).*(1./(exp(h.*c./(k.*T.*w))-1));

%%%Integrating Planck's law between specified limits
q = integral(fun,3.8e-7,7.5e-7)

%%%Energy contained in those limits from two lines ago
q/integral(fun,0,inf)

