% Calculate irradiance in equivalent 460 nm for melanospin 

% Load melanopsin adsorption (nomogram from Mike's mimicking_sunlight.pxp,
% for a mixture of the R and E states)
absorp = readtable('melanopsin RE.txt');
absorp.wavelength = (300:(300 + length(absorp.nomo) - 1))';

% Load the measured stimulus spectrum
%spec = readtable('460 LED staircase step 6.xlsx');
%spec = readtable('step 7_00018.txt');
latestfile = getlatestfile('C:\Users\Franklin\Desktop\pupil\20220615 opn4 het and ko persistence\calibration\');
spec = readtable(latestfile);
spec.Properties.VariableNames = {'nm', 'uW_per_cm2'};


% Interpolate each function so they are defined at the same wavelength
% values 
wavelength = 400:0.1:780;
absorption = interp1(absorp.wavelength, absorp.nomo, wavelength); 
emission = interp1(spec.nm, spec.uW_per_cm2, wavelength); 


% Plot original and interpolated spectra
a = figure;
plot(absorp.wavelength, absorp.nomo, 'k')
hold on 
plot(wavelength, absorption, 'r');
xlabel('Wavelength (nm)')
ylabel('Normalized Sensitivity');
box off 
set(gca, 'TickDir', 'out');
title('Melanopsin (R/E) Absorption');
legend({'Original', 'Interp.'});
xlim([min(wavelength), max(wavelength)]);

s = figure;
plot(spec.nm, spec.uW_per_cm2, 'k')
hold on 
plot(wavelength, emission, 'r');
xlabel('Wavelength (nm)')
ylabel('Power (uW*cm^{-2})');
box off 
set(gca, 'TickDir', 'out');
title('LED Emission');
legend({'Original', 'Interp.'});
xlim([min(wavelength), max(wavelength)]);


% Convolve emission with normalized absorption  
scaled_emission = absorption .* emission; 


% Convert power to irradiance
single_photon_energy = 6.626*10^(-34)*2.9979*10^8 ./(wavelength*10^(-9)); % Joules(wavelength) 
irradiance = (scaled_emission * 10^-6) ./ single_photon_energy * 10^-8; % at each wavelength photons / um^2 / s


% Integrate to obtain equivalent monochromatic irradiance 
integrated_irradiance = trapz(wavelength, irradiance); 
[~, idx] = max(absorp.nomo);
disp(['melanopsin-equivalent ', num2str(absorp.wavelength(idx)), ' nm photons:']);
disp([num2str(integrated_irradiance), ' photons * um^-2 * s^-1']);
disp([num2str(log10(integrated_irradiance)), ' Log units']);


