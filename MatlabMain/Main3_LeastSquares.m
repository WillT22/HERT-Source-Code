% Least Squares and Count Rate Analysis for HERT
% Last modified: 5/29/2024

% finding initial flux guess from bowtie analysis
initial_flux = j_nom;

j_lsqr = lsqr(geo_EC,hits_whole_EC',[],100);


% Load data
count_rate = hits_whole_EC;
geo_factor = geo_EC;
bowtie_flux = j_nom;

%% Least Squares Functions
% Define least-squares function (anonymous function)
func = @(flux, count_rate, geo_factor) trapz(M_energy_midpoints, flux .* geo_factor');

% Define chi-square function
chisq = @(flux, count_rate, geo_factor) count_rate - func(flux, count_rate, geo_factor);

% Use 'lsqnonlin' for least-squares minimization (similar to 'lm' in scipy)
channel_flux = lsqnonlin(@(flux) chisq(flux,count_rate,geo_factor), bowtie_flux);


%% Plots calculated flux
lsqr_flux = load('E:\HERT_Drive\Matlab Main\Bow Tie\lsqr_flux.txt');
f = figure;
f.Position = [100 100 1200 720];
hold on

plot(M_energy_midpoints,M_energy_bin/(4 * (pi^2) * (r_source^2)),'*', 'Color', 'black');
plot(E_eff,j_nom,'o', 'Color', '#0072BD');
plot(E_eff,lsqr_flux,'x', 'Color', '#D95319');
plot(E_eff,channel_flux,'square', 'Color', '#77AC30');
plot(M_energy_midpoints,j_lsqr,'square', 'Color', '#77AC30');

legend({['Incident Particle Measurements'],['Bowtie Method'],['Python Least Squares'], ...
                 ['MATLAB Least Squares']}, 'Location', 'northeast');

xlim([1 7])
xticks((1:1:7))
set(gca, 'YScale', 'log')
%ylim([10^5, 10^5.5])

ylabel('Flux')
%ylabel('I_{nom} (#/(cm^2 sr s MeV)')
xlabel('Nominal Effective Energy (MeV)')
hold off
