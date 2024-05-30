% Least Squares and Count Rate Analysis for HERT
% Last modified: 5/29/2024


%% Least Squares Function for Energy Channels

% Define least-squares function (anonymous function)
func = @(flux, count_rate, geo_factor) trapz(M_energy_midpoints, flux' .* geo_factor');

% Define chi-square function
chisq = @(flux, count_rate, geo_factor) count_rate - func(flux, count_rate, geo_factor);

% Calculating approximation from count_rate/geo_factor
flux_approx = zeros(1,length(geo_total));
for bin = 1:length(geo_total)
    if geo_total(bin) > 0
        flux_approx(bin) = hits_log_total(bin)/geo_total(bin)/bin_width;
    else
        flux_approx(bin) = 10^-20;
    end
end

% Use 'lsqnonlin' for least-squares minimization (similar to 'lm' in scipy)
channel_flux = lsqnonlin(@(flux) chisq(flux,hits_whole_EC,geo_EC), flux_approx);


%% Plots calculated flux
f = figure;
f.Position = [100 100 1200 720];
hold on

plot(M_energy_midpoints,M_energy_bin/(4 * (pi^2) * (r_source^2))/bin_width, 'Color', 'black','LineWidth',2);
plot(E_eff,j_nom,'o', 'Color', '#0072BD');
plot(M_energy_midpoints(channel_flux>1),...
    channel_flux(channel_flux>1),'.', 'Color', 'r', 'MarkerSize',6);

legend({['Incident Particle Measurements'],['Bowtie Method'], ...
                 ['Least Squares Method']}, 'Location', 'northeast');

xlim([0 8])
xticks((0:1:8))
set(gca, 'YScale', 'log')
ylim([10^5.8, 10^6.2])

ylabel('Flux')
%ylabel('I_{nom} (#/(cm^2 sr s MeV)')
xlabel('Nominal Effective Energy (MeV)')
hold off
