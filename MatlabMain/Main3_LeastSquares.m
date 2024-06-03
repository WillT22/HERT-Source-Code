% Least Squares and Count Rate Analysis for HERT
% Last modified: 5/29/2024

%% Inputs
geo_EC = readmatrix('E:\HERT_Drive\Matlab Main\Result\geometric_factor_EC.txt');

%% Least Squares Function for Energy Channels

% Define least-squares function (anonymous function)
func = @(flux, count_rate, geo_factor) trapz(energy_midpoints, flux' .* geo_factor');

% Define chi-square function
chisq = @(flux, count_rate, geo_factor) count_rate - func(flux, count_rate, geo_factor);

flux_approx = lsqr(geo_EC,hits_whole_EC',1e-6,20)'./bin_width;

% Use 'lsqnonlin' for least-squares minimization (similar to 'lm' in scipy)
channel_flux = lsqnonlin(@(flux) chisq(flux,hits_whole_EC,geo_EC), flux_approx);

%% Error Analysis
flux_actual = M_energy_bin./(4 * (pi^2) * (r_source^2))./bin_width;

% Choose appropriate values
app_ind = energy_midpoints>=1;

% Find errors of summation versus integration techniques
approx_error = abs(flux_approx(app_ind)-flux_actual(app_ind))./flux_actual(app_ind);
avg_approx_error = mean(approx_error);
lsqr_error = abs(channel_flux(app_ind)-flux_actual(app_ind))./flux_actual(app_ind);
avg_lsqr_error = mean(lsqr_error);

error_count = zeros(1,sum(app_ind));
for bin = 1:sum(app_ind)
    if  lsqr_error(bin) > approx_error(bin)
        error_count(bin) = 1;
    end
end
int_sum_error = sum(error_count);


%% Plots calculated flux
f = figure;
f.Position = [0 0 1700 900];
hold on

% Plot simulated flux
plot(energy_midpoints,M_energy_bin/(4 * (pi^2) * (r_source^2))/bin_width, 'Color', 'black','LineWidth',2);

%
% Plot Bowtie points and line of best fit
plot(E_eff,j_nom,'o', 'Color', '#0072BD','MarkerSize',10);
fit_result = fit(E_eff',j_nom','exp1');
fitted_curve = @(x) feval(fit_result, x);
plot(E_eff, fitted_curve(E_eff), 'LineStyle', '--','LineWidth',2,'Color', '#0072BD');
%

% Plot LSQ Integration Method
plot(energy_midpoints(channel_flux>1),...
    channel_flux(channel_flux>1),'.', 'Color', 'r', 'MarkerSize',12);

% Plot LSQ Summation Method and line of best fit
plot(energy_midpoints(channel_flux>1),flux_approx(channel_flux>1),'x','Color','#77AC30','MarkerSize',10)
fit_result = fit(energy_midpoints(channel_flux>1)',flux_approx(channel_flux>1)','exp1');
fitted_curve = @(x) feval(fit_result, x);
plot(energy_midpoints(channel_flux>1), fitted_curve(energy_midpoints(channel_flux>1)),...
    'LineStyle', '--','LineWidth',2,'Color', '#77AC30');

%
legend({['Incident Particle Measurements'],['Bowtie Method'],['Bowtie Method (Best Fit)'], ...
                 ['Least Squares Method (Integral)'], ['Least Squares Method (Summation)'],['Least Squares Method (Summation Best Fit)']},...
                 'Location', 'northeast','FontSize',18);
%{
legend({['Incident Particle Measurements'], ...
                 ['Least Squares Method (Integral)'], ['Least Squares Method (Summation)'],['Least Squares Method (Summation Best Fit)']},...
                 'Location', 'northeast','FontSize',18);
%}

set(gca, 'FontSize', textsize)
xlim([0 8])
xticks((0:1:8))
set(gca, 'YScale', 'log')
%ylim([10^5.5, 10^6.5])

%ylabel('Flux')
ylabel('I (#/(cm^2 sr s MeV)','FontSize',textsize)
xlabel('Energy (MeV)','FontSize',textsize)
hold off
