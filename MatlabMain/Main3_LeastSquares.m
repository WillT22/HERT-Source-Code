% Least Squares and Count Rate Analysis for HERT
% Last modified: 5/29/2024

%% Inputs
geo_EC = readmatrix('E:\HERT_Drive\Matlab Main\Result\geometric_factor_EC.txt');
dt = 1;

%% Least Squares Function for Energy Channels (initial attempt)
flux_approx = lsqr(geo_EC,hits_whole_EC',1e-6,20)'./bin_width;

%% Least Squares Function for Energy Channels (Selesnick/Khoo)
% Initialize variables
m0 = log(flux_approx');
m0_int = zeros(size(energy_channels,1),1);
d0 = zeros(size(energy_channels,1),1);
for channel = 1:size(energy_channels,1)
    m0_int(channel) = sum(trapz(energy_midpoints,geo_EC(channel,:)'.*exp(m0 + log(energy_midpoints'))));
    d0(channel) = log(m0_int(channel));
end
g_n = 1./hits_whole_EC' .* m0_int;
g_m0 = log(m0_int);

C_d = zeros(size(energy_channels,1),size(energy_channels,1));
for channel = 1:size(energy_channels,1)
    C_d(channel,channel) = (hits_whole_EC(channel)*dt+1)/dt^2;
end
inv_C_d = inv(C_d);

C_m = 2;
inv_C_m = inv(C_m);

% Initialize iteration variables
mn = m0;
mn_int = m0_int;
dn = d0;
error = mean(abs(log(hits_whole_EC')-dn));
iterations = 1;

% Iterate to desired tolerance or iterations
while error >= 1e-8 && iterations <= 10000
    mn = m0 +  inv(g_n' * inv_C_d * g_n + inv_C_m) * g_n' * inv_C_d * (log(hits_whole_EC') - g_m0);
    for channel = 1:size(energy_channels,1)
        mn_int(channel) = sum(trapz(energy_midpoints,geo_EC(channel,:)'.*exp(mn + log(energy_midpoints'))));
        dn(channel) = log(mn_int(channel));
    end
    error = mean(abs(log(hits_whole_EC')-dn));
    iterations = iterations + 1;
end
fprintf("Iteration Number: %.0d \n",iterations-1)
fprintf("Tolerance: %.6e \n",error)

flux_lsqr = exp(mn);

%% Error Analysis
flux_actual = M_energy_bin./(4 * (pi^2) * (r_source^2))./bin_width;

% Choose appropriate values
app_ind = energy_midpoints>=1;

% Find errors of summation versus integration techniques
approx_error = abs(flux_approx(app_ind)-flux_actual(app_ind))./flux_actual(app_ind);
avg_approx_error = mean(approx_error);
lsqr_error = abs(flux_lsqr(app_ind)-flux_actual(app_ind))./flux_actual(app_ind);
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
plot(energy_midpoints,M_energy_bin/(4*pi^2*r_source^2)/bin_width, 'Color', 'black','LineWidth',2);

%{
% Plot Bowtie points and line of best fit
plot(E_eff,j_nom,'o', 'Color', '#0072BD','MarkerSize',10);
fit_result = fit(E_eff',j_nom','exp1');
fitted_curve = @(x) feval(fit_result, x);
plot(E_eff, fitted_curve(E_eff), 'LineStyle', '--','LineWidth',2,'Color', '#0072BD');
%}

% Plot LSQ Integration Method
plot(energy_midpoints,flux_lsqr,'.', 'Color', 'r', 'MarkerSize',12);

% Plot LSQ Summation Method and line of best fit
plot(energy_midpoints,flux_approx,'x','Color','#77AC30','MarkerSize',10)
fit_result = fit(energy_midpoints',flux_approx','exp1');
fitted_curve = @(x) feval(fit_result, x);
plot(energy_midpoints, fitted_curve(energy_midpoints),...
    'LineStyle', '--','LineWidth',2,'Color', '#77AC30');

%{
legend({['Incident Particle Measurements'],['Bowtie Method'],['Bowtie Method (Best Fit)'], ...
                 ['Least Squares Method (Integral)'], ['Least Squares Method (Summation)'],['Least Squares Method (Summation Best Fit)']},...
                 'Location', 'northeast','FontSize',18);
%}
legend({['Incident Particle Measurements'], ...
                 ['Least Squares Method (Selesnick/Khoo)'], ['Least Squares Method (Summation)'],['Least Squares Method (Summation Best Fit)']},...
                 'Location', 'northeast','FontSize',18);
%

set(gca, 'FontSize', textsize)
xlim([0 8])
xticks((0:1:8))
set(gca, 'YScale', 'log')
%ylim([10^5.5, 10^6.5])

%ylabel('Flux')
ylabel('I (#/(cm^2 sr s MeV)','FontSize',textsize)
xlabel('Energy (MeV)','FontSize',textsize)
hold off
