% Least Squares and Count Rate Analysis for HERT
% Last modified: 5/29/2024

%% Inputs
geo_EC = readmatrix('E:\HERT_Drive\Matlab Main\Result\geometric_factor_EC.txt');
dt = 1;

%% Least Squares Function for Energy Channels (initial attempt)
flux_approx = lsqr(geo_EC.*bin_width,hits_whole_EC',0.025,20)';

% using idicies for E>0.5
indices = energy_midpoints > 1;
fit_result = fit(energy_midpoints(indices)',flux_approx(indices)','exp2');

%% Least Squares Function for Energy Channels (Selesnick/Khoo)
% Initialize variables
    it_max = 10000;   % maximum number of possible iterations 

% Setting up constant matrices
    G = geo_EC.*bin_width;
    
    C_d = zeros(size(energy_channels,1));
    sigma_d = sqrt(hits_whole_EC(channel)*dt+1)/dt;
    for channel = 1:size(energy_channels,1)
        C_d(channel,channel) = sigma_d^2;
    end
    inv_C_d = inv(C_d);

    C_d_new = zeros(bins);
    for bin = 1:bins
        C_d_new(bin,bin) = sigma_d^2;
    end

    % Calculate delta variable
    fwhm = zeros(1,length(energy_channels));
    for u = 1:length(energy_channels)
        fwhm(u) = findFWHM(energy_midpoints,geo_EC(u,:));
    end
    [~, max_geo_EC] = max(geo_EC, [], 2);
    channel_peak = energy_midpoints(max_geo_EC);
    delta =2;

    % Initialize variance parameter
    sigma_flux = 105;
    %sigma_flux_array = logspace(0,3,2000);
    %sigma_flux_array = linspace(0,20,2001);

    % Initialize arrays to store results
    %chi_sq_array = zeros(length(sigma_flux_array),1);
    %norm_array = zeros(length(sigma_flux_array),1);

% Loop over variance parameter
%for sf_i = 1:length(sigma_flux_array)
    % Current sigma_flux value
    %sigma_flux = sigma_flux_array(sf_i);
    
    % Create C_m covariance matrix
    C_m = zeros(bins);
    for i = 1:bins
        for j = 1:bins
            C_m(i,j) = sigma_flux^2 * exp(-((energy_midpoints(i)-energy_midpoints(j))^2)./(2*delta.^2));
        end
    end
    inv_C_m = inv(C_m);
    
    mat_mult = inv(G'* inv_C_d * G + inv_C_m) * G' * inv_C_d;

    
% Defining constants over each iteration
    d_obs = log(hits_whole_EC');

% Defining initial values 
    mn = zeros(bins,it_max);
    %mn(:,1) = ones(1,bins).*log(10^3);
    % starting from best fit line of numeric lsqr
    %mn(:,1) = log(fit_result.a * exp(fit_result.b * energy_midpoints));
    mn(:,1) = log(fit_result.a * exp(fit_result.b * energy_midpoints)...
        +fit_result.c * exp(fit_result.d * energy_midpoints));

    
    x_edges = log(energy_edges);
    x = log(energy_midpoints');
    dx = x_edges(2:end)-x_edges(1:end-1);
    dx(dx>100) = 100;
    
    integ = geo_EC .* dx .* exp(mn(:,1) + x)';
    G_n = 1./hits_whole_EC' .* integ;
    
    g_mn = zeros(size(energy_channels,1),it_max);
    g_mn(:,1) = log(sum(integ,2));
      
    error = zeros(it_max,1);
    error_max = zeros(it_max,1);
    error_temp = (d_obs-g_mn(:,1)).^2 ./ sigma_d;
    error_max(1) = max(error_temp);
    error(1) = sum((d_obs-g_mn(:,1)).^2 ./ sigma_d);

    iteration = 2;
    error_cond = false;
% Begin iterations
while error_cond == false && iteration <= it_max
    mn(:,iteration) = mn(:,iteration-1) + mat_mult * (d_obs - g_mn(:,iteration-1) + G_n*(mn(:,iteration-1)-mn(:,1)));
    
    integ = geo_EC .* dx .* exp(mn(:,iteration) + x)';
    G_n = 1./hits_whole_EC' .* integ;
    g_mn(:,iteration) = log(sum(integ,2));
    error_temp = (d_obs-g_mn(:,iteration)).^2 ./ sigma_d;
    error_max(iteration) = max(error_temp);
    error(iteration) = sum((d_obs-g_mn(:,iteration)).^2 ./ sigma_d);
    %{
    if error(iteration) > error(iteration-1)
        error_cond = true;
    else
    %}
        iteration = iteration + 1;
    %end
end

[~,min_error_ind] = min(error);
fprintf("Iteration Number: %.0d \n",min_error_ind)
fprintf("Tolerance: %.6e \n",error(min_error_ind))
fprintf("Max Error: %.6e \n",error_max(min_error_ind))

flux_lsqr = exp(mn(:,min_error_ind));
%{
chi_sq_array(sf_i) = error(min_error_ind);

dx1 = (dx(2:end)-dx(1:end-1));
flux_deriv = diff(flux_lsqr) ./ dx1';

% Calculate second derivative
dx2 = (dx1(2:end)-dx1(1:end-1));
flux_d2 = diff(flux_deriv) ./ dx2';

% Calculate norm (sum of squared second derivatives)
norm_array(sf_i) = sum(flux_d2.^2);
end
%}

%% Plots calculated flux
f = figure;
f.Position = [0 0 1700 900];
hold on

% Plot simulated flux
plot(energy_midpoints,M_energy_bin/(4*pi^2*r_source^2)./bin_width, 'Color', 'black','LineWidth',2.5);
fit_flux_result = fit(energy_midpoints',(M_energy_bin/(4*pi^2*r_source^2)./bin_width)','exp2');
fitted_flux_curve = @(x) feval(fit_flux_result, x);
plot(energy_midpoints, fitted_flux_curve(energy_midpoints),...
    'LineStyle', '--','LineWidth',2,'Color', 'black');

%
% Plot Bowtie points
%plot(E_eff,j_nom,'o', 'Color', '#0072BD','MarkerSize',10);

% Plot LSQR Selesnick Method
plot(energy_midpoints(indices),flux_lsqr(indices),'.', 'Color', 'r', 'MarkerSize',12);

% Plot LSQR Summation Method
plot(energy_midpoints(indices),flux_approx(indices),'x','Color','#77AC30','MarkerSize',10)
% and fitted curve
fitted_curve = @(x) feval(fit_result, x);
plot(energy_midpoints, fitted_curve(energy_midpoints),...
    'LineStyle', '--','LineWidth',2,'Color', '#77AC30');

%{
legend({['Incident Particle Measurements'],['Bowtie Method'], ...
                 ['Least Squares Method (Selesnick)'], ['Least Squares Method (Simple)']},...
                 'Location', 'northeast','FontSize',18);
%}
%
legend({['Incident Particle Measurements'],['Least Squares Method (Selesnick)'],['Least Squares Method (Simple)']},...
                 'Location', 'northeast','FontSize',18);
%

set(gca, 'FontSize', textsize)
xlim([0 8])
xticks((0:1:8))
set(gca, 'YScale', 'log')
ylim([10^3, 10^4])

%ylabel('Flux')
ylabel('I  #/(cm^2 sr s MeV)','FontSize',textsize)
xlabel('Energy (MeV)','FontSize',textsize)
hold off

%{
figure
% Create the log-log plot
loglog(chi_sq_array, norm_array,'.','MarkerSize',12);

% Label the axes
xlabel('Data Misfit (chi-squared)');
ylabel('Model Norm');

% Add title (optional)
title('L-curve analysis for flux spectrum estimation');

figure
plot(sigma_flux_array,chi_sq_array,'.','MarkerSize',12);
xlabel('Sigma');
ylabel('Data Misfit (chi-squared)');
set(gca, 'YScale', 'log')

figure
plot(sigma_flux_array,norm_array,'.','MarkerSize',12);
xlabel('Sigma')
ylabel('Model Norm');
set(gca, 'YScale', 'log')
%}
