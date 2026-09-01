% Least Squares and Count Rate Analysis for HERT
% Last modified: 5/5/2026
% Added Weibull SEP Spectra
close all
cd 'C:\Users\wzt0020\Box\HERT_Box\Inverse Problem\FS_figures'
%warning('off','last')

%% Inputs
% Initialize Variables
r_source = 8.5;
parttype = 0; % 0 = electron, 1 = proton
geo_EC = readmatrix('D:\HERT_Drive\Matlab Main\Result\Electron_FS\electron_FS_31000_v2_GFbyEC.txt');
geo_EC = geo_EC(:,2:end)';
%geo_EC = readmatrix('C:\Users\William Teague\Box\HERT_Box\Matlab Main\Result\geometric_factor_EC.txt');
% Finds energy edges and midpoints based off of bins
if parttype == 0
    bins = 400;
    energy_min = 0.01;
    energy_max = 10;
    energy_edges = logspace(log10(energy_min),log10(energy_max),bins+1);
    energy_edges(end) = energy_max + 0.01;
elseif parttype == 1
    bins = 400;
    energy_min = 10;
    energy_max = 1000;
    energy_edges = logspace(log10(energy_min),log10(energy_max),bins+1);
end

energy_midpoints = (energy_edges(2:end) + energy_edges(1:end-1))/2;
dt = 1;

energy_channels = readmatrix('D:\HERT_Drive\Matlab Main\Result\channel_select\electron_channels_v2.txt');
bin_width = diff(energy_edges);
%energy_channels = readmatrix('C:\Users\William Teague\Box\HERT_Box\Matlab Main\Result\channel_select\electron_channels_v1.txt');

%% Creating Test Fluxes %%
% --- Select Flux Model Here ---
% Options: 'Main1', 'Linear', 'Exponential', 'BOT1', 'BOT2', 'Power Law', 'Gaussian'
flux_model = 'Power Law'; 
rng(42);
%rng('shuffle');

% --- Flux Calculation Tree ---
switch flux_model  
    case 'Main1'
        % Calculating Flux from Main1
        flux = M_energy_bin ./ (4*pi^2*r_source^2) ./ bin_width_red;
        sigma_m = 1500;
        delta = 50;      
    case 'Linear'
        flux = ones(1, bins) * 10^3;
        sigma_m = 2000;
        delta = 50;      
    case 'Exponential'
        % Vary coefficient from 0.2 to 2
        coeff = 1; 
        flux = 10^5 * exp(-(energy_midpoints) / coeff);
        sigma_m = 2000;
        delta = 50;    
    case 'BOT1'
        % BOT/Inverse Option 1
        flux = (1/0.01) * exp(log(energy_midpoints.^-0.69)) + ...
               (1/0.001) .* exp(-(log(energy_midpoints) - log(2.365)).^2 ./ (2*0.14));
        sigma_m = 2000;
        delta = 6;         
    case 'BOT2'
        % BOT/Inverse Option 2
        flux = (1/0.001) * exp(log(energy_midpoints.^-1.2)) + ...
               (1/0.001) .* exp(-(log(energy_midpoints) - log(4)).^2 ./ (2*0.08));
        sigma_m = 2000;
        delta = 6;     
    case 'Power Law'
        % Power Law: alpha can be between 2 and 6
        % Note: use j0=10^2 and alpha=3.14 or j0=10^1 and alpha=4.73 from Hong's paper
        j0 = 10^3;
        alpha = 4;
        flux = j0 .* energy_midpoints.^-alpha;
        sigma_m = 2000;
        delta = 40;        
    case 'Gaussian'
        flux = (1/0.000001) .* exp(-(log(energy_midpoints) - log(2)).^2 ./ (2*0.004));
        sigma_m = 2000;
        delta = 50;    
    otherwise
        error('Invalid flux model selected. Check your spelling in the flux_model variable.');
end

hits_whole_EC = sum(geo_EC * (flux' .* bin_width'),2)*dt; 
% add poisson noise to counts
poisson_err = sqrt(hits_whole_EC);
random_scatter = randn(size(hits_whole_EC));
hits_whole_EC_mod = hits_whole_EC + (poisson_err .* random_scatter);
idx_valid = hits_whole_EC_mod>=1;
hits_whole_EC_red = hits_whole_EC_mod(idx_valid);
hits_whole_EC_red = round(hits_whole_EC_red);

% One Count %
%hits_whole_EC= ones(size(geo_EC,1),1);
%

%% Generate Weibull Flux Spectra from Python Export
%{
% 1. Read the exported parameters from Python
param_file = 'D:\HERT_Drive\Matlab Main\Result\Proton_FS\Weibull_Parameters_Export.txt';
weibull_data = readtable(param_file);

% 2. Extract the parameters for 25%, 50%, and 75% dynamically
% Find the row index where the CI column matches the desired confidence interval
idx_25 = find(weibull_data.CI == 25);
k_25   = weibull_data.k(idx_25);
E0_25  = weibull_data.E_0(idx_25);
b_25   = weibull_data.b(idx_25);

idx_50 = find(weibull_data.CI == 50);
k_50   = weibull_data.k(idx_50);
E0_50  = weibull_data.E_0(idx_50);
b_50   = weibull_data.b(idx_50);

idx_75 = find(weibull_data.CI == 75);
k_75   = weibull_data.k(idx_75);
E0_75  = weibull_data.E_0(idx_75);
b_75   = weibull_data.b(idx_75);

% 3. Calculate the Differential Flux Spectra
% Equation: Phi(E) = k * E^(b-1) * exp(-(E/E_0)^b)
flux_25 = k_25 .* (E.^(b_25 - 1)) .* exp(-(E ./ E0_25).^b_25);
flux_50 = k_50 .* (E.^(b_50 - 1)) .* exp(-(E ./ E0_50).^b_50);
flux_75 = k_75 .* (E.^(b_75 - 1)) .* exp(-(E ./ E0_75).^b_75);

% 4. Plot to verify
%{
figure('Position', [100, 100, 800, 600]);
loglog(E, flux_25, 'b-', 'LineWidth', 2, 'DisplayName', '25% CI'); hold on;
loglog(E, flux_50, 'g-', 'LineWidth', 2, 'DisplayName', '50% CI');
loglog(E, flux_75, 'r-', 'LineWidth', 2, 'DisplayName', '75% CI');

grid on;
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on');
xlabel('Energy (MeV)', 'FontSize', 14);
ylabel('Average Flux (cm^{-2} sr^{-1} s^{-1} MeV^{-1})', 'FontSize', 14);
title('Weibull Flux Spectra by Confidence Interval', 'FontSize', 16);
legend('Location', 'southwest', 'FontSize', 12);
%}

% Find hits in each energy channel from simulated flux %
flux = flux_75;

% add noise to flux spectrum
noise_level = 0.15;
random_scatter = randn(size(flux));
flux_noisy = flux .* exp(noise_level .* random_scatter);
flux = flux_noisy;

dt = 600;
hits_whole_EC= sum(geo_EC * (flux' .* bin_width'),2)*dt; 
%}

%% Reducing equations to remove zero hit counts
energy_channels_red = energy_channels(hits_whole_EC_red>=1,:);
energy_edges = energy_edges(energy_edges<energy_channels_red(end,2)+1);
bounds = energy_midpoints<energy_edges(end);
energy_midpoints_red = energy_midpoints(bounds);
geo_EC_red = geo_EC(hits_whole_EC_red>=1,bounds);
bin_width_red = bin_width(bounds);
flux = flux(bounds);
energy_max = max(energy_midpoints_red);
%bound_plot = energy_midpoints_red >= 0.5 & energy_midpoints_red<=7;
%indices = find(bound_plot);

%% Simple Multiple Linear Regression Method %%
%
A = geo_EC_red .* bin_width_red;            % Calculate known/"independent variable"
inv_A = pinv(A);                    % take pseudo inverse (not square matrix)
flux_lin = inv_A * hits_whole_EC_red;  % find flux from linear algebra

% Plot
%{
f = figure;
f.Position = [0 0 1200 900];
hold on

% bound to HERT's acceptable energy range
plot(energy_midpoints_red,flux, 'Color', 'black','LineWidth',4);
plot(energy_midpoints_red,flux_lin,'.', 'Color', 'r','MarkerSize',12);
xline(0.6, '--','color', [.5 .5 .5],'LineWidth',2);

% Plot Bowtie points
plot(E_eff,j_nom,'o', 'Color', '#0072BD','MarkerSize',10,'LineWidth',2);

legend({['Acutal Flux'],['Basic Linear Regression'],['Low Energy Threshold'],['Bowtie']},...
                 'Location', 'northeast','FontSize',18);

textsize = 24;
set(gca, 'FontSize', textsize)
xlim([0 7])
%ylim([10^0 10^5])
xticks((0:1:8))
%set(gca, 'XScale', 'log')
%set(gca, 'YScale', 'log')

ylabel('Flux  (# cm^{-2} sr^{-1} s^{-1} MeV^{-1})','FontSize',textsize)
xlabel('Energy (MeV)','FontSize',textsize)
hold off
%}

%% Least Squares Function for Energy Channels (Selesnick/Khoo) %%
if length(hits_whole_EC_red)>2
    % Initialize variables
    flux_lsqr = 0;
    jsig = 0;
    it_max = 100;   % maximum number of possible iterations 

    % FIX 1: Convert simulated total counts back to an observed count rate
    rate_obs = hits_whole_EC_red / dt; 

    % Setting up constant matrices
    Cd = zeros(size(energy_channels_red,1));
    sigma_d = zeros(size(energy_channels_red,1),1);
    for channel = 1:size(energy_channels_red,1)
        sigma_d(channel) = sqrt(rate_obs(channel)*dt+1)/dt; 
        Cd(channel,channel) = (sigma_d(channel)/(rate_obs(channel))).^2;
    end
    inv_Cd = inv(Cd); 

    Cm = sigma_m.^2 .* exp(-((energy_midpoints_red' - energy_midpoints_red).^2) ./ (2 * delta.^2));
    
    % Adds a tiny mathematical floor to the diagonal to prevent matrix singularity
    Cm = Cm + eye(length(energy_midpoints_red)) * 1e-4;
    
    inv_Cm = pinv(Cm);
    
    d_obs = log(rate_obs);

    % Find initial exponential guess via pre-fitting
    mn = zeros(it_max, length(energy_midpoints_red));

    cost_function = @(x) sum( (log(rate_obs) - log( sum(geo_EC_red .* bin_width_red .* (exp(x(1)) .* exp(-energy_midpoints_red ./ x(2))), 2) )) .^ 2 );
    
    % Initial guess for the optimizer [ln(j0), E0]
    % Based on your documentation, j0 ~ 10^6 and E0 ranges from 0.2 to 2.0 MeV
    guess_params = [log(1e6), 1.0];
    
    % Run the Nelder-Mead simplex optimizer silently
    options = optimset('Display', 'off');
    best_fit = fminsearch(cost_function, guess_params, options);
    
    % Extract the mathematically optimal parameters
    j0_opt = exp(best_fit(1));
    E0_opt = best_fit(2);
    mn(1,:) = log(j0_opt) - (energy_midpoints_red ./ E0_opt);
    %mn(1,:) = real(log(k_50) + (b_50 - 1).*log(energy_midpoints_red) - (energy_midpoints_red ./ E0_50).^b_50);

    x_edges = log(energy_edges);
    x = log(energy_midpoints_red);
    dx = x_edges(2:end)-x_edges(1:end-1);
    dx(dx>100) = 100;
    
    % initial count rate guess
    g_mn = zeros(it_max,size(energy_channels_red,1));
    g_mn(1,:) = log(sum(geo_EC_red .* dx .* exp(mn(1,:)+x),2)); 
    
    Gn = 1./exp(g_mn(1,:))' .* geo_EC_red .* dx .* exp(mn(1,:)+x);
   
    iteration = 2;
    convergence = false;
    Cmm = zeros(length(energy_midpoints_red));

%% Begin iterations %%
while convergence == false && iteration <= it_max
    Hessian = Gn'* inv_Cd' * Gn + inv_Cm';
    Gradient = Gn' * inv_Cd * (d_obs - g_mn(iteration-1,:)' + Gn * (mn(iteration-1,:)' - mn(1,:)'));
    
    mn(iteration,:) = (mn(1,:)' + pinv(Hessian) * Gradient)';
    
    g_mn(iteration,:) = log(sum(geo_EC_red .* dx .* exp(mn(iteration,:)+x),2));
    Gn = 1./exp(g_mn(iteration,:))' .* geo_EC_red .* dx .* exp(mn(iteration,:)+x);
    
    Cmm = pinv(Hessian);
    
    if max((mn(iteration,:)- mn(iteration-1,:)).^2) < 0.01
        convergence = true;
        disp("Converges")
    else
        iteration = iteration + 1;
    end
end

if convergence == true
    flux_lsqr = exp(mn(iteration,:));
    jsig = flux_lsqr.*sqrt(diag(Cmm))';
    fprintf("Iteration Number: %.0d \n",iteration)
end

else
    error('More than 2 energy channels must have counts!')
end

%% Plots calculated flux
% 0. Setup Publication Typography & Physical Dimensions
textsize = 14;              
fig_width_in  = 7.0;        
fig_height_in = 2.8;        

f = figure('Units', 'inches', ...
           'Position', [1, 1, fig_width_in, fig_height_in], ...
           'PaperUnits', 'inches', ...
           'PaperPosition', [0, 0, fig_width_in, fig_height_in], ...
           'PaperPositionMode', 'auto'); 
ax = axes(f);
hold(ax, 'on');

% 1. Plot simulated flux
p1 = plot(ax, energy_midpoints_red, flux, 'LineWidth', 2, 'Color', 'black', 'DisplayName', 'Theoretical Flux');

% 2. Plot Bowtie points and error
ypos = j_nom_err;
yneg = j_nom_err;
lower_bound = j_nom - yneg;
below_zero_idx = lower_bound < 10^-4;
yneg(below_zero_idx) = j_nom(below_zero_idx) - 10^-4;
idx_valid(end) = false;
p2 = errorbar(ax, E_eff(idx_valid), j_nom(idx_valid), yneg(idx_valid), ypos(idx_valid), 'o', 'Color', '#0072BD', ...
    'MarkerSize', 5, 'MarkerFaceColor', '#0072BD', 'LineWidth', 1.5, 'CapSize', 10, 'DisplayName', 'Bowtie Analysis');

% 3. Plot LSQR Selesnick Method Fit
p3 = plot(ax, energy_midpoints_red, flux_lsqr, 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'GNLSM Fit');

% 4. Plot standard deviation from fit (Hidden from legend)
plot(ax, energy_midpoints_red, flux_lsqr + jsig, 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(ax, energy_midpoints_red, flux_lsqr - jsig, 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Apply Legend
legend(ax, [p1, p2, p3], 'Location', 'southwest', 'FontSize', textsize - 2);

% --- SHARED FORMATTING ---
set(ax, 'FontSize', textsize, 'Layer', 'top');
set(ax, 'XScale', 'log', 'YScale', 'log');

set(ax, 'XTick', [0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
grid(ax, 'on'); 
ax.YMinorGrid = 'off';

ylabel(ax, 'Flux (cm^{-2} sr^{-1} s^{-1} MeV^{-1})', 'FontSize', textsize);
xlabel(ax, 'Energy (MeV)', 'FontSize', textsize);

% Dynamic Axis Limits
switch flux_model
    case 'BOT2'
        xlim(ax, [0.5, 10]); 
        ylim(ax, [40, max(flux_lsqr(energy_midpoints_red > 0.5) + jsig(energy_midpoints_red > 0.5)) * 2]);
    case 'Power Law'
        xlim(ax, [0.5, min(max(energy_midpoints_red), 10)]); 
        % 1. Calculate the highest power of 10 needed for the upper limit
        max_flux = max(flux_lsqr(energy_midpoints_red > 0.5) + jsig(energy_midpoints_red > 0.5));
        max_power = ceil(log10(max_flux));
        % 2. Set the Y-axis limits
        ylim(ax, [10^0, 10^max_power]);
        % 3. Force the Y-ticks to only appear on even powers (10^0, 10^2, 10^4...)
        set(ax, 'YTick', 10.^(0:2:max_power));
    otherwise
        xlim(ax, [0.5, 10]); 
        ylim(ax, [10^0, max(flux_lsqr(energy_midpoints_red > 0.5) + jsig(energy_midpoints_red > 0.5)) * 2]);
end

hold(ax, 'off');

cd("C:\Users\wzt0020\Box\HERT_Box\Paper\Figures\GNLSM examples");

% Determine File Name based on Selected Model
switch flux_model
    case 'Main1'
        base_filename = 'main1_noise_50_Emax';      
    case 'Linear'
        base_filename = 'linear_103_noise_50_Emax';       
    case 'Exponential'
        base_filename = 'exp_j0_105_E0_1_noise_50_Emax'; 
        title('Exponential Electron Flux Spectrum')
    case 'BOT1'
        base_filename = 'bot1_noise_50_Emax';
        title('Bump-on-Tail Electron Flux Spectrum')
    case 'BOT2'
        base_filename = 'bot2_noise_50_Emax';
        title('Bump-on-Tail Electron Flux Spectrum')
    case 'Power Law'
        base_filename = 'pow_j0_103_alpha_4_noise_50_Emax_reduced';
        title('Power Law Electron Flux Spectrum')
    case 'Gaussian'
        base_filename = 'gauss_noise_50_Emax';        
    otherwise
        base_filename = 'custom_flux_model_export';
end

% 1. Save the interactive MATLAB figure file (.fig)
savefig(f, base_filename);

% 2. Export High-Resolution PNG and Vector PDF
print(f, [base_filename, '.png'], '-dpng', '-r300');
exportgraphics(f, [base_filename, '.pdf'], 'ContentType', 'vector');

%% Chi Squared Goodness of Fit Test
%{
null = chi2inv(0.05,size(energy_channels,1)-2); %2 DOF stripped from sigma and delta?
chi_squared = sum(((d_obs - g_mn(iteration,:)')./sigma_d).^2);
chi_squared_bowtie = sum(((d_obs - log(j_nom))./sigma_d).^2);

if chi_squared < null
    fprintf("Chi Squared %.6f < Null %.6f \n",chi_squared,null)
    fprintf("Do Not Reject Null Hypotheses (fit is good) \n")
else
    fprintf("Chi Squared %.6f > Null %.6f \n",chi_squared,null)
    fprintf("Reject Null Hypotheses (fit is poor) \n")
end
if chi_squared_bowtie < null
    fprintf("Bowtie Chi Squared %.6f < Null %.6f \n",chi_squared_bowtie,null)
    fprintf("Do Not Reject Null Hypotheses (fit is good) \n")
else
    fprintf("Bowtie Chi Squared %.6f > Null %.6f \n",chi_squared_bowtie,null)
    fprintf("Reject Null Hypotheses (fit is poor) \n")
end
%}