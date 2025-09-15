% Least Squares and Count Rate Analysis for HERT
% Last modified: 5/29/2024
close all
%warning('off','last')

%% Inputs
% Initialize Variables
r_source = 8.5;
geo_EC = readmatrix('E:\HERT_Drive\Matlab Main\Result\geofactor_EC_FS.txt');
%geo_EC = readmatrix('C:\Users\William Teague\Box\HERT_Box\Matlab Main\Result\geometric_factor_EC.txt');
bins = size(geo_EC,2);
energy_edges = logspace(log10(0.01),log10(8),bins+1);
bin_width = diff(energy_edges);
energy_midpoints = (energy_edges(2:end) + energy_edges(1:end-1))/2;

energy_channels = readmatrix('E:\HERT_Drive\Matlab Main\Result\channel_select\electron_channels_v1.txt');
%energy_channels = readmatrix('C:\Users\William Teague\Box\HERT_Box\Matlab Main\Result\channel_select\electron_channels_v1.txt');

%% Creating Test Fluxes %%
% Calculating Flux from Main1 %
%flux = M_energy_bin/(4*pi^2*r_source^2)./bin_width;

% Linear %
%flux = ones(1,bins)*10^3;

% Exponential % Vary coefficient from 0.2 to 2
%flux = 10^6 * exp(-(energy_midpoints)/1);

% BOT/Inverse %
%flux = 1/0.01 * exp(log(energy_midpoints.^-0.69))+ 1/0.001 .* exp(-(log(energy_midpoints)-log(2.365)).^2./(2*0.14));
%flux = 1/0.01 * exp(log(energy_midpoints.^-1.2))+ 1/0.001 .* exp(-(log(energy_midpoints)-log(4)).^2./(2*0.08));

% Power Law % alpha can be between 2 and 6
flux = 10^4 .* energy_midpoints.^-5;

% Gaussian %
%flux = 1/0.000001 .* exp(-(log(energy_midpoints)-log(2)).^2./(2*0.004));

% Find hits in each energy channel from simulated flux %
hits_whole_EC= sum(geo_EC * (flux' .* bin_width'),2); 

%% Reducing equations to remove zero hit counts
energy_channels = energy_channels(hits_whole_EC~=0,:);
energy_edges = energy_edges(energy_edges<energy_channels(end,2)+1);
bounds = energy_midpoints<energy_edges(end);
energy_midpoints = energy_midpoints(bounds);
geo_EC = geo_EC(hits_whole_EC~=0,bounds);
hits_whole_EC = hits_whole_EC(hits_whole_EC~=0);
bin_width = bin_width(bounds);
flux = flux(bounds);
bound_plot = energy_midpoints >= 0.5 & energy_midpoints<=7;
indices = find(bound_plot);

dt = 10;

%% Simple Multiple Linear Regression Method %%
%
A = geo_EC .* bin_width;            % Calculate known/"independent variable"
inv_A = pinv(A);                    % take pseudo inverse (not square matrix)
flux_lin = inv_A * hits_whole_EC;  % find flux from linear algebra

% Plot
%
f = figure;
f.Position = [0 0 1200 900];
hold on

% bound to HERT's acceptable energy range
plot(energy_midpoints,flux, 'Color', 'black','LineWidth',4);
plot(energy_midpoints,flux_lin,'.', 'Color', 'r','MarkerSize',12);
xline(0.6, '--','color', [.5 .5 .5],'LineWidth',2);

% Plot Bowtie points
%plot(E_eff,j_nom,'o', 'Color', '#0072BD','MarkerSize',10,'LineWidth',2);

%legend({['Acutal Flux'],['Basic Linear Regression'],['Low Energy Threshold'],['Bowtie']},...
%                 'Location', 'northeast','FontSize',18);

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
%

%% Least Squares Function for Energy Channels (Selesnick/Khoo) %%
% Initialize variables
    flux_lsqr = 0;
    jsig = 0;
    it_max = 100;   % maximum number of possible iterations 

% Setting up constant matrices
    % Cd is the covariance matrix of the data
    Cd = zeros(size(energy_channels,1));
    sigma_d = zeros(size(energy_channels,1),1);
    for channel = 1:size(energy_channels,1)
        sigma_d(channel) = sqrt(hits_whole_EC(channel)*dt+1)/dt; % this is the data std dev
        Cd(channel,channel) = (sigma_d(channel)/(hits_whole_EC(channel))).^2;
    end
    inv_Cd = inv(Cd); % finding the inverse for later use

    % Initialize variance parameter %
    %sigma_m = 1600; % Exp = 16000   BOT = 700,   POW = 270 100?

    % Initialize smoothness parameter
    %delta = 100; % Exp = 100   BOT = 2,   POW = 27 60?

    % Define the range of parameters to search
    % These ranges should be wide enough to find the optimal values.
    sigma_m_range = [100, 500, 1000, 2000, 5000, 10000];
    delta_range = [10, 50, 100, 200, 500, 1000];
    
    % Initialize a matrix to store the results
    results = zeros(length(sigma_m_range), length(delta_range)); % 3 for chi2, conv_status, flux_val

    % Main loop for parameter search
    for j = 1:length(sigma_m_range)
        sigma_m = sigma_m_range(j);
        for i = 1:length(delta_range)
            delta = delta_range(i);

            % Create C_m covariance matrix, covariance of model/guess
            Cm = sigma_m.^2 .* exp(-((energy_midpoints' - energy_midpoints).^2) ./ (2 * delta.^2));
            inv_Cm = inv(Cm);
        
            % Defining constants over each iteration
            d_obs = log(hits_whole_EC);
        
            % Defining initial values 
            % Create initial estimate
            mn = zeros(it_max,length(energy_midpoints));
        
            % finding x and dx for integrals from existing edges
            x_edges = log(energy_edges);
            x = log(energy_midpoints);
            dx = x_edges(2:end)-x_edges(1:end-1);
            dx(dx>100) = 100;
            
            % Set up first iteration of model
            
            % initial count rate guess based on initial flux estimate
            g_mn = zeros(it_max,size(energy_channels,1));
            g_mn(1,:) = log(sum(geo_EC .* dx .* exp(mn(1,:)+x),2)); 
            Gn = 1./exp(g_mn(1,:))' .* geo_EC .* exp(mn(1,:)+x);
           
            % initialize loop variables
            iteration = 2;
            convergence = false;
            Cmm = zeros(length(energy_midpoints));

%% Begin iterations %%
            while convergence == false && iteration <= it_max
                
                % find the log(flux) for this iteration
                %mn(iteration,:) = (Gn'* inv_Cd' * Gn + inv_Cm') \ (Gn'* inv_Cd' * (d_obs - g_mn(iteration-1,:)' - Gn * mn(iteration-1,:)') + inv_Cm' * mn(iteration-1,:)');
                 mn(iteration,:) = mn(1,:)' + (Gn'* inv_Cd' * Gn + inv_Cm') \ Gn' * inv_Cd * (d_obs - g_mn(iteration-1,:)' + Gn * (mn(iteration-1,:)' - mn(1,:)'));
                
                % Calculate the new matrices for the iteration
                g_mn(iteration,:) = log(sum(geo_EC .* dx .* exp(mn(iteration,:)+x),2));
                Gn = 1./exp(g_mn(iteration,:))' .* geo_EC .* dx .* exp(mn(iteration,:)+x);
            
                % Calculate the new model covariance matrix
                Cmm = inv(Gn'* inv_Cd' * Gn + inv_Cm');
            
                % determine if the flux converges
                if max((mn(iteration,:)- mn(iteration-1,:)).^2)<0.01
                    convergence = true;
                    disp("Converges")
                else
                    iteration = iteration + 1;
                end
            end

            % If convergence is true, find the flux and the error from the last iteration
            if convergence == true
                flux_lsqr = exp(mn(iteration,:));
                jsig = flux_lsqr.*sqrt(diag(Cmm))';
                fprintf("Iteration Number: %.0d \n",iteration)
                %ind_act_error_avg(d,sf_i) = actual_error_avg(iteration);
                results(i, j) = sum(((d_obs - g_mn(iteration,:)')./sigma_d).^2);
            end
        end
    end

%% Plots calculated flux
f = figure;
f.Position = [0 0 1200 900];
hold on

% Plot simulated flux
plot(energy_midpoints(bound_plot),flux(bound_plot), 'Color', 'black','LineWidth',4);

%plot(energy_midpoints(bounds),M_energy_bin(bounds)./(4*pi^2*r_source^2)./bin_width(bounds),'.', 'Color', 'black','MarkerSize',8);

% Plot Bowtie points
%plot(E_eff,j_nom,'o', 'Color', '#0072BD','MarkerSize',10,'LineWidth',2);

% Plot LSQR Selesnick Method
% Plot calculated fit
plot(energy_midpoints(bound_plot),flux_lsqr(bound_plot), 'Color', 'r', 'LineWidth',2);
% plot standard deviation from fit
plot(energy_midpoints(bound_plot),flux_lsqr(bound_plot)+jsig(bound_plot),'r--','LineWidth',2);
plot(energy_midpoints(bound_plot),flux_lsqr(bound_plot)-jsig(bound_plot),'r--','LineWidth',2);

legend({['Acutal Flux'],['LSQR'],['Standard Deviation']},...
                 'Location', 'northeast','FontSize',18);

%legend({['Theoretical Flux'],['Incident Particle Measurement'],['Bowtie Analysis'],['LSQR'],['Standard Deviation']},...
 %                'Location', 'northeast','FontSize',18);

textsize = 24;
set(gca, 'FontSize', textsize)
xlim([0 7])
ylim([10^0 10^5])
xticks((0:1:8))
%set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
ylabel('Flux  (# cm^{-2} sr^{-1} s^{-1} MeV^{-1})','FontSize',textsize)
xlabel('Energy (MeV)','FontSize',textsize)
hold off

%% Chi Squared Goodness of Fit Test
null = chi2inv(0.05,size(energy_channels,1)-2); %2 DOF stripped from sigma and delta?
chi_squared = sum(((d_obs - g_mn(iteration,:)')./sigma_d).^2);

if chi_squared < null
    fprintf("Chi Squared %.6f < Null %.6f \n",chi_squared,null)
    fprintf("Do Not Reject Null Hypotheses (fit is good) \n")
else
    fprintf("Chi Squared %.6f > Null %.6f \n",chi_squared,null)
    fprintf("Reject Null Hypotheses (fit is poor) \n")
end