% Least Squares and Count Rate Analysis for HERT
% Last modified: 5/29/2024

%warning('off','last')

%% Inputs
% Initialize Variables
r_source = 8.5;
geo_EC = readmatrix('E:\HERT_Drive\Matlab Main\Result\geometric_factor_EC.txt');
bins = size(geo_EC,2);
energy_edges = linspace(0,8,bins+1);
bin_width = energy_edges(2:end)-energy_edges(1:end-1);
energy_midpoints = (energy_edges(2:end) + energy_edges(1:end-1))/2;

energy_channels = readmatrix('E:\HERT_Drive\Matlab Main\Result\channel_select\electron_channels_v1.txt');

% Creating Test Fluxes
%flux = ones(1,bins)*10^3;
%flux = 10^6 * exp(-(energy_midpoints)/1) ./ (4*pi^2*r_source^2) ./bin_width;
flux = 1/0.01 * exp(log(energy_midpoints.^-0.69))+ 1/0.001 .* exp(-(log(energy_midpoints)-log(2.365)).^2./(2*0.14));
%flux = 1/0.01 * exp(log(energy_midpoints.^-1.2))+ 1/0.001 .* exp(-(log(energy_midpoints)-log(4)).^2./(2*0.08));
%flux = 10^5 .* energy_midpoints.^-4.73;
%flux = 10^5 .* energy_midpoints.^-3.14;

%
hits_whole_EC = zeros(1,size(geo_EC,1));
for channel = 1:size(energy_channels,1)
    hits_whole_EC(channel) = sum(geo_EC(channel,:) .* flux .* bin_width); 
end
%

% Reducing equations to remove zero hit counts
energy_channels = energy_channels(hits_whole_EC~=0,:);
energy_edges = energy_edges(energy_edges<energy_channels(end,2)+1);
bounds = energy_midpoints>0.6 & energy_midpoints<7 & energy_midpoints<energy_edges(end);
energy_midpoints = energy_midpoints(bounds);
energy_edges = energy_edges(energy_edges>=0.6 & energy_edges<=7);
geo_EC = geo_EC(hits_whole_EC~=0,bounds);
hits_whole_EC = hits_whole_EC(hits_whole_EC~=0);
%flux = M_energy_bin/(4*pi^2*r_source^2)./bin_width;
bin_width = bin_width(bounds);
flux = flux(bounds);

dt = 10;

%% Least Squares Function for Energy Channels (Selesnick/Khoo)
% Initialize variables
    flux_lsqr = 0;
    jsig = 0;
    it_max = 25;   % maximum number of possible iterations 

% Setting up constant matrices
    % Cd is the covariance matrix of the data
    Cd = zeros(size(energy_channels,1));
    sigma_d = zeros(size(energy_channels,1),1);
    for channel = 1:size(energy_channels,1)
        sigma_d(channel) = sqrt(hits_whole_EC(channel)*dt+1)/dt; % this is the data std dev
        Cd(channel,channel) = (sigma_d(channel)/(hits_whole_EC(channel))).^2;
    end
    inv_Cd = inv(Cd); % finding the inverse for later use

    % Initialize variance parameter
    sigma_m = 10^7;
    %sigma_m_array = logspace(0,8,100);

    % Initialize smoothness parameter
    delta = 4; 
    %delta_array = logspace(-2,2,40);


% Loop over variance and smoothness parameters
%{
ind_act_error_avg = zeros(length(delta_array),length(sigma_m_array)); 

%
for d = 1:length(delta_array)
    fprintf("Delta Iteration: %.0d \n",d);
    delta = delta_array(d);
%
for sf_i = 1:length(sigma_m_array)
    % Current sigma_m value
    sigma_m = sigma_m_array(sf_i);
%}
    % Create C_m covariance matrix, covariance of model/guess
    Cm = zeros(length(energy_midpoints));
    Cmm = zeros(length(energy_midpoints));
    for bin = 1:length(energy_midpoints)
        Cm(:,bin) = sigma_m.^2 .* exp(-((energy_midpoints-energy_midpoints(bin)).^2)./(2*delta.^2));
    end
    inv_Cm = inv(Cm);
    
% Defining constants over each iteration
    d_obs = log(hits_whole_EC');

% Defining initial values 
    % Create initial estimate
    mn = zeros(it_max,length(energy_midpoints));
    mn(1,:) = ones(1,length(energy_midpoints));
    
    x_edges = log(energy_edges);
    x = log(energy_midpoints);
    dx = x_edges(2:end)-x_edges(1:end-1);
    dx(dx>100) = 100;
    
    % Set up first iteration of model
    integ = geo_EC .* dx .* exp(mn(1,:)+x);
    % initial count rate guess based on initial flux estimate
    g_mn = zeros(it_max,size(energy_channels,1));
    g_mn(1,:) = log(sum(integ,2)); 
    % matrix relation between energy bins and count rates
    Gn = 1./exp(g_mn(1,:)') .* integ;       
    mat_mult = inv(Gn'* inv_Cd * Gn + inv_Cm) * Gn' * inv_Cd;
   
    % initialize loop variables
    iteration = 2;
    convergence = false;

    % actual error initialization & calculation (REMOVE BEFORE ACTUAL USE)
    actual_error = zeros(it_max,length(energy_midpoints));
    actual_error_avg = zeros(it_max,1);
    actual_error_max = zeros(it_max,1);
    actual_error(1,:) = (exp(mn(1,:))-flux)./flux;
    actual_error_max(1) = max(abs(actual_error(1,:)));
    actual_error_avg(1) = mean(abs(actual_error(1,:)));

% Begin iterations
while convergence == false && iteration <= it_max
    
    mn(iteration,:) = mn(1,:)' + mat_mult * (d_obs - g_mn(iteration-1,:)' + Gn*(mn(iteration-1,:)-mn(1,:))');
    
    integ = geo_EC .* dx .* exp(mn(iteration,:)+x);
    g_mn(iteration,:) = log(sum(integ,2));
    Gn = 1./exp(g_mn(iteration,:))' .* integ;
    mat_mult = (Gn'* inv_Cd * Gn + inv_Cm) \ Gn' * inv_Cd;

    Cmm = inv(Gn'* inv_Cd * Gn + inv_Cm);

    % actual error calculation (REMOVE BEFORE ACTUAL USE)
    actual_error(iteration,:) = (exp(mn(iteration,:))-flux)./flux;
    actual_error_max(iteration) = max(abs(actual_error(iteration,actual_error(iteration,:)~=Inf)));
    actual_error_avg(iteration) = mean(abs(actual_error(iteration,actual_error(iteration,:)~=Inf)));

    if max(abs(mn(iteration,:)- mn(iteration-1,:)))<0.1
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
    %ind_act_error_avg(d,sf_i) = actual_error_avg(iteration);
end

%end
%end

%% Plots calculated flux
f = figure;
f.Position = [0 0 1700 900];
hold on

% Plot simulated flux
plot(energy_midpoints,flux, 'Color', 'black','LineWidth',4);

%plot(energy_midpoints,M_energy_bin(bounds)./(4*pi^2*r_source^2)./bin_width,'x', 'Color', 'black','MarkerSize',10);

% Plot Bowtie points
plot(E_eff,j_nom,'o', 'Color', '#0072BD','MarkerSize',10);

% Plot LSQR Selesnick Method
plot(energy_midpoints,flux_lsqr, 'Color', 'r', 'LineWidth',2);
plot(energy_midpoints,flux_lsqr+jsig,'r--','LineWidth',2);
plot(energy_midpoints,flux_lsqr-jsig,'r--','LineWidth',2);

legend({['Acutal Flux'],['Bowtie Analysis'],['LSQR'],['Standard Deviation']},...
                 'Location', 'northeast','FontSize',18);

%legend({['Theoretical Flux'],['Incident Particle Measurement'],['Bowtie Analysis'],['LSQR'],['Standard Deviation']},...
 %                'Location', 'northeast','FontSize',18);

textsize = 24;
set(gca, 'FontSize', textsize)
%xlim([0 8])
%ylim([10^0 10^4])
xticks((0:1:8))
%set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

ylabel('I  #/(cm^2 sr s MeV)','FontSize',textsize)
xlabel('Energy (MeV)','FontSize',textsize)
hold off
