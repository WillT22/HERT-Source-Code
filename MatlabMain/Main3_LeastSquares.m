% Least Squares and Count Rate Analysis for HERT
% Last modified: 5/29/2024

%warning('off','last')

%% Inputs
geo_EC = readmatrix('E:\HERT_Drive\Matlab Main\Result\geometric_factor_EC.txt');
bins = size(geo_EC,2);
energy_edges = linspace(0,8,bins+1);
bin_width = energy_edges(2:end)-energy_edges(1:end-1);
energy_midpoints = (energy_edges(2:end) + energy_edges(1:end-1))/2;

energy_channels = readmatrix('E:\HERT_Drive\Matlab Main\Result\channel_select\electron_channels_v1.txt');

r_source = 8.5;
flux_fov = 10^4*exp(-energy_midpoints/1.5)./(4*pi^2*r_source^2)./bin_width;
hits_whole_EC = zeros(1,size(geo_EC,1));
for bin = 1:bins
    for channel = 1:size(energy_channels,1)
        if energy_midpoints(bin)>=energy_channels(channel,1) && energy_midpoints(bin)<energy_channels(channel,2)
            hits_whole_EC(channel) = hits_whole_EC(channel) + flux_fov(bin);
        end
    end
end
 
hits_whole_EC(hits_whole_EC==0)=1;

flux = flux_fov.*(4*pi^2*r_source^2).*bin_width;

dt = 1;

%% Least Squares Function for Energy Channels (Selesnick/Khoo)
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

    % Initialize variance parameter
    sigma_m = 80;
    %sigma_m_array = logspace(0,5,50);
    %sigma_m_array = linspace(70,90,81);

    % Initialize smoothness parameter
    delta = 39;
    %delta_array = logspace(0,5,50);
    %delta_array = linspace(35,45,41);


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
    Cm = zeros(bins);
    Cmm = zeros(bins);
    for bin = 1:bins
        Cm(:,bin) = sigma_m^2 * exp(-((energy_midpoints-energy_midpoints(bin)).^2)./(2*delta^2));
    end
    inv_Cm = inv(Cm);
    
% Defining constants over each iteration
    d_obs = log(hits_whole_EC');

% Defining initial values 
    % Create initial estimate
    mn = zeros(it_max,bins);
    mn(1,:) = ones(1,bins);
    
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
    actual_flux = log(M_energy_bin/(4*pi^2*r_source^2)./bin_width);
    actual_error = zeros(it_max,bins);
    actual_error_avg = zeros(it_max,1);
    actual_error_max = zeros(it_max,1);
    actual_error(1,:) = mn(1,:)-actual_flux;
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
    actual_error(iteration,:) = mn(iteration,:)-actual_flux;
    actual_error_max(iteration) = max(abs(actual_error(iteration,:)));
    actual_error_avg(iteration) = mean(abs(actual_error(iteration,:)));

    if max(abs(mn(iteration,:)- mn(iteration-1,:)))<0.08
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
    %ind_act_error_avg(sf_i) = actual_error_avg(iteration);
end

%end
%end

%% Plots calculated flux
f = figure;
f.Position = [0 0 1700 900];
hold on

% Plot simulated flux
plot(energy_midpoints,flux, 'Color', 'black','LineWidth',2.5);

%plot(energy_midpoints,M_energy_bin/(4*pi^2*r_source^2)./bin_width,'x', 'Color', 'black','MarkerSize',8);
%plot(energy_midpoints,10^6*0.5*exp(-(energy_midpoints)/0.5), 'Color', 'black','LineWidth',2.5);

% Plot Bowtie points
%plot(E_eff,j_nom,'o', 'Color', '#0072BD','MarkerSize',10);

% Plot LSQR Selesnick Method
plot(energy_midpoints,flux_lsqr, 'Color', 'r', 'LineWidth',3);
plot(energy_midpoints,flux_lsqr+jsig,'r--','LineWidth',2);
plot(energy_midpoints,flux_lsqr-jsig,'r--','LineWidth',2);

legend({['Theoretical Flux'],['LSQR'],['Standard Deviation']},...
                 'Location', 'northeast','FontSize',18);

set(gca, 'FontSize', textsize)
xlim([0 8])
ylim([10^0 10^6])
xticks((0:1:8))
set(gca, 'YScale', 'log')

ylabel('I  #/(cm^2 sr s MeV)','FontSize',textsize)
xlabel('Energy (MeV)','FontSize',textsize)
hold off

%{
logic_array_act = zeros(size(ind_act_error_avg));
for i = 1:length(delta_array)
    for j = 1:length(sigma_m_array)
        if ind_act_error_avg(i,j) < 0.025 && ind_act_error_avg(i,j) > 0
            logic_array_act(i,j) = 1;
        end
    end
end

value= min(ind_act_error_avg(ind_act_error_avg>0),[],'all');
[index1,index2] = find(ind_act_error_avg == value);

f = figure;
f.Position = [0 0 1700 900];
hold on
plot(sigma_m_array./delta_array',ind_act_error_avg,'.', 'Color', 'black', 'MarkerSize',12);
ylabel('Actual Error','FontSize',textsize)
xlabel('Sigma/Delta','FontSize',textsize)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

f = figure;
f.Position = [0 0 1700 900];
hold on
plot(sigma_m_array,ind_act_error_avg,'.', 'Color', 'black', 'MarkerSize',12);
ylabel('Actual Error','FontSize',textsize)
xlabel('Sigma','FontSize',textsize)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

f = figure;
f.Position = [0 0 1700 900];
hold on
plot(delta_array,ind_act_error_avg,'.', 'Color', 'black', 'MarkerSize',12);
ylabel('Actual Error','FontSize',textsize)
xlabel('Delta','FontSize',textsize)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off
%}