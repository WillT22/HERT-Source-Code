% Least Squares and Count Rate Analysis for HERT
% Last modified: 5/29/2024

%warning('off','last')

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
    it_max = 100;   % maximum number of possible iterations 

% Setting up constant matrices
    Cd = zeros(size(energy_channels,1));
    sigma_d = zeros(size(energy_channels,1),1);
    for channel = 1:size(energy_channels,1)
        sigma_d(channel) = sqrt(hits_whole_EC(channel)*dt+1)/dt;
        Cd(channel,channel) = (sigma_d(channel)/(hits_whole_EC(channel))).^2;;
    end
    inv_Cd = inv(Cd);

    % Initialize variance parameter
    sigma_m = 80;
    %sigma_m_array = logspace(0,5,50);
    %sigma_m_array = linspace(20,100,41);

    % Initialize smoothness parameter
    delta = 38;
    %delta_array = logspace(0,5,50);
    %delta_array = linspace(30,60,31);


% Loop over variance parameter
%{
ind_error_avg = zeros(length(delta_array),length(sigma_m_array)); 
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
    % Create C_m covariance matrix
    Cm = zeros(bins);
    Cmm = zeros(bins);
    for bin = 1:bins
        Cm(:,bin) = sigma_m^2 * exp(-((energy_midpoints-energy_midpoints(bin)).^2)./(2*delta^2));
    end
    inv_Cm = inv(Cm);
    
% Defining constants over each iteration
    d_obs = log(hits_whole_EC');

% Defining initial values 
    mn = zeros(it_max,bins);
    mn(1,:) = ones(1,bins);%.*log(10^3.5);
    % starting from best fit line of numeric lsqr
    %mn(1,:) = log(fit_result.a * exp(fit_result.b * energy_midpoints));
    %mn(1,:) = log(fit_result.a * exp(fit_result.b * energy_midpoints)...
    %    +fit_result.c * exp(fit_result.d * energy_midpoints));
    
    x_edges = log(energy_edges);
    x = log(energy_midpoints);
    dx = x_edges(2:end)-x_edges(1:end-1);
    dx(dx>100) = 100;
    
    integ = geo_EC .* dx .* exp(mn(1,:)+x);
    g_mn = zeros(it_max,size(energy_channels,1));
    g_mn(1,:) = log(sum(integ,2));
    Gn = 1./exp(g_mn(1,:)') .* integ;
    mat_mult = inv(Gn'* inv_Cd * Gn + inv_Cm) * Gn' * inv_Cd;
   

    iteration = 2;
    convergence = false;

    % error matrix initialization
    error = zeros(it_max,size(energy_channels,1));
    error_avg = zeros(it_max,1);
    error_max = zeros(it_max,1);
    
    % error inital calculations
    error(1,:) = (d_obs-g_mn(1,:)').^2 ./ sigma_d;
    error_max(1) = max(error(1,:));
    error_avg(1) = mean(error(1,:));

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

    error(iteration,:) = (d_obs-g_mn(iteration,:)').^2 ./ sigma_d;
    error_max(iteration) = max(abs(error(iteration,:)));
    error_avg(iteration) = mean(abs(error(iteration,:)));

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
    %ind_error_avg(d,sf_i) = error_avg(iteration);
    %ind_error_max(d,sf_i) = error_max(iteration);
    %ind_act_error_avg(d,sf_i) = actual_error_avg(iteration);
    %ind_error_avg(sf_i) = error_avg(iteration);
    %ind_error_max(sf_i) = error_max(iteration);
    %ind_act_error_avg(sf_i) = actual_error_avg(iteration);
end

%end
%end

%% Plots calculated flux
f = figure;
f.Position = [0 0 1700 900];
hold on

% Plot simulated flux
plot(energy_midpoints,M_energy_bin/(4*pi^2*r_source^2)./bin_width, 'Color', 'black','LineWidth',2.5);
%fit_flux_result = fit(energy_midpoints',(M_energy_bin/(4*pi^2*r_source^2)./bin_width)','exp2');
%fitted_flux_curve = @(x) feval(fit_flux_result, x);
%plot(energy_midpoints, fitted_flux_curve(energy_midpoints),...
%    'LineStyle', '--','LineWidth',2,'Color', 'black');

%
% Plot Bowtie points
%plot(E_eff,j_nom,'o', 'Color', '#0072BD','MarkerSize',10);

% Plot LSQR Selesnick Method
plot(energy_midpoints(indices),flux_lsqr(indices),'.', 'Color', 'r', 'MarkerSize',12);
plot(energy_midpoints(indices),flux_lsqr(indices)+jsig(indices),'blue--','LineWidth',2);
plot(energy_midpoints(indices),flux_lsqr(indices)-jsig(indices),'blue--','LineWidth',2);

% Plot LSQR Summation Method
%plot(energy_midpoints(indices),flux_approx(indices),'x','Color','#77AC30','MarkerSize',10)
% and fitted curve
%fitted_curve = @(x) feval(fit_result, x);
%plot(energy_midpoints, fitted_curve(energy_midpoints),...
%    'LineStyle', '--','LineWidth',2,'Color', '#77AC30');

%{
legend({['Incident Particle Measurements'],['Bowtie Method'], ...
                 ['Least Squares Method (Selesnick)'], ['Least Squares Method (Simple)']},...
                 'Location', 'northeast','FontSize',18);
%}
%
legend({['Incident Particle Measurements'],['LSQR'],['Standard Deviation']},...
                 'Location', 'northeast','FontSize',18);
%

set(gca, 'FontSize', textsize)
xlim([0 8])
xticks((0:1:8))
set(gca, 'YScale', 'log')
%ylim([10^3, 10^4])

%ylabel('Flux')
ylabel('I  #/(cm^2 sr s MeV)','FontSize',textsize)
xlabel('Energy (MeV)','FontSize',textsize)
hold off

%{
logic_array = zeros(size(ind_error_avg));
for i = 1:length(delta_array)
    for j = 1:length(sigma_m_array)
        if ind_error_avg(i,j) < 0.00125 && ind_error_avg(i,j) > 0
            logic_array(i,j) = 1;
        end
    end
end

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