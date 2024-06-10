% Least Squares and Count Rate Analysis for HERT
% Last modified: 5/29/2024

%% Inputs
geo_EC = readmatrix('E:\HERT_Drive\Matlab Main\Result\geometric_factor_EC.txt');
dt = 1;

%% Least Squares Function for Energy Channels (initial attempt)
flux_approx = lsqr(geo_EC.*bin_width,hits_whole_EC',1e-6,20)';

% using idicies for E>0.5
indices = energy_midpoints > 1;
fit_result = fit(energy_midpoints(indices)',flux_approx(indices)','exp1');

%% Least Squares Function for Energy Channels (Selesnick/Khoo)
% Initialize variables
    it_max = 10000;   % maximum number of possible iterations 

% Setting up constant matrices
    G = geo_EC.*bin_width;
    
    C_d = zeros(size(energy_channels,1));
    for channel = 1:size(energy_channels,1)
        C_d(channel,channel) = (hits_whole_EC(channel)*dt+1)/dt^2;
    end
    inv_C_d = inv(C_d);
    
    C_m = zeros(bins);
    for bin= 1:bins
        C_m(bin,bin) = 5000;
    end
    inv_C_m = inv(C_m);
    
    mat_mult = inv(G'* inv_C_d * G + inv_C_m) * G' * inv_C_d;

% Defining constants over each iteration
    d_obs = log(hits_whole_EC');

% Defining initial values 
    mn = zeros(bins,it_max);
    %mn(:,1) = ones(1,bins).*log(10^3);
    % starting from best fit line of numeric lsqr
    mn(:,1) = log(fit_result.a * exp(fit_result.b * energy_midpoints));
    %mn(:,1) = log(fit_result.a * exp(fit_result.b * energy_midpoints)...
    %    +fit_result.c * exp(fit_result.d * energy_midpoints));

    
    x_edges = log(energy_edges);
    x = log(energy_midpoints');
    dx = x_edges(2:end)-x_edges(1:end-1);
    dx(dx>100) = 100;
    
    integ = geo_EC .* dx .* exp(mn(:,1) + x)';
    G_n = 1./hits_whole_EC' .* integ;
    
    g_mn = zeros(size(energy_channels,1),it_max);
    g_mn(:,1) = log(sum(integ,2));
      
    error = zeros(it_max,1);
    error(1) = mean(abs((d_obs-g_mn(:,1))./d_obs));

    iteration = 2;
    error_cond = false;
% Begin iterations
while error_cond == false && iteration <= it_max
    mn(:,iteration) = mn(:,iteration-1) + mat_mult * (d_obs - g_mn(:,iteration-1) + G_n*(mn(:,iteration-1)-mn(:,1)));
    
    integ = geo_EC .* dx .* exp(mn(:,iteration) + x)';
    G_n = 1./hits_whole_EC' .* integ;
    g_mn(:,iteration) = log(sum(integ,2));
    error(iteration) = mean(abs((d_obs-g_mn(:,iteration))./d_obs));
    if error(iteration) > error(iteration-1)
        error_cond = true;
    else
        iteration = iteration + 1;
    end
end

fprintf("Iteration Number: %.0d \n",iteration-1)
fprintf("Tolerance: %.6e \n",error(iteration-1))

flux_lsqr = exp(mn(:,iteration-1));

%% Plots calculated flux
f = figure;
f.Position = [0 0 1700 900];
hold on

% Plot simulated flux
plot(energy_midpoints,M_energy_bin/(4*pi^2*r_source^2)./bin_width, 'Color', 'black','LineWidth',2.5);

%
% Plot Bowtie points
plot(E_eff,j_nom,'o', 'Color', '#0072BD','MarkerSize',10);

% Plot LSQR Selesnick Method
plot(energy_midpoints(indices),flux_lsqr(indices),'.', 'Color', 'r', 'MarkerSize',12);

% Plot LSQR Summation Method
plot(energy_midpoints(indices),flux_approx(indices),'x','Color','#77AC30','MarkerSize',10)
% and fitted curve
fitted_curve = @(x) feval(fit_result, x);
plot(energy_midpoints, fitted_curve(energy_midpoints),...
    'LineStyle', '--','LineWidth',2,'Color', '#77AC30');

%
legend({['Incident Particle Measurements'],['Bowtie Method'], ...
                 ['Least Squares Method (Selesnick)'], ['Least Squares Method (Summation)']},...
                 'Location', 'northeast','FontSize',18);

%legend({['Incident Particle Measurements'],['Least Squares Method (Selesnick)'],['Least Squares Method (Summation)']},...
%                 'Location', 'northeast','FontSize',18);


set(gca, 'FontSize', textsize)
xlim([0 8])
xticks((0:1:8))
set(gca, 'YScale', 'log')
ylim([10^2.5, 10^4])

%ylabel('Flux')
ylabel('I  #/(cm^2 sr s MeV)','FontSize',textsize)
xlabel('Energy (MeV)','FontSize',textsize)
hold off
