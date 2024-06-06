% Least Squares and Count Rate Analysis for HERT
% Last modified: 5/29/2024

%% Inputs
geo_EC = readmatrix('E:\HERT_Drive\Matlab Main\Result\geometric_factor_EC.txt');
dt = 1;

%% Least Squares Function for Energy Channels (initial attempt)
C_d = zeros(size(energy_channels,1));
for channel = 1:size(energy_channels,1)
    C_d(channel,channel) = (hits_whole_EC(channel)*dt+1)/dt^2;
end
inv_C_d = inv(C_d);

C_m = zeros(bins);
for bin= 1:bins
    C_m(bin,bin) = 0.85;
end
inv_C_m = inv(C_m);

flux_approx = lsqr(geo_EC,hits_whole_EC',1e-6,20)'/bin_width;

% Define chi-square function
%chisq = @(m, d, G) (G * m - d)'/(G * C_m * G' + C_d)*(G * m - d);

% Use 'lsqnonlin' for least-squares minimization (similar to 'lm' in scipy)
%m_output = lsqnonlin(@(flux) chisq(flux,log(hits_whole_EC'),geo_EC), zeros(bins,1));
%int_flux = exp(m_output);


%% Least Squares Function for Energy Channels (Selesnick/Khoo)
% Initialize variables
C_d = zeros(size(energy_channels,1));
for channel = 1:size(energy_channels,1)
    C_d(channel,channel) = (hits_whole_EC(channel)*dt+1)/dt^2;
end
inv_C_d = inv(C_d);

C_m = zeros(bins);
for bin= 1:bins
    C_m(bin,bin) = 8;
end
inv_C_m = inv(C_m);



mn(iteration+1,:) = mn(iteration,:) + inv(G'*inv_C_d*G+inv_C_m)...
    *G'*inv_C_d*(d_obs-g_mn+G_n*(mn(iteration,:)-mn(1,:)));

% Initialize iteration variables
it_max = 1000;
mn = zeros(bins,it_max);
mn(:,1) = ones(bins,1).*log(10^3);
G = geo_EC.*bin_width;
d_0 = log(sum(geo_EC*exp(mn(:,1) + log(energy_midpoints')),2));
G_n = 1./hits_whole_EC' .* geo_EC.*exp(mn(:,1) + log(energy_midpoints'))';
d_obs = log(hits_whole_EC');
error = zeros(it_max,1);
error(1) = mean(abs(log(hits_whole_EC')-d_n));

iteration = 1;
% Iterate to desired tolerance or iterations
while error(iteration) >= 1e-8 && iteration <= it_max
    iteration = iteration + 1;
    G_n = 1./hits_whole_EC' .* geo_EC.*exp(mn(:,iteration) + log(energy_midpoints'))';
    d_n = d_n + G_n * (mn(:,iteration)-mn(:,1));
    mn(:,iteration+1) = mn(:,iteration) +  inv(G' * inv_C_d * G + inv_C_m)...
        * G' * inv_C_d * (d_obs - d_n + G_n*(mn(:,iteration)-mn(:,1)));
    error(iteration) = mean(abs(d_obs-d_n));
end

fprintf("Iteration Number: %.0d \n",iteration-1)
fprintf("Tolerance: %.6e \n",error(iteration))

flux_lsqr = exp(mn(:,iteration));

%% Plots calculated flux
f = figure;
f.Position = [0 0 1700 900];
hold on

% Plot simulated flux
plot(energy_midpoints,M_energy_bin/(4*pi^2*r_source^2)/bin_width, 'Color', 'black','LineWidth',2.5);

%{
% Plot Bowtie points and line of best fit
plot(E_eff,j_nom,'o', 'Color', '#0072BD','MarkerSize',10);
fit_result = fit(E_eff',j_nom','exp1');
fitted_curve = @(x) feval(fit_result, x);
plot(E_eff, fitted_curve(E_eff), 'LineStyle', '--','LineWidth',2,'Color', '#0072BD');
%}

% using idicies for E>0.5
indices = energy_midpoints > 1;

% Plot LSQR Selesnick Method
plot(energy_midpoints(indices),flux_lsqr(indices),'.', 'Color', 'r', 'MarkerSize',12);

% Plot LSQR Khoo Method
%plot(energy_midpoints(indices),int_flux(indices),'o', 'Color', 'b', 'MarkerSize',10);

% Plot LSQR Summation Method
plot(energy_midpoints(indices),flux_approx(indices),'x','Color','#77AC30','MarkerSize',10)


%{
legend({['Incident Particle Measurements'],['Bowtie Method'],['Bowtie Method (Best Fit)'], ...
                 ['Least Squares Method (Integral)'], ['Least Squares Method (Summation)'],['Least Squares Method (Summation Best Fit)']},...
                 'Location', 'northeast','FontSize',18);
%
legend({['Incident Particle Measurements'], ['Least Squares Method (Selesnick)'],...
    ['Least Squares Method (Khoo)'],['Least Squares Method (Will)']},...
                 'Location', 'northeast','FontSize',18);
%}
legend({['Incident Particle Measurements'],['Least Squares Method (Selesnick)'],['Least Squares Method (Summation)']},...
                 'Location', 'northeast','FontSize',18);

set(gca, 'FontSize', textsize)
xlim([0 8])
xticks((0:1:8))
set(gca, 'YScale', 'log')
%ylim([10^5.5, 10^6.5])

%ylabel('Flux')
ylabel('I (#/(cm^2 sr s MeV)','FontSize',textsize)
xlabel('Energy (MeV)','FontSize',textsize)
hold off
