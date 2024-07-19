%% Inputs
geo_EC = readmatrix('E:\HERT_Drive\Matlab Main\Result\geometric_factor_EC.txt');
dt = 1;

% Creating Initialized Matrices
sigma_d = sqrt(hits_whole_EC*dt+1)/dt;
for channel = 1:size(energy_channels,1)
   Cd(channel,channel) = (sigma_d(channel)/(hits_whole_EC(channel))).^2;
end

sigma_m = 18;
delta = 12;
for bin = 1:bins
   Cm2(:,bin) = sigma_m^2 * exp(-((energy_midpoints-energy_midpoints(bin)).^2)./(2*delta.^2));
end

x_edges = log(energy_edges);
x = log(energy_midpoints);
dx = x_edges(2:end)-x_edges(1:end-1);
%dx(dx>100) = 100;

mn = zeros(it_max,bins);
%mn(1,:) = ones(1,bins).*log(10^3.5);
mn(1,:) = log(fit_result.a * exp(fit_result.b * energy_midpoints)...
        +fit_result.c * exp(fit_result.d * energy_midpoints));

g_mn = zeros(it_max,size(energy_channels,1));
g_mn(1,:) = log(sum(geo_EC.*exp(mn(1,:)+x).*dx,2));

Gn2 = 1./exp(g_mn(1,:)') .* geo_EC.*exp(mn(1,:)+x).*dx;

Sn = Cd + Gn2*Cm2*Gn2';
B = Cm2*Gn2'*Sn^-1;

d_obs = log(hits_whole_EC');

iteration = 2;
it_max = 1000;
convergence = false;

while convergence == false && iteration <= it_max
    mn(iteration,:) = mn(1,:)' + B*(d_obs-g_mn(iteration-1,:)'+Gn2*(mn(iteration-1,:)-mn(1,:))');
    
    g_mn(iteration,:) = log(sum(geo_EC.*exp(mn(iteration,:)+x).*dx,2));
    Gn2 = 1./exp(g_mn(iteration,:)') .* geo_EC.*exp(mn(iteration,:)+x).*dx;
    Sn = Cd + Gn2*Cm2*Gn2';
    B = Cm2*Gn2'*Sn^-1;

    Cmm2 =  Cm2 - Cm2*Gn2'*(Sn^-1)*Gn2*Cm2;
    jn = exp(mn(iteration,:));
    jsig = jn.*sqrt(diag(Cmm2))';

    if max(abs(mn(iteration,:)- mn(iteration-1,:)))<0.1
        convergence = true;
        disp("Converges")
    else
        iteration = iteration+1;
    end
end

flux_lsqr = exp(mn(iteration,:));

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
plot(energy_midpoints(indices),flux_lsqr(indices)+jsig(indices),'c--','LineWidth',2);
plot(energy_midpoints(indices),flux_lsqr(indices)-jsig(indices),'c--','LineWidth',2);

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
%ylim([10^3, 10^4])

%ylabel('Flux')
ylabel('I  #/(cm^2 sr s MeV)','FontSize',textsize)
xlabel('Energy (MeV)','FontSize',textsize)
hold off