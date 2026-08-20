%%BowTie.m
% BowTie and Count Rate Analysis for HERT
% Last modified: 5/29/2024

clc
close all

%% Bow Tie Analysis-Selesnick/Blake
%Changes directory for png of graphs
cd 'D:\HERT_Drive\Matlab Main\Bow Tie'
%geo_EC = readmatrix('D:\HERT_Drive\Matlab Main\Result\geofactor_EC_FS.txt');

if parttype == 0    % --- ELECTRONS ---
    energy_max_select = 8; % MeV max for electrons
    energy_range = energy_midpoints(energy_midpoints < energy_max_select);
    
    % Sets Ei and Range of Eo
    % Ei is the incident energy from the GEANT4 results, Eo set by user.
    Eo = 0.2:0.2:2.0; 
    valid_geo_EC = geo_EC(:, 1:size(energy_range,2));
    
elseif parttype == 1    % --- PROTONS ---
    energy_max_select = 1000; % MeV max for protons (D9 triggering @ 52.29 MeV)
    energy_range = energy_midpoints(energy_midpoints < energy_max_select);
    
    calc_type = 'pen'; % select calculation type: 'full', 'range', 'pen'
    
    % Sets Ei and Range of Eo
    % Ei is the incident energy from the GEANT4 results, Eo set by user.
    if strcmp(calc_type, 'full')
        Eo = 19:0.5:38; % full,  (x > 14.15) & (x < 1000)
        valid_geo_EC = geo_EC(:, 1:size(energy_range,2));
    elseif strcmp(calc_type, 'range')
        Eo = 8:0.5:16; % range, (x > 14.15) & (x < 52.29)
        valid_geo_EC = geo_EC(:, 1:size(energy_range,2));
    elseif strcmp(calc_type, 'pen')
        Eo = 19:0.5:38; % pen,   (x > 52.29) & (x < 1000)
        valid_geo_EC = geo_back_EC(:, 1:size(energy_range,2));
    end
end

% Apply the selected bounds
energy_range = energy_midpoints(energy_midpoints < energy_max_select);

%Sets up color vectors for plotting the different Eo curves
Eo_color = magma(length(Eo)+1);

%Preallocates all variables prior to For Loops
J_e = zeros(length(energy_range),length(Eo));
G_int = zeros(length(energy_range),length(Eo),length(energy_channels));
G_term = zeros(length(energy_range),length(Eo),length(energy_channels));
G_E_eff = zeros(length(energy_range),length(Eo),length(energy_channels));
xi = zeros(sum(1:length(Eo)-1),length(energy_channels));
yi = zeros(sum(1:length(Eo)-1),length(energy_channels));
E_eff = zeros(1,length(energy_channels));
G_eff_dE= zeros(1,length(energy_channels));
BowTieLegend = strings([1,length(Eo)]);
%Count_Rate =  zeros(length(energy_channels),length(Eo));
BinWidth = zeros(1, length(energy_channels));
Geff = zeros(1, length(energy_channels));

%Creates J(e) and creates String Array for Plot Legends
for i = 1:length(Eo)
    J_e(:,i) = exp(-energy_range/Eo(i));
    BowTieLegend(i) = num2str(Eo(i));
end

%Adds 'Average Intersection Point' for plotting
BowTieLegend = [BowTieLegend,'Intersection Point','Average Intersection Point'];

%Creates J(e)^-1
J_e_inv = 1./J_e;

% Finds FWHM for each energy channel
fprintf('\nFull Width at Half Max Values:\n')

% Preallocates FWHM vectors
fwhm = zeros(1,length(energy_channels));
E_max = zeros(1,length(energy_channels));

for u = 1:length(energy_channels)
    if parttype == 0
        [fwhm(u),E_max(u)] = findFWHM_limited(energy_range,valid_geo_EC(u,:),parttype);
    elseif parttype == 1
        if strcmp(calc_type, 'range')
            [fwhm(u),E_max(u)] = findFWHM_limited(energy_range,valid_geo_EC(u,:),parttype);
        elseif strcmp(calc_type, 'pen')
            [fwhm(u),E_max(u)] = findFWHM_max(energy_range,valid_geo_EC(u,:),parttype);
        end
    end
    % Print full width half max values into command window
    fprintf('%.2f - %.2f MeV: %.4f\n',energy_channels(u,1),energy_channels(u,2),fwhm(u))
end

%For Loop for calculating a line for each Eo and finding the average intersection point
fprintf('Energy Channel Processing: ')
next_threshold = 10; 
num_channels = length(energy_channels);

for c = 1:num_channels
    % Calculate current completion percentage
    percentage = (c / num_channels) * 100;
    
    % Check if we crossed one or more 10% thresholds
    while percentage >= next_threshold
        fprintf('%.0f%% ', next_threshold);
        next_threshold = next_threshold + 10;
    end
        
    %Calculates line for each Eo
    for i = 1:length(Eo)
        valid_geo = ~isnan(valid_geo_EC(c,:));
        G_int(valid_geo,i,c) = valid_geo_EC(c,valid_geo)'.*J_e(valid_geo,i);
        G_term(:,i,c) = trapz(energy_range,G_int(:,i,c));
        G_E_eff(:,i,c)= G_term(:,i,c).*J_e_inv(:,i);  
    end
    
    %Find Intersections of Eo Lines
    num = 0;
    for j = 1:(length(Eo)-1)
        for k = 1:(length(Eo)-j)
            [xi(num+k,c),yi(num+k,c)] = polyxpoly(energy_range,G_E_eff(:,j,c),energy_range,G_E_eff(:,j+k,c),'unique');   
        end
        num = num + k;
    end
    
    %Finds average intersection point
    E_eff(c) = mean(xi(:,c));
    G_eff_dE(c) = mean(yi(:,c));
    % Calculates the Bow Tie systematic error (Standard Deviation)
    E_eff_err(c) = std(xi(:,c));
    G_eff_dE_err(c) = std(yi(:,c));
end

%Calculate Bin Characteristics
Geff = G_eff_dE./fwhm;
BinWidth = fwhm;
Energy_Resolution = 100 * BinWidth ./ E_max;

% Calculate Flux
N_valid = min(length(hits_whole_EC),40);
j_nom = hits_whole_EC(1:N_valid)./G_eff_dE(1:N_valid)';

% Assuming C is your count rate matrix
relative_error_counts = sqrt(hits_whole_EC(1:N_valid)) ./ hits_whole_EC(1:N_valid);
relative_error_bowtie = (G_eff_dE_err(1:N_valid) ./ G_eff_dE(1:N_valid))';

% Total fractional uncertainty
total_relative_error = sqrt(relative_error_counts.^2 + relative_error_bowtie.^2);
j_nom_err = j_nom .* total_relative_error;
fprintf('\n');


%% Variable Archiving by Calculation Type (Protons Only)
% Dynamically copies all calculated variables and appends _range or _pen to their names
if parttype == 1
    var_list = {'J_e', 'J_e_inv', 'G_int', 'G_term', 'G_E_eff', 'xi', 'yi', ...
                'E_eff', 'G_eff_dE', 'E_eff_err', 'G_eff_dE_err', 'fwhm', ...
                'E_max', 'Geff', 'BinWidth', 'j_nom', 'j_nom_err', ...
                'Energy_Resolution'};
            
    for v = 1:length(var_list)
        % e.g., creates E_max_range = E_max;
        eval([var_list{v}, '_', calc_type, ' = ', var_list{v}, ';']);
    end
end


%% Energy Resolution Plot (Dual Colormap for Protons)
% 0. Setup and Typography
textsize = 14;              
idx = 3:(length(energy_channels)-4);
tick_spacing = 3; 

f = figure('Units', 'inches', 'Position', [1, 1, 7.0, 3.3], ...
           'PaperUnits', 'inches', 'PaperPosition', [0, 0, 7.0, 3.3], 'PaperPositionMode', 'auto'); 
ax1 = axes(f); hold(ax1, 'on');

if parttype == 0
    % ==========================================
    % ELECTRON PLOT 
    % ==========================================
    % Changed starting X from 0 to 0.5 to prevent log(0) error
    plot(ax1, [0.5, 520], [15, 15], 'k--', 'LineWidth', 1.5, 'DisplayName', 'Energy Resolution Requirement');
    scatter(ax1, E_max(idx), Energy_Resolution(idx), 40, idx, 'filled', 'HandleVisibility', 'off');
    colormap(ax1, plasma(length(energy_channels)));
    
    cb1 = colorbar(ax1);
    clim(ax1, [min(idx) - 0.5, max(idx) + 0.5]);
    cb1.Ticks = min(idx):tick_spacing:max(idx); 
    cb1.Label.String = 'Channel Number'; cb1.Label.FontSize = textsize; cb1.FontSize = textsize - 2; 
    
    title(ax1, 'Electron Energy Channel Resolution');
    xlim(ax1, [0.5, 6]);
    ylim(ax1, [0, 30]);
    
    % Targeted ax1 directly and removed '0' from the tick array
    set(ax1, 'XTick', [0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
    legend(ax1, 'Location', 'northeast', 'FontSize', textsize - 2);
    plot_name = 'electron_FS_31000_v2_publishable_energyres_Epeak';

elseif parttype == 1
    % ==========================================
    % PROTON PLOT (Combined)
    % ==========================================
    scatter(ax1, E_max_range(idx), Energy_Resolution_range(idx), 40, idx, 'filled');
    colormap(ax1, plasma(length(energy_channels)));
    
    cb1 = colorbar(ax1);
    clim(ax1, [min(idx) - 0.5, max(idx) + 0.5]);
    cb1.Ticks = min(idx):tick_spacing:max(idx); 
    cb1.Label.String = 'Range Channel Number'; cb1.Label.FontSize = textsize; cb1.FontSize = textsize - 2; 
    
    % Invisible Overlay Axes for Penetrating Protons
    ax2 = axes('Position', [0.10, 0.18, 0.68, 0.72], 'Color', 'none', 'XColor', 'none', 'YColor', 'none'); 
    hold(ax2, 'on');
    
    % Corrected scatter marker syntax to 's'
    scatter(ax2, E_max_pen(idx), Energy_Resolution_pen(idx), 40, idx, 's', 'filled');
    colormap(ax2, winter(size(energy_channels, 1)));
    
    cb2 = colorbar(ax2);
    clim(ax2, [min(idx) - 0.5, max(idx) + 0.5]);
    cb2.Ticks = min(idx):tick_spacing:max(idx); 
    cb2.Label.String = 'Penetrating Channel Number'; cb2.Label.FontSize = textsize; cb2.FontSize = textsize - 2; 
    
    % Lock Layout and Scales
    set(cb2, 'Position', [0.9, 0.18, 0.025, 0.72]); 
    linkaxes([ax1, ax2], 'xy');    
    
    % Plot invisible points with NaN coordinates just to generate legend graphics
    h_range = plot(ax1, NaN, NaN, 'o', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'DisplayName', 'Range Channels');
    h_pen   = plot(ax1, NaN, NaN, 's', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'DisplayName', 'Penetrating Channels');
    legend(ax1, [h_range, h_pen], 'Location', 'northwest', 'FontSize', textsize - 2);
    
    xlim(ax1, [10, 550]);
    ylim(ax1, [0, 120]);

    set(ax1, 'XTick', [10, 20, 30, 50, 100, 200, 300, 500, 1000]);
    t = title(ax1, 'Proton Energy Channel Resolution');
    plot_name = 'proton_FS_17000_v6_combined_publishable_energyres_Epeak';
end

% ==========================================
% SHARED FORMATTING
% ==========================================
set(ax1, 'FontSize', textsize, 'Layer', 'top');
set(ax1, 'Position', [0.10, 0.18, 0.68, 0.72]);
set(cb1, 'Position', [0.8, 0.18, 0.025, 0.72]); 

% Apply Log Scale to the X-Axis
set(ax1, 'XScale', 'log');
if exist('ax2', 'var')
    set(ax2, 'XScale', 'log'); 
end

ylabel(ax1, 'Energy Resolution dE/E (%)', 'FontSize', textsize);
xlabel(ax1, 'Peak Energy (MeV)', 'FontSize', textsize);
grid(ax1, 'on'); ax1.YMinorGrid = 'off';

hold(ax1, 'off'); if exist('ax2', 'var'), hold(ax2, 'off'); end

% High-Resolution Output
print(f, [plot_name, '.png'], '-dpng', '-r300');
exportgraphics(f, [plot_name, '.pdf'], 'ContentType', 'vector');

%% Write to Text File
%{
% Replace '.txt' with nothing
clean_name = strrep(Selected_Channel_name, '.txt', '');

% Append the new text
filename = [clean_name, '_resolution_',calc_type,'.txt'];
fileID = fopen(filename, 'w');
for i = 1:length(E_max)
    % When accessing variables here after a run, you can now use E_max_range or E_max_pen
    fprintf(fileID,'%.6f %.6f\n',E_max(i),Energy_Resolution(i));
end
fclose(fileID);
%}

%{
% Range Channels
data_to_save = [E_max_range', geo_back_EC];
writematrix(data_to_save, 'proton_FS_14000_v6_range_GFbyEC.txt');

% Penetrating Channels
data_to_save = [E_max_pen', geo_back_EC];
writematrix(data_to_save, 'proton_FS_14000_v6_pen_GFbyEC.txt');
%}