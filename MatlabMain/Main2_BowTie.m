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
    energy_max_select = 10; % MeV max for electrons
    energy_range = energy_midpoints(energy_midpoints < energy_max_select);
    
    % Sets Ei and Range of Eo
    % Ei is the incident energy from the GEANT4 results, Eo set by user.
    Eo = 0.2:0.2:2.0; 
    energy_channels_use = energy_channels;
    valid_geo_EC = geo_EC(:, 1:size(energy_range,2));
    
elseif parttype == 1    % --- PROTONS ---
    energy_max_select = 1000; % MeV max for protons (D9 triggering @ 52.29 MeV)
    energy_range = energy_midpoints(energy_midpoints < energy_max_select);
    
    calc_type = 'pen'; % select calculation type: 'full', 'range', 'pen'
    % Sets Ei and Range of Eo
    % Ei is the incident energy from the GEANT4 results, Eo set by user.
    if strcmp(calc_type, 'full')
        Eo = 19:0.5:38; % full,  (x > 14.15) & (x < 1000)
        energy_channels_use = energy_channels;
        valid_geo_EC = geo_EC(:, 1:size(energy_range,2));
    elseif strcmp(calc_type, 'range')
        energy_channels_use = energy_channels;
        Eo = 8:0.5:16; % range, (x > 14.15) & (x < 52.29)
        valid_geo_EC = geo_EC(:, 1:size(energy_range,2));
    elseif strcmp(calc_type, 'pen')
        energy_channels_use = energy_channels_combined;
        Eo = 19:0.5:38; % pen,   (x > 52.29) & (x < 1000)
        valid_geo_EC = geo_back_EC_combined(:, 1:size(energy_range,2));
    end
end

% Apply the selected bounds
energy_range = energy_midpoints(energy_midpoints < energy_max_select);

%Sets up color vectors for plotting the different Eo curves
Eo_color = magma(length(Eo)+1);

%Preallocates all variables prior to For Loops
J_e = zeros(length(energy_range),length(Eo));

G_int = zeros(length(energy_range),length(Eo),length(energy_channels_use));
G_term = zeros(length(energy_range),length(Eo),length(energy_channels_use));
G_E_eff = zeros(length(energy_range),length(Eo),length(energy_channels_use));

xi = zeros(sum(1:length(Eo)-1),length(energy_channels_use));
yi = zeros(sum(1:length(Eo)-1),length(energy_channels_use));

E_eff = zeros(1,length(energy_channels_use));
G_eff_dE= zeros(1,length(energy_channels_use));
BowTieLegend = strings([1,length(Eo)]);

%Count_Rate =  zeros(length(energy_channels_use),length(Eo));

BinWidth = zeros(1, length(energy_channels_use));
Geff = zeros(1, length(energy_channels_use));

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
fwhm = zeros(1,length(energy_channels_use));
E_max = zeros(1,length(energy_channels_use));

for u = 1:length(energy_channels_use)
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
    fprintf('%.2f - %.2f MeV: %.4f\n',energy_channels_use(u,1),energy_channels_use(u,2),fwhm(u))
end

%For Loop for calculating a line for each Eo and finding the average intersection point
fprintf('Energy Channel Processing: ')

next_threshold = 10; 
num_channels = length(energy_channels_use);
for c = 1:num_channels
    % [Your existing processing code goes here]
    
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
end
%Calculate Bin Characteristics
Geff = G_eff_dE./fwhm;
BinWidth = fwhm;

% Calculate Flux
j_nom = hits_whole_EC./G_eff_dE';

fprintf('\n');

%% Plots for Config
%{
textsize = 28;
%Plots Graph for each energy channel
ymax = round(max(G_eff_dE)*1.1,3);
%
for c=1:height(energy_channels_use)
    
    f = figure;
    f.Position = [100 100 1000 720];
    hold on
    for i = 1:length(Eo)
        plot(energy_range,G_E_eff(:,i,c),'Color',Eo_color(i,:),'DisplayName',BowTieLegend(1,i),'Linewidth',2);
    end
    
    plot(xi(:,c),yi(:,c),'*b','DisplayName',BowTieLegend(1,end),'MarkerSize',12)
    plot(E_eff(c),G_eff_dE(c),'o','Linewidth',2,'DisplayName',BowTieLegend(1,end-1),'Color','green','MarkerFaceColor', 'green','MarkerSize',12)
    y_min = min(yi(:,c));
    y_max = max(yi(:,c));
    y_span = y_max - y_min;
    % Fallback in case all data points are exactly the same (flat line)
    if y_span == 0
        if y_max == 0
            y_span = 1; % Arbitrary span if all data is exactly 0
        else
            y_span = abs(y_max); % Scale buffer to the data magnitude
        end
    end
    % Calculate limits using a 5% buffer based on the total range
    ylim_l = round(y_min - (0.05 * y_span), 3);
    ylim_u = round(y_max + (0.05 * y_span), 3);
    % Safety check in case extreme rounding forces limits to be identical
    if ylim_l >= ylim_u
        ylim_u = ylim_l + 0.001; 
    end
    
    ylim([ylim_l ylim_u])
    plot([0 E_eff(c)],[G_eff_dE(c) G_eff_dE(c)],'--g')
    plot([E_eff(c) E_eff(c)],[0 G_eff_dE(c)],'--g')
    set(gca,'FontSize',18)
    text(E_eff(c), G_eff_dE(c), append('  E_Eff = ', num2str(E_eff(c)), ' MeV'), ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left');
    
    legend(BowTieLegend,'Location', 'southoutside','NumColumns',round(length(BowTieLegend)/2))
    xlim_l = round(min(xi(:,c))*0.95,3);
    xlim_u = round(max(xi(:,c))*1.05,3);
    xlim([xlim_l xlim_u])
    %xticks((1.5:0.05:2.25))
    title(append('Bow Tie Analysis Energy Channel ',num2str(c),' Deposited Energy Range: ',num2str(energy_channels_use(c,1)),' to ',num2str(energy_channels_use(c,2)),' MeV'))
    ylabel('G_{eff} * \Delta E ','FontSize',textsize)
    xlabel('Peak Energy (MeV)','FontSize',textsize)
    hold off
    
    effsave = append(date(),'Bow Tie Energy Channel ',num2str(c),' Eo ',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV','_pen_',addin,num2str(length(energy_channels_use)),'.png');
    saveas(gcf,effsave)  
end
%

%Plots graph of all the average intersection points of each energy channel
f = figure;
f.Position = [100 100 1200 720];
hold on
for c=1:length(energy_channels_use)
    plot(E_eff(c),G_eff_dE(c),'o','Color',Effplotcolor(c,:))
end
legend(EngLegend,'Location', 'southoutside','NumColumns',8)
%ylim([0 ymax])
%yticks((0:0.005:ymax))
xlim([1 7])
xticks((1:1:7))

title(append('Bow Tie Analysis Configuration All Energy Channels ',' Eo Range: ',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV'))
ylabel('G_{eff} * \Delta E ')
xlabel('Nominal Effective Energy (MeV)')
hold off

effsave = append(date(),' Bow Tie  All Energy Channels',' Eo',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV','_pen_',addin,num2str(length(energy_channels_use)),'.png');
saveas(gcf,effsave)
%}

% Plot nominal effective energies with counts
%{
f = figure;
f.Position = [100 100 1200 720];
hold on
for c=1:length(energy_channels_use)
    plot(E_eff(c),hits_whole_EC(c),'o','Color',Effplotcolor(c,:))
end
legend(EngLegend,'Location', 'southoutside','NumColumns',8)

xlim([1 7])
xticks((1:1:7))

title(append('Bow Tie Analysis Configuration All Energy Channels ',' Eo Range: ',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV'))
ylabel('Counts')
xlabel('Nominal Effective Energy (MeV)')
hold off
%}

%Plots each energy channel FWHM value in the same color as the Eff. Curve
f = figure;
f.Position = [100 100 1600 720];
hold on
for c = 1:length(energy_channels_use)
    bar(E_eff(c),G_eff_dE(c)/fwhm(c),fwhm(c),'EdgeColor','k','FaceColor',Effplotcolor(c,:))
    
end
legend(EngLegend,'Location', 'southoutside','NumColumns',8)
%title(append('Energy Channel Bins-',' Eo Range: ',num2str(Cen(Eo)),' to ',num2str(max(Eo)),' MeV'),'FontSize',20)
ylabel('Effective Geometric Factor','FontSize',28)
xlabel('Nominal Effective Energy (MeV)','FontSize',28)
hold off

effsave = append(date(),' Energy Channel Bins',' Eo ',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV','_pen_',addin,num2str(length(energy_channels_use)),'.png');
saveas(gcf,effsave)
%

%% Energy Resolution Plot
Energy_Resolution = 100 * BinWidth ./ E_max;
% Create a matrix with the rounded values
data_to_export = [energy_channels_use(Energy_Resolution~=0,:), E_max(Energy_Resolution~=0)', Energy_Resolution(Energy_Resolution~=0)'];
% Write the matrix to a text file
dlmwrite('output.txt', data_to_export, 'delimiter', '\t', 'precision', '%.4f');

figure
plot(Energy_Resolution,'xb')
title('Energy Resolution per Energy Channel')
xlabel('Energy Channel Number')
ylabel('Energy Resolution (%)')

energy_channel_list = 1:1:length(energy_channels_use);
Effplotcolor = plasma(length(energy_channel_list));

f = figure;
f.Position = [100 100 1600 720];
textsize = 28;

hold on
plot([0,max(max(M_energy_beam))],[12,12],'k--','LineWidth',2)

% Scatter plot applies the sliced colormap matrix to the sliced data points
scatter(E_max(3:end-1), Energy_Resolution(3:end-1), 64, Effplotcolor(3:end-1, :), 'filled');

set(gca,'FontSize',textsize)
ylabel('Energy Resolution dE/E(%)','FontSize',textsize)
xlabel('Peak Energy (MeV)','FontSize',textsize)

if parttype == 0
    xlim([0,8])
elseif parttype == 1
    if strcmp(calc_type, 'full') || strcmp(calc_type, 'range')
        xlim([0,120])
    elseif strcmp(calc_type, 'pen')
        xlim([0,520])
    end
end
ylim([0,30])

legend('Energy Resolution Requirement','Location', 'northeast')
hold off

effsave = append(date(),' Energy Resolution',' Eo ',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV_',calc_type,'_',addin,num2str(length(energy_channels_use)),'.png');
saveas(gcf,effsave)
close all
%

%% Write to Text File
%{
% Replace '.txt' with nothing
clean_name = strrep(Selected_Channel_name, '.txt', '');
% Append the new text
filename = [clean_name, '_resolution_',calc_type,'.txt'];
fileID = fopen(filename, 'w');
for i = 1:length(E_max)
fprintf(fileID,'%.6f %.6f\n',E_max(i),Energy_Resolution(i));
end
fclose(fileID);
%}
%{
# Ragne Channels
data_to_save = [E_max', geo_back_EC];
writematrix(data_to_save, 'proton_FS_14000_v6_range_GFbyEC.txt');

# Penetrating Channels
data_to_save = [E_max', geo_back_EC_combined];
writematrix(data_to_save, 'proton_FS_14000_v6_pen_GFbyEC.txt');
%}


