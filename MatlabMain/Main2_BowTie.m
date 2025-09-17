%%BowTie.m
% BowTie and Count Rate Analysis for HERT
% Last modified: 5/29/2024

clc
close all

%% Bow Tie Analysis-Selesnick/Blake
%Changes directory for png of graphs
cd 'E:\HERT_Drive\Matlab Main\Bow Tie'

geo_EC = readmatrix('E:\HERT_Drive\Matlab Main\Result\geofactor_EC_FS.txt');
valid_geo_EC = geo_EC(:,1:size(energy_midpoints,2));
%hits_whole_EC= ones(size(geo_EC,1),1);

%Sets Ei and Range of Eo
%Ei is the incident energy from the GEANT4 results
%Eo set by user.
Eo = 0.2:0.2:2.0;

%Sets up color vectors for plotting the different Eo curves
Eo_color = magma(length(Eo)+1);

%Preallocates all variables prior to For Loops
J_e = zeros(length(energy_midpoints),length(Eo));

G_int = zeros(length(energy_midpoints),length(Eo),length(energy_channels));
G_term = zeros(length(energy_midpoints),length(Eo),length(energy_channels));
G_E_eff = zeros(length(energy_midpoints),length(Eo),length(energy_channels));

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
    J_e(:,i) = exp(-energy_midpoints/Eo(i));
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

for u = 1:length(energy_channels)
    fwhm(u) = findFWHM(energy_midpoints,valid_geo_EC(u,:));
    % Print full width half max values into command window
    fprintf('%.2f - %.2f MeV: %.4f\n',energy_channels(u,1),energy_channels(u,2),fwhm(u))
end

%For Loop for calculating a line for each Eo and finding the average intersection point
fprintf('Energy Channel Processing: ')
for c=1:length(energy_channels)
    % Display the percentage complete (every 10% for example)
    percentage = floor(c / length(energy_channels) * 100);  
    if (mod(c / length(energy_channels) * 100, 10) == 0)
        fprintf('%.0f ',percentage);
    end
        
    %Calculates line for each Eo
    for i = 1:length(Eo)
        valid_geo = ~isnan(valid_geo_EC(c,:));
        G_int(valid_geo,i,c) = valid_geo_EC(c,valid_geo)'.*J_e(valid_geo,i);
        G_term(:,i,c) = trapz(energy_midpoints,G_int(:,i,c));
        G_E_eff(:,i,c)= G_term(:,i,c).*J_e_inv(:,i);  
    end
    
    %Find Intersections of Eo Lines
    num = 0;
    for j = 1:(length(Eo)-1)
        for k = 1:(length(Eo)-j)
            [xi(num+k,c),yi(num+k,c)] = polyxpoly(energy_midpoints,G_E_eff(:,j,c),energy_midpoints,G_E_eff(:,j+k,c),'unique');   
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
%{
for c=1:height(energy_channels)
    
    f = figure;
    f.Position = [100 100 1000 720];
    hold on
    for i = 1:length(Eo)
        plot(Ei,G_E_eff(:,i,c),'Color',Eo_color(i,:),'DisplayName',BowTieLegend(1,i),'Linewidth',2);
    end
    
    plot(xi(:,c),yi(:,c),'*b','DisplayName',BowTieLegend(1,end),'MarkerSize',12)
    plot(E_eff(c),G_eff_dE(c),'o','Linewidth',2,'DisplayName',BowTieLegend(1,end-1),'Color','green','MarkerFaceColor', 'green','MarkerSize',12)
    ylim_l = round(min(yi(:,c))*0.95,3);
    ylim_u = round(max(yi(:,c))*1.05,3);
    ylim([ylim_l ylim_u])
    plot([0 E_eff(c)],[G_eff_dE(c) G_eff_dE(c)],'--g')
    plot([E_eff(c) E_eff(c)],[0 G_eff_dE(c)],'--g')
    set(gca,'FontSize',18)
    labelpoints(E_eff(c),G_eff_dE(c),append('E_Eff = ',num2str(E_eff(c)),' MeV'),'SE',0.15)
    %labelpoints(E_eff(c),G_eff_dE(c),append('G_eff_dE = ',num2str(G_eff_dE(c)),' cm^2 sr MeV'),'SE',0.65)
    
    legend(BowTieLegend,'Location', 'southoutside','NumColumns',round(length(BowTieLegend)/2))
    xlim_l = round(min(xi(:,c))*0.95,3);
    xlim_u = round(max(xi(:,c))*1.05,3);
    xlim([xlim_l xlim_u])
    %xticks((1.5:0.05:2.25))
    title(append('Bow Tie Analysis Energy Channel ',num2str(c),' Deposited Energy Range: ',num2str(energy_channels(c,1)),' to ',num2str(energy_channels(c,2)),' MeV'))
    ylabel('G_{eff} * \Delta E ','FontSize',textsize)
    xlabel('Nominal Energy (MeV)','FontSize',textsize)
    hold off
    
    effsave = append(date(),'Bow Tie Energy Channel ',num2str(c),' Eo ',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV','_',addin,num2str(length(energy_channels)),'.png');
    saveas(gcf,effsave)  
end
%}

%Plots graph of all the average intersection points of each energy channel
f = figure;
f.Position = [100 100 1200 720];
hold on
for c=1:length(energy_channels)
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

effsave = append(date(),' Bow Tie  All Energy Channels',' Eo',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV','_',addin,num2str(length(energy_channels)),'.png');
saveas(gcf,effsave)

% Plot nominal effective energies with counts
%{
f = figure;
f.Position = [100 100 1200 720];
hold on
for c=1:length(energy_channels)
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
for c = 1:length(energy_channels)
    bar(E_eff(c),G_eff_dE(c)/fwhm(c),fwhm(c),'EdgeColor','k','FaceColor',Effplotcolor(c,:))
    
end
legend(EngLegend,'Location', 'southoutside','NumColumns',8)
%title(append('Energy Channel Bins-',' Eo Range: ',num2str(Cen(Eo)),' to ',num2str(max(Eo)),' MeV'),'FontSize',20)
ylabel('Effective Geometric Factor','FontSize',28)
xlabel('Nominal Effective Energy (MeV)','FontSize',28)
hold off

effsave = append(date(),' Energy Channel Bins',' Eo ',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV','_',addin,num2str(length(energy_channels)),'.png');
saveas(gcf,effsave)
%}

%% Energy Resolution Plot
Energy_Resolution= 100*BinWidth./E_eff;

% Create a matrix with the rounded values
data_to_export = [energy_channels(Energy_Resolution~=0,:), E_eff(Energy_Resolution~=0)', Energy_Resolution(Energy_Resolution~=0)'];

% Write the matrix to a text file
dlmwrite('output.txt', data_to_export, 'delimiter', '\t', 'precision', '%.4f');

%{
figure
plot(Energy_Resolution,'xb')
title('Energy Resolution per Energy Channel')
xlabel('Energy Channel Number')
ylabel('Energy Resolution (%)')

energy_channel_list = 1:1:length(energy_channels);
Effplotcolor = plasma(length(energy_channel_list));
f = figure;
f.Position = [100 100 1600 720];
textsize = 28;

%Plots each energy channel FWHM value
hold on
plot(E_eff(Energy_Resolution~=0),Energy_Resolution(Energy_Resolution~=0),'o','MarkerSize',8,...
        'MarkerEdgeColor','b','MarkerFaceColor','b')
set(gca,'FontSize',textsize)
%title('HERT Energy Resolution','FontSize',textsize)
ylabel('Spectral Resolution dE/E(%)','FontSize',textsize)
xlabel('Nominal Energy (MeV)','FontSize',textsize)
ylim([0,40])
plot([0,max(max(M_energy_beam))],[12,12],'k--','LineWidth',2)%,'DisplayName','Energy Resolution Requirement')

%legend([EngLegend,'Energy Resolution Requirement'],'Location', 'southoutside','NumColumns',8)
legend(['Energy Resolution Requirement'],'Location', 'northeast','NumColumns',8)

hold off
effsave = append(date(),' Energy Resolution',' Eo ',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV','_',addin,num2str(length(energy_channels)),'.png');
saveas(gcf,effsave)
%}

%% Write to Text File
%{
fileID = fopen('effective_energies_DARTBe.txt','w');
for i = 1:length(E_eff)
fprintf(fileID,'%.6f \n',E_eff(i));
end
fclose(fileID);
%}


