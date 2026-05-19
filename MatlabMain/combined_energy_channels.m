%% Load Data into Structures to Prevent Variable Overwriting
% Replace these filenames with your actual .mat file paths
data1 = load("D:\HERT_Drive\Matlab Main\Result\Proton_FS\proton_FS_14000_v6_range.mat"); 
data2 = load("D:\HERT_Drive\Matlab Main\Result\Proton_FS\proton_FS_14000_v6_pen.mat");

% Calculate Energy Resolution for Both Sets
ER_1 = 100 .* data1.BinWidth ./ data1.E_max;
ER_2 = 100 .* data2.BinWidth ./ data2.E_max;

%% Generate Master Colormaps
% 1. Get the total number of channels from the raw data
N_range = length(data1.E_max);
N_pen   = length(data2.E_max);

% 2. Define the specific slices you want to plot 
idx_range = 1:(N_range - 4); 
idx_pen   = 2:(N_pen - 4);

% 3. Create the exact same master color palettes used in your geometric factor scripts
c_plasma = plasma(N_range);

% Shift the winter colormap exactly how you did for penetrating channels
full_winter = winter(N_range); 
c_winter_pen = full_winter((end - N_pen + 1):end, :);

%% Setup Figure
f = figure;
f.Position = [100 100 1600 720];
textsize = 28;
hold on

% 1. Plot the Energy Resolution Requirement Line
max_peak_energy = max([max(data1.E_max), max(data2.E_max)]);
h_req = plot([0, max_peak_energy], [12, 12], 'k--', 'LineWidth', 2, 'DisplayName', 'Energy Resolution Requirement');

% 2. Plot Dataset 1 (Range - Plasma)
% We add 'HandleVisibility', 'off' so the legend strictly ignores this scatter plot
scatter(data1.E_max(idx_range), ER_1(idx_range), 100, c_plasma(idx_range, :), 'o', 'filled', 'HandleVisibility', 'off');

% 3. Plot Dataset 2 (Penetrating - Winter)
% Also ignored by the legend
scatter(data2.E_max(idx_pen), ER_2(idx_pen), 100, c_winter_pen(idx_pen, :), 's', 'filled', 'HandleVisibility', 'off');

% 4. Dummy Plots for a Clean Legend
% We capture the outputs as handles (h_range, h_pen) and set the faces to black ('k')
h_range = plot(NaN, NaN, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Range Channels');

h_pen = plot(NaN, NaN, 's', 'MarkerSize', 10, 'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Penetrating Channels');

% Apply Formatting
set(gca, 'FontSize', textsize)
ylabel('Energy Resolution dE/E (%)', 'FontSize', textsize)
xlabel('Peak Energy (MeV)', 'FontSize', textsize)

% X-Axis Limits
xlim([0, 520]) 

% Add the dynamic legend by EXPLICITLY passing all three handles
legend([h_req, h_range, h_pen], 'Location', 'northwest', 'NumColumns', 1, 'FontSize', 20)
grid on
hold off

%% Save Figure
effsave = append(string(datetime("today")), '_Dual_Energy_Resolution.png');
saveas(gcf, effsave)

%% --- LaTeX Table Generator to File ---
% Define the output filename
filename = 'Energy_Channels_Table.txt';

% Open the file for writing
fileID = fopen(filename, 'w');

if fileID == -1
    error('Cannot open file for writing.');
end

% ==========================================
% Subtable 1: Electron Channels
% ==========================================
fprintf(fileID, '%% --- Subtable 1: Electrons ---\n');
fprintf(fileID, '\\begin{table}[H]\n');
fprintf(fileID, '    \\centering\n');
fprintf(fileID, '    \\caption{Energy Channel Parameters for Electrons and Protons}\n');
fprintf(fileID, '    \\label{tab:energy_channels}\n');
fprintf(fileID, '    \\begin{subtable}{\\textwidth}\n');
fprintf(fileID, '        \\centering\n');
fprintf(fileID, '        \\caption{Electron Energy Channel Parameters: Core channels are highlighted in light gray. Energy resolution for supplemental and integral channels is not physically meaningful and is therefore not shown for these channels.}\n');
fprintf(fileID, '        \\label{tab:electron_channels}\n');
fprintf(fileID, '        \\begin{tabular}{cccccc}\n');
fprintf(fileID, '            \\toprule\n');
fprintf(fileID, '            Color & Min $E_{dep}$ (MeV) & Max $E_{dep}$ (MeV) & Peak Energy (MeV) & FWHM (MeV) & $\\Delta E/E$ (\\%%) \\\\\n');
fprintf(fileID, '            \\midrule\n');

% Loop through electron data
N_e = length(E_max);
for i = 1:N_e
    hex_color = sprintf('%02X%02X%02X', round(Effplotcolor(i, :) * 255));
    
    % Blank out first two (1, 2) and last channel (N_e)
    if i == 1 || i == 2 || i == N_e
        fprintf(fileID, '            \\textcolor[HTML]{%s}{\\LARGE $\\bullet$} & %.3f & %.3f & %.3f & --- & --- \\\\\n', ...
            hex_color, energy_channels(i,1), energy_channels(i,2), E_max(i));
    else
        % Added \rowcolor{gray!10} for core channels
        fprintf(fileID, '            \\rowcolor{gray!10} \\textcolor[HTML]{%s}{\\LARGE $\\bullet$} & %.3f & %.3f & %.3f & %.3f & %.3f \\\\\n', ...
            hex_color, energy_channels(i,1), energy_channels(i,2), E_max(i), fwhm(i), Energy_Resolution(i));
    end
end

fprintf(fileID, '            \\bottomrule\n');
fprintf(fileID, '        \\end{tabular}\n');
fprintf(fileID, '    \\end{subtable}\n');
fprintf(fileID, '\\end{table}\n\n');

% ==========================================
% Subtable 2: Proton Range Channels
% ==========================================
fprintf(fileID, '%% --- Subtable 2: Proton Range ---\n');
fprintf(fileID, '\\begin{table}[htbp]\n');
fprintf(fileID, '    \\ContinuedFloat %% Tells LaTeX this is a continuation of the previous table\n');
fprintf(fileID, '    \\centering\n');
fprintf(fileID, '    \\begin{subtable}{\\textwidth}\n');
fprintf(fileID, '        \\centering\n');
fprintf(fileID, '        \\caption{Range Proton Energy Channel Parameters: Range protons are classified as depositing ${>}0.1$ MeV onto detector 1 while not depositing ${>}0.1$ MeV onto detector 9. Core channels are highlighted in light gray. Energy resolution for supplemental and integral channels is not physically meaningful and is therefore not shown for these channels.}\n');
fprintf(fileID, '        \\label{tab:proton_range}\n');
fprintf(fileID, '        \\begin{tabular}{cccccc}\n');
fprintf(fileID, '            \\toprule\n');
fprintf(fileID, '            Color & Min $E_{dep}$ (MeV) & Max $E_{dep}$ (MeV) & Peak Energy (MeV) & FWHM (MeV) & $\\Delta E/E$ (\\%%) \\\\\n');
fprintf(fileID, '            \\midrule\n');

% Loop through Proton Range data
N_pr = length(data1.E_max);
for i = 1:N_pr
    hex_color = sprintf('%02X%02X%02X', round(c_plasma(i, :) * 255)); 
    
    % Blank out the last four channels
    if i >= (N_pr - 3)
        fprintf(fileID, '            \\textcolor[HTML]{%s}{\\LARGE $\\bullet$} & %.3f & %.3f & %.3f & --- & --- \\\\\n', ...
            hex_color, data1.energy_channels(i,1), data1.energy_channels(i,2), data1.E_max(i));
    else
        % Added \rowcolor{gray!10} for core channels
        fprintf(fileID, '            \\rowcolor{gray!10} \\textcolor[HTML]{%s}{\\LARGE $\\bullet$} & %.3f & %.3f & %.3f & %.3f & %.3f \\\\\n', ...
            hex_color, data1.energy_channels(i,1), data1.energy_channels(i,2), data1.E_max(i), data1.fwhm(i), data1.Energy_Resolution(i));
    end
end

fprintf(fileID, '            \\bottomrule\n');
fprintf(fileID, '        \\end{tabular}\n');
fprintf(fileID, '    \\end{subtable}\n');
fprintf(fileID, '\\end{table}\n\n');

% ==========================================
% Subtable 3: Proton Penetrating Channels
% ==========================================
fprintf(fileID, '%% --- Subtable 3: Proton Penetrating ---\n');
fprintf(fileID, '\\begin{table}[htbp]\n');
fprintf(fileID, '    \\ContinuedFloat\n');
fprintf(fileID, '    \\centering\n');
fprintf(fileID, '    \\begin{subtable}{\\textwidth}\n');
fprintf(fileID, '        \\centering\n');
fprintf(fileID, '        \\caption{Penetrating Proton Energy Channel Parameters: Penetrating protons are classified as depositing ${>}0.1$ MeV onto each detector 1 and detector 9. Core channels are highlighted in light gray. Energy resolution for supplemental and integral channels is not physically meaningful and is therefore not shown for these channels.}\n');
fprintf(fileID, '        \\label{tab:proton_pen}\n');
fprintf(fileID, '        \\begin{tabular}{cccccc}\n');
fprintf(fileID, '            \\toprule\n');
fprintf(fileID, '            Color & Min $E_{dep}$ (MeV) & Max $E_{dep}$ (MeV) & Peak Energy (MeV) & FWHM (MeV) & $\\Delta E/E$ (\\%%) \\\\\n');
fprintf(fileID, '            \\midrule\n');

% Loop through Proton Penetrating data
N_pp = length(data2.E_max);
for i = 1:N_pp
    hex_color = sprintf('%02X%02X%02X', round(c_winter_pen(i, :) * 255));
    
    % Blank out the first one (1) and last four channels
    if i == 1 || i >= (N_pp - 3)
        fprintf(fileID, '            \\textcolor[HTML]{%s}{\\LARGE $\\bullet$} & %.3f & %.3f & %.3f & --- & --- \\\\\n', ...
            hex_color, data2.energy_channels_combined(i,1), data2.energy_channels_combined(i,2), data2.E_max(i));
    else
        % Added \rowcolor{gray!10} for core channels
        fprintf(fileID, '            \\rowcolor{gray!10} \\textcolor[HTML]{%s}{\\LARGE $\\bullet$} & %.3f & %.3f & %.3f & %.3f & %.3f \\\\\n', ...
            hex_color, data2.energy_channels_combined(i,1), data2.energy_channels_combined(i,2), data2.E_max(i), data2.fwhm(i), data2.Energy_Resolution(i));
    end
end

fprintf(fileID, '            \\bottomrule\n');
fprintf(fileID, '        \\end{tabular}\n');
fprintf(fileID, '    \\end{subtable}\n');
fprintf(fileID, '\\end{table}\n');

% Close the file
fclose(fileID);

disp(['LaTeX table successfully written to ', filename]);