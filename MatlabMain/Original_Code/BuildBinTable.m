function BuildBinTable(energy_channels,E_eff,Bin_Width,G_Eff,ER,Apo25,Apo50,Apo75,exportname)

cd 'D:\HERT\MATLAB Main\Bow Tie'


% Inputs

tlength = length(energy_channels);
titles = {'Channel Number','Deposited Energy Range','Ei (MeV)','Bin Width (MeV)','G_Eff (cm^2 sr)','Energy Resolution','25th at Apo','50th at Apo','75th at Apo'};
twidth = length(titles);

% Memeory Prelocation
table = cell(tlength,twidth);

% Build BinWidth Table
for j = 1:length(energy_channels)
    
    %channel = sprintf('%.2f - %.2f',energy_channels(j,1),energy_channels(j,2));
    table{j,1} = j;
    table{j,2} = sprintf('%.2f - %.2f',energy_channels(j,1),energy_channels(j,2));
    table{j,3} = num2str(round(E_eff(j),4));
    table{j,4} = round(Bin_Width(j),4);
    table{j,5} = round(G_Eff(j),2);
    table{j,6} = round(ER(j),2);
    table{j,7} = round(Apo25(j),2);
    table{j,8} = round(Apo50(j),2);
    table{j,9} = round(Apo75(j),2);
    
    
end

FinalTable = [titles;table];

% Write the FWHM tables to their Respective Folder

exportfile = append(exportname,'.xls');
writecell(FinalTable,exportfile)


cd ..