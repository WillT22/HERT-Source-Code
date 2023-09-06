function BuildCountTable(energy_channels,count_rates,exportname)

cd 'Bow Tie'


% Inputs
twidth = 4;
tlength = length(energy_channels);
titles = {'Channel Number','Deposited Energy Range','Count Rate Lower Limit','Count Rate Upper Limit'};


% Memeory Prelocation
table = cell(tlength,twidth);

% Build BinWidth Table
for j = 1:length(energy_channels)
    
    table{j,1} = j;
    table{j,2} = sprintf('%.2f - %.2f',energy_channels(j,1),energy_channels(j,2));
    table{j,3} = num2str(round(count_rates(1,j),4));
    table{j,4} = round(count_rates(2,j),4);
    
end

FinalTable = [titles;table];

% Write the FWHM tables to their Respective Folder

exportfile = append(exportname,'.xls');
writecell(FinalTable,exportfile)


cd ..