function buildtable(energy_channels,fwhm,exportname)

cd Plots/Tables

% Inputs
twidth = 1;
uniquec = 3;
tlength = length(energy_channels);
titles = {'Energy Channel','FWHM','FWHM/ Geometric Center'};

% Duplicate titles to match parameters
t = 1;
titlestr = cell(1,twidth*uniquec);
for i = 1:twidth
    titlestr(t:t+uniquec-1) = titles;
    t = t + uniquec;
end

% Calculate Geometric Centers of Energy Channels
geo_center = zeros(1,length(energy_channels));
for i = 1:length(energy_channels)
    ind1 = energy_channels(i,1);
    ind2 = energy_channels(i,2);
    geo_center(i) = 10^((log10(ind2)+log10(ind1))/2);
end

% Memeory Prelocation
table = cell(3,length(energy_channels));

% Build FWHM table
for i = 1:length(energy_channels)
    channel = sprintf('%.2f - %.2f',energy_channels(i,1),energy_channels(i,2));
    table{i,1} = channel;
    table{i,2} = round(fwhm(i),6);
    table{i,3} = append(num2str(100*round(fwhm(i)/geo_center(i),3)),'%');
end

% Format table to appropriate size
c = 1;
cc = 1;
t = 1;
final_table = cell(tlength,uniquec*twidth);
for i = 1:twidth
    final_table(1:tlength,c:c+uniquec-1) = table(t:t+tlength-1,cc:cc+uniquec-1);
    t = t + tlength;
    c = c + uniquec;
end

final_table = [titlestr;final_table];
% Write the FWHM tables to their Respective Folder

exportfile = append(exportname,'.xls');
writecell(final_table,exportfile)

cd ..
cd ..

