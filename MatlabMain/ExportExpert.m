% Explot Least-Squares necessary variables to Python

fileID = fopen('count_rate.txt','w');
for i = 1:length(energy_channels)
    fprintf(fileID,'%.0f\n',hits_whole_EC(i));
end
fclose(fileID);

fileID = fopen('geo_factor.txt','w');
for i = 1:length(energy_channels)
    for j = 1:bins
        fprintf(fileID,'%.16f ',geo_EC(i,j));
    end
    fprintf(fileID,'\n');
end
fclose(fileID);

fileID = fopen('energy_bins.txt','w');
for i = 1:length(M_energy_midpoints)
    fprintf(fileID,'%.2f\n',M_energy_midpoints(i));
end
fclose(fileID);

fileID = fopen('bowtie_flux.txt','w');
for i = 1:length(j_nom)
    fprintf(fileID,'%.16f\n',j_nom(i));
end
fclose(fileID);