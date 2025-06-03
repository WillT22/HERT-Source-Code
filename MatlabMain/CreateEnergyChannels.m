
clc
clear all
%% Electrons
%{
ExLow_Limit = 0.1;
Low_Limit = 0.5; %0.5
Med_Limit = 1.0; %1.0
High_Limit = 2.5; %2.5
Max_Limit = 8; %8
%}

%% Protons
%
ExLow_Limit = 0.1;
Low_Limit = 0.5;
Med_Limit = 100;
%High_Limit = 50;
%Max_Limit = 100;
%

% Number of channels in each range:
ExLownum = 1; %2
Lownum = 39; %12
%Mednum = 16; %16
%Highnum = 10; %10

% Excluded energy channels (below threshold)
x = 1:(ExLownum+1);

DE = (log10(Low_Limit(1))-log10(ExLow_Limit(1)))/ExLownum; %energy resolution

for i = 1:(length(x)-1)
    ExLow_Limit(i+1) = 10^(DE+log10(ExLow_Limit(i)));
    
    ExLowChannels(i,1) = ExLow_Limit(i);
    ExLowChannels(i,2) = ExLow_Limit(i+1);
end

% Low energy channels
x = 1:(Lownum+1);

DE = (log10(Med_Limit(1))-log10(Low_Limit(1)))/Lownum;


for i = 1:(length(x)-1)
    Low_Limit(i+1) = 10^(DE+log10(Low_Limit(i)));
    
    LowChannels(i,1) = Low_Limit(i);
    LowChannels(i,2) = Low_Limit(i+1);
end
%{
% Mid-range energy channels
x_med = 1:(Mednum+1);

DE_med = (log10(High_Limit(1))-log10(Med_Limit(1)))/Mednum;

for i = 1:(length(x_med)-1)
    Med_Limit(i+1) = 10^(DE_med+log10(Med_Limit(i)));
    
    MedChannels(i,1) = Med_Limit(i);
    MedChannels(i,2) = Med_Limit(i+1);
end

% High energy channels
x_high = 1:(Highnum+1);

DE_high = (log10(Max_Limit)-log10(High_Limit(1)))/Highnum;

for i = 1:(length(x_high)-1)
    High_Limit(i+1) = 10^(DE_high+log10(High_Limit(i)));
    
    HighChannels(i,1) = High_Limit(i);
    HighChannels(i,2) = High_Limit(i+1);
end
%}
% Combining energy channels into one variable
%Channels = [ExLowChannels;LowChannels;MedChannels;HighChannels;MaxChannels]
Channels = [ExLowChannels;LowChannels]
NumChannels = length(Channels)

% Combining enery resolution of each channel
for i = 1:length(Channels)
   Resolution(i,1) = 100*(-Channels(i,1)+Channels(i,2))/((Channels(i,1)+Channels(i,2))/2);
end
Resolution

%% Write to Text File

fileID = fopen('channel_select\proton_channels_v1.txt','w');
for i = 1:size(Channels,1)
fprintf(fileID,'%6.3f,%6.3f \n',Channels(i,:));
end
fclose(fileID);
