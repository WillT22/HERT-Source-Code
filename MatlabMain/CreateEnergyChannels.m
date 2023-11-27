
clc
clear all
%% Electrons
%{
ExLimit = 0.1;
Limit = 0.5;
MedLimit = 1.0;
High_Limit = 2.5;
Max_Limit = 7;
%}

%% Protons
%
ExLimit = 15;
Limit = 20;
MedLimit = 30.0;
High_Limit = 70.0;
Max_Limit = 120;
%

% Number of channels in each range:
ExLownum = 2;
Lownum = 12; %12
Mednum = 16; %16
Highnum = 10; %10

% Excluded energy channels (below threshold)
x = 1:(ExLownum+1);

DE = (log10(Limit(1))-log10(ExLimit(1)))/ExLownum; %energy resolution

for i = 1:(length(x)-1)
    ExLimit(i+1) = 10^(DE+log10(ExLimit(i)));
    
    ExLowChannels(i,1) = ExLimit(i);
    ExLowChannels(i,2) = ExLimit(i+1);
end

% Low energy channels
x = 1:(Lownum+1);

DE = (log10(MedLimit(1))-log10(Limit(1)))/Lownum;


for i = 1:(length(x)-1)
    Limit(i+1) = 10^(DE+log10(Limit(i)));
    
    LowChannels(i,1) = Limit(i);
    LowChannels(i,2) = Limit(i+1);
end

% Mid-range energy channels
x_med = 1:(Mednum+1);

DE_med = (log10(High_Limit(1))-log10(MedLimit(1)))/Mednum;

for i = 1:(length(x_med)-1)
    MedLimit(i+1) = 10^(DE_med+log10(MedLimit(i)));
    
    MedChannels(i,1) = MedLimit(i);
    MedChannels(i,2) = MedLimit(i+1);
end

% High energy channels
x_high = 1:(Highnum+1);

DE_high = (log10(Max_Limit)-log10(High_Limit(1)))/Highnum;

for i = 1:(length(x_high)-1)
    High_Limit(i+1) = 10^(DE_high+log10(High_Limit(i)));
    
    HighChannels(i,1) = High_Limit(i);
    HighChannels(i,2) = High_Limit(i+1);
end

% Combining energy channels into one variable
Channels = [ExLowChannels;LowChannels;MedChannels;HighChannels]
NumChannels = length(Channels)

% Combining enery resolution of each channel
for i = 1:length(Channels)
   Resolution(i,1) = 100*(-Channels(i,1)+Channels(i,2))/((Channels(i,1)+Channels(i,2))/2);
end
Resolution

%% Write to Text File

fileID = fopen('proton_channels_v1.txt','w');
for i = 1:size(Channels,1)
fprintf(fileID,'%6.3f,%6.3f \n',Channels(i,:));
end
fclose(fileID);
