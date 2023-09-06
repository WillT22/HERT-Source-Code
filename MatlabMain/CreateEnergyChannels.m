
clc
clear all
ExLimits(1) = 0.1;
Limits(1) = 0.5;
MedLimits(1) = 1.0;
High_Limits(1) = 2.5;

ExLownum = 2;
x = 1:(ExLownum+1);

DE = (log10(Limits(1))-log10(ExLimits(1)))/ExLownum;


for i = 1:(length(x)-1)
    ExLimits(i+1) = 10^(DE+log10(ExLimits(i)));
    
    ExLowChannels(i,1) = ExLimits(i);
    ExLowChannels(i,2) = ExLimits(i+1);
    
    
    
end

Lownum = 12;
x = 1:(Lownum+1);

DE = (log10(MedLimits(1))-log10(Limits(1)))/Lownum;


for i = 1:(length(x)-1)
    Limits(i+1) = 10^(DE+log10(Limits(i)));
    
    LowChannels(i,1) = Limits(i);
    LowChannels(i,2) = Limits(i+1);
    
    
    
end

Mednum = 16;
x_med = 1:(Mednum+1);

DE_med = (log10(High_Limits(1))-log10(MedLimits(1)))/Mednum;


for i = 1:(length(x_med)-1)
    MedLimits(i+1) = 10^(DE_med+log10(MedLimits(i)));
    
    MedChannels(i,1) = MedLimits(i);
    MedChannels(i,2) = MedLimits(i+1);
    
    
    
end

Highnum = 10;
x_high = 1:(Highnum+1);

DE_high = (log10(7)-log10(High_Limits(1)))/Highnum;
High_Limits(1) = 2.5;

for i = 1:(length(x_high)-1)
    High_Limits(i+1) = 10^(DE_high+log10(High_Limits(i)));
    
    HighChannels(i,1) = High_Limits(i);
    HighChannels(i,2) = High_Limits(i+1);
    
    
end



Channels = [ExLowChannels;LowChannels;MedChannels;HighChannels]
NumChannels = length(Channels)

for i = 1:length(Channels)
    
   Resolution(i,1) = 100*(-Channels(i,1)+Channels(i,2))/((Channels(i,1)+Channels(i,2))/2);
    
end

Resolution


%% Write to Text File

fileID = fopen('electron_channels_v1.txt','w');
for i = 1:size(Channels,1)
fprintf(fileID,'%6.3f,%6.3f \n',Channels(i,:));
end
fclose(fileID);
