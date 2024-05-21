%%BowTie.m
% BowTie and Count Rate Analysis for HERT
% Last modified: 4/15/2023

clc
close all
%% Create expontential fits to AE9AP9 data for Omnidirectional https://www.vdl.afrl.af.mil/programs/ae9ap9/architecture.php
CF = 4*pi;
M_energy_beam = [0.5,0.75,1:0.5:7];
%Values from AE9
%25th Percentile
MaxFlux25_AE9_Omni= [8111977.47	1573149.04	570832.259	145269.2	29218.1317	6903.77267	2063.2227	742.263646	266.679309	109.348528	49.3155003	24.100156	13.0634938	7.06531883	3.28402254]./CF;
CenFlux25_AE9Omni= [119336.699	8225.82878	2706.75572	1120.70882	441.600604	172.554702	70.32912	28.4551312	12.1755939	5.72643673	3.11791022	1.9170406	1.23915444	0.804382395	0.466120109]./CF;
ApoFlux25_AE9Omni= [1329043.75	252562.735	54297.3124	12851.4779	2957.10699	756.164159	234.966797	90.8712388	47.1895398	22.0176569	10.0054712	4.57520635	2.71797272	1.77125091	0.887628471]./CF;

%50th Percentile
MaxFlux50_AE9Omni= [19893090.7	5097010.91	2097907.42	624573.54	154628.957	40909.9701	13350.1255	5049.95988	2026.19355	857.250336	375.862092	170.394972	85.9078482	50.603342	26.4107426]./CF;
CenFlux50_AE9Omni= [831159.082	75215.3761	24933.8101	10050.2331	4041.88421	1583.18182	688.277682	298.083931	133.99576	67.5783659	38.6913855	24.7525623	16.1691211	10.7614667	6.38568385]./CF;
ApoFlux50_AE9Omni= [3560799.93	801096.429	209891.941	55885.6218	13987.2765	3724.35075	1213.94357	471.837509	231.23887	106.863573	50.3731833	25.199803	15.7532854	10.6828812	5.7595443]./CF;

%75th Percentile
MaxFlux75_AE9Omni= [40883035.7	12922351.4	6049174.9	2077265.16	590046.286	166606.194	58924.2062	23885.7512	10125.2594	4511.16945	1978.85991	914.656732	473.531305	291.004977	163.06214]./CF;
CenFlux75_AE9Omni= [3846756.52	434054.751	144124.757	56744.2736	23170.7812	9094.12836	4159.27749	1901.0335	888.89288	473.988988	282.537685	186.68609	123.24219	84.0403946	51.0708825]./CF;
ApoFlux75_AE9Omni= [7745852.95	1990748.82	609664.598	178140.465	47641.6501	13094.5472	4431.71299	1730.04776	810.168494	371.511692	180.126638	96.8010134	63.2327947	44.4305583	25.490612]./CF;

%95th Percentile
MaxFlux95_AE9Omni= [91558289.4	41140003.2	21777989.6	8588339.04	2658679.36	821612.364	310389.614	138246.272	62029.605	28853.2118	13085.9523	6466.13666	3657.18038	2383.31583	1430.15343]./CF;
CenFlux95_AE9Omni= [21177994.6	3071319.81	1017724.32	389492.874	161615.176	63576.8611	30758.273	14929.137	7296.80974	4142.65193	2584.25855	1771.65054	1185.81287	834.132484	522.86174]./CF;
ApoFlux95_AE9Omni= [18385065.5	5479128.5	1996072.17	646850.622	186225.504	53014.2791	18706.5562	7341.4144	3269.1016	1485.84951	742.754015	432.57455	297.754015	218.423896	134.799101]./CF;

%Mean
MaxFluxMean_AE9Omni= [29719534	10294550.9	5145356.75	1972260.08	591967.488	180958.825	67634.9545	29785.0655	13309.3245	6179.96504	2801.18165	1385.17797	791.346898	521.590481	314.831095]./CF;
CenFluxMean_AE9Omni= [4590162.32	658963.392	218115.238	83384.8825	34597.477	13610.6814	6596.79335	3214.55575	1577.85667	902.813172	567.118013	391.331833	262.762424	186.186808	117.636928]./CF;
ApoFluxMean_AE9Omni= [5680714.31	1525170.91	504193.498	156031.64	43687.6255	12289.7226	4280.69018	1677.55596	758.781532	345.866514	171.16136	97.4543412	66.3050963	48.2490523	29.4378704]./CF;


%% Create expontential fits to AE9AP9 data for 90 degree pitch angle
%Values from AE9
%25th Percentile
MaxFlux25_AE9_90 = [1299128.18	159571.67	58982.02	15129.9237	3041.56181	718.45857	219.159343	77.6956967	26.2584103	10.3683737	4.73543853	2.43125066	1.31421209	0.694297103	0.328897155];
CenFlux25_AE9_90 = [236773.987	102478.867	42665.4286	12486.5777	2898.92566	683.399265	201.636631	69.2078129	24.9863608	10.0010258	4.63892914	2.28897503	1.19219726	0.648284687	0.312877121];
ApoFlux25_AE9_90 = [127986.22	24967.9964	5444.69973	1319.99226	312.878706	81.8187931	25.3755998	9.72523144	4.90561484	2.21073182	0.957248098	0.416308498	0.239633178	0.152650546	0.07342178];

%50th Percentile
MaxFlux50_AE9_90 = [3043556.52	506402.542	211356.954	62748.4502	15499.6257	4167.08963	1362.59182	508.053908	196.611127	82.051654	37.154501	17.4685642	8.97551526	5.28536699	2.81916956];
CenFlux50_AE9_90 = [876087.841	397179.856	181014.549	58780.2084	15120.1328	3942.52241	1290.06064	487.952961	187.386107	76.4888646	33.3951236	16.0558237	8.53704862	5.05846258	2.61591484];
ApoFlux50_AE9_90 = [341776.994	78702.5833	20845.5994	5693.35574	1472.38512	399.944599	129.090197	49.2612793	23.5157668	10.5373774	4.81524394	2.36038762	1.4794959	1.00421865	0.532351405];

%75th Percentile
MaxFlux75_AE9_90 = [5954877.83	1258557.25	590380.762	202629.404	58159.668	16656.5776	5906.61956	2397.14709	991.163556	436.275356	196.258873	94.3908837	49.8709683	30.3881283	16.5763025];
CenFlux75_AE9_90 =[2469165.81	1160706.51	567296.774	199717.547	55612.4246	15693.0662	5571.19522	2275.12434	917.33847	380.385147	158.469032	74.7436904	40.4205575	25.6070785	13.9684838];
ApoFlux75_AE9_90 = [741554.833	194643.441	60101.779	18034.1354	4996.24183	1398.42223	465.736857	177.111453	80.9439822	36.1004976	17.2120895	9.2718457	6.21568293	4.43625712	2.53942722];

%95th Percentile
MaxFlux95_AE9_90 = [12559884.1	3917196.38	2094895.39	826481.191	256040.313	79982.653	30444.5749	13671.416	6070.77165	2795.68906	1292.31374	654.009611	366.498171	235.290299	139.450819];
CenFlux95_AE9_90 = [7854048.94	3842396.1	2026654.85	779563.225	236687.829	72888.9963	28329.6745	12599.2993	5363.30435	2263.84066	895.985629	414.265436	228.444014	155.761882	90.0617119];
ApoFlux95_AE9_90 = [1753727.44	532323.177	194919.825	64926.6752	19416.2462	5617.60949	1937.2571	734.009612	319.68837	141.834443	70.8966219	42.4119355	30.6420444	23.1319331	14.4237781];

%Mean
MaxFluxMean_AE9_90 = [2015489.32	966624.012	492501.967	183035.525	54030.6479	16263.9459	6207.6321	2727.12992	1154.56719	486.691215	193.553499	89.7218939	49.3789625	33.4564409	19.2888158];
CenFluxMean_AE9_90 = [2015489.32	966624.012	492501.967	183035.525	54030.6479	16263.9459	6207.6321	2727.12992	1154.56719	486.691215	193.553499	89.7218939	49.3789625	33.4564409	19.2888158];
ApoFluxMean_AE9_90 = [543134.812	148679.823	49426.7817	15701.1899	4559.39628	1304.11571	444.805326	168.743676	74.6672731	33.1902862	16.3350572	9.49303049	6.74878269	5.04578025	3.11202649];


%% AE9 Flux Model using Omni/4pi
MaxFlux_25_Omni = zeros(size(M_energy_beam));
CenFlux_25_Omni = zeros(size(M_energy_beam));
ApoFlux_25_Omni = zeros(size(M_energy_beam));
MaxFlux_50_Omni = zeros(size(M_energy_beam));
CenFlux_50_Omni = zeros(size(M_energy_beam));
ApoFlux_50_Omni = zeros(size(M_energy_beam));
MaxFlux_75_Omni = zeros(size(M_energy_beam));
CenFlux_75_Omni = zeros(size(M_energy_beam));
ApoFlux_75_Omni = zeros(size(M_energy_beam));
MaxFlux_95_Omni = zeros(size(M_energy_beam));
CenFlux_95_Omni = zeros(size(M_energy_beam));
ApoFlux_95_Omni = zeros(size(M_energy_beam));
MaxFlux_Mean_Omni = zeros(size(M_energy_beam));
CenFlux_Mean_Omni = zeros(size(M_energy_beam));
ApoFlux_Mean_Omni = zeros(size(M_energy_beam));


for i = 1:length(M_energy_beam(1,:))
    
    [MaxFlux_25_Omni(i),CenFlux_25_Omni(i),ApoFlux_25_Omni(i),MaxFlux_50_Omni(i),CenFlux_50_Omni(i),ApoFlux_50_Omni(i),MaxFlux_75_Omni(i),CenFlux_75_Omni(i),ApoFlux_75_Omni(i),MaxFlux_95_Omni(i),CenFlux_95_Omni(i),ApoFlux_95_Omni(i),MaxFlux_Mean_Omni(i),CenFlux_Mean_Omni(i),ApoFlux_Mean_Omni(i)] = FluxRateOmni(M_energy_beam(1,i));

end

figure
semilogy(M_energy_beam(1,:),MaxFlux_25_Omni,M_energy_beam(1,:),CenFlux_25_Omni,M_energy_beam(1,:),ApoFlux_25_Omni,'Linewidth',2)
Title = append('OMNI/4*pi 25th Percentile Flux(E) Average',' E Range: ',num2str(min(M_energy_beam(1,:))),' to ',num2str(max(M_energy_beam(1,:))),' MeV');
title(Title)
ylim ([0 10^8])
ylabel('Flux (#/cm^2/s/sr/Mev)')
xlabel('Incident Energy (MeV)')
legend('Max','Cen','Apo','Max_Test','Cen_Test','Apo_Test')
saveas(gcf,append('Omni 25th Flux','.jpg'))

figure
semilogy(M_energy_beam(1,:),MaxFlux_50_Omni,M_energy_beam(1,:),CenFlux_50_Omni,M_energy_beam(1,:),ApoFlux_50_Omni,'Linewidth',2)
Title = append('OMNI/4*pi 50th Percentile Flux(E) Average',' E Range: ',num2str(min(M_energy_beam(1,:))),' to ',num2str(max(M_energy_beam(1,:))),' MeV');
title(Title)
ylim ([0 10^8])
ylabel('Flux (#/cm^2/s/sr/Mev)')
xlabel('Incident Energy (MeV)')
legend('Max','Cen','Apo')
saveas(gcf,append('Omni 50th Flux','.jpg'))

figure
semilogy(M_energy_beam(1,:),MaxFlux_75_Omni,M_energy_beam(1,:),CenFlux_75_Omni,M_energy_beam(1,:),ApoFlux_75_Omni,'Linewidth',2)
Title = append('OMNI/4*pi 75th Percentile Flux(E) Average',' E Range: ',num2str(min(M_energy_beam(1,:))),' to ',num2str(max(M_energy_beam(1,:))),' MeV');
title(Title)
ylim ([0 10^8])
ylabel('Flux (#/cm^2/s/sr/Mev)')
xlabel('Incident Energy (MeV)')
legend('Max','Cen','Apo')
saveas(gcf,append('Omni 75th Flux','.jpg'))

figure
semilogy(M_energy_beam(1,:),MaxFlux_95_Omni,M_energy_beam(1,:),CenFlux_95_Omni,M_energy_beam(1,:),ApoFlux_95_Omni,'Linewidth',2)
Title = append('OMNI/4*pi 95th Percentile Flux(E) Average',' E Range: ',num2str(min(M_energy_beam(1,:))),' to ',num2str(max(M_energy_beam(1,:))),' MeV');
title(Title)
ylim ([0 10^8])
ylabel('Flux (#/cm^2/s/sr/Mev)')
xlabel('Incident Energy (MeV)')
legend('Max','Cen','Apo')
saveas(gcf,append('Omni 95th Flux','.jpg'))

figure
semilogy(M_energy_beam(1,:),MaxFlux_Mean_Omni,M_energy_beam(1,:),CenFlux_Mean_Omni,M_energy_beam(1,:),ApoFlux_Mean_Omni,'Linewidth',2)
Title = append('OMNI/4*pi Mean Flux(E) Average',' E Range: ',num2str(min(M_energy_beam(1,:))),' to ',num2str(max(M_energy_beam(1,:))),' MeV');
title(Title)
ylim ([0 10^8])
ylabel('Flux (#/cm^2/s/sr/Mev)')
xlabel('Incident Energy (MeV)')
legend('Max','Cen','Apo')
saveas(gcf,append('Omni Mean Flux','.jpg'))

%% AE9 Flux Model using 90 deg pitch angle
MaxFlux25_90 = zeros(size(M_energy_beam));
CenFlux25_90 = zeros(size(M_energy_beam));
ApoFlux25_90 = zeros(size(M_energy_beam));
MaxFlux50_90 = zeros(size(M_energy_beam));
CenFlux50_90 = zeros(size(M_energy_beam));
ApoFlux50_90 = zeros(size(M_energy_beam));
MaxFlux75_90 = zeros(size(M_energy_beam));
CenFlux75_90 = zeros(size(M_energy_beam));
ApoFlux75_90 = zeros(size(M_energy_beam));
MaxFlux95_90 = zeros(size(M_energy_beam));
CenFlux95_90 = zeros(size(M_energy_beam));
ApoFlux95_90 = zeros(size(M_energy_beam));
MaxFluxMean_90 = zeros(size(M_energy_beam));
CenFluxMean_90 = zeros(size(M_energy_beam));
ApoFluxMean_90 = zeros(size(M_energy_beam));


for i = 1:length(M_energy_beam(1,:))
    
    [MaxFlux25_90(i),CenFlux25_90(i),ApoFlux25_90(i),MaxFlux50_90(i),CenFlux50_90(i),ApoFlux50_90(i),MaxFlux75_90(i),CenFlux75_90(i),ApoFlux75_90(i),MaxFlux95_90(i),CenFlux95_90(i),ApoFlux95_90(i),MaxFluxMean_90(i),CenFluxMean_90(i),ApoFluxMean_90(i)] = FluxRate90(M_energy_beam(1,i));
    
end

figure
semilogy(M_energy_beam(1,:),CenFlux25_90,M_energy_beam(1,:),ApoFlux25_90,'Linewidth',2)
Title = append('90 deg pitch 25th Percentile Flux(E) Average',' E Range: ',num2str(min(M_energy_beam(1,:))),' to ',num2str(max(M_energy_beam(1,:))),' MeV');
title(Title)
ylim ([0 10^8])
ylabel('Flux (#/cm^2/s/sr/Mev)')
xlabel('Incident Energy (MeV)')
legend('Cen','Apo')
saveas(gcf,append('90 deg Pitch Angle 25th Flux','.jpg'))

figure
semilogy(M_energy_beam(1,:),CenFlux50_90,M_energy_beam(1,:),ApoFlux50_90,'Linewidth',2)
Title = append('90 deg pitch 50th Percentile Flux(E) Average',' E Range: ',num2str(min(M_energy_beam(1,:))),' to ',num2str(max(M_energy_beam(1,:))),' MeV');
title(Title)
ylim ([0 10^8])
ylabel('Flux (#/cm^2/s/sr/Mev)')
xlabel('Incident Energy (MeV)')
legend('Cen','Apo')
saveas(gcf,append('90 deg Pitch Angle 50th Flux','.jpg'))

figure
semilogy(M_energy_beam(1,:),CenFlux75_90,M_energy_beam(1,:),ApoFlux75_90,'Linewidth',2)
Title = append('90 deg pitch 75th Percentile Flux(E) Average',' E Range: ',num2str(min(M_energy_beam(1,:))),' to ',num2str(max(M_energy_beam(1,:))),' MeV');
title(Title)
ylim ([0 10^8])
ylabel('Flux (#/cm^2/s/sr/Mev)')
xlabel('Incident Energy (MeV)')
legend('Cen','Apo')
saveas(gcf,append('90 deg Pitch Angle 75th Flux','.jpg'))

figure
semilogy(M_energy_beam(1,:),CenFlux95_90,M_energy_beam(1,:),ApoFlux95_90,'Linewidth',2)
Title = append('90 deg pitch 95th Percentile Flux(E) Average',' E Range: ',num2str(min(M_energy_beam(1,:))),' to ',num2str(max(M_energy_beam(1,:))),' MeV');
title(Title)
ylim ([0 10^8])
ylabel('Flux (#/cm^2/s/sr/Mev)')
xlabel('Incident Energy (MeV)')
legend('Cen','Apo')
saveas(gcf,append('90 deg Pitch Angle 95th Flux','.jpg'))

figure
semilogy(M_energy_beam(1,:),CenFluxMean_90,M_energy_beam(1,:),ApoFluxMean_90,'Linewidth',2)
Title = append('90 deg pitch Mean Flux(E) Average',' E Range: ',num2str(min(M_energy_beam(1,:))),' to ',num2str(max(M_energy_beam(1,:))),' MeV');
title(Title)
ylim ([0 10^8])
ylabel('Flux (#/cm^2/s/sr/Mev)')
xlabel('Incident Energy (MeV)')
legend('Cen','Apo')
saveas(gcf,append('90 deg Pitch Angle Mean Flux','.jpg'))


%% Bow Tie Analysis-Selesnick/Blake
%Changes directory for jpg of graphs
cd '..\Bow Tie'

%Sets Ei and Range of Eo
%Ei is the incident energy from the GEANT4 results
Ei = M_energy_beam(1,:);
%Eo set by user.
Eo = 0.2:0.2:2.0;


%Sets up color vectors for plotting the different Eo curves
Eo_color = magma(length(Eo)+1);

%Preallocates all variables prior to For Loops
J_e = zeros(length(M_energy_beam),length(Eo));

G_int = zeros(length(M_energy_beam),length(Eo));
G_term = zeros(length(M_energy_beam),length(Eo));
G_E_eff = zeros(length(M_energy_beam),length(Eo));

xi = zeros(sum(1:length(Eo)-1),length(energy_channels));
yi = xi;

E_eff = zeros(1,length(energy_channels));
G_eff_dE= E_eff;
BowTieLegend = strings([1,length(Eo)]);

%Count_Rate =  zeros(length(energy_channels),length(Eo));

BinWidth = zeros(1, length(energy_channels));
Geff = zeros(1, length(energy_channels));

%Creates J(e) and creates String Array for Plot Legends
for i = 1:length(Eo)
    J_e(:,i) = exp(-Ei/Eo(i));
    BowTieLegend(i) = num2str(Eo(i));
end

%Adds 'Average Intersection Point' for plotting
BowTieLegend = [BowTieLegend,'Intersection Point','Average Intersection Point'];

%Creates J(e)^-1
J_e_inv = 1./J_e;

%Preallocation
BowCount_Rate_Max25_Omni = zeros(length(energy_channels),1);
BowCount_Rate_Cen25_Omni = zeros(length(energy_channels),1);
BowCount_Rate_Apo25_Omni = zeros(length(energy_channels),1);

%50thPercentile
BowCount_Rate_Max50_Omni = zeros(length(energy_channels),1);
BowCount_Rate_Cen50_Omni = zeros(length(energy_channels),1);
BowCount_Rate_Apo50_Omni = zeros(length(energy_channels),1);

%75thPercentile
BowCount_Rate_Max75_Omni =zeros(length(energy_channels),1);
BowCount_Rate_Cen75_Omni = zeros(length(energy_channels),1);
BowCount_Rate_Apo75_Omni = zeros(length(energy_channels),1);

%95th Percentile
BowCount_Rate_Max95_Omni = zeros(length(energy_channels),1);
BowCount_Rate_Cen95_Omni = zeros(length(energy_channels),1);
BowCount_Rate_Apo95_Omni = zeros(length(energy_channels),1);

%Average
BowCount_Rate_MaxMean_Omni = zeros(length(energy_channels),1);
BowCount_Rate_CenMean_Omni = zeros(length(energy_channels),1);
BowCount_Rate_ApoMean_Omni = zeros(length(energy_channels),1);

BowCount_Rate_Max25_90 = zeros(length(energy_channels),1);
BowCount_Rate_Cen25_90 = zeros(length(energy_channels),1);
BowCount_Rate_Apo25_90 = zeros(length(energy_channels),1);

%50thPercentile
BowCount_Rate_Max50_90 = zeros(length(energy_channels),1);
BowCount_Rate_Cen50_90 = zeros(length(energy_channels),1);
BowCount_Rate_Apo50_90 = zeros(length(energy_channels),1);

%75thPercentile
BowCount_Rate_Max75_90 = zeros(length(energy_channels),1);
BowCount_Rate_Cen75_90 = zeros(length(energy_channels),1);
BowCount_Rate_Apo75_90 = zeros(length(energy_channels),1);

%95th Percentile
BowCount_Rate_Max95_90 = zeros(length(energy_channels),1);
BowCount_Rate_Cen95_90 = zeros(length(energy_channels),1);
BowCount_Rate_Apo95_90 = zeros(length(energy_channels),1);

%Average
BowCount_Rate_MaxMean_90 = zeros(length(energy_channels),1);
BowCount_Rate_CenMean_90 = zeros(length(energy_channels),1);
BowCount_Rate_ApoMean_90 = zeros(length(energy_channels),1);


%For Loop for calculating a line for each Eo and finding the average intersection point
for c=1:length(energy_channels)
    %Prints energy channel to show progress
    fprintf('Energy Channel-# %.2f \n',c)
    
    %Calculates line for each Eo
    for i = 1:length(Eo)
        G_int(:,i,c) = geo_EC(:,c).*J_e(:,i);
        G_term(i,c) = trapz(Ei,G_int(:,i,c));
        G_E_eff(:,i,c)= G_term(i,c)*J_e_inv(:,i);
        
    end
    
    %Find Intersections of Eo Lines
    num = 0;
    for j = 1:(length(Eo)-1)
        
        for k = 1:(length(Eo)-j)
            [xi(num+k,c),yi(num+k,c)] = polyxpoly(M_energy_beam(1,:),G_E_eff(:,j,c),M_energy_beam(1,:),G_E_eff(:,j+k,c),'unique');
            
        end
        num = num + k;
    end
    
    %Finds average intersection point
    E_eff(c) = mean(xi(:,c));
    G_eff_dE(c) = mean(yi(:,c));
    
    %Calculate Bin Characteristics
    Geff(c) = G_eff_dE(c)/fwhm(c);
    BinWidth(c) = fwhm(c);
    
end


%% Plots for Config
textsize = 28;
%Plots Graph for each energy channel
ymax = round(max(G_eff_dE)*1.1,3);
for c=1:height(energy_channels)
    
    f = figure;
    f.Position = [100 100 1000 720];
    hold on
    for i = 1:length(Eo)
        plot(M_energy_beam(1,:),G_E_eff(:,i,c),'Color',Eo_color(i,:),'DisplayName',BowTieLegend(1,i),'Linewidth',2);
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
    %title(append('Bow Tie Analysis Energy Channel ',num2str(c),' Deposited Energy Range: ',num2str(energy_channels(c,1)),' to ',num2str(energy_channels(c,2)),' MeV'))
    ylabel('G_{eff} * \Delta E ','FontSize',textsize)
    xlabel('Nominal Energy (MeV)','FontSize',textsize)
    hold off
    
    effsave = append(datetime("today"),'Bow Tie Energy Channel ',num2str(c),' Eo ',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
    saveas(gcf,effsave)
    
end

%Plots graph of all the average intersection points of each energy channel
f = figure;
f.Position = [100 100 1200 720];
hold on
for c=1:length(energy_channels)
    plot(E_eff(c),G_eff_dE(c),'o','Color',Effplotcolor(c,:))
end
legend(EngLegend,'Location', 'southoutside','NumColumns',8)
ylim([0 ymax])
yticks((0:0.005:ymax))

xlim([1 7])
xticks((1:1:7))

title(append('Bow Tie Analysis Configuration All Energy Channels ',' Eo Range: ',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV'))
ylabel('G_{eff} * \Delta E ')
xlabel('Nominal Effective Energy (MeV)')
hold off

effsave = append(datetime("today"),' Bow Tie  All Energy Channels',' Eo',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
saveas(gcf,effsave)

f = figure;
f.Position = [100 100 1600 720];
hold on
%Plots each energy channel FWHM value in the same color as
%the Eff. Curve
for c = 1:length(energy_channels)
    bar(E_eff(c),G_eff_dE(c)/fwhm(c),fwhm(c),'EdgeColor','k','FaceColor',Effplotcolor(c,:))
    
end
legend(EngLegend,'Location', 'southoutside','NumColumns',8)
%title(append('Energy Channel Bins-',' Eo Range: ',num2str(Cen(Eo)),' to ',num2str(max(Eo)),' MeV'),'FontSize',20)
ylabel('Effective Geometric Factor','FontSize',28)
xlabel('Nominal Effective Energy (MeV)','FontSize',28)
hold off

effsave = append(datetime("today"),' Energy Channel Bins',' Eo ',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
saveas(gcf,effsave)


% Count Rate Plot
%{
energy_channel_list = 1:1:length(energy_channels);
%count_rate_colors = magma(length(FluxRange));
count_rate_colors = magma(1);

f = figure;
f.Position = [100 100 1200 720];
%hold on
%Plots each energy channel FWHM value in the same color as
%the Eff. Curve

semilogy(energy_channel_list,BowCount_Rate_Max50_Omni,'o',energy_channel_list,BowCount_Rate_Max75_Omni,'o',energy_channel_list,BowCount_Rate_Max95_Omni,'o')

xlim([0 length(energy_channels)+1])
xticks((0:1:length(energy_channels)+1))

title(append('Count Rates ',' Eo Range: ',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV'))
ylabel('Count Rates')
xlabel('Energy Channels')
%hold off

effsave = append(datetime("today"),' Count Rate',' Eo ',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
saveas(gcf,effsave)


%Returns to main directory
cd ..
cd ..
Energy_Resolution= 100*BinWidth./E_eff;
%Exports Bin Info to Excel File
exportname = append(datetime("today"),' Center- Bin Width Table ',' Eo ',num2str(min(Eo)),' to ',num2str(max(Eo)),'_',addin,num2str(length(energy_channels)));
BuildBinTable(energy_channels,E_eff,BinWidth,Geff,Energy_Resolution,IntEnergyChannel_Cen25_90,IntEnergyChannel_Cen50_90,IntEnergyChannel_Cen75_90,exportname)

exportname = append(datetime("today"),' Apo- Bin Width Table ',' Eo ',num2str(min(Eo)),' to ',num2str(max(Eo)),'_',addin,num2str(length(energy_channels)));
BuildBinTable(energy_channels,E_eff,BinWidth_,Geff_,Energy_Resolution,IntEnergyChannel_Apo25_90,IntEnergyChannel_Apo50_90,IntEnergyChannel_Apo75_90,exportname)


exportname = append(datetime("today"),' Configuration- Count Rate Table ',' Eo ',num2str(min(Eo)),' to ',num2str(max(Eo)),'_',addin,num2str(length(energy_channels)));
%BuildCountTable(energy_channels,Count_Rate_,exportname)
%}
%{
figure
plot(Energy_Resolution,'xb')
title('Energy Resolution per Energy Channel')
xlabel('Energy Channel Number')
ylabel('Energy Resolution (%)')
%}

%% Energy Resolution Plot
%{
Effplotcolor = plasma(length(energy_channel_list));
f = figure;
f.Position = [100 100 1600 720];
textsize = 28;

%Plots each energy channel FWHM value in the same color as
%the Eff. Curve

hold on
for c = 1:length(energy_channels)
%     plot(E_eff(c),Energy_Resolution(c),'o','MarkerSize',8,...
%         'MarkerEdgeColor',Effplotcolor(c,:),...
%         'MarkerFaceColor',Effplotcolor(c,:))
    
    plot(E_eff(c),Energy_Resolution(c),'o','MarkerSize',8,...
        'MarkerEdgeColor','b',...
        'MarkerFaceColor','b')
end
set(gca,'FontSize',textsize)
%title('HERT Energy Resolution','FontSize',textsize)
ylabel('Spectral Resolution dE/E(%)','FontSize',textsize)
xlabel('Nominal Energy (MeV)','FontSize',textsize)
ylim([0,40])
plot([min(min(M_energy_beam)),max(max(M_energy_beam))],[12,12],'k--','LineWidth',2)%,'DisplayName','Energy Resolution Requirement')

%legend([EngLegend,'Energy Resolution Requirement'],'Location', 'southoutside','NumColumns',8)
legend([strings(1,length(energy_channels)),'Energy Resolution Requirement'],'Location', 'northeast','NumColumns',8)

hold off
effsave = append(datetime("today"),' Energy Resolution',' Eo ',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
saveas(gcf,effsave)

% %% Intersection Plot
% hold on
% num = 0;
% for j = 1:(length(Eo)-1)
%
%     for k = 1:(length(Eo)-j)
%         plot(xi(num+k,16),yi(num+k,16),'o','Color',Eo_color(j,:))
%
%     end
%     num = num + k;
% end
% plot(E_eff(16),G_eff_dE(16),'or')
% hold off
%}

%% Count Rate Comparision
%{
cd 'Bow Tie'

f = figure;
f.Position = [100 100 1200 720];

%Plots each energy channel FWHM value in the same color as the Eff. Curve
semilogy(energy_channel_list,BowCount_Rate_Max50_Omni,'o')

xlim([0 length(energy_channels)+1])
xticks((0:1:length(energy_channels)+1))

legend({'Count Rate',' Count Rate'},'Location', 'southoutside')
title(append('Count Rates Comparision 50th ',' Eo Range: ',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV'))
ylabel('Count Rates')
xlabel('Energy Channels')

effsave = append(datetime("today"),' Count Rate Comparision 50th',' Eo ',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
saveas(gcf,effsave)

f = figure;
f.Position = [100 100 1200 720];
semilogy(energy_channel_list,BowCount_Rate_Max75_Omni,'o',energy_channel_list,BowCount_Rate_InnerMax75_Omni,'og')


xlim([0 length(energy_channels)+1])
xticks((0:1:length(energy_channels)+1))

legend({'Count Rate',' Count Rate'},'Location', 'southoutside')
title(append('Count Rates Comparision 75th',' Eo Range: ',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV'))
ylabel('Count Rates')
xlabel('Energy Channels')


effsave = append(datetime("today"),' Count Rate Comparision 75th',' Eo ',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
saveas(gcf,effsave)
 
f = figure;
f.Position = [100 100 1200 720];
semilogy(energy_channel_list,BowCount_Rate_Max95_Omni,'o',energy_channel_list,BowCount_Rate_InnerMax95_Omni,'og')
 
xlim([0 length(energy_channels)+1])
xticks((0:1:length(energy_channels)+1))
 
legend({'Count Rate',' Count Rate'},'Location', 'southoutside')
title(append('Count Rates Comparision 95th',' Eo Range: ',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV'))
ylabel('Count Rates')
xlabel('Energy Channels')
 
effsave = append(datetime("today"),' Count Rate Comparision 95th',' Eo ',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
saveas(gcf,effsave)

cd ..
%}

%% Count Rate Comparision-per Energy
%{
%All Channels
    channel_select = 1:length(energy_channels);
    %channel_select = [2,10,20,30,35];
                
f = figure;
f.Position = [100 100 1200 720];
color_iter = 1;

%Plots each energy channel FWHM value in the same color as the Eff. Curve
EC_plot_color = plasma(length(channel_select));

hold on
for i = 1:width(M_energy_beam)
    if max(i == channel_select)
        plot(E_eff(i),BowCount_Rate_Apo50_90(i,1),'o','MarkerEdgeColor','b',...
        'MarkerFaceColor','b');
    %EC_plot_color(i,:)
        EngLegend_EC(i) = append(sprintf('#%.0f: ',i),EngLegend(i));
        color_iter = color_iter+1;
        
    else
        plot(E_eff(i),BowCount_Rate_Apo50_90(i),'Color',[0.75, 0.75, 0.75],'LineWidth',line_width);
    end
    
end

for i = 1:width(M_energy_beam)
    if max(i == channel_select)
        plot(E_eff(i),BowCount_Rate_Cen50_90(i,1),'o','MarkerEdgeColor','g',...
        'MarkerFaceColor','g');
    %EC_plot_color(i,:)
        EngLegend_EC(i) = append(sprintf('#%.0f: ',i),EngLegend(i));
        color_iter = color_iter+1;
        
    else
        plot(E_eff(i),BowCount_Rate_Cen50_90(i),'Color',[0.75, 0.75, 0.75],'LineWidth',line_width);
    end
    
end

set(gca, 'YScale', 'log')
set(gca, 'Fontsize', 20)
hold off

%xlim([0 8])
%xticks((0:0.5:8))

%legend(EngLegend_EC,'Location', 'southoutside','NumColumns',6)
%title(append('Count Rates Comparision 50th ',' Eo Range: ',num2str(Cen(Eo)),' to ',num2str(max(Eo)),' MeV'))
ylabel('Count Rates (#/sec)')
xlabel('Nominal Energy (MeV)')

effsave = append(datetime("today"),' Count Rate Comparision 50th',' Eo ',num2str(min(Eo)),' to ',num2str(max(Eo)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
saveas(gcf,effsave)

% f = figure;
% f.Position = [100 100 1200 720];
% semilogy(energy_channel_list,BowCount_Rate_Max75_Omni,'o',energy_channel_list,BowCount_Rate_InnerMax75_Omni,'og')
% 
% 
% xlim([0 length(energy_channels)+1])
% xticks((0:1:length(energy_channels)+1))
% 
% legend({'Count Rate',' Count Rate'},'Location', 'southoutside')
% title(append('Count Rates Comparision 75th',' Eo Range: ',num2str(Cen(Eo)),' to ',num2str(max(Eo)),' MeV'))
% ylabel('Count Rates')
% xlabel('Energy Channels')
% 
% 
% effsave = append(datetime("today"),' Count Rate Comparision 75th',' Eo ',num2str(Cen(Eo)),' to ',num2str(max(Eo)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
% saveas(gcf,effsave)
% 
% f = figure;
% f.Position = [100 100 1200 720];
% semilogy(energy_channel_list,BowCount_Rate_Max95_Omni,'o',energy_channel_list,BowCount_Rate_InnerMax95_Omni,'og')
% 
% 
% xlim([0 length(energy_channels)+1])
% xticks((0:1:length(energy_channels)+1))
% 
% legend({'Count Rate',' Count Rate'},'Location', 'southoutside')
% title(append('Count Rates Comparision 95th',' Eo Range: ',num2str(Cen(Eo)),' to ',num2str(max(Eo)),' MeV'))
% ylabel('Count Rates')
% xlabel('Energy Channels')
% 
% 
% effsave = append(datetime("today"),' Count Rate Comparision 95th',' Eo ',num2str(Cen(Eo)),' to ',num2str(max(Eo)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
% saveas(gcf,effsave)
%}


%% Energy per Second for Each Detector
%{
detector_energy_S = FluxRangeEL(:,1:width(detector_energy)).*detector_energy;

inner_detector_energy_S = FluxRangeEL(:,1:width(inner_detector_energy)).*inner_detector_energy;

total_rate = sum(detector_energy_S,1)'
Inner_total_rate = sum(inner_detector_energy_S,1)'

f = figure;
f.Position = [100 100 1200 720];

semilogy(M_energy_beam(1,:),detector_energy_S(:,1),M_energy_beam(1,:),detector_energy_S(:,2),M_energy_beam(1,:),detector_energy_S(:,3),M_energy_beam(1,:),detector_energy_S(:,4),M_energy_beam(1,:),detector_energy_S(:,5),M_energy_beam(1,:),detector_energy_S(:,6),M_energy_beam(1,:),detector_energy_S(:,7),M_energy_beam(1,:),detector_energy_S(:,8),M_energy_beam(1,:),detector_energy_S(:,9))
title('Config: Energy Deposited per Second (MeV/s) vs. Incident Energy (MeV)')
ylabel('Energy Deposited per Second (MeV/s)')
xlabel('Incident Energy (MeV)')

effsave = append(datetime("today"),' Energy Rate-',addin,num2str(length(energy_channels)),'.jpg');
saveas(gcf,effsave)

f = figure;
f.Position = [100 100 1200 720];

semilogy(M_energy_beam(1,:),inner_detector_energy_S(:,1),M_energy_beam(1,:),inner_detector_energy_S(:,2),M_energy_beam(1,:),inner_detector_energy_S(:,3),M_energy_beam(1,:),inner_detector_energy_S(:,4),M_energy_beam(1,:),inner_detector_energy_S(:,5),M_energy_beam(1,:),inner_detector_energy_S(:,6),M_energy_beam(1,:),inner_detector_energy_S(:,7),M_energy_beam(1,:),inner_detector_energy_S(:,8),M_energy_beam(1,:),inner_detector_energy_S(:,9),M_energy_beam(1,:),inner_detector_energy_S(:,10),M_energy_beam(1,:),inner_detector_energy_S(:,11),M_energy_beam(1,:),inner_detector_energy_S(:,12),M_energy_beam(1,:),inner_detector_energy_S(:,13),M_energy_beam(1,:),inner_detector_energy_S(:,14),M_energy_beam(1,:),inner_detector_energy_S(:,15),M_energy_beam(1,:),inner_detector_energy_S(:,16),M_energy_beam(1,:),inner_detector_energy_S(:,17))
title('Config: Energy Deposited per Second (MeV/s) vs. Incident Energy (MeV)')
ylabel('Energy Deposited per Second (MeV/s)')
xlabel('Incident Energy (MeV)')

effsave = append(datetime("today"),' Energy Rate ',addin,num2str(length(energy_channels)),'.jpg');
saveas(gcf,effsave)
%}

%% Back Pen Calculations Electron
x_inter = [0.5,0.75,1:0.5:10];
ApoFlux50_AE9_90 = [341776.9940	78702.5833	20845.5994	5693.3557	1472.3851	399.9446	129.0902	49.2613	23.5158	10.5374	4.8152	2.3604	1.4795	1.0042	0.5324	0.4305	0.3286	0.2267	0.1948	0.1629	0.1310];
MaxFlux95_AE9_90 = [12559884.1	3917196.38	2094895.39	826481.191	256040.313	79982.653	30444.5749	13671.416	6070.77165	2795.68906	1292.31374	654.009611	366.498171	235.290299	139.450819	119.126969	98.8031184	78.9011666	71.2501658	63.5991651	56.0650802];

function [MaxFlux25_90,CenFlux25_90,ApoFlux25_90,MaxFlux50_90,CenFlux50_90,ApoFlux50_90,MaxFlux75_90,CenFlux75_90,ApoFlux75_90,MaxFlux95_90,CenFlux95_90,ApoFlux95_90,MaxFluxMean_90,CenFluxMean_90,ApoFluxMean_90] = FluxRate90(r)
%FluxRate Summary of this function goes here
%   Inputs:
%     M_energy_beam - energy at flux
%     Outputs:
%     Flux50:Expected average flux at energy M_energy_beam
%     Flux75:Expected 75th percentile flux at energy M_energy_beam
%     Flux95:Expected 95th percentile flux at energy M_energy_beam

M_energy_beam = [0.5,0.75,1:0.5:8];

%% 90 Degrees
%Values from AE9
%25th Percentile
MaxFlux25_AE9_90 = [1299128.18	159571.67	58982.02	15129.9237	3041.56181	718.45857	219.159343	77.6956967	26.2584103	10.3683737	4.73543853	2.43125066	1.31421209	0.694297103	0.328897155 0.265913281	0.202929408];
CenFlux25_AE9_90 = [236773.987	102478.867	42665.4286	12486.5777	2898.92566	683.399265	201.636631	69.2078129	24.9863608	10.0010258	4.63892914	2.28897503	1.19219726	0.648284687	0.312877121	0.253244987	0.193612852];
ApoFlux25_AE9_90 = [127986.22	24967.9964	5444.69973	1319.99226	312.878706	81.8187931	25.3755998	9.72523144	4.90561484	2.21073182	0.957248098	0.416308498	0.239633178	0.152650546	0.07342178 0.057816453 0.042211125];

%50th Percentile
MaxFlux50_AE9_90 = [3043556.52	506402.542	211356.954	62748.4502	15499.6257	4167.08963	1362.59182	508.053908	196.611127	82.051654	37.154501	17.4685642	8.97551526	5.28536699	2.81916956 2.33	1.84];
CenFlux50_AE9_90 = [876087.841	397179.856	181014.549	58780.2084	15120.1328	3942.52241	1290.06064	487.952961	187.386107	76.4888646	33.3951236	16.0558237	8.53704862	5.05846258	2.61591484 2.14087987	1.66584491];
ApoFlux50_AE9_90 = [341776.994	78702.5833	20845.5994	5693.35574	1472.38512	399.944599	129.090197	49.2612793	23.5157668	10.5373774	4.81524394	2.36038762	1.4794959	1.00421865	0.532351405 0.430483052	0.328614699];

%75th Percentile
MaxFlux75_AE9_90 = [5954877.83	1258557.25	590380.762	202629.404	58159.668	16656.5776	5906.61956	2397.14709	991.163556	436.275356	196.258873	94.3908837	49.8709683	30.3881283	16.5763025 13.7416615 10.9904901];
CenFlux75_AE9_90 = [2469165.81	1160706.51	567296.774	199717.547	55612.4246	15693.0662	5571.19522	2275.12434	917.33847	380.385147	158.469032	74.7436904	40.4205575	25.6070785	13.9684838	11.5349933	9.10150277];
ApoFlux75_AE9_90 = [741554.833	194643.441	60101.779	18034.1354	4996.24183	1398.42223	465.736857	177.111453	80.9439822	36.1004976	17.2120895	9.2718457	6.21568293	4.43625712	2.53942722 2.10258632 1.66562056];

%95th Percentile
MaxFlux95_AE9_90 = [12559884.1	3917196.38	2094895.39	826481.191	256040.313	79982.653	30444.5749	13671.416	6070.77165	2795.68906	1292.31374	654.009611	366.498171	235.290299	139.450819 119.126969 98.8031184];
CenFlux95_AE9_90 = [7854048.94	3842396.1	2026654.85	779563.225	236687.829	72888.9963	28329.6745	12599.2993	5363.30435	2263.84066	895.985629	414.265436	228.444014	155.761882	90.0617119	75.1379942	60.2142765];
ApoFlux95_AE9_90 = [1753727.44	532323.177	194919.825	64926.6752	19416.2462	5617.60949	1937.2571	734.009612	319.68837	141.834443	70.8966219	42.4119355	30.6420444	23.1319331	14.4237781 12.2959433 10.1681086];

%Mean
MaxFluxMean_AE9_90 = [4284294.4	987847.371	497075.967	190568.703	57386.6082	17693.1617	6656.09705	2948.71527	1302.47442	598.629251	276.547916	140.044867	79.0806048	50.969725	31.2524332	26.9827706	22.713108];
CenFluxMean_AE9_90 = [2015489.32	966624.012	492501.967	183035.525	54030.6479	16263.9459	6207.6321	2727.12992	1154.56719	486.691215	193.553499	89.7218939	49.3789625	33.4564409	19.2888158	16.0897458	12.8906759];
ApoFluxMean_AE9_90 = [543134.812	148679.823	49426.7817	15701.1899	4559.39628	1304.11571	444.805326	168.743676	74.6672731	33.1902862	16.3350572	9.49303049	6.74878269	5.04578025	3.11202649 2.64802017 2.18401386];

%Calculate flux value
% 25th Percetile
CenFlux25_90 = interp1(M_energy_beam,CenFlux25_AE9_90,r);
MaxFlux25_90 = interp1(M_energy_beam,MaxFlux25_AE9_90,r);
ApoFlux25_90 = interp1(M_energy_beam,ApoFlux25_AE9_90,r);

%50th Percetile
CenFlux50_90 = interp1(M_energy_beam,CenFlux50_AE9_90,r);
MaxFlux50_90 = interp1(M_energy_beam,MaxFlux50_AE9_90,r);
ApoFlux50_90 = interp1(M_energy_beam,ApoFlux50_AE9_90,r);

%75th Percetile
CenFlux75_90 = interp1(M_energy_beam,CenFlux75_AE9_90,r);
MaxFlux75_90 = interp1(M_energy_beam,MaxFlux75_AE9_90,r);
ApoFlux75_90 = interp1(M_energy_beam,ApoFlux75_AE9_90,r);

%95th Percetile
CenFlux95_90 = interp1(M_energy_beam,CenFlux95_AE9_90,r);
MaxFlux95_90 = interp1(M_energy_beam,MaxFlux95_AE9_90,r);
ApoFlux95_90 = interp1(M_energy_beam,ApoFlux95_AE9_90,r);

%Average
CenFluxMean_90 = interp1(M_energy_beam,CenFluxMean_AE9_90,r);
MaxFluxMean_90 = interp1(M_energy_beam,MaxFluxMean_AE9_90,r);
ApoFluxMean_90 = interp1(M_energy_beam,ApoFluxMean_AE9_90,r);

end

function [MaxFlux25_Omni,CenFlux25_Omni,ApoFlux25_Omni,MaxFlux50_Omni,CenFlux50_Omni,ApoFlux50_Omni,MaxFlux75_Omni,CenFlux75_Omni,ApoFlux75_Omni,MaxFlux95_Omni,CenFlux95_Omni,ApoFlux95_Omni,MaxFluxMean_Omni,CenFluxMean_Omni,ApoFluxMean_Omni] = FluxRateOmni(r)
%FluxRate Summary of this function goes here
%   Inputs:
%     M_energy_beam- energy at flux
%     Outputs:
%     Flux50:Expected average flux at energy M_energy_beam
%     Flux75:Expected 75th percentile flux at energy M_energy_beam
%     Flux95:Expected 95th percentile flux at energy M_energy_beam

M_energy_beam = [0.5,0.75,1:0.5:8];

% Conversion Factor to /sr from Omni
CF = 4*pi;

%% Omnidirectional
%Values from AE9
%25th Percentile
MaxFlux25_AE9_Omni= [8111977.47	1573149.04	570832.259	145269.2	29218.1317	6903.77267	2063.2227	742.263646	266.679309	109.348528	49.3155003	24.100156	13.0634938	7.06531883	3.28402254 2.64043001	1.99701516]./CF;
CenFlux25_AE9_Omni= [2363391.15	1010236.28	420136.43	123138.153	28252.2425	6593.3688	1950.68566	681.678876	259.525518	106.342774	48.4952391	23.0311398	11.8035589	6.45582063	3.0726934	2.47457182	1.87645024]./CF;
ApoFlux25_AE9_Omni= [1329508.46	252674.19	54324.3703	12858.4732	2958.76327	756.557986	235.070404	90.9015458	47.2026265	22.0241695	10.0087299	4.57676091	2.71873186	1.77166374	0.887809304	0.706650242	0.525491181]./CF;

%50th Percentile
MaxFlux50_AE9_Omni= [19893090.7	5097010.91	2097907.42	624573.54	154628.957	40909.9701	13350.1255	5049.95988	2026.19355	857.250336	375.862092	170.394972	85.9078482	50.603342	26.4107426	21.6188826	16.864596]./CF;
CenFlux50_AE9_Omni= [8902106.27	3997011.55	1816783.24	591581.383	150674.077	38761.4475	12597.5379	4815.77357	1932.24052	811.280178	348.745134	160.521736	82.9799251	49.2802646	25.5049115	20.8261552	16.1473989]./CF;
ApoFlux50_AE9_Omni= [3561954.92	801437.742	209995.397	55915.559	13994.7966	3726.22034	1214.48544	472.01705	231.316141	106.900138	50.3908547	25.2086036	15.7579626	10.6856645	5.76092933	4.68958599	3.61824264]./CF;

%75th Percentile
MaxFlux75_AE9_Omni= [40883035.7	12922351.4	6049174.9	2077265.16	590046.286	166606.194	58924.2062	23885.7512	10125.2594	4511.16945	1978.85991	914.656732	473.531305	291.004977	163.06214	136.183404	109.948726]./CF;
CenFlux75_AE9_Omni= [25466597	11884475.2	5783791.43	2043482.19	564388.867	156691.424	54819.331	22491.5776	9408.0567	4029.48609	1654.74889	743.967198	387.630602	245.668857	135.648925	112.031465	88.4140052]./CF;
ApoFlux75_AE9_Omni= [7748204.15	1991570.86	609961.739	178234.465	47666.3592	13100.9096	4433.70582	1730.77206	810.478391	371.653477	180.193512	96.8353773	63.2519282	44.4424177	25.4969211	21.1858925	16.8748639]./CF;

%95th Percentile
MaxFlux95_AE9_Omni= [91558289.4	41140003.2	21777989.6	8588339.04	2658679.36	821612.364	310389.614	138246.272	62029.605	28853.2118	13085.9523	6466.13666	3657.18038	2383.31583	1430.15343	1220.62378	1011.13351]./CF;
CenFlux95_AE9_Omni= [82438977.2	40152604.7	21041509.5	8128343.18	2453176.17	740869.067	281208.686	124800.997	54691.1672	23963.8502	9361.59773	4105.59847	2160.44605	1471.82903	872.166313	729.849019	587.531725]./CF;
ApoFlux95_AE9_Omni= [18390202.3	5481304.18	1997029.96	647185.244	186317.874	53039.0122	18715.0165	7344.80766	3270.54056	1486.48666	743.047183	432.730168	297.843794	218.480884	134.831044	114.986834	95.1426237]./CF;

%Mean
MaxFluxMean_AE9_Omni= [29719534   10294550.9 5145356.75	1972260.08	591967.488	180958.825	67634.9545	29785.0655	13309.3245	6179.96504	2801.18165	1385.17797	791.346898	521.590481	314.831095	269.160117	223.703666]./CF;
CenFluxMean_AE9_Omni= [20975174.9 10008126   5076315.94	1896363.99	556804.147	164682.291	61532.9153	27007.8912	11780.2207	5153.02796	2022.54965	889.781998	467.824031	316.607283	186.893584	156.369112	125.84464]./CF;
ApoFluxMean_AE9_Omni= [5682405.63 1525791.21 504437.068	156112.947	43709.6229	12295.5287	4282.62169	1678.31053	759.102006	346.00969	171.227699	97.4892416	66.3249877	48.2615504	29.4447683	25.0441188	20.6434693]./CF;

%Calculate flux value
% 25th Percetile
CenFlux25_Omni = interp1(M_energy_beam,CenFlux25_AE9_Omni,r);
MaxFlux25_Omni = interp1(M_energy_beam,MaxFlux25_AE9_Omni,r);
ApoFlux25_Omni = interp1(M_energy_beam,ApoFlux25_AE9_Omni,r);

%50th Percetile
CenFlux50_Omni = interp1(M_energy_beam,CenFlux50_AE9_Omni,r);
MaxFlux50_Omni = interp1(M_energy_beam,MaxFlux50_AE9_Omni,r);
ApoFlux50_Omni = interp1(M_energy_beam,ApoFlux50_AE9_Omni,r);

%75th Percetile
CenFlux75_Omni = interp1(M_energy_beam,CenFlux75_AE9_Omni,r);
MaxFlux75_Omni = interp1(M_energy_beam,MaxFlux75_AE9_Omni,r);
ApoFlux75_Omni = interp1(M_energy_beam,ApoFlux75_AE9_Omni,r);

%95th Percetile
CenFlux95_Omni = interp1(M_energy_beam,CenFlux95_AE9_Omni,r);
MaxFlux95_Omni = interp1(M_energy_beam,MaxFlux95_AE9_Omni,r);
ApoFlux95_Omni = interp1(M_energy_beam,ApoFlux95_AE9_Omni,r);

%Average
CenFluxMean_Omni = interp1(M_energy_beam,CenFluxMean_AE9_Omni,r);
MaxFluxMean_Omni = interp1(M_energy_beam,MaxFluxMean_AE9_Omni,r);
ApoFluxMean_Omni = interp1(M_energy_beam,ApoFluxMean_AE9_Omni,r);

end


