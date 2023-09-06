%%BowTie .m
% BowTie and Count Rate Analysis for HERT
% Last modified: 4/15/2023

clc
close all
%% Create expontential fits to AE9AP9 data for Omnidirectional
CF = 4*pi;
x_fit = [0.5,0.75,1:0.5:7];
%Values from AE9
%25th Percentile
MaxFlux25_AE9Omni= [8111977.47	1573149.04	570832.259	145269.2	29218.1317	6903.77267	2063.2227	742.263646	266.679309	109.348528	49.3155003	24.100156	13.0634938	7.06531883	3.28402254]./CF;
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

%Calculate flux value
% 25th Percetile
f_25_MaxOmni = fit(x_fit',MaxFlux25_AE9Omni','exp1');
f_25_CenOmni = fit(x_fit',CenFlux25_AE9Omni','exp1');
f_25_ApoOmni = fit(x_fit',ApoFlux25_AE9Omni','exp1');

%50th Percetile
f_50_MaxOmni = fit(x_fit',MaxFlux50_AE9Omni','exp1');
f_50_CenOmni = fit(x_fit',CenFlux50_AE9Omni','exp1');
f_50_ApoOmni = fit(x_fit',ApoFlux50_AE9Omni','exp1');

%75th Percetile
f_75_MaxOmni = fit(x_fit',MaxFlux75_AE9Omni','exp1');
f_75_CenOmni = fit(x_fit',CenFlux75_AE9Omni','exp1');
f_75_ApoOmni = fit(x_fit',ApoFlux75_AE9Omni','exp1');

%95th Percetile
f_95_MaxOmni = fit(x_fit',MaxFlux95_AE9Omni','exp1');
f_95_CenOmni = fit(x_fit',CenFlux95_AE9Omni','exp1');
f_95_ApoOmni = fit(x_fit',ApoFlux95_AE9Omni','exp1');

%95th Percetile
f_Mean_MaxOmni = fit(x_fit',MaxFluxMean_AE9Omni','exp1');
f_Mean_CenOmni = fit(x_fit',CenFluxMean_AE9Omni','exp1');
f_Mean_ApoOmni = fit(x_fit',ApoFluxMean_AE9Omni','exp1');

%% Create expontential fits to AE9AP9 data for 90 degree pitch angle
%Values from AE9
%25th Percentile
MaxFlux25_AE9_90= [1299128.18	159571.67	58982.02	15129.9237	3041.56181	718.45857	219.159343	77.6956967	26.2584103	10.3683737	4.73543853	2.43125066	1.31421209	0.694297103	0.328897155];
CenFlux25_AE9_90= [236773.987	102478.867	42665.4286	12486.5777	2898.92566	683.399265	201.636631	69.2078129	24.9863608	10.0010258	4.63892914	2.28897503	1.19219726	0.648284687	0.312877121];%	0.253244987	0.193612852
%[10912.9057	825.198408	330.449704	169.54976	77.7162659	31.8123459	13.0916374	4.78518407	1.8178498	0.769833692	0.39312206	0.230325859	0.140320451	0.085666892	0.046729895];
ApoFlux25_AE9_90= [127986.22	24967.9964	5444.69973	1319.99226	312.878706	81.8187931	25.3755998	9.72523144	4.90561484	2.21073182	0.957248098	0.416308498	0.239633178	0.152650546	0.07342178];

%50th Percentile
MaxFlux50_AE9_90= [3043556.52	506402.542	211356.954	62748.4502	15499.6257	4167.08963	1362.59182	508.053908	196.611127	82.051654	37.154501	17.4685642	8.97551526	5.28536699	2.81916956];
CenFlux50_AE9_90= [876087.841	397179.856	181014.549	58780.2084	15120.1328	3942.52241	1290.06064	487.952961	187.386107	76.4888646	33.3951236	16.0558237	8.53704862	5.05846258	2.61591484];	%2.14087987	1.66584491
%87065.4841	9000.52661	3466.00034	1653.26674	737.498429	299.675391	128.223723	51.527098	21.5263988	10.2522321	5.63969532	3.47792809	2.17690499	1.39891194	0.80894955];
ApoFlux50_AE9_90= [341776.994	78702.5833	20845.5994	5693.35574	1472.38512	399.944599	129.090197	49.2612793	23.5157668	10.5373774	4.81524394	2.36038762	1.4794959	1.00421865	0.532351405];

%75th Percentile
MaxFlux75_AE9_90= [5954877.83	1258557.25	590380.762	202629.404	58159.668	16656.5776	5906.61956	2397.14709	991.163556	436.275356	196.258873	94.3908837	49.8709683	30.3881283	16.5763025];
CenFlux75_AE9_90=[2469165.81	1160706.51	567296.774	199717.547	55612.4246	15693.0662	5571.19522	2275.12434	917.33847	380.385147	158.469032	74.7436904	40.4205575	25.6070785	13.9684838];%	11.5349933	9.10150277]
%[430398.486	44763.2346	14269.359	5545.31911	2183.6765	847.127292	380.657981	172.490309	81.783271	44.5794151	27.1835619	18.1506675	12.0351114	8.21203002	5.02883858];
ApoFlux75_AE9_90= [741554.833	194643.441	60101.779	18034.1354	4996.24183	1398.42223	465.736857	177.111453	80.9439822	36.1004976	17.2120895	9.2718457	6.21568293	4.43625712	2.53942722];

%95th Percentile
MaxFlux95_AE9_90= [12559884.1	3917196.38	2094895.39	826481.191	256040.313	79982.653	30444.5749	13671.416	6070.77165	2795.68906	1292.31374	654.009611	366.498171	235.290299	139.450819];
CenFlux95_AE9_90= [7854048.94	3842396.1	2026654.85	779563.225	236687.829	72888.9963	28329.6745	12599.2993	5363.30435	2263.84066	895.985629	414.265436	228.444014	155.761882	90.0617119];%	75.1379942	60.2142765
    %2404183.7	286581.807	80278.5175	27311.1463	10245.9702	3822.18999	1828.73641	863.673911	425.782876	241.123197	154.292021	105.154032	71.3507142	50.6479176	32.2673942];
ApoFlux95_AE9_90= [1753727.44	532323.177	194919.825	64926.6752	19416.2462	5617.60949	1937.2571	734.009612	319.68837	141.834443	70.8966219	42.4119355	30.6420444	23.1319331	14.4237781];

%Mean
MaxFluxMean_AE9_90= [2015489.32	966624.012	492501.967	183035.525	54030.6479	16263.9459	6207.6321	2727.12992	1154.56719	486.691215	193.553499	89.7218939	49.3789625	33.4564409	19.2888158];%	16.0897458	12.8906759]
CenFluxMean_AE9_90= [2015489.32	966624.012	492501.967	183035.525	54030.6479	16263.9459	6207.6321	2727.12992	1154.56719	486.691215	193.553499	89.7218939	49.3789625	33.4564409	19.2888158];%	16.0897458	12.8906759]
ApoFluxMean_AE9_90= [543134.812	148679.823	49426.7817	15701.1899	4559.39628	1304.11571	444.805326	168.743676	74.6672731	33.1902862	16.3350572	9.49303049	6.74878269	5.04578025	3.11202649];

%Calculate flux value
% 25th Percetile
f_25_Max_90 = fit(x_fit',MaxFlux25_AE9_90','exp1');
f_25_Cen_90 = fit(x_fit',CenFlux25_AE9_90','exp1');
f_25_Apo_90 = fit(x_fit',ApoFlux25_AE9_90','exp1');

%50th Percetile
f_50_Max_90 = fit(x_fit',MaxFlux50_AE9_90','exp1');
f_50_Cen_90 = fit(x_fit',CenFlux50_AE9_90','exp1');
f_50_Apo_90 = fit(x_fit',ApoFlux50_AE9_90','exp1');

%75th Percetile
f_75_Max_90 = fit(x_fit',MaxFlux75_AE9_90','exp1');
f_75_Cen_90 = fit(x_fit',CenFlux75_AE9_90','exp1');
f_75_Apo_90 = fit(x_fit',ApoFlux75_AE9_90','exp1');

%95th Percetile
f_95_Max_90 = fit(x_fit',MaxFlux95_AE9_90','exp1');
f_95_Cen_90 = fit(x_fit',CenFlux95_AE9_90','exp1');
f_95_Apo_90 = fit(x_fit',ApoFlux95_AE9_90','exp1');

%95th Percetile
f_Mean_Max_90 = fit(x_fit',MaxFluxMean_AE9_90','exp1');
f_Mean_Cen_90 = fit(x_fit',CenFluxMean_AE9_90','exp1');
f_Mean_Apo_90 = fit(x_fit',ApoFluxMean_AE9_90','exp1');


%% AE9 Flux Model using Omni/4pi
MaxFlux25Omni_ = zeros(length(x(:,1)),1);
CenFlux25Omni_ = zeros(length(x(:,1)),1);
ApoFlux25Omni_ = zeros(length(x(:,1)),1);
MaxFlux50Omni_ = zeros(length(x(:,1)),1);
CenFlux50Omni_ = zeros(length(x(:,1)),1);
ApoFlux50Omni_ = zeros(length(x(:,1)),1);
MaxFlux75Omni_ = zeros(length(x(:,1)),1);
CenFlux75Omni_ = zeros(length(x(:,1)),1);
ApoFlux75Omni_ = zeros(length(x(:,1)),1);
MaxFlux95Omni_ = zeros(length(x(:,1)),1);
CenFlux95Omni_ = zeros(length(x(:,1)),1);
ApoFlux95Omni_ = zeros(length(x(:,1)),1);
MaxFluxMeanOmni_ = zeros(length(x(:,1)),1);
CenFluxMeanOmni_ = zeros(length(x(:,1)),1);
ApoFluxMeanOmni_ = zeros(length(x(:,1)),1);


for i = 1:length(x(:,1))
    
    [MaxFlux25Omni_(i),CenFlux25Omni_(i),ApoFlux25Omni_(i),MaxFlux50Omni_(i),CenFlux50Omni_(i),ApoFlux50Omni_(i),MaxFlux75Omni_(i),CenFlux75Omni_(i),ApoFlux75Omni_(i),MaxFlux95Omni_(i),CenFlux95Omni_(i),ApoFlux95Omni_(i),MaxFluxMeanOmni_(i),CenFluxMeanOmni_(i),ApoFluxMeanOmni_(i)] = FluxRateOmni(x(i,1));
    MaxFlux25Omni_test(i) = f_25_MaxOmni(x(i,1));
    CenFlux25Omni_test(i) = f_25_CenOmni(x(i,1));
    ApoFlux25Omni_test(i) = f_25_ApoOmni(x(i,1));
    
    MaxFlux50Omni_test(i) = f_50_MaxOmni(x(i,1));
    CenFlux50Omni_test(i) = f_50_CenOmni(x(i,1));
    ApoFlux50Omni_test(i) = f_50_ApoOmni(x(i,1));
    
    MaxFlux75Omni_test(i) = f_75_MaxOmni(x(i,1));
    CenFlux75Omni_test(i) = f_75_CenOmni(x(i,1));
    ApoFlux75Omni_test(i) = f_75_ApoOmni(x(i,1));
    
    MaxFlux95Omni_test(i) = f_95_MaxOmni(x(i,1));
    CenFlux95Omni_test(i) = f_95_CenOmni(x(i,1));
    ApoFlux95Omni_test(i) = f_95_ApoOmni(x(i,1));
    
    MaxFluxMeanOmni_test(i) = f_Mean_MaxOmni(x(i,1));
    CenFluxMeanOmni_test(i) = f_Mean_CenOmni(x(i,1));
    ApoFluxMeanOmni_test(i) = f_Mean_ApoOmni(x(i,1));
    
end

figure
semilogy(x(:,1),MaxFlux25Omni_,x(:,1),CenFlux25Omni_,x(:,1),ApoFlux25Omni_,x(:,1),MaxFlux25Omni_test,x(:,1),CenFlux25Omni_test,x(:,1),ApoFlux25Omni_test,'Linewidth',2)
Title = append('OMNI/4*pi 25th Percentile Flux(E) Average',' E Range: ',num2str(min(x(:,1))),' to ',num2str(max(x(:,1))),' MeV');
title(Title)
ylim ([0 10^8])
ylabel('Flux (#/cm^2/s/sr/Mev)')
xlabel('Incident Energy (MeV)')
legend('Max','Cen','Apo','Max_Test','Cen_Test','Apo_Test')
saveas(gcf,append('Omni 25th Flux','.jpg'))

figure
semilogy(x(:,1),MaxFlux50Omni_,x(:,1),CenFlux50Omni_,x(:,1),ApoFlux50Omni_,'Linewidth',2)
Title = append('OMNI/4*pi 50th Percentile Flux(E) Average',' E Range: ',num2str(min(x(:,1))),' to ',num2str(max(x(:,1))),' MeV');
title(Title)
ylim ([0 10^8])
ylabel('Flux (#/cm^2/s/sr/Mev)')
xlabel('Incident Energy (MeV)')
legend('Max','Cen','Apo')
saveas(gcf,append('Omni 50th Flux','.jpg'))

figure
semilogy(x(:,1),MaxFlux75Omni_,x(:,1),CenFlux75Omni_,x(:,1),ApoFlux75Omni_,'Linewidth',2)
Title = append('OMNI/4*pi 75th Percentile Flux(E) Average',' E Range: ',num2str(min(x(:,1))),' to ',num2str(max(x(:,1))),' MeV');
title(Title)
ylim ([0 10^8])
ylabel('Flux (#/cm^2/s/sr/Mev)')
xlabel('Incident Energy (MeV)')
legend('Max','Cen','Apo')
saveas(gcf,append('Omni 75th Flux','.jpg'))

figure
semilogy(x(:,1),MaxFlux95Omni_,x(:,1),CenFlux95Omni_,x(:,1),ApoFlux95Omni_,'Linewidth',2)
Title = append('OMNI/4*pi 95th Percentile Flux(E) Average',' E Range: ',num2str(min(x(:,1))),' to ',num2str(max(x(:,1))),' MeV');
title(Title)
ylim ([0 10^8])
ylabel('Flux (#/cm^2/s/sr/Mev)')
xlabel('Incident Energy (MeV)')
legend('Max','Cen','Apo')
saveas(gcf,append('Omni 95th Flux','.jpg'))

figure
semilogy(x(:,1),MaxFluxMeanOmni_,x(:,1),CenFluxMeanOmni_,x(:,1),ApoFluxMeanOmni_,'Linewidth',2)
Title = append('OMNI/4*pi Mean Flux(E) Average',' E Range: ',num2str(min(x(:,1))),' to ',num2str(max(x(:,1))),' MeV');
title(Title)
ylim ([0 10^8])
ylabel('Flux (#/cm^2/s/sr/Mev)')
xlabel('Incident Energy (MeV)')
legend('Max','Cen','Apo')
saveas(gcf,append('Omni Mean Flux','.jpg'))

%% AE9 Flux Model using 90 deg pitch angle
MaxFlux25_90_ = zeros(length(x(:,1)),1);
CenFlux25_90_ = zeros(length(x(:,1)),1);
ApoFlux25_90_ = zeros(length(x(:,1)),1);
MaxFlux50_90_ = zeros(length(x(:,1)),1);
CenFlux50_90_ = zeros(length(x(:,1)),1);
ApoFlux50_90_ = zeros(length(x(:,1)),1);
MaxFlux75_90_ = zeros(length(x(:,1)),1);
CenFlux75_90_ = zeros(length(x(:,1)),1);
ApoFlux75_90_ = zeros(length(x(:,1)),1);
MaxFlux95_90_ = zeros(length(x(:,1)),1);
CenFlux95_90_ = zeros(length(x(:,1)),1);
ApoFlux95_90_ = zeros(length(x(:,1)),1);
MaxFluxMean_90_ = zeros(length(x(:,1)),1);
CenFluxMean_90_ = zeros(length(x(:,1)),1);
ApoFluxMean_90_ = zeros(length(x(:,1)),1);


for i = 1:length(x(:,1))
    [MaxFlux25_90_(i),CenFlux25_90_(i),ApoFlux25_90_(i),MaxFlux50_90_(i),CenFlux50_90_(i),ApoFlux50_90_(i),MaxFlux75_90_(i),CenFlux75_90_(i),ApoFlux75_90_(i),MaxFlux95_90_(i),CenFlux95_90_(i),ApoFlux95_90_(i),MaxFluxMean_90_(i),CenFluxMean_90_(i),ApoFluxMean_90_(i)] = FluxRate90(x(i,1));
    
    MaxFlux25_90_test(i) = f_25_Max_90(x(i,1));
    CenFlux25_90_test(i) = f_25_Cen_90(x(i,1));
    ApoFlux25_90_test(i) = f_25_Apo_90(x(i,1));
    
    MaxFlux50_90_test(i) = f_50_Max_90(x(i,1));
    CenFlux50_90_test(i) = f_50_Cen_90(x(i,1));
    ApoFlux50_90_test(i) = f_50_Apo_90(x(i,1));
    
    MaxFlux75_90_test(i) = f_75_Max_90(x(i,1));
    CenFlux75_90_test(i) = f_75_Cen_90(x(i,1));
    ApoFlux75_90_test(i) = f_75_Apo_90(x(i,1));
    
    MaxFlux95_90_test(i) = f_95_Max_90(x(i,1));
    CenFlux95_90_test(i) = f_95_Cen_90(x(i,1));
    ApoFlux95_90_test(i) = f_95_Apo_90(x(i,1));
    
    MaxFluxMean_90_test(i) = f_Mean_Max_90(x(i,1));
    CenFluxMean_90_test(i) = f_Mean_Cen_90(x(i,1));
    ApoFluxMean_90_test(i) = f_Mean_Apo_90(x(i,1));
end

figure
semilogy(x(:,1),CenFlux25_90_,x(:,1),ApoFlux25_90_,'Linewidth',2)
Title = append('90 deg pitch 25th Percentile Flux(E) Average',' E Range: ',num2str(min(x(:,1))),' to ',num2str(max(x(:,1))),' MeV');
title(Title)
ylim ([0 10^8])
ylabel('Flux (#/cm^2/s/sr/Mev)')
xlabel('Incident Energy (MeV)')
legend('Cen','Apo')
saveas(gcf,append('90 deg Pitch Angle 25th Flux','.jpg'))

figure
semilogy(x(:,1),CenFlux50_90_,x(:,1),ApoFlux50_90_,'Linewidth',2)
Title = append('90 deg pitch 50th Percentile Flux(E) Average',' E Range: ',num2str(min(x(:,1))),' to ',num2str(max(x(:,1))),' MeV');
title(Title)
ylim ([0 10^8])
ylabel('Flux (#/cm^2/s/sr/Mev)')
xlabel('Incident Energy (MeV)')
legend('Cen','Apo')
saveas(gcf,append('90 deg Pitch Angle 50th Flux','.jpg'))

figure
semilogy(x(:,1),CenFlux75_90_,x(:,1),ApoFlux75_90_,'Linewidth',2)
Title = append('90 deg pitch 75th Percentile Flux(E) Average',' E Range: ',num2str(min(x(:,1))),' to ',num2str(max(x(:,1))),' MeV');
title(Title)
ylim ([0 10^8])
ylabel('Flux (#/cm^2/s/sr/Mev)')
xlabel('Incident Energy (MeV)')
legend('Cen','Apo')
saveas(gcf,append('90 deg Pitch Angle 75th Flux','.jpg'))

figure
semilogy(x(:,1),CenFlux95_90_,x(:,1),ApoFlux95_90_,'Linewidth',2)
Title = append('90 deg pitch 95th Percentile Flux(E) Average',' E Range: ',num2str(min(x(:,1))),' to ',num2str(max(x(:,1))),' MeV');
title(Title)
ylim ([0 10^8])
ylabel('Flux (#/cm^2/s/sr/Mev)')
xlabel('Incident Energy (MeV)')
legend('Cen','Apo')
saveas(gcf,append('90 deg Pitch Angle 95th Flux','.jpg'))

figure
semilogy(x(:,1),CenFluxMean_90_,x(:,1),ApoFluxMean_90_,'Linewidth',2)
Title = append('90 deg pitch Mean Flux(E) Average',' E Range: ',num2str(min(x(:,1))),' to ',num2str(max(x(:,1))),' MeV');
title(Title)
ylim ([0 10^8])
ylabel('Flux (#/cm^2/s/sr/Mev)')
xlabel('Incident Energy (MeV)')
legend('Cen','Apo')
saveas(gcf,append('90 deg Pitch Angle Mean Flux','.jpg'))



%% Bow Tie Analysis-Selesnick/Blake-Inner


if filetype==0
    % Changes directory for jpg of graphs
    cd 'Bow Tie\Inner'
    
    %Sets Ei and Range of Eo
    %Ei_inner is the same as X from the GEANT4 results
    Ei_inner = x(:,1);
    %Eo_inner set by user.
    Eo_inner = 0.2:0.2:2.0;
    
    
    %Sets up color vectors for plotting the different Eo curves
    Eo_color_inner = magma(length(Eo_inner)+1);
    
    %Preallocates all variables prior to For Loops
    J_e_inner = zeros(length(x),length(Eo_inner));
    
    G_int_inner = zeros(length(x),length(Eo_inner),length(energy_channels));
    G_term_inner = zeros(length(energy_channels),length(Eo_inner));
    G_E_eff_inner = zeros(length(x),length(Eo_inner),length(energy_channels));
    
    xi_inner = zeros(sum(1:length(Eo_inner)-1),length(energy_channels));
    yi_inner = xi_inner;
    
    E_eff_inner = zeros(1,length(energy_channels));
    G_eff_dE_inner= E_eff_inner;
    BowTieLegend_inner = strings([1,length(Eo_inner)]);
    
    BowCount_Rate_InnerMax50_ = zeros(length(energy_channels),1);
    BowCount_Rate_InnerMax75_ = zeros(length(energy_channels),1);
    BowCount_Rate_InnerMax95_ = zeros(length(energy_channels),1);
    
    %Count_Rate_Inner =  zeros(length(energy_channels),length(Eo_inner));
    
    BinWidth_Inner = zeros(1, length(energy_channels));
    Geff_Inner = zeros(1, length(energy_channels));
    
    %Creates J(e) and creates String Array for Plot Legends
    for i = 1:length(Eo_inner)
        J_e_inner(:,i) = exp(-Ei_inner/Eo_inner(i));
        BowTieLegend_inner(i) = num2str(Eo_inner(i));
        
    end
    
    %Adds 'Average Intersection Point' for plotting
    BowTieLegend_inner = [BowTieLegend_inner, 'Average Intersection Point'];
    
    %Creates J(e)^-1
    J_e_inv_inner = 1./J_e_inner;
    
    %% Preallocation
    BowCount_Rate_InnerMax25_Omni = zeros(length(energy_channels),1);
    BowCount_Rate_InnerCen25_Omni = zeros(length(energy_channels),1);
    BowCount_Rate_InnerApo25_Omni = zeros(length(energy_channels),1);
    
    %50thPercentile
    BowCount_Rate_InnerMax50_Omni = zeros(length(energy_channels),1);
    BowCount_Rate_InnerCen50_Omni = zeros(length(energy_channels),1);
    BowCount_Rate_InnerApo50_Omni = zeros(length(energy_channels),1);
    
    %75thPercentile
    BowCount_Rate_InnerMax75_Omni =zeros(length(energy_channels),1);
    BowCount_Rate_InnerCen75_Omni = zeros(length(energy_channels),1);
    BowCount_Rate_InnerApo75_Omni = zeros(length(energy_channels),1);
    
    %95th Percentile
    BowCount_Rate_InnerMax95_Omni = zeros(length(energy_channels),1);
    BowCount_Rate_InnerCen95_Omni = zeros(length(energy_channels),1);
    BowCount_Rate_InnerApo95_Omni = zeros(length(energy_channels),1);
    
    %Average
    BowCount_Rate_InnerMaxMean_Omni = zeros(length(energy_channels),1);
    BowCount_Rate_InnerCenMean_Omni = zeros(length(energy_channels),1);
    BowCount_Rate_InnerApoMean_Omni = zeros(length(energy_channels),1);
    
    BowCount_Rate_InnerMax25__90 = zeros(length(energy_channels),1);
    BowCount_Rate_InnerCen25__90 = zeros(length(energy_channels),1);
    BowCount_Rate_InnerApo25__90 = zeros(length(energy_channels),1);
    
    %50thPercentile
    BowCount_Rate_InnerMax50__90 = zeros(length(energy_channels),1);
    BowCount_Rate_InnerCen50__90 = zeros(length(energy_channels),1);
    BowCount_Rate_InnerApo50__90 = zeros(length(energy_channels),1);
    
    %75thPercentile
    BowCount_Rate_InnerMax75__90 = zeros(length(energy_channels),1);
    BowCount_Rate_InnerCen75__90 = zeros(length(energy_channels),1);
    BowCount_Rate_InnerApo75__90 = zeros(length(energy_channels),1);
    
    %95th Percentile
    BowCount_Rate_InnerMax95__90 = zeros(length(energy_channels),1);
    BowCount_Rate_InnerCen95__90 = zeros(length(energy_channels),1);
    BowCount_Rate_InnerApo95__90 = zeros(length(energy_channels),1);
    
    %Average
    BowCount_Rate_InnerMaxMean__90 = zeros(length(energy_channels),1);
    BowCount_Rate_InnerCenMean__90 = zeros(length(energy_channels),1);
    BowCount_Rate_InnerApoMean__90 = zeros(length(energy_channels),1);
    
    %25th Percentile
    IntDetector_AllRate_InnerMax25_Omni = zeros(width(inner_detector_GAllCounts),1);
    IntDetector_AllRate_InnerCen25_Omni = zeros(width(inner_detector_GAllCounts),1);
    IntDetector_AllRate_InnerApo25_Omni = zeros(width(inner_detector_GAllCounts),1);
    
    %50th Percentile
    IntDetector_AllRate_InnerMax50_Omni = zeros(width(inner_detector_GAllCounts),1);
    IntDetector_AllRate_InnerCen50_Omni = zeros(width(inner_detector_GAllCounts),1);
    IntDetector_AllRate_InnerApo50_Omni = zeros(width(inner_detector_GAllCounts),1);
    
    %75th Percentile
    IntDetector_AllRate_InnerMax75_Omni = zeros(width(inner_detector_GAllCounts),1);
    IntDetector_AllRate_InnerCen75_Omni = zeros(width(inner_detector_GAllCounts),1);
    IntDetector_AllRate_InnerApo75_Omni = zeros(width(inner_detector_GAllCounts),1);
    
    %95th Percentile
    IntDetector_AllRate_InnerMax95_Omni = zeros(width(inner_detector_GAllCounts),1);
    IntDetector_AllRate_InnerCen95_Omni =zeros(width(inner_detector_GAllCounts),1);
    IntDetector_AllRate_InnerApo95_Omni =zeros(width(inner_detector_GAllCounts),1);
    
    %Mean
    IntDetector_AllRate_InnerMaxMean_Omni = zeros(width(inner_detector_GAllCounts),1);
    IntDetector_AllRate_InnerCenMean_Omni = zeros(width(inner_detector_GAllCounts),1);
    IntDetector_AllRate_InnerApoMean_Omni = zeros(width(inner_detector_GAllCounts),1);
    
    %%
    %For Loop for calculating a line for each Eo and finding the average intersection point
    
    
    for c=1:length(energy_channels)
        %Prints energy channel to show progress
        fprintf('Energy Channel-Inner # %.2f \n',c)
        
        %Calculates line for each Eo
        for i = 1:length(Eo_inner)
            G_int_inner(:,i,c) = geo_EC_inner(:,c).*J_e_inner(:,i);
            G_term_inner(i,c) = trapz(Ei_inner,G_int_inner(:,i,c));
            G_E_eff_inner(:,i,c)= G_term_inner(i,c)*J_e_inv_inner(:,i);
            
        end
        
        %Find Intersections of Eo Lines
        num = 0;
        for j = 1:(length(Eo_inner)-1)
            
            for k = 1:(length(Eo_inner)-j)
                [xi_inner(num+k,c),yi_inner(num+k,c)] = polyxpoly(x(:,1),G_E_eff_inner(:,j,c),x(:,1),G_E_eff_inner(:,j+k,c),'unique');
                
            end
            num = num + k;
        end
        
        %Finds average intersection point
        E_eff_inner(c) = mean(xi_inner(:,c));
        G_eff_dE_inner(c) = mean(yi_inner(:,c));
        
        %Flux Rate from AE9 Model
        [MaxFlux25Omni,CenFlux25Omni,ApoFlux25Omni,MaxFlux50Omni,CenFlux50Omni,ApoFlux50Omni,MaxFlux75Omni,CenFlux75Omni,ApoFlux75Omni,MaxFlux95Omni,CenFlux95Omni,ApoFlux95Omni,MaxFluxMeanOmni,CenFluxMeanOmni,ApoFluxMeanOmni] = FluxRateOmni(E_eff_inner(c));
        
        %Calculates Estimated Count Rate for each Eo
        %25thPercentile
        BowCount_Rate_InnerMax25_Omni(c) = G_eff_dE_inner(c)*MaxFlux25Omni;
        BowCount_Rate_InnerCen25_Omni(c) = G_eff_dE_inner(c)*CenFlux25Omni;
        BowCount_Rate_InnerApo25_Omni(c) = G_eff_dE_inner(c)*ApoFlux25Omni;
        
        %50thPercentile
        BowCount_Rate_InnerMax50_Omni(c) = G_eff_dE_inner(c)*MaxFlux50Omni;
        BowCount_Rate_InnerCen50_Omni(c) = G_eff_dE_inner(c)*CenFlux50Omni;
        BowCount_Rate_InnerApo50_Omni(c) = G_eff_dE_inner(c)*ApoFlux50Omni;
        
        %75thPercentile
        BowCount_Rate_InnerMax75_Omni(c) = G_eff_dE_inner(c)*MaxFlux75Omni;
        BowCount_Rate_InnerCen75_Omni(c) = G_eff_dE_inner(c)*CenFlux75Omni;
        BowCount_Rate_InnerApo75_Omni(c) = G_eff_dE_inner(c)*ApoFlux75Omni;
        
        %95th Percentile
        BowCount_Rate_InnerMax95_Omni(c) = G_eff_dE_inner(c)*MaxFlux95Omni;
        BowCount_Rate_InnerCen95_Omni(c) = G_eff_dE_inner(c)*CenFlux95Omni;
        BowCount_Rate_InnerApo95_Omni(c) = G_eff_dE_inner(c)*ApoFlux95Omni;
        
        %Average
        BowCount_Rate_InnerMaxMean_Omni(c) = G_eff_dE_inner(c)*MaxFluxMeanOmni;
        BowCount_Rate_InnerCenMean_Omni(c) = G_eff_dE_inner(c)*CenFluxMeanOmni;
        BowCount_Rate_InnerApoMean_Omni(c) = G_eff_dE_inner(c)*ApoFluxMeanOmni;
        
        %Flux Rate from AE9 Model
        [MaxFlux25_90,CenFlux25_90,ApoFlux25_90,MaxFlux50_90,CenFlux50_90,ApoFlux50_90,MaxFlux75_90,CenFlux75_90,ApoFlux75_90,MaxFlux95_90,CenFlux95_90,ApoFlux95_90,MaxFluxMean_90,CenFluxMean_90,ApoFluxMean_90] = FluxRate90(E_eff_inner(c));
        
        %Calculates Estimated Count Rate for each Eo
        %25thPercentile
        BowCount_Rate_InnerMax25__90(c) = G_eff_dE_inner(c)*MaxFlux25_90;
        BowCount_Rate_InnerCen25__90(c) = G_eff_dE_inner(c)*CenFlux25_90;
        BowCount_Rate_InnerApo25__90(c) = G_eff_dE_inner(c)*ApoFlux25_90;
        
        %50thPercentile
        BowCount_Rate_InnerMax50__90(c) = G_eff_dE_inner(c)*MaxFlux50_90;
        BowCount_Rate_InnerCen50__90(c) = G_eff_dE_inner(c)*CenFlux50_90;
        BowCount_Rate_InnerApo50__90(c) = G_eff_dE_inner(c)*ApoFlux50_90;
        
        %75thPercentile
        BowCount_Rate_InnerMax75__90(c) = G_eff_dE_inner(c)*MaxFlux75_90;
        BowCount_Rate_InnerCen75__90(c) = G_eff_dE_inner(c)*CenFlux75_90;
        BowCount_Rate_InnerApo75__90(c) = G_eff_dE_inner(c)*ApoFlux75_90;
        
        %95th Percentile
        BowCount_Rate_InnerMax95__90(c) = G_eff_dE_inner(c)*MaxFlux95_90;
        BowCount_Rate_InnerCen95__90(c) = G_eff_dE_inner(c)*CenFlux95_90;
        BowCount_Rate_InnerApo95__90(c) = G_eff_dE_inner(c)*ApoFlux95_90;
        
        %Average
        BowCount_Rate_InnerMaxMean__90(c) = G_eff_dE_inner(c)*MaxFluxMean_90;
        BowCount_Rate_InnerCenMean__90(c) = G_eff_dE_inner(c)*CenFluxMean_90;
        BowCount_Rate_InnerApoMean__90(c) = G_eff_dE_inner(c)*ApoFluxMean_90;
        
        %Calcaulte Bin Characteristics
        Geff_Inner(c) = G_eff_dE_inner(c)/fwhm_inner(c);
        BinWidth_Inner(c) = fwhm_inner(c);
        
    end
    
    
    
    %% Count rates using Omnidirectinola divided by 4*pi
    %BowTie Method
    %25thPercentile
    BowDetector1_Rate_InnerMax25_Omni = sum(BowCount_Rate_InnerMax25_Omni)
    BowDetector1_Rate_InnerCen25_Omni = sum(BowCount_Rate_InnerCen25_Omni)
    BowDetector1_Rate_InnerApo25_Omni = sum(BowCount_Rate_InnerApo25_Omni)
    
    %50thPercentile
    BowDetector1_Rate_InnerMax50_Omni = sum(BowCount_Rate_InnerMax50_Omni)
    BowDetector1_Rate_InnerCen50_Omni = sum(BowCount_Rate_InnerCen50_Omni)
    BowDetector1_Rate_InnerApo50_Omni = sum(BowCount_Rate_InnerApo50_Omni)
    
    %75thPercentile
    BowDetector1_Rate_InnerMax75_Omni = sum(BowCount_Rate_InnerMax75_Omni)
    BowDetector1_Rate_InnerCen75_Omni = sum(BowCount_Rate_InnerCen75_Omni)
    BowDetector1_Rate_InnerApo75_Omni = sum(BowCount_Rate_InnerApo75_Omni)
    
    %95thPercentile
    BowDetector1_Rate_InnerMax95_Omni = sum(BowCount_Rate_InnerMax95_Omni)
    BowDetector1_Rate_InnerCen95_Omni = sum(BowCount_Rate_InnerCen95_Omni)
    BowDetector1_Rate_InnerApo95_Omni = sum(BowCount_Rate_InnerApo95_Omni)
    
    %Mean
    BowDetector1_Rate_InnerMaxMean_Omni = sum(BowCount_Rate_InnerMaxMean_Omni)
    BowDetector1_Rate_InnerCenMean_Omni = sum(BowCount_Rate_InnerCenMean_Omni)
    BowDetector1_Rate_InnerApoMean_Omni = sum(BowCount_Rate_InnerApoMean_Omni)
    
    
    %Integral Method
    %25th Percentile
    IntDetector1_BinRate_InnerMax25_Omni = trapz(x(:,1),geo_EL_inner(:).*MaxFlux25Omni_)
    IntDetector1_BinRate_InnerCen25_Omni = trapz(x(:,1),geo_EL_inner(:).*CenFlux25Omni_)
    IntDetector1_BinRate_InnerApo25_Omni = trapz(x(:,1),geo_EL_inner(:).*ApoFlux25Omni_)
    
    %50th Percentile
    IntDetector1_BinRate_InnerMax50_Omni = trapz(x(:,1),geo_EL_inner(:).*MaxFlux50Omni_)
    IntDetector1_BinRate_InnerCen50_Omni = trapz(x(:,1),geo_EL_inner(:).*CenFlux50Omni_)
    IntDetector1_BinRate_InnerApo50_Omni = trapz(x(:,1),geo_EL_inner(:).*ApoFlux50Omni_)
    
    %75th Percentile
    IntDetector1_BinRate_InnerMax75_Omni = trapz(x(:,1),geo_EL_inner(:).*MaxFlux75Omni_)
    IntDetector1_BinRate_InnerCen75_Omni = trapz(x(:,1),geo_EL_inner(:).*CenFlux75Omni_)
    IntDetector1_BinRate_InnerApo75_Omni = trapz(x(:,1),geo_EL_inner(:).*ApoFlux75Omni_)
    
    %95th Percentile
    IntDetector1_BinRate_InnerMax95_Omni = trapz(x(:,1),geo_EL_inner(:).*MaxFlux95Omni_)
    IntDetector1_BinRate_InnerCen95_Omni = trapz(x(:,1),geo_EL_inner(:).*CenFlux95Omni_)
    IntDetector1_BinRate_InnerApo95_Omni = trapz(x(:,1),geo_EL_inner(:).*ApoFlux95Omni_)
    
    %Mean
    IntDetector1_BinRate_InnerMaxMean_Omni = trapz(x(:,1),geo_EL_inner(:).*MaxFluxMeanOmni_)
    IntDetector1_BinRate_InnerCenMean_Omni = trapz(x(:,1),geo_EL_inner(:).*CenFluxMeanOmni_)
    IntDetector1_BinRate_InnerApoMean_Omni = trapz(x(:,1),geo_EL_inner(:).*ApoFluxMeanOmni_)
    
    for i = 1:width(inner_detector_GAllCounts)
        %25th Percentile
        IntDetector_AllRate_InnerMax25_Omni(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*MaxFlux25Omni_);
        IntDetector_AllRate_InnerCen25_Omni(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*CenFlux25Omni_);
        IntDetector_AllRate_InnerApo25_Omni(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*ApoFlux25Omni_);
        
        %50th Percentile
        IntDetector_AllRate_InnerMax50_Omni(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*MaxFlux50Omni_);
        IntDetector_AllRate_InnerCen50_Omni(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*CenFlux50Omni_);
        IntDetector_AllRate_InnerApo50_Omni(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*ApoFlux50Omni_);
        
        %75th Percentile
        IntDetector_AllRate_InnerMax75_Omni(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*MaxFlux75Omni_);
        IntDetector_AllRate_InnerCen75_Omni(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*CenFlux75Omni_);
        IntDetector_AllRate_InnerApo75_Omni(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*ApoFlux75Omni_);
        
        %95th Percentile
        IntDetector_AllRate_InnerMax95_Omni(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*MaxFlux95Omni_);
        IntDetector_AllRate_InnerCen95_Omni(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*CenFlux95Omni_);
        IntDetector_AllRate_InnerApo95_Omni(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*ApoFlux95Omni_);
        
        %Mean
        IntDetector_AllRate_InnerMaxMean_Omni(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*MaxFluxMeanOmni_);
        IntDetector_AllRate_InnerCenMean_Omni(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*CenFluxMeanOmni_);
        IntDetector_AllRate_InnerApoMean_Omni(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*ApoFluxMeanOmni_);
        
    end
    
    for i = 1:length(energy_channels)
        %25th Percentile
        IntEnergyChannel_InnerMax25_Omni(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*MaxFlux25Omni_);
        IntEnergyChannel_InnerCen25_Omni(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*CenFlux25Omni_);
        IntEnergyChannel_InnerApo25_Omni(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*ApoFlux25Omni_);
        
        %50th Percentile
        IntEnergyChannel_InnerMax50_Omni(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*MaxFlux50Omni_);
        IntEnergyChannel_InnerCen50_Omni(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*CenFlux50Omni_);
        IntEnergyChannel_InnerApo50_Omni(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*ApoFlux50Omni_);
        
        %75th Percentile
        IntEnergyChannel_InnerMax75_Omni(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*MaxFlux75Omni_);
        IntEnergyChannel_InnerCen75_Omni(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*CenFlux75Omni_);
        IntEnergyChannel_InnerApo75_Omni(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*ApoFlux75Omni_);
        
        %95th Percentile
        IntEnergyChannel_InnerMax95_Omni(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*MaxFlux95Omni_);
        IntEnergyChannel_InnerCen95_Omni(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*CenFlux95Omni_);
        IntEnergyChannel_InnerApo95_Omni(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*ApoFlux95Omni_);
        
        %Mean
        IntEnergyChannel_InnerMaxMean_Omni(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*MaxFluxMeanOmni_);
        IntEnergyChannel_InnerCenMean_Omni(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*CenFluxMeanOmni_);
        IntEnergyChannel_InnerApoMean_Omni(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*ApoFluxMeanOmni_);
        
    end
    
    
    for i = 1:width(inner_detector_GEnergy)
        
        %25th Percentile
        IntDetector_ERate_InnerMax25_Omni(i) = trapz(x(:,1),inner_detector_GEnergy(:,i).*MaxFlux25Omni_);
        IntDetector_ERate_InnerCen25_Omni(i) = trapz(x(:,1),inner_detector_GEnergy(:,i).*CenFlux25Omni_);
        IntDetector_ERate_InnerApo25_Omni(i) = trapz(x(:,1),inner_detector_GEnergy(:,i).*ApoFlux25Omni_);
        
        %50th Percentile
        IntDetector_ERate_InnerMax50_Omni(i) = trapz(x(:,1),inner_detector_GEnergy(:,i).*MaxFlux50Omni_);
        IntDetector_ERate_InnerCen50_Omni(i) = trapz(x(:,1),inner_detector_GEnergy(:,i).*CenFlux50Omni_);
        IntDetector_ERate_InnerApo50_Omni(i) = trapz(x(:,1),inner_detector_GEnergy(:,i).*ApoFlux50Omni_);
        
        %75th Percentile
        IntDetector_ERate_InnerMax75_Omni(i) = trapz(x(:,1),inner_detector_GEnergy(:,i).*MaxFlux75Omni_);
        IntDetector_ERate_InnerCen75_Omni(i) = trapz(x(:,1),inner_detector_GEnergy(:,i).*CenFlux75Omni_);
        IntDetector_ERate_InnerApo75_Omni(i) = trapz(x(:,1),inner_detector_GEnergy(:,i).*ApoFlux75Omni_);
        
        %95th Percentile
        IntDetector_ERate_InnerMax95_Omni(i) = trapz(x(:,1),inner_detector_GEnergy(:,i).*MaxFlux95Omni_);
        IntDetector_ERate_InnerCen95_Omni(i) = trapz(x(:,1),inner_detector_GEnergy(:,i).*CenFlux95Omni_);
        IntDetector_ERate_InnerApo95_Omni(i) = trapz(x(:,1),inner_detector_GEnergy(:,i).*ApoFlux95Omni_);
        
        %Mean
        IntDetector_ERate_InnerMaxMean_Omni(i) = trapz(x(:,1),inner_detector_GEnergy(:,i).*MaxFluxMeanOmni_);
        IntDetector_ERate_InnerCenMean_Omni(i) = trapz(x(:,1),inner_detector_GEnergy(:,i).*CenFluxMeanOmni_);
        IntDetector_ERate_InnerApoMean_Omni(i) = trapz(x(:,1),inner_detector_GEnergy(:,i).*ApoFluxMeanOmni_);
        
        
    end
    
    %% Count Rates Using 90 degree pitch angle
    %BowTie Method
    %25thPercentile
    BowDetector1_Rate_InnerMax25__90 = sum(BowCount_Rate_InnerMax25__90)
    BowDetector1_Rate_InnerCen25__90 = sum(BowCount_Rate_InnerCen25__90)
    BowDetector1_Rate_InnerApo25__90 = sum(BowCount_Rate_InnerApo25__90)
    
    %50thPercentile
    BowDetector1_Rate_InnerMax50__90 = sum(BowCount_Rate_InnerMax50__90)
    BowDetector1_Rate_InnerCen50__90 = sum(BowCount_Rate_InnerCen50__90)
    BowDetector1_Rate_InnerApo50__90 = sum(BowCount_Rate_InnerApo50__90)
    
    %75thPercentile
    BowDetector1_Rate_InnerMax75__90 = sum(BowCount_Rate_InnerMax75__90)
    BowDetector1_Rate_InnerCen75__90 = sum(BowCount_Rate_InnerCen75__90)
    BowDetector1_Rate_InnerApo75__90 = sum(BowCount_Rate_InnerApo75__90)
    
    %95thPercentile
    BowDetector1_Rate_InnerMax95__90 = sum(BowCount_Rate_InnerMax95__90)
    BowDetector1_Rate_InnerCen95__90 = sum(BowCount_Rate_InnerCen95__90)
    BowDetector1_Rate_InnerApo95__90 = sum(BowCount_Rate_InnerApo95__90)
    
    %Mean
    BowDetector1_Rate_InnerMaxMean__90 = sum(BowCount_Rate_InnerMaxMean__90)
    BowDetector1_Rate_InnerCenMean__90 = sum(BowCount_Rate_InnerCenMean__90)
    BowDetector1_Rate_InnerApoMean__90 = sum(BowCount_Rate_InnerApoMean__90)
    
    
    %Integral Method
    %25th Percentile
    IntDetector1_BinRate_InnerMax25__90 = trapz(x(:,1),geo_EL_inner(:).*MaxFlux25_90_)
    IntDetector1_BinRate_InnerCen25__90 = trapz(x(:,1),geo_EL_inner(:).*CenFlux25_90_)
    IntDetector1_BinRate_InnerApo25__90 = trapz(x(:,1),geo_EL_inner(:).*ApoFlux25_90_)
    
    %50th Percentile
    IntDetector1_BinRate_InnerMax50__90 = trapz(x(:,1),geo_EL_inner(:).*MaxFlux50_90_)
    IntDetector1_BinRate_InnerCen50__90 = trapz(x(:,1),geo_EL_inner(:).*CenFlux50_90_)
    IntDetector1_BinRate_InnerApo50__90 = trapz(x(:,1),geo_EL_inner(:).*ApoFlux50_90_)
    
    %75th Percentile
    IntDetector1_BinRate_InnerMax75__90 = trapz(x(:,1),geo_EL_inner(:).*MaxFlux75_90_)
    IntDetector1_BinRate_InnerCen75__90 = trapz(x(:,1),geo_EL_inner(:).*CenFlux75_90_)
    IntDetector1_BinRate_InnerApo75__90 = trapz(x(:,1),geo_EL_inner(:).*ApoFlux75_90_)
    
    %95th Percentile
    IntDetector1_BinRate_InnerMax95__90 = trapz(x(:,1),geo_EL_inner(:).*MaxFlux95_90_)
    IntDetector1_BinRate_InnerCen95__90 = trapz(x(:,1),geo_EL_inner(:).*CenFlux95_90_)
    IntDetector1_BinRate_InnerApo95__90 = trapz(x(:,1),geo_EL_inner(:).*ApoFlux95_90_)
    
    %Mean
    IntDetector1_BinRate_InnerMaxMean__90 = trapz(x(:,1),geo_EL_inner(:).*MaxFluxMean_90_)
    IntDetector1_BinRate_InnerCenMean__90 = trapz(x(:,1),geo_EL_inner(:).*CenFluxMean_90_)
    IntDetector1_BinRate_InnerApoMean__90 = trapz(x(:,1),geo_EL_inner(:).*ApoFluxMean_90_)
    
    for i = 1:width(inner_detector_GAllCounts)
        %25th Percentile
        IntDetector_AllRate_InnerMax25__90(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*MaxFlux25_90_)';
        IntDetector_AllRate_InnerCen25__90(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*CenFlux25_90_)';
        IntDetector_AllRate_InnerApo25__90(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*ApoFlux25_90_)';
        
        %50th Percentile
        IntDetector_AllRate_InnerMax50__90(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*MaxFlux50_90_)';
        IntDetector_AllRate_InnerCen50__90(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*CenFlux50_90_)';
        IntDetector_AllRate_InnerApo50__90(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*ApoFlux50_90_)';
        
        %75th Percentile
        IntDetector_AllRate_InnerMax75__90(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*MaxFlux75_90_)';
        IntDetector_AllRate_InnerCen75__90(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*CenFlux75_90_)';
        IntDetector_AllRate_InnerApo75__90(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*ApoFlux75_90_)';
        
        %95th Percentile
        IntDetector_AllRate_InnerMax95__90(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*MaxFlux95_90_)';
        IntDetector_AllRate_InnerCen95__90(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*CenFlux95_90_)';
        IntDetector_AllRate_InnerApo95__90(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*ApoFlux95_90_)';
        
        %Mean
        IntDetector_AllRate_InnerMaxMean__90(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*MaxFluxMean_90_);
        IntDetector_AllRate_InnerCenMean__90(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*CenFluxMean_90_);
        IntDetector_AllRate_InnerApoMean__90(i,1) = trapz(x(:,1),inner_detector_GAllCounts(:,i).*ApoFluxMean_90_);
        
    end
    
    for i = 1:length(energy_channels)
        %25th Percentile
        IntEnergyChannel_InnerMax25__90(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*MaxFlux25_90_);
        IntEnergyChannel_InnerCen25__90(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*CenFlux25_90_);
        IntEnergyChannel_InnerApo25__90(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*ApoFlux25_90_);
        
        %50th Percentile
        IntEnergyChannel_InnerMax50__90(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*MaxFlux50_90_);
        IntEnergyChannel_InnerCen50__90(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*CenFlux50_90_);
        IntEnergyChannel_InnerApo50__90(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*ApoFlux50_90_);
        
        %75th Percentile
        IntEnergyChannel_InnerMax75__90(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*MaxFlux75_90_);
        IntEnergyChannel_InnerCen75__90(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*CenFlux75_90_);
        IntEnergyChannel_InnerApo75__90(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*ApoFlux75_90_);
        
        %95th Percentile
        IntEnergyChannel_InnerMax95__90(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*MaxFlux95_90_);
        IntEnergyChannel_InnerCen95__90(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*CenFlux95_90_);
        IntEnergyChannel_InnerApo95__90(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*ApoFlux95_90_);
        
        %Mean
        IntEnergyChannel_InnerMaxMean__90(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*MaxFluxMean_90_);
        IntEnergyChannel_InnerCenMean__90(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*CenFluxMean_90_);
        IntEnergyChannel_InnerApoMean__90(i,1)=trapz(x(:,1),geo_EC_inner(:,i).*ApoFluxMean_90_);
        
    end
    
    
    for i = 1:width(inner_detector_GEnergy)
        
        %25th Percentile
        IntDetector_ERate_InnerMax25__90(i,1) = trapz(x(:,1),inner_detector_GEnergy(:,i).*MaxFlux25_90_);
        IntDetector_ERate_InnerCen25__90(i,1) = trapz(x(:,1),inner_detector_GEnergy(:,i).*CenFlux25_90_);
        IntDetector_ERate_InnerApo25__90(i,1) = trapz(x(:,1),inner_detector_GEnergy(:,i).*ApoFlux25_90_);
        
        %50th Percentile
        IntDetector_ERate_InnerMax50__90(i,1) = trapz(x(:,1),inner_detector_GEnergy(:,i).*MaxFlux50_90_);
        IntDetector_ERate_InnerCen50__90(i,1) = trapz(x(:,1),inner_detector_GEnergy(:,i).*CenFlux50_90_);
        IntDetector_ERate_InnerApo50__90(i,1) = trapz(x(:,1),inner_detector_GEnergy(:,i).*ApoFlux50_90_);
        
        %75th Percentile
        IntDetector_ERate_InnerMax75__90(i,1) = trapz(x(:,1),inner_detector_GEnergy(:,i).*MaxFlux75_90_);
        IntDetector_ERate_InnerCen75__90(i,1) = trapz(x(:,1),inner_detector_GEnergy(:,i).*CenFlux75_90_);
        IntDetector_ERate_InnerApo75__90(i,1) = trapz(x(:,1),inner_detector_GEnergy(:,i).*ApoFlux75_90_);
        
        %95th Percentile
        IntDetector_ERate_InnerMax95__90(i,1) = trapz(x(:,1),inner_detector_GEnergy(:,i).*MaxFlux95_90_);
        IntDetector_ERate_InnerCen95__90(i,1) = trapz(x(:,1),inner_detector_GEnergy(:,i).*CenFlux95_90_);
        IntDetector_ERate_InnerApo95__90(i,1) = trapz(x(:,1),inner_detector_GEnergy(:,i).*ApoFlux95_90_);
        
        %Mean
        IntDetector_ERate_InnerMaxMean__90(i,1) = trapz(x(:,1),inner_detector_GEnergy(:,i).*MaxFluxMean_90_);
        IntDetector_ERate_InnerCenMean__90(i,1) = trapz(x(:,1),inner_detector_GEnergy(:,i).*CenFluxMean_90_);
        IntDetector_ERate_InnerApoMean__90(i,1) = trapz(x(:,1),inner_detector_GEnergy(:,i).*ApoFluxMean_90_);
        
        
    end
    
    
    % semilogy(1:1:17,IntDetector_ERate_InnerMax50_Omni,'o',1:1:17,IntDetector_ERate_Inner75_Omni,'o',1:1:17,IntDetector_ERate_Inner95_Omni,'o')
    % xlim([0 18])
    % xticks((0:1:18))
    % title(append('Energy Deposit Rates Inner ',' Eo Range: ',num2str(Cen(Eo_inner)),' to ',num2str(max(Eo_inner)),' MeV'))
    % ylabel('Energy Deposit Rates')
    % xlabel('Detectors')
    % %hold off
    % effsave = append(date(),' Energy Deposit Rate',' Eo ',num2str(Cen(Eo_inner)),' to ',num2str(max(Eo_inner)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
    % saveas(gcf,effsave)
    
    
    %Plots Graph for each energy channel
    ymax = round(max(G_eff_dE_inner)*1.1,3);
    for c=1:length(energy_channels)
        
        f = figure;
        f.Position = [100 100 720 720];
        hold on
        for i = 1:length(Eo_inner)
            plot(x(:,1),G_E_eff_inner(:,i,c),'Color',Eo_color_inner(i,:),'DisplayName',BowTieLegend_inner(1,i));
        end
        
        
        plot(E_eff_inner(c),G_eff_dE_inner(c),'o','Linewidth',2,'DisplayName',BowTieLegend_inner(1,end),'Color','b')
        ylim([0 ymax])
        yticks((0:0.005:ymax))
        plot([0 E_eff_inner(c)],[G_eff_dE_inner(c) G_eff_dE_inner(c)],'--k')
        plot([E_eff_inner(c) E_eff_inner(c)],[0 G_eff_dE_inner(c)],'--k','MarkerSize',2)
        labelpoints(E_eff_inner(c),G_eff_dE_inner(c),append('E_Eff = ',num2str(E_eff_inner(c)),' MeV'),'SE',0.15)
        %labelpoints(E_eff_inner(c),G_eff_dE_inner(c),append('G_eff_dE = ',num2str(G_eff_dE_inner(c)),' cm^2 sr MeV'),'SE',0.65)
        
        legend(BowTieLegend_inner,'Location', 'southoutside','NumColumns',round(length(BowTieLegend_inner)/2))
        xlim([1 7])
        xticks((1:1:7))
        title(append('Bow Tie Analysis Inner Energy Channel ',num2str(c),' Deposited Energy Range: ',num2str(energy_channels(c,1)),' to ',num2str(energy_channels(c,2)),' MeV'))
        ylabel('G_{eff} * \Delta E ')
        xlabel('Nominal Energy (MeV)')
        hold off
        
        effsave = append(date(),'Bow Tie Inner Energy Channel ',num2str(c),' Eo ',num2str(min(Eo_inner)),' to ',num2str(max(Eo_inner)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
        saveas(gcf,effsave)
        
    end
    
    %Plots graph of all the average intersection points of each energy channel
    f = figure;
    f.Position = [100 100 1200 720];
    hold on
    for c=1:length(energy_channels)
        plot(E_eff_inner(c),G_eff_dE_inner(c),'o','Color',Effplotcolor(c,:))
    end
    legend(EngLegend,'Location', 'southoutside','NumColumns',8)
    ylim([0 ymax])
    yticks((0:0.005:ymax))
    
    xlim([1 7])
    xticks((1:1:7))
    title(append('Bow Tie Analysis Inner Configuration All Energy Channels ',' Eo Range: ',num2str(min(Eo_inner)),' to ',num2str(max(Eo_inner)),' MeV'))
    ylabel('G_{eff} * \Delta E ')
    xlabel('Nominal Effective Energy (MeV)')
    hold off
    
    effsave = append(date(),' Bow Tie Inner  All Energy Channels',' Eo',num2str(min(Eo_inner)),' to ',num2str(max(Eo_inner)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
    saveas(gcf,effsave)
    
    f = figure;
    f.Position = [100 100 1600 720];
    hold on
    %Plots each energy channel FWHM value in the same color as
    %the Eff. Curve
    for c = 1:length(energy_channels)
        bar(E_eff_inner(c),G_eff_dE_inner(c)/fwhm_inner(c),fwhm_inner(c),'EdgeColor','k','FaceColor',Effplotcolor(c,:))
        
    end
    legend(EngLegend,'Location', 'southoutside','NumColumns',8)
    title(append('Energy Channel Bins-Inner ',' Eo Range: ',num2str(min(Eo_inner)),' to ',num2str(max(Eo_inner)),' MeV'),'FontSize',20)
    ylabel('Effective Geometric Factor','FontSize',20)
    xlabel('Nominal Effective Energy (MeV)','FontSize',20)
    hold off
    
    effsave = append(date(),' Energy Channel Bins-Inner',' Eo ',num2str(min(Eo_inner)),' to ',num2str(max(Eo_inner)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
    saveas(gcf,effsave)
    
    
    % Count Rate Plot
    energy_channel_list = 1:1:length(energy_channels);
    %count_rate_colors = magma(length(FluxRange));
    count_rate_colors = magma(1);
    
    f = figure;
    f.Position = [100 100 1200 720];
    %hold on
    %Plots each energy channel FWHM value in the same color as
    %the Eff. Curve
    
    semilogy(energy_channel_list,BowCount_Rate_InnerMax50_Omni,'o',energy_channel_list,BowCount_Rate_InnerMax75_Omni,'o',energy_channel_list,BowCount_Rate_InnerMax95_Omni,'o')
    
    
    xlim([0 length(energy_channels)+1])
    xticks((0:1:length(energy_channels)+1))
    
    
    title(append('Count Rates Inner ',' Eo Range: ',num2str(min(Eo_inner)),' to ',num2str(max(Eo_inner)),' MeV'))
    ylabel('Count Rates')
    xlabel('Energy Channels')
    %hold off
    
    
    
    effsave = append(date(),' Count Rate -Inner',' Eo ',num2str(min(Eo_inner)),' to ',num2str(max(Eo_inner)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
    saveas(gcf,effsave)
    
    %Returns to main directory
    cd ..
    cd ..
    
    %Exports Bin Info to Excel File
    exportname = append(date(),' Inner- Center Bin Width Table ',' Eo ',num2str(min(Eo_inner)),' to ',num2str(max(Eo_inner)),'_',addin,num2str(length(energy_channels)));
    Energy_Resolution_inner= 100*BinWidth_Inner./E_eff_inner;
    BuildBinTable(energy_channels,E_eff_inner,BinWidth_Inner,Geff_Inner,Energy_Resolution_inner,IntEnergyChannel_InnerCen25__90,IntEnergyChannel_InnerCen50__90,IntEnergyChannel_InnerCen75__90,exportname)

    
    exportname = append(date(),' Inner- Apo Count Rate Table ',' Eo ',num2str(min(Eo_inner)),' to ',num2str(max(Eo_inner)),'_',addin,num2str(length(energy_channels)));
    BuildBinTable(energy_channels,E_eff_inner,BinWidth_Inner,Geff_Inner,Energy_Resolution_inner,IntEnergyChannel_InnerApo25__90,IntEnergyChannel_InnerApo50__90,IntEnergyChannel_InnerApo75__90,exportname)

    %BuildCountTable(energy_channels,Count_Rate_Inner,exportname)
    
    % %% Intersection Plot
    % hold on
    % num = 0;
    % for j = 1:(length(Eo_inner)-1)
    %
    %     for k = 1:(length(Eo_inner)-j)
    %         plot(xi_inner(num+k,16),yi_inner(num+k,16),'o','Color',Eo_color_inner(j,:))
    %
    %     end
    %     num = num + k;
    % end
    % plot(E_eff_inner(16),G_eff_dE_inner(16),'or')
    % hold off
    
    
    
    
    % End of Inner Configuration
end



%% Bow Tie Analysis-Selesnick/Blake-Whole
%Changes directory for jpg of graphs
cd 'Bow Tie\Whole'

%Sets Ei and Range of Eo
%Ei_Whole is the same as X from the GEANT4 results
Ei_Whole = x(:,1);
%Eo_Whole set by user.
Eo_Whole = 0.2:0.2:2.0;


%Sets up color vectors for plotting the different Eo curves
Eo_color_Whole = magma(length(Eo_Whole)+1);

%Preallocates all variables prior to For Loops
J_e_Whole = zeros(length(x),length(Eo_Whole));

G_int_Whole = zeros(length(x),length(Eo_Whole),length(energy_channels));
G_term_Whole = zeros(length(energy_channels),length(Eo_Whole));
G_E_eff_Whole = zeros(length(x),length(Eo_Whole),length(energy_channels));

xi_Whole = zeros(sum(1:length(Eo_Whole)-1),length(energy_channels));
yi_Whole = xi_Whole;

E_eff_Whole = zeros(1,length(energy_channels));
G_eff_dE_Whole= E_eff_Whole;
BowTieLegend_Whole = strings([1,length(Eo_Whole)]);

BowCount_Rate_WholeMax50_ = zeros(length(energy_channels),1);
BowCount_Rate_WholeMax75_ = zeros(length(energy_channels),1);
BowCount_Rate_WholeMax95_ = zeros(length(energy_channels),1);

%Count_Rate_Whole =  zeros(length(energy_channels),length(Eo_Whole));

BinWidth_Whole = zeros(1, length(energy_channels));
Geff_Whole = zeros(1, length(energy_channels));

%Creates J(e) and creates String Array for Plot Legends
for i = 1:length(Eo_Whole)
    J_e_Whole(:,i) = exp(-Ei_Whole/Eo_Whole(i));
    BowTieLegend_Whole(i) = num2str(Eo_Whole(i));
    
end

%Adds 'Average Intersection Point' for plotting
BowTieLegend_Whole = [BowTieLegend_Whole,'Intersection Point','Average Intersection Point'];

%Creates J(e)^-1
J_e_inv_Whole = 1./J_e_Whole;

%Preallocation
BowCount_Rate_WholeMax25_Omni = zeros(length(energy_channels),1);
BowCount_Rate_WholeCen25_Omni = zeros(length(energy_channels),1);
BowCount_Rate_WholeApo25_Omni = zeros(length(energy_channels),1);

%50thPercentile
BowCount_Rate_WholeMax50_Omni = zeros(length(energy_channels),1);
BowCount_Rate_WholeCen50_Omni = zeros(length(energy_channels),1);
BowCount_Rate_WholeApo50_Omni = zeros(length(energy_channels),1);

%75thPercentile
BowCount_Rate_WholeMax75_Omni =zeros(length(energy_channels),1);
BowCount_Rate_WholeCen75_Omni = zeros(length(energy_channels),1);
BowCount_Rate_WholeApo75_Omni = zeros(length(energy_channels),1);

%95th Percentile
BowCount_Rate_WholeMax95_Omni = zeros(length(energy_channels),1);
BowCount_Rate_WholeCen95_Omni = zeros(length(energy_channels),1);
BowCount_Rate_WholeApo95_Omni = zeros(length(energy_channels),1);

%Average
BowCount_Rate_WholeMaxMean_Omni = zeros(length(energy_channels),1);
BowCount_Rate_WholeCenMean_Omni = zeros(length(energy_channels),1);
BowCount_Rate_WholeApoMean_Omni = zeros(length(energy_channels),1);

BowCount_Rate_WholeMax25__90 = zeros(length(energy_channels),1);
BowCount_Rate_WholeCen25__90 = zeros(length(energy_channels),1);
BowCount_Rate_WholeApo25__90 = zeros(length(energy_channels),1);

%50thPercentile
BowCount_Rate_WholeMax50__90 = zeros(length(energy_channels),1);
BowCount_Rate_WholeCen50__90 = zeros(length(energy_channels),1);
BowCount_Rate_WholeApo50__90 = zeros(length(energy_channels),1);

%75thPercentile
BowCount_Rate_WholeMax75__90 = zeros(length(energy_channels),1);
BowCount_Rate_WholeCen75__90 = zeros(length(energy_channels),1);
BowCount_Rate_WholeApo75__90 = zeros(length(energy_channels),1);

%95th Percentile
BowCount_Rate_WholeMax95__90 = zeros(length(energy_channels),1);
BowCount_Rate_WholeCen95__90 = zeros(length(energy_channels),1);
BowCount_Rate_WholeApo95__90 = zeros(length(energy_channels),1);

%Average
BowCount_Rate_WholeMaxMean__90 = zeros(length(energy_channels),1);
BowCount_Rate_WholeCenMean__90 = zeros(length(energy_channels),1);
BowCount_Rate_WholeApoMean__90 = zeros(length(energy_channels),1);


%For Loop for calculating a line for each Eo and finding the average intersection point
for c=1:length(energy_channels)
    %Prints energy channel to show progress
    fprintf('Energy Channel-Whole # %.2f \n',c)
    
    %Calculates line for each Eo
    for i = 1:length(Eo_Whole)
        G_int_Whole(:,i,c) = geo_EC(:,c).*J_e_Whole(:,i);
        G_term_Whole(i,c) = trapz(Ei_Whole,G_int_Whole(:,i,c));
        G_E_eff_Whole(:,i,c)= G_term_Whole(i,c)*J_e_inv_Whole(:,i);
        
    end
    
    %Find Intersections of Eo Lines
    num = 0;
    for j = 1:(length(Eo_Whole)-1)
        
        for k = 1:(length(Eo_Whole)-j)
            [xi_Whole(num+k,c),yi_Whole(num+k,c)] = polyxpoly(x(:,1),G_E_eff_Whole(:,j,c),x(:,1),G_E_eff_Whole(:,j+k,c),'unique');
            
        end
        num = num + k;
    end
    
    %Finds average intersection point
    E_eff_Whole(c) = mean(xi_Whole(:,c));
    G_eff_dE_Whole(c) = mean(yi_Whole(:,c));
    
    %Flux Rate from AE9 Model
    [MaxFlux25Omni,CenFlux25Omni,ApoFlux25Omni,MaxFlux50Omni,CenFlux50Omni,ApoFlux50Omni,MaxFlux75Omni,CenFlux75Omni,ApoFlux75Omni,MaxFlux95Omni,CenFlux95Omni,ApoFlux95Omni,MaxFluxMeanOmni,CenFluxMeanOmni,ApoFluxMeanOmni] = FluxRateOmni(E_eff_Whole(c));
    
    %Calculates Estimated Count Rate for each Eo
    %25thPercentile
    BowCount_Rate_WholeMax25_Omni(c) = G_eff_dE_Whole(c)*MaxFlux25Omni;
    BowCount_Rate_WholeCen25_Omni(c) = G_eff_dE_Whole(c)*CenFlux25Omni;
    BowCount_Rate_WholeApo25_Omni(c) = G_eff_dE_Whole(c)*ApoFlux25Omni;
    
    %50thPercentile
    BowCount_Rate_WholeMax50_Omni(c) = G_eff_dE_Whole(c)*MaxFlux50Omni;
    BowCount_Rate_WholeCen50_Omni(c) = G_eff_dE_Whole(c)*CenFlux50Omni;
    BowCount_Rate_WholeApo50_Omni(c) = G_eff_dE_Whole(c)*ApoFlux50Omni;
    
    %75thPercentile
    BowCount_Rate_WholeMax75_Omni(c) = G_eff_dE_Whole(c)*MaxFlux75Omni;
    BowCount_Rate_WholeCen75_Omni(c) = G_eff_dE_Whole(c)*CenFlux75Omni;
    BowCount_Rate_WholeApo75_Omni(c) = G_eff_dE_Whole(c)*ApoFlux75Omni;
    
    %95th Percentile
    BowCount_Rate_WholeMax95_Omni(c) = G_eff_dE_Whole(c)*MaxFlux95Omni;
    BowCount_Rate_WholeCen95_Omni(c) = G_eff_dE_Whole(c)*CenFlux95Omni;
    BowCount_Rate_WholeApo95_Omni(c) = G_eff_dE_Whole(c)*ApoFlux95Omni;
    
    %Average
    BowCount_Rate_WholeMaxMean_Omni(c) = G_eff_dE_Whole(c)*MaxFluxMeanOmni;
    BowCount_Rate_WholeCenMean_Omni(c) = G_eff_dE_Whole(c)*CenFluxMeanOmni;
    BowCount_Rate_WholeApoMean_Omni(c) = G_eff_dE_Whole(c)*ApoFluxMeanOmni;
    
    %Flux Rate from AE9 Model
    [MaxFlux25_90,CenFlux25_90,ApoFlux25_90,MaxFlux50_90,CenFlux50_90,ApoFlux50_90,MaxFlux75_90,CenFlux75_90,ApoFlux75_90,MaxFlux95_90,CenFlux95_90,ApoFlux95_90,MaxFluxMean_90,CenFluxMean_90,ApoFluxMean_90] = FluxRate90(E_eff_Whole(c));
    
    %Calculates Estimated Count Rate for each Eo
    %25thPercentile
    BowCount_Rate_WholeMax25__90(c) = G_eff_dE_Whole(c)*MaxFlux25_90;
    BowCount_Rate_WholeCen25__90(c) = G_eff_dE_Whole(c)*CenFlux25_90;
    BowCount_Rate_WholeApo25__90(c) = G_eff_dE_Whole(c)*ApoFlux25_90;
    
    %50thPercentile
    BowCount_Rate_WholeMax50__90(c) = G_eff_dE_Whole(c)*MaxFlux50_90;
    BowCount_Rate_WholeCen50__90(c) = G_eff_dE_Whole(c)*CenFlux50_90;
    BowCount_Rate_WholeApo50__90(c) = G_eff_dE_Whole(c)*ApoFlux50_90;
    
    %75thPercentile
    BowCount_Rate_WholeMax75__90(c) = G_eff_dE_Whole(c)*MaxFlux75_90;
    BowCount_Rate_WholeCen75__90(c) = G_eff_dE_Whole(c)*CenFlux75_90;
    BowCount_Rate_WholeApo75__90(c) = G_eff_dE_Whole(c)*ApoFlux75_90;
    
    %95th Percentile
    BowCount_Rate_WholeMax95__90(c) = G_eff_dE_Whole(c)*MaxFlux95_90;
    BowCount_Rate_WholeCen95__90(c) = G_eff_dE_Whole(c)*CenFlux95_90;
    BowCount_Rate_WholeApo95__90(c) = G_eff_dE_Whole(c)*ApoFlux95_90;
    
    %Average
    BowCount_Rate_WholeMaxMean__90(c) = G_eff_dE_Whole(c)*MaxFluxMean_90;
    BowCount_Rate_WholeCenMean__90(c) = G_eff_dE_Whole(c)*CenFluxMean_90;
    BowCount_Rate_WholeApoMean__90(c) = G_eff_dE_Whole(c)*ApoFluxMean_90;
    
    
    %Calculate Bin Characteristics
    Geff_Whole(c) = G_eff_dE_Whole(c)/fwhm_whole(c);
    BinWidth_Whole(c) = fwhm_whole(c);
    
end

%% Count Rates using Omnidirectional divied by 4*pi
%25thPercentile
BowDetector1_Rate_WholeMax25_Omni = sum(BowCount_Rate_WholeMax25_Omni)
BowDetector1_Rate_WholeCen25_Omni = sum(BowCount_Rate_WholeCen25_Omni)
BowDetector1_Rate_WholeApo25_Omni = sum(BowCount_Rate_WholeApo25_Omni)

%50thPercentile
BowDetector1_Rate_WholeMax50_Omni = sum(BowCount_Rate_WholeMax50_Omni)
BowDetector1_Rate_WholeCen50_Omni = sum(BowCount_Rate_WholeCen50_Omni)
BowDetector1_Rate_WholeApo50_Omni = sum(BowCount_Rate_WholeApo50_Omni)

%75thPercentile
BowDetector1_Rate_WholeMax75_Omni = sum(BowCount_Rate_WholeMax75_Omni)
BowDetector1_Rate_WholeCen75_Omni = sum(BowCount_Rate_WholeCen75_Omni)
BowDetector1_Rate_WholeApo75_Omni = sum(BowCount_Rate_WholeApo75_Omni)

%95thPercentile
BowDetector1_Rate_WholeMax95_Omni = sum(BowCount_Rate_WholeMax95_Omni)
BowDetector1_Rate_WholeCen95_Omni = sum(BowCount_Rate_WholeCen95_Omni)
BowDetector1_Rate_WholeApo95_Omni = sum(BowCount_Rate_WholeApo95_Omni)

%Mean
BowDetector1_Rate_WholeMaxMean_Omni = sum(BowCount_Rate_WholeMaxMean_Omni)
BowDetector1_Rate_WholeCenMean_Omni = sum(BowCount_Rate_WholeCenMean_Omni)
BowDetector1_Rate_WholeApoMean_Omni = sum(BowCount_Rate_WholeApoMean_Omni)


%Integral Method
%25th Percentile
IntDetector1_BinRate_WholeMax25_Omni = trapz(x(:,1),geo_EL(:).*MaxFlux25Omni_)
IntDetector1_BinRate_WholeCen25_Omni = trapz(x(:,1),geo_EL(:).*CenFlux25Omni_)
IntDetector1_BinRate_WholeApo25_Omni = trapz(x(:,1),geo_EL(:).*ApoFlux25Omni_)

%50th Percentile
IntDetector1_BinRate_WholeMax50_Omni = trapz(x(:,1),geo_EL(:).*MaxFlux50Omni_)
IntDetector1_BinRate_WholeCen50_Omni = trapz(x(:,1),geo_EL(:).*CenFlux50Omni_)
IntDetector1_BinRate_WholeApo50_Omni = trapz(x(:,1),geo_EL(:).*ApoFlux50Omni_)

%75th Percentile
IntDetector1_BinRate_WholeMax75_Omni = trapz(x(:,1),geo_EL(:).*MaxFlux75Omni_)
IntDetector1_BinRate_WholeCen75_Omni = trapz(x(:,1),geo_EL(:).*CenFlux75Omni_)
IntDetector1_BinRate_WholeApo75_Omni = trapz(x(:,1),geo_EL(:).*ApoFlux75Omni_)

%95th Percentile
IntDetector1_BinRate_WholeMax95_Omni = trapz(x(:,1),geo_EL(:).*MaxFlux95Omni_)
IntDetector1_BinRate_WholeCen95_Omni = trapz(x(:,1),geo_EL(:).*CenFlux95Omni_)
IntDetector1_BinRate_WholeApo95_Omni = trapz(x(:,1),geo_EL(:).*ApoFlux95Omni_)

%Mean
IntDetector1_BinRate_WholeMaxMean_Omni = trapz(x(:,1),geo_EL(:).*MaxFluxMeanOmni_)
IntDetector1_BinRate_WholeCenMean_Omni = trapz(x(:,1),geo_EL(:).*CenFluxMeanOmni_)
IntDetector1_BinRate_WholeApoMean_Omni = trapz(x(:,1),geo_EL(:).*ApoFluxMeanOmni_)


%Preallocate Variables

%25th Percentile
    IntDetector_AllRate_WholeMax25_Omni = zeros(width(whole_detector_GAllCounts),1);
    IntDetector_AllRate_WholeCen25_Omni = zeros(width(whole_detector_GAllCounts),1);
    IntDetector_AllRate_WholeApo25_Omni = zeros(width(whole_detector_GAllCounts),1);
    
    %50th Percentile
    IntDetector_AllRate_WholeMax50_Omni = zeros(width(whole_detector_GAllCounts),1);
    IntDetector_AllRate_WholeCen50_Omni = zeros(width(whole_detector_GAllCounts),1);
    IntDetector_AllRate_WholeApo50_Omni = zeros(width(whole_detector_GAllCounts),1);
    
    %75th Percentile
    IntDetector_AllRate_WholeMax75_Omni =zeros(width(whole_detector_GAllCounts),1);
    IntDetector_AllRate_WholeCen75_Omni = zeros(width(whole_detector_GAllCounts),1);
    IntDetector_AllRate_WholeApo75_Omni = zeros(width(whole_detector_GAllCounts),1);
    
    %95th Percentile
    IntDetector_AllRate_WholeMax95_Omni = zeros(width(whole_detector_GAllCounts),1);
    IntDetector_AllRate_WholeCen95_Omni = zeros(width(whole_detector_GAllCounts),1);
    IntDetector_AllRate_WholeApo95_Omni = zeros(width(whole_detector_GAllCounts),1);
    
    %Mean
    IntDetector_AllRate_WholeMaxMean_Omni = zeros(width(whole_detector_GAllCounts),1);
    IntDetector_AllRate_WholeCenMean_Omni =zeros(width(whole_detector_GAllCounts),1);
    IntDetector_AllRate_WholeApoMean_Omni =zeros(width(whole_detector_GAllCounts),1);
    

for i = 1:width(whole_detector_GAllCounts)
    
    %25th Percentile
    IntDetector_AllRate_WholeMax25_Omni(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*MaxFlux25Omni_)';
    IntDetector_AllRate_WholeCen25_Omni(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*CenFlux25Omni_)';
    IntDetector_AllRate_WholeApo25_Omni(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*ApoFlux25Omni_)';
    
    %50th Percentile
    IntDetector_AllRate_WholeMax50_Omni(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*MaxFlux50Omni_)';
    IntDetector_AllRate_WholeCen50_Omni(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*CenFlux50Omni_)';
    IntDetector_AllRate_WholeApo50_Omni(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*ApoFlux50Omni_)';
    
    %75th Percentile
    IntDetector_AllRate_WholeMax75_Omni(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*MaxFlux75Omni_)';
    IntDetector_AllRate_WholeCen75_Omni(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*CenFlux75Omni_)';
    IntDetector_AllRate_WholeApo75_Omni(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*ApoFlux75Omni_)';
    
    %95th Percentile
    IntDetector_AllRate_WholeMax95_Omni(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*MaxFlux95Omni_)';
    IntDetector_AllRate_WholeCen95_Omni(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*CenFlux95Omni_)';
    IntDetector_AllRate_WholeApo95_Omni(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*ApoFlux95Omni_)';
    
    %Mean
    IntDetector_AllRate_WholeMaxMean_Omni(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*MaxFluxMeanOmni_)';
    IntDetector_AllRate_WholeCenMean_Omni(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*CenFluxMeanOmni_)';
    IntDetector_AllRate_WholeApoMean_Omni(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*ApoFluxMeanOmni_)';
    
end

for i = 1:length(energy_channels)
    %25th Percentile
    IntEnergyChannel_WholeMax25_Omni(i,1)=trapz(x(:,1),geo_EC(:,i).*MaxFlux25Omni_);
    IntEnergyChannel_WholeCen25_Omni(i,1)=trapz(x(:,1),geo_EC(:,i).*CenFlux25Omni_);
    IntEnergyChannel_WholeApo25_Omni(i,1)=trapz(x(:,1),geo_EC(:,i).*ApoFlux25Omni_);
    
    %50th Percentile
    IntEnergyChannel_WholeMax50_Omni(i,1)=trapz(x(:,1),geo_EC(:,i).*MaxFlux50Omni_);
    IntEnergyChannel_WholeCen50_Omni(i,1)=trapz(x(:,1),geo_EC(:,i).*CenFlux50Omni_);
    IntEnergyChannel_WholeApo50_Omni(i,1)=trapz(x(:,1),geo_EC(:,i).*ApoFlux50Omni_);
    
    %75th Percentile
    IntEnergyChannel_WholeMax75_Omni(i,1)=trapz(x(:,1),geo_EC(:,i).*MaxFlux75Omni_);
    IntEnergyChannel_WholeCen75_Omni(i,1)=trapz(x(:,1),geo_EC(:,i).*CenFlux75Omni_);
    IntEnergyChannel_WholeApo75_Omni(i,1)=trapz(x(:,1),geo_EC(:,i).*ApoFlux75Omni_);
    
    %95th Percentile
    IntEnergyChannel_WholeMax95_Omni(i,1)=trapz(x(:,1),geo_EC(:,i).*MaxFlux95Omni_);
    IntEnergyChannel_WholeCen95_Omni(i,1)=trapz(x(:,1),geo_EC(:,i).*CenFlux95Omni_);
    IntEnergyChannel_WholeApo95_Omni(i,1)=trapz(x(:,1),geo_EC(:,i).*ApoFlux95Omni_);
    
    %Mean
    IntEnergyChannel_WholeMaxMean_Omni(i,1)=trapz(x(:,1),geo_EC(:,i).*MaxFluxMeanOmni_);
    IntEnergyChannel_WholeCenMean_Omni(i,1)=trapz(x(:,1),geo_EC(:,i).*CenFluxMeanOmni_);
    IntEnergyChannel_WholeApoMean_Omni(i,1)=trapz(x(:,1),geo_EC(:,i).*ApoFluxMeanOmni_);
    
end

IntDetector_ERate_Whole50_ = zeros(width(whole_detector_GEnergy),1);
IntDetector_ERate_Whole75_ = zeros(width(whole_detector_GEnergy),1);
IntDetector_ERate_Whole95_ = zeros(width(whole_detector_GEnergy),1);

IntDetector_ERate_WholeCen50_ = zeros(width(whole_detector_GEnergy),1);

for i = 1:width(whole_detector_GEnergy)
    %25th Percentile
    IntDetector_ERate_WholeMax25_Omni(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*MaxFlux25Omni_);
    IntDetector_ERate_WholeCen25_Omni(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*CenFlux25Omni_);
    IntDetector_ERate_WholeApo25_Omni(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*ApoFlux25Omni_);
    
    %50th Percentile
    IntDetector_ERate_WholeMax50_Omni(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*MaxFlux50Omni_);
    IntDetector_ERate_WholeCen50_Omni(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*CenFlux50Omni_);
    IntDetector_ERate_WholeApo50_Omni(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*ApoFlux50Omni_);
    
    %75th Percentile
    IntDetector_ERate_WholeMax75_Omni(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*MaxFlux75Omni_);
    IntDetector_ERate_WholeCen75_Omni(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*CenFlux75Omni_);
    IntDetector_ERate_WholeApo75_Omni(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*ApoFlux75Omni_);
    
    %95th Percentile
    IntDetector_ERate_WholeMax95_Omni(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*MaxFlux95Omni_);
    IntDetector_ERate_WholeCen95_Omni(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*CenFlux95Omni_);
    IntDetector_ERate_WholeApo95_Omni(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*ApoFlux95Omni_);
    
    %Mean
    IntDetector_ERate_WholeMaxMean_Omni(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*MaxFluxMeanOmni_);
    IntDetector_ERate_WholeCenMean_Omni(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*CenFluxMeanOmni_);
    IntDetector_ERate_WholeApoMean_Omni(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*ApoFluxMeanOmni_);
    
    
end

%% Count Rates using 90 deg pitch angle
%25thPercentile
BowDetector1_Rate_WholeMax25__90 = sum(BowCount_Rate_WholeMax25__90)
BowDetector1_Rate_WholeCen25__90 = sum(BowCount_Rate_WholeCen25__90)
BowDetector1_Rate_WholeApo25__90 = sum(BowCount_Rate_WholeApo25__90)

%50thPercentile
BowDetector1_Rate_WholeMax50__90 = sum(BowCount_Rate_WholeMax50__90)
BowDetector1_Rate_WholeCen50__90 = sum(BowCount_Rate_WholeCen50__90)
BowDetector1_Rate_WholeApo50__90 = sum(BowCount_Rate_WholeApo50__90)

%75thPercentile
BowDetector1_Rate_WholeMax75__90 = sum(BowCount_Rate_WholeMax75__90)
BowDetector1_Rate_WholeCen75__90 = sum(BowCount_Rate_WholeCen75__90)
BowDetector1_Rate_WholeApo75__90 = sum(BowCount_Rate_WholeApo75__90)

%95thPercentile
BowDetector1_Rate_WholeMax95__90 = sum(BowCount_Rate_WholeMax95__90)
BowDetector1_Rate_WholeCen95__90 = sum(BowCount_Rate_WholeCen95__90)
BowDetector1_Rate_WholeApo95__90 = sum(BowCount_Rate_WholeApo95__90)

%Mean
BowDetector1_Rate_WholeMaxMean__90 = sum(BowCount_Rate_WholeMaxMean__90)
BowDetector1_Rate_WholeCenMean__90 = sum(BowCount_Rate_WholeCenMean__90)
BowDetector1_Rate_WholeApoMean__90 = sum(BowCount_Rate_WholeApoMean__90)


%Integral Method
%25th Percentile
IntDetector1_BinRate_WholeMax25__90 = trapz(x(:,1),geo_EL(:).*MaxFlux25_90_)
IntDetector1_BinRate_WholeCen25__90 = trapz(x(:,1),geo_EL(:).*CenFlux25_90_)
IntDetector1_BinRate_WholeApo25__90 = trapz(x(:,1),geo_EL(:).*ApoFlux25_90_)

%50th Percentile
IntDetector1_BinRate_WholeMax50__90 = trapz(x(:,1),geo_EL(:).*MaxFlux50_90_)
IntDetector1_BinRate_WholeCen50__90 = trapz(x(:,1),geo_EL(:).*CenFlux50_90_)
IntDetector1_BinRate_WholeApo50__90 = trapz(x(:,1),geo_EL(:).*ApoFlux50_90_)

%75th Percentile
IntDetector1_BinRate_WholeMax75__90 = trapz(x(:,1),geo_EL(:).*MaxFlux75_90_)
IntDetector1_BinRate_WholeCen75__90 = trapz(x(:,1),geo_EL(:).*CenFlux75_90_)
IntDetector1_BinRate_WholeApo75__90 = trapz(x(:,1),geo_EL(:).*ApoFlux75_90_)

%95th Percentile
IntDetector1_BinRate_WholeMax95__90 = trapz(x(:,1),geo_EL(:).*MaxFlux95_90_)
IntDetector1_BinRate_WholeCen95__90 = trapz(x(:,1),geo_EL(:).*CenFlux95_90_)
IntDetector1_BinRate_WholeApo95__90 = trapz(x(:,1),geo_EL(:).*ApoFlux95_90_)

%Mean
IntDetector1_BinRate_WholeMaxMean__90 = trapz(x(:,1),geo_EL(:).*MaxFluxMean_90_)
IntDetector1_BinRate_WholeCenMean__90 = trapz(x(:,1),geo_EL(:).*CenFluxMean_90_)
IntDetector1_BinRate_WholeApoMean__90 = trapz(x(:,1),geo_EL(:).*ApoFluxMean_90_)


%Preallocate Variables
IntDetector_AllRate_Whole50 = zeros(width(whole_detector_GAllCounts),1);
IntDetector_AllRate_Whole75 = zeros(width(whole_detector_GAllCounts),1);
IntDetector_AllRate_Whole95 = zeros(width(whole_detector_GAllCounts),1);
IntDetector_AllRate_WholeCen50 = zeros(width(whole_detector_GAllCounts),1);
IntDetector_AllRate_WholeApo50 = zeros(width(whole_detector_GAllCounts),1);

for i = 1:width(whole_detector_GAllCounts)
    
    %25th Percentile
    IntDetector_AllRate_WholeMax25__90(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*MaxFlux25_90_);
    IntDetector_AllRate_WholeCen25__90(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*CenFlux25_90_);
    IntDetector_AllRate_WholeApo25__90(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*ApoFlux25_90_);
    
    %50th Percentile
    IntDetector_AllRate_WholeMax50__90(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*MaxFlux50_90_);
    IntDetector_AllRate_WholeCen50__90(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*CenFlux50_90_);
    IntDetector_AllRate_WholeApo50__90(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*ApoFlux50_90_);
    
    %75th Percentile
    IntDetector_AllRate_WholeMax75__90(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*MaxFlux75_90_);
    IntDetector_AllRate_WholeCen75__90(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*CenFlux75_90_);
    IntDetector_AllRate_WholeApo75__90(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*ApoFlux75_90_);
    
    %95th Percentile
    IntDetector_AllRate_WholeMax95__90(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*MaxFlux95_90_);
    IntDetector_AllRate_WholeCen95__90(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*CenFlux95_90_);
    IntDetector_AllRate_WholeApo95__90(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*ApoFlux95_90_);
    
    %Mean
    IntDetector_AllRate_WholeMaxMean__90(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*MaxFluxMean_90_);
    IntDetector_AllRate_WholeCenMean__90(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*CenFluxMean_90_);
    IntDetector_AllRate_WholeApoMean__90(i,1) = trapz(x(:,1),whole_detector_GAllCounts(:,i).*ApoFluxMean_90_);
    
end

for i = 1:length(energy_channels)
    %25th Percentile
    IntEnergyChannel_WholeMax25__90(i,1)=trapz(x(:,1),geo_EC(:,i).*MaxFlux25_90_);
    IntEnergyChannel_WholeCen25__90(i,1)=trapz(x(:,1),geo_EC(:,i).*CenFlux25_90_);
    IntEnergyChannel_WholeApo25__90(i,1)=trapz(x(:,1),geo_EC(:,i).*ApoFlux25_90_);
    
    %50th Percentile
    IntEnergyChannel_WholeMax50__90(i,1)=trapz(x(:,1),geo_EC(:,i).*MaxFlux50_90_);
    IntEnergyChannel_WholeCen50__90(i,1)=trapz(x(:,1),geo_EC(:,i).*CenFlux50_90_);
    IntEnergyChannel_WholeApo50__90(i,1)=trapz(x(:,1),geo_EC(:,i).*ApoFlux50_90_);
    
    %75th Percentile
    IntEnergyChannel_WholeMax75__90(i,1)=trapz(x(:,1),geo_EC(:,i).*MaxFlux75_90_);
    IntEnergyChannel_WholeCen75__90(i,1)=trapz(x(:,1),geo_EC(:,i).*CenFlux75_90_);
    IntEnergyChannel_WholeApo75__90(i,1)=trapz(x(:,1),geo_EC(:,i).*ApoFlux75_90_);
    
    %95th Percentile
    IntEnergyChannel_WholeMax95__90(i,1)=trapz(x(:,1),geo_EC(:,i).*MaxFlux95_90_);
    IntEnergyChannel_WholeCen95__90(i,1)=trapz(x(:,1),geo_EC(:,i).*CenFlux95_90_);
    IntEnergyChannel_WholeApo95__90(i,1)=trapz(x(:,1),geo_EC(:,i).*ApoFlux95_90_);
    
    %Mean
    IntEnergyChannel_WholeMaxMean__90(i,1)=trapz(x(:,1),geo_EC(:,i).*MaxFluxMean_90_);
    IntEnergyChannel_WholeCenMean__90(i,1)=trapz(x(:,1),geo_EC(:,i).*CenFluxMean_90_);
    IntEnergyChannel_WholeApoMean__90(i,1)=trapz(x(:,1),geo_EC(:,i).*ApoFluxMean_90_);
    
end

IntDetector_ERate_Whole50_ = zeros(width(whole_detector_GEnergy),1);
IntDetector_ERate_Whole75_ = zeros(width(whole_detector_GEnergy),1);
IntDetector_ERate_Whole95_ = zeros(width(whole_detector_GEnergy),1);

IntDetector_ERate_WholeCen50_ = zeros(width(whole_detector_GEnergy),1);

for i = 1:width(whole_detector_GEnergy)
    %25th Percentile
    IntDetector_ERate_WholeMax25__90(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*MaxFlux25_90_);
    IntDetector_ERate_WholeCen25__90(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*CenFlux25_90_);
    IntDetector_ERate_WholeApo25__90(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*ApoFlux25_90_);
    
    %50th Percentile
    IntDetector_ERate_WholeMax50__90(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*MaxFlux50_90_);
    IntDetector_ERate_WholeCen50__90(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*CenFlux50_90_);
    IntDetector_ERate_WholeApo50__90(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*ApoFlux50_90_);
    
    %75th Percentile
    IntDetector_ERate_WholeMax75__90(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*MaxFlux75_90_);
    IntDetector_ERate_WholeCen75__90(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*CenFlux75_90_);
    IntDetector_ERate_WholeApo75__90(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*ApoFlux75_90_);
    
    %95th Percentile
    IntDetector_ERate_WholeMax95__90(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*MaxFlux95_90_);
    IntDetector_ERate_WholeCen95__90(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*CenFlux95_90_);
    IntDetector_ERate_WholeApo95__90(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*ApoFlux95_90_);
    
    %Mean
    IntDetector_ERate_WholeMaxMean__90(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*MaxFluxMean_90_);
    IntDetector_ERate_WholeCenMean__90(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*CenFluxMean_90_);
    IntDetector_ERate_WholeApoMean__90(i) = trapz(x(:,1),whole_detector_GEnergy(:,i).*ApoFluxMean_90_);
    
    
end


%% Plots for Whole Config


% semilogy(1:1:17,IntDetector_ERate_Whole50_,'o',1:1:17,IntDetector_ERate_Whole75_,'o',1:1:17,IntDetector_ERate_Whole95_,'o')
% xlim([0 18])
% xticks((0:1:18))
% title(append('Energy Deposit Rates Whole ',' Eo Range: ',num2str(Cen(Eo_Whole)),' to ',num2str(max(Eo_Whole)),' MeV'))
% ylabel('Energy Deposit Rates')
% xlabel('Detectors')
% %hold off
% effsave = append(date(),' Energy Deposit Rate',' Eo ',num2str(Cen(Eo_Whole)),' to ',num2str(max(Eo_Whole)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
% saveas(gcf,effsave)

textsize =28;
%Plots Graph for each energy channel
ymax = round(max(G_eff_dE_Whole)*1.1,3);
for c=1:height(energy_channels)
    
    f = figure;
    f.Position = [100 100 1000 720];
    hold on
    for i = 1:length(Eo_Whole)
        plot(x(:,1),G_E_eff_Whole(:,i,c),'Color',Eo_color_Whole(i,:),'DisplayName',BowTieLegend_Whole(1,i),'Linewidth',2);
    end
    
    plot(xi_Whole(:,c),yi_Whole(:,c),'*b','DisplayName',BowTieLegend_Whole(1,end),'MarkerSize',12)
    plot(E_eff_Whole(c),G_eff_dE_Whole(c),'o','Linewidth',2,'DisplayName',BowTieLegend_Whole(1,end-1),'Color','green','MarkerFaceColor', 'green','MarkerSize',12)
    ylim_l = round(min(yi_Whole(:,c))*0.95,3);
    ylim_u = round(max(yi_Whole(:,c))*1.05,3);
    ylim([ylim_l ylim_u])
    plot([0 E_eff_Whole(c)],[G_eff_dE_Whole(c) G_eff_dE_Whole(c)],'--g')
    plot([E_eff_Whole(c) E_eff_Whole(c)],[0 G_eff_dE_Whole(c)],'--g')
    set(gca,'FontSize',18)
    labelpoints(E_eff_Whole(c),G_eff_dE_Whole(c),append('E_Eff = ',num2str(E_eff_Whole(c)),' MeV'),'SE',0.15)
    %labelpoints(E_eff_Whole(c),G_eff_dE_Whole(c),append('G_eff_dE = ',num2str(G_eff_dE_Whole(c)),' cm^2 sr MeV'),'SE',0.65)
    
    legend(BowTieLegend_Whole,'Location', 'southoutside','NumColumns',round(length(BowTieLegend_Whole)/2))
    xlim_l = round(min(xi_Whole(:,c))*0.95,3);
    xlim_u = round(max(xi_Whole(:,c))*1.05,3);
    xlim([xlim_l xlim_u])
    %xticks((1.5:0.05:2.25))
    %title(append('Bow Tie Analysis Whole Energy Channel ',num2str(c),' Deposited Energy Range: ',num2str(energy_channels(c,1)),' to ',num2str(energy_channels(c,2)),' MeV'))
    ylabel('G_{eff} * \Delta E ','FontSize',textsize)
    xlabel('Nominal Energy (MeV)','FontSize',textsize)
    hold off
    
    effsave = append(date(),'Bow Tie Whole Energy Channel ',num2str(c),' Eo ',num2str(min(Eo_Whole)),' to ',num2str(max(Eo_Whole)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
    saveas(gcf,effsave)
    
end
%%
%Plots graph of all the average intersection points of each energy channel
f = figure;
f.Position = [100 100 1200 720];
hold on
for c=1:length(energy_channels)
    plot(E_eff_Whole(c),G_eff_dE_Whole(c),'o','Color',Effplotcolor(c,:))
end
legend(EngLegend,'Location', 'southoutside','NumColumns',8)
ylim([0 ymax])
yticks((0:0.005:ymax))

xlim([1 7])
xticks((1:1:7))

title(append('Bow Tie Analysis Whole Configuration All Energy Channels ',' Eo Range: ',num2str(min(Eo_Whole)),' to ',num2str(max(Eo_Whole)),' MeV'))
ylabel('G_{eff} * \Delta E ')
xlabel('Nominal Effective Energy (MeV)')
hold off

effsave = append(date(),' Bow Tie Whole  All Energy Channels',' Eo',num2str(min(Eo_Whole)),' to ',num2str(max(Eo_Whole)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
saveas(gcf,effsave)

f = figure;
f.Position = [100 100 1600 720];
hold on
%Plots each energy channel FWHM value in the same color as
%the Eff. Curve
for c = 1:length(energy_channels)
    bar(E_eff_Whole(c),G_eff_dE_Whole(c)/fwhm_whole(c),fwhm_whole(c),'EdgeColor','k','FaceColor',Effplotcolor(c,:))
    
end
legend(EngLegend,'Location', 'southoutside','NumColumns',8)
%title(append('Energy Channel Bins-Whole ',' Eo Range: ',num2str(Cen(Eo_Whole)),' to ',num2str(max(Eo_Whole)),' MeV'),'FontSize',20)
ylabel('Effective Geometric Factor','FontSize',28)
xlabel('Nominal Effective Energy (MeV)','FontSize',28)
hold off

effsave = append(date(),' Energy Channel Bins-Whole',' Eo ',num2str(min(Eo_Whole)),' to ',num2str(max(Eo_Whole)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
saveas(gcf,effsave)


% Count Rate Plot
energy_channel_list = 1:1:length(energy_channels);
%count_rate_colors = magma(length(FluxRange));
count_rate_colors = magma(1);

f = figure;
f.Position = [100 100 1200 720];
%hold on
%Plots each energy channel FWHM value in the same color as
%the Eff. Curve

semilogy(energy_channel_list,BowCount_Rate_WholeMax50_Omni,'o',energy_channel_list,BowCount_Rate_WholeMax75_Omni,'o',energy_channel_list,BowCount_Rate_WholeMax95_Omni,'o')


xlim([0 length(energy_channels)+1])
xticks((0:1:length(energy_channels)+1))


title(append('Count Rates Whole ',' Eo Range: ',num2str(min(Eo_Whole)),' to ',num2str(max(Eo_Whole)),' MeV'))
ylabel('Count Rates')
xlabel('Energy Channels')
%hold off



effsave = append(date(),' Count Rate -Whole',' Eo ',num2str(min(Eo_Whole)),' to ',num2str(max(Eo_Whole)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
saveas(gcf,effsave)

%Returns to main directory
cd ..
cd ..
Energy_Resolution= 100*BinWidth_Whole./E_eff_Whole;
%Exports Bin Info to Excel File
exportname = append(date(),' Whole Center- Bin Width Table ',' Eo ',num2str(min(Eo_Whole)),' to ',num2str(max(Eo_Whole)),'_',addin,num2str(length(energy_channels)));
BuildBinTable(energy_channels,E_eff_Whole,BinWidth_Whole,Geff_Whole,Energy_Resolution,IntEnergyChannel_WholeCen25__90,IntEnergyChannel_WholeCen50__90,IntEnergyChannel_WholeCen75__90,exportname)

exportname = append(date(),' Wholen Apo- Bin Width Table ',' Eo ',num2str(min(Eo_Whole)),' to ',num2str(max(Eo_Whole)),'_',addin,num2str(length(energy_channels)));
BuildBinTable(energy_channels,E_eff_Whole,BinWidth_Whole,Geff_Whole,Energy_Resolution,IntEnergyChannel_WholeApo25__90,IntEnergyChannel_WholeApo50__90,IntEnergyChannel_WholeApo75__90,exportname)


exportname = append(date(),' Whole Configuration- Count Rate Table ',' Eo ',num2str(min(Eo_Whole)),' to ',num2str(max(Eo_Whole)),'_',addin,num2str(length(energy_channels)));
%BuildCountTable(energy_channels,Count_Rate_Whole,exportname)


figure
plot(Energy_Resolution,'xb')
title('Energy Resolution per Energy Channel')
xlabel('Energy Channel Number')
ylabel('Energy Resolution (%)')

%% Energy Resolution Plot
Effplotcolor = plasma(length(energy_channel_list));
f = figure;
f.Position = [100 100 1600 720];
textsize = 28;

%Plots each energy channel FWHM value in the same color as
%the Eff. Curve

hold on
for c = 1:length(energy_channels)
%     plot(E_eff_Whole(c),Energy_Resolution(c),'o','MarkerSize',8,...
%         'MarkerEdgeColor',Effplotcolor(c,:),...
%         'MarkerFaceColor',Effplotcolor(c,:))
    
    plot(E_eff_Whole(c),Energy_Resolution(c),'o','MarkerSize',8,...
        'MarkerEdgeColor','b',...
        'MarkerFaceColor','b')
end
set(gca,'FontSize',textsize)
%title('HERT Energy Resolution','FontSize',textsize)
ylabel('Spectral Resolution dE/E(%)','FontSize',textsize)
xlabel('Nominal Energy (MeV)','FontSize',textsize)
ylim([0,40])
plot([min(min(x)),max(max(x))],[12,12],'k--','LineWidth',2)%,'DisplayName','Energy Resolution Requirement')

%legend([EngLegend,'Energy Resolution Requirement'],'Location', 'southoutside','NumColumns',8)
legend([strings(1,length(energy_channels)),'Energy Resolution Requirement'],'Location', 'northeast','NumColumns',8)

hold off
effsave = append(date(),' Energy Resolution',' Eo ',num2str(min(Eo_Whole)),' to ',num2str(max(Eo_Whole)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
saveas(gcf,effsave)

% %% Intersection Plot
% hold on
% num = 0;
% for j = 1:(length(Eo_Whole)-1)
%
%     for k = 1:(length(Eo_Whole)-j)
%         plot(xi_Whole(num+k,16),yi_Whole(num+k,16),'o','Color',Eo_color_Whole(j,:))
%
%     end
%     num = num + k;
% end
% plot(E_eff_Whole(16),G_eff_dE_Whole(16),'or')
% hold off

%% Count Rate Comparision
% cd 'Bow Tie'
% 
% f = figure;
% f.Position = [100 100 1200 720];
% 
% %Plots each energy channel FWHM value in the same color as
% %the Eff. Curve
% 
% semilogy(energy_channel_list,BowCount_Rate_WholeMax50_Omni,'o')
% 
% 
% xlim([0 length(energy_channels)+1])
% xticks((0:1:length(energy_channels)+1))
% 
% legend({'Whole Count Rate','Inner Count Rate'},'Location', 'southoutside')
% title(append('Count Rates Comparision 50th ',' Eo Range: ',num2str(min(Eo_Whole)),' to ',num2str(max(Eo_Whole)),' MeV'))
% ylabel('Count Rates')
% xlabel('Energy Channels')
% 
% 
% effsave = append(date(),' Count Rate Comparision 50th',' Eo ',num2str(min(Eo_Whole)),' to ',num2str(max(Eo_Whole)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
% saveas(gcf,effsave)
% 
% f = figure;
% f.Position = [100 100 1200 720];
% semilogy(energy_channel_list,BowCount_Rate_WholeMax75_Omni,'o',energy_channel_list,BowCount_Rate_InnerMax75_Omni,'og')
% 
% 
% xlim([0 length(energy_channels)+1])
% xticks((0:1:length(energy_channels)+1))
% 
% legend({'Whole Count Rate','Inner Count Rate'},'Location', 'southoutside')
% title(append('Count Rates Comparision 75th',' Eo Range: ',num2str(min(Eo_Whole)),' to ',num2str(max(Eo_Whole)),' MeV'))
% ylabel('Count Rates')
% xlabel('Energy Channels')
% 
% 
% effsave = append(date(),' Count Rate Comparision 75th',' Eo ',num2str(min(Eo_Whole)),' to ',num2str(max(Eo_Whole)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
% saveas(gcf,effsave)
% 
% f = figure;
% f.Position = [100 100 1200 720];
% semilogy(energy_channel_list,BowCount_Rate_WholeMax95_Omni,'o',energy_channel_list,BowCount_Rate_InnerMax95_Omni,'og')
% 
% 
% xlim([0 length(energy_channels)+1])
% xticks((0:1:length(energy_channels)+1))
% 
% legend({'Whole Count Rate','Inner Count Rate'},'Location', 'southoutside')
% title(append('Count Rates Comparision 95th',' Eo Range: ',num2str(min(Eo_Whole)),' to ',num2str(max(Eo_Whole)),' MeV'))
% ylabel('Count Rates')
% xlabel('Energy Channels')
% 
% 
% effsave = append(date(),' Count Rate Comparision 95th',' Eo ',num2str(min(Eo_Whole)),' to ',num2str(max(Eo_Whole)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
% saveas(gcf,effsave)

% cd ..

%% Count Rate Comparision-per Energy

%All Channels
                channel_select = [1:length(energy_channels)];
                %channel_select = [2,10,20,30,35];
                
f = figure;
f.Position = [100 100 1200 720];
color_iter =1;
EC_plot_color = plasma(length(channel_select));
%Plots each energy channel FWHM value in the same color as
%the Eff. Curve
hold on
for i = 1:width(x)
    if max(i == channel_select)
        plot(E_eff_Whole(i),BowCount_Rate_WholeApo50__90(i,1),'o','MarkerEdgeColor','b',...
        'MarkerFaceColor','b');
    %EC_plot_color(i,:)
        EngLegend_EC(i) = append(sprintf('#%.0f: ',i),EngLegend(i));
        color_iter = color_iter+1;
        
    else
        plot(E_eff_Whole(i),BowCount_Rate_WholeApo50__90(i),'Color',[0.75, 0.75, 0.75],'LineWidth',line_width);
    end
    
end

for i = 1:width(x)
    if max(i == channel_select)
        plot(E_eff_Whole(i),BowCount_Rate_WholeCen50__90(i,1),'o','MarkerEdgeColor','g',...
        'MarkerFaceColor','g');
    %EC_plot_color(i,:)
        EngLegend_EC(i) = append(sprintf('#%.0f: ',i),EngLegend(i));
        color_iter = color_iter+1;
        
    else
        plot(E_eff_Whole(i),BowCount_Rate_WholeCen50__90(i),'Color',[0.75, 0.75, 0.75],'LineWidth',line_width);
    end
    
end

set(gca, 'YScale', 'log')
set(gca, 'Fontsize', 20)
hold off



%xlim([0 8])
%xticks((0:0.5:8))

%legend(EngLegend_EC,'Location', 'southoutside','NumColumns',6)
%title(append('Count Rates Comparision 50th ',' Eo Range: ',num2str(Cen(Eo_Whole)),' to ',num2str(max(Eo_Whole)),' MeV'))
ylabel('Count Rates (#/sec)')
xlabel('Nominal Energy (MeV)')


effsave = append(date(),' Count Rate Comparision 50th',' Eo ',num2str(min(Eo_Whole)),' to ',num2str(max(Eo_Whole)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
saveas(gcf,effsave)

%%
% f = figure;
% f.Position = [100 100 1200 720];
% semilogy(energy_channel_list,BowCount_Rate_WholeMax75_Omni,'o',energy_channel_list,BowCount_Rate_InnerMax75_Omni,'og')
% 
% 
% xlim([0 length(energy_channels)+1])
% xticks((0:1:length(energy_channels)+1))
% 
% legend({'Whole Count Rate','Inner Count Rate'},'Location', 'southoutside')
% title(append('Count Rates Comparision 75th',' Eo Range: ',num2str(Cen(Eo_Whole)),' to ',num2str(max(Eo_Whole)),' MeV'))
% ylabel('Count Rates')
% xlabel('Energy Channels')
% 
% 
% effsave = append(date(),' Count Rate Comparision 75th',' Eo ',num2str(Cen(Eo_Whole)),' to ',num2str(max(Eo_Whole)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
% saveas(gcf,effsave)
% 
% f = figure;
% f.Position = [100 100 1200 720];
% semilogy(energy_channel_list,BowCount_Rate_WholeMax95_Omni,'o',energy_channel_list,BowCount_Rate_InnerMax95_Omni,'og')
% 
% 
% xlim([0 length(energy_channels)+1])
% xticks((0:1:length(energy_channels)+1))
% 
% legend({'Whole Count Rate','Inner Count Rate'},'Location', 'southoutside')
% title(append('Count Rates Comparision 95th',' Eo Range: ',num2str(Cen(Eo_Whole)),' to ',num2str(max(Eo_Whole)),' MeV'))
% ylabel('Count Rates')
% xlabel('Energy Channels')
% 
% 
% effsave = append(date(),' Count Rate Comparision 95th',' Eo ',num2str(Cen(Eo_Whole)),' to ',num2str(max(Eo_Whole)),' MeV','_',addin,num2str(length(energy_channels)),'.jpg');
% saveas(gcf,effsave)



%% Energy per Second for Each Detector

%whole_detector_energy_S = FluxRangeEL(:,1:width(whole_detector_energy)).*whole_detector_energy;
%
% inner_detector_energy_S = FluxRangeEL(:,1:width(inner_detector_energy)).*inner_detector_energy;
%
% Whole_total_rate = sum(whole_detector_energy_S,1)'
% Inner_total_rate = sum(inner_detector_energy_S,1)'
%
% f = figure;
% f.Position = [100 100 1200 720];
%
%
% semilogy(x(:,1),whole_detector_energy_S(:,1),x(:,1),whole_detector_energy_S(:,2),x(:,1),whole_detector_energy_S(:,3),x(:,1),whole_detector_energy_S(:,4),x(:,1),whole_detector_energy_S(:,5),x(:,1),whole_detector_energy_S(:,6),x(:,1),whole_detector_energy_S(:,7),x(:,1),whole_detector_energy_S(:,8),x(:,1),whole_detector_energy_S(:,9))
% title('Whole Config: Energy Deposited per Second (MeV/s) vs. Incident Energy (MeV)')
% ylabel('Energy Deposited per Second (MeV/s)')
% xlabel('Incident Energy (MeV)')
%
% effsave = append(date(),' Energy Rate-Whole ',addin,num2str(length(energy_channels)),'.jpg');
% saveas(gcf,effsave)
%
% f = figure;
% f.Position = [100 100 1200 720];
%
%
% semilogy(x(:,1),inner_detector_energy_S(:,1),x(:,1),inner_detector_energy_S(:,2),x(:,1),inner_detector_energy_S(:,3),x(:,1),inner_detector_energy_S(:,4),x(:,1),inner_detector_energy_S(:,5),x(:,1),inner_detector_energy_S(:,6),x(:,1),inner_detector_energy_S(:,7),x(:,1),inner_detector_energy_S(:,8),x(:,1),inner_detector_energy_S(:,9),x(:,1),inner_detector_energy_S(:,10),x(:,1),inner_detector_energy_S(:,11),x(:,1),inner_detector_energy_S(:,12),x(:,1),inner_detector_energy_S(:,13),x(:,1),inner_detector_energy_S(:,14),x(:,1),inner_detector_energy_S(:,15),x(:,1),inner_detector_energy_S(:,16),x(:,1),inner_detector_energy_S(:,17))
% title('Inner Config: Energy Deposited per Second (MeV/s) vs. Incident Energy (MeV)')
% ylabel('Energy Deposited per Second (MeV/s)')
% xlabel('Incident Energy (MeV)')
%
% effsave = append(date(),' Energy Rate-Inner ',addin,num2str(length(energy_channels)),'.jpg');
% saveas(gcf,effsave)

%% Back Pen Calculations Electron
x_inter = [0.5,0.75,1:0.5:10];
ApoFlux50_AE9_90= [341776.9940	78702.5833	20845.5994	5693.3557	1472.3851	399.9446	129.0902	49.2613	23.5158	10.5374	4.8152	2.3604	1.4795	1.0042	0.5324	0.4305	0.3286	0.2267	0.1948	0.1629	0.1310];
MaxFlux95_AE9_90= [12559884.1	3917196.38	2094895.39	826481.191	256040.313	79982.653	30444.5749	13671.416	6070.77165	2795.68906	1292.31374	654.009611	366.498171	235.290299	139.450819	119.126969	98.8031184	78.9011666	71.2501658	63.5991651	56.0650802];

%ApoFlux50_90_(17:end) = interp1(x_inter,ApoFlux50_AE9_90,[8.5:0.5:10],'linear');
%MaxFlux95_90_(17:end) = interp1(x_inter,MaxFlux95_AE9_90,[8.5:0.5:10],'linear');

%ApoFlux50_90_(end) = 0.01;


Apo50BackCountRate = trapz(x(:,1),(Back_Hits_Whole./output_number)*(pi*(pi*2^2)).*ApoFlux50_AE9_90)
Max95BackCountRate = trapz(x(:,1),(Back_Hits_Whole./output_number)*(pi*(pi*2^2)).*MaxFlux95_AE9_90)

function [MaxFlux25_90,CenFlux25_90,ApoFlux25_90,MaxFlux50_90,CenFlux50_90,ApoFlux50_90,MaxFlux75_90,CenFlux75_90,ApoFlux75_90,MaxFlux95_90,CenFlux95_90,ApoFlux95_90,MaxFluxMean_90,CenFluxMean_90,ApoFluxMean_90] = FluxRate90(r)
%FluxRate Summary of this function goes here
%   Inputs:
%     X- energy at flux
%     Outputs:
%     Flux50:Expected average flux at energy X
%     Flux75:Expected 75th percentile flux at energy X
%     Flux95:Expected 95th percentile flux at energy X

x = [0.5,0.75,1:0.5:8];

%% Omnidirectional
%Values from AE9
%25th Percentile
MaxFlux25_AE9_90= [1299128.18	159571.67	58982.02	15129.9237	3041.56181	718.45857	219.159343	77.6956967	26.2584103	10.3683737	4.73543853	2.43125066	1.31421209	0.694297103	0.328897155 0.265913281	0.202929408];
CenFlux25_AE9_90= [236773.987	102478.867	42665.4286	12486.5777	2898.92566	683.399265	201.636631	69.2078129	24.9863608	10.0010258	4.63892914	2.28897503	1.19219726	0.648284687	0.312877121	0.253244987	0.193612852];
ApoFlux25_AE9_90= [127986.22	24967.9964	5444.69973	1319.99226	312.878706	81.8187931	25.3755998	9.72523144	4.90561484	2.21073182	0.957248098	0.416308498	0.239633178	0.152650546	0.07342178 0.057816453 0.042211125];

%50th Percentile
MaxFlux50_AE9_90= [3043556.52	506402.542	211356.954	62748.4502	15499.6257	4167.08963	1362.59182	508.053908	196.611127	82.051654	37.154501	17.4685642	8.97551526	5.28536699	2.81916956 2.33	1.84];
CenFlux50_AE9_90= [876087.841	397179.856	181014.549	58780.2084	15120.1328	3942.52241	1290.06064	487.952961	187.386107	76.4888646	33.3951236	16.0558237	8.53704862	5.05846258	2.61591484 2.14087987	1.66584491];
ApoFlux50_AE9_90= [341776.994	78702.5833	20845.5994	5693.35574	1472.38512	399.944599	129.090197	49.2612793	23.5157668	10.5373774	4.81524394	2.36038762	1.4794959	1.00421865	0.532351405 0.430483052	0.328614699];

%75th Percentile
MaxFlux75_AE9_90= [5954877.83	1258557.25	590380.762	202629.404	58159.668	16656.5776	5906.61956	2397.14709	991.163556	436.275356	196.258873	94.3908837	49.8709683	30.3881283	16.5763025 13.7416615 10.9904901];
CenFlux75_AE9_90= [2469165.81	1160706.51	567296.774	199717.547	55612.4246	15693.0662	5571.19522	2275.12434	917.33847	380.385147	158.469032	74.7436904	40.4205575	25.6070785	13.9684838	11.5349933	9.10150277];
ApoFlux75_AE9_90= [741554.833	194643.441	60101.779	18034.1354	4996.24183	1398.42223	465.736857	177.111453	80.9439822	36.1004976	17.2120895	9.2718457	6.21568293	4.43625712	2.53942722 2.10258632 1.66562056];

%95th Percentile
MaxFlux95_AE9_90= [12559884.1	3917196.38	2094895.39	826481.191	256040.313	79982.653	30444.5749	13671.416	6070.77165	2795.68906	1292.31374	654.009611	366.498171	235.290299	139.450819 119.126969 98.8031184];
CenFlux95_AE9_90= [7854048.94	3842396.1	2026654.85	779563.225	236687.829	72888.9963	28329.6745	12599.2993	5363.30435	2263.84066	895.985629	414.265436	228.444014	155.761882	90.0617119	75.1379942	60.2142765];
ApoFlux95_AE9_90= [1753727.44	532323.177	194919.825	64926.6752	19416.2462	5617.60949	1937.2571	734.009612	319.68837	141.834443	70.8966219	42.4119355	30.6420444	23.1319331	14.4237781 12.2959433 10.1681086];

%Mean
MaxFluxMean_AE9_90= [4284294.4	987847.371	497075.967	190568.703	57386.6082	17693.1617	6656.09705	2948.71527	1302.47442	598.629251	276.547916	140.044867	79.0806048	50.969725	31.2524332	26.9827706	22.713108];
CenFluxMean_AE9_90= [2015489.32	966624.012	492501.967	183035.525	54030.6479	16263.9459	6207.6321	2727.12992	1154.56719	486.691215	193.553499	89.7218939	49.3789625	33.4564409	19.2888158	16.0897458	12.8906759];
ApoFluxMean_AE9_90= [543134.812	148679.823	49426.7817	15701.1899	4559.39628	1304.11571	444.805326	168.743676	74.6672731	33.1902862	16.3350572	9.49303049	6.74878269	5.04578025	3.11202649 2.64802017 2.18401386];

%Values from AE9
%25th Percentile

%50th Percentile

%75th Percentile

%95th Percentile
    %2404183.7	286581.807	80278.5175	27311.1463	10245.9702	3822.18999	1828.73641	863.673911	425.782876	241.123197	154.292021	105.154032	71.3507142	50.6479176	32.2673942];


%Calculate flux value
% 25th Percetile
CenFlux25_90 = interp1(x,CenFlux25_AE9_90,r);
MaxFlux25_90 = interp1(x,MaxFlux25_AE9_90,r);
ApoFlux25_90 = interp1(x,ApoFlux25_AE9_90,r);

%50th Percetile
CenFlux50_90 = interp1(x,CenFlux50_AE9_90,r);
MaxFlux50_90 = interp1(x,MaxFlux50_AE9_90,r);
ApoFlux50_90 = interp1(x,ApoFlux50_AE9_90,r);

%75th Percetile
CenFlux75_90 = interp1(x,CenFlux75_AE9_90,r);
MaxFlux75_90 = interp1(x,MaxFlux75_AE9_90,r);
ApoFlux75_90 = interp1(x,ApoFlux75_AE9_90,r);

%95th Percetile
CenFlux95_90 = interp1(x,CenFlux95_AE9_90,r);
MaxFlux95_90 = interp1(x,MaxFlux95_AE9_90,r);
ApoFlux95_90 = interp1(x,ApoFlux95_AE9_90,r);

%Average
CenFluxMean_90 = interp1(x,CenFluxMean_AE9_90,r);
MaxFluxMean_90= interp1(x,MaxFluxMean_AE9_90,r);
ApoFluxMean_90 = interp1(x,ApoFluxMean_AE9_90,r);






end

function [MaxFlux25Omni,CenFlux25Omni,ApoFlux25Omni,MaxFlux50Omni,CenFlux50Omni,ApoFlux50Omni,MaxFlux75Omni,CenFlux75Omni,ApoFlux75Omni,MaxFlux95Omni,CenFlux95Omni,ApoFlux95Omni,MaxFluxMeanOmni,CenFluxMeanOmni,ApoFluxMeanOmni] = FluxRateOmni(r)
%FluxRate Summary of this function goes here
%   Inputs:
%     X- energy at flux
%     Outputs:
%     Flux50:Expected average flux at energy X
%     Flux75:Expected 75th percentile flux at energy X
%     Flux95:Expected 95th percentile flux at energy X

x = [0.5,0.75,1:0.5:8];

% Conversion Factor to /sr from Omni
CF = 4*pi;

%% Omnidirectional
%Values from AE9
%25th Percentile
MaxFlux25_AE9Omni= [8111977.47	1573149.04	570832.259	145269.2	29218.1317	6903.77267	2063.2227	742.263646	266.679309	109.348528	49.3155003	24.100156	13.0634938	7.06531883	3.28402254 2.64043001	1.99701516]./CF;
CenFlux25_AE9Omni= [2363391.15	1010236.28	420136.43	123138.153	28252.2425	6593.3688	1950.68566	681.678876	259.525518	106.342774	48.4952391	23.0311398	11.8035589	6.45582063	3.0726934	2.47457182	1.87645024]./CF;
ApoFlux25_AE9Omni= [1329508.46	252674.19	54324.3703	12858.4732	2958.76327	756.557986	235.070404	90.9015458	47.2026265	22.0241695	10.0087299	4.57676091	2.71873186	1.77166374	0.887809304	0.706650242	0.525491181]./CF;

%50th Percentile
MaxFlux50_AE9Omni= [19893090.7	5097010.91	2097907.42	624573.54	154628.957	40909.9701	13350.1255	5049.95988	2026.19355	857.250336	375.862092	170.394972	85.9078482	50.603342	26.4107426	21.6188826	16.864596]./CF;
CenFlux50_AE9Omni= [8902106.27	3997011.55	1816783.24	591581.383	150674.077	38761.4475	12597.5379	4815.77357	1932.24052	811.280178	348.745134	160.521736	82.9799251	49.2802646	25.5049115	20.8261552	16.1473989]./CF;
ApoFlux50_AE9Omni= [3561954.92	801437.742	209995.397	55915.559	13994.7966	3726.22034	1214.48544	472.01705	231.316141	106.900138	50.3908547	25.2086036	15.7579626	10.6856645	5.76092933	4.68958599	3.61824264]./CF;

%75th Percentile
MaxFlux75_AE9Omni= [40883035.7	12922351.4	6049174.9	2077265.16	590046.286	166606.194	58924.2062	23885.7512	10125.2594	4511.16945	1978.85991	914.656732	473.531305	291.004977	163.06214	136.183404	109.948726]./CF;
CenFlux75_AE9Omni= [25466597	11884475.2	5783791.43	2043482.19	564388.867	156691.424	54819.331	22491.5776	9408.0567	4029.48609	1654.74889	743.967198	387.630602	245.668857	135.648925	112.031465	88.4140052]./CF;
ApoFlux75_AE9Omni= [7748204.15	1991570.86	609961.739	178234.465	47666.3592	13100.9096	4433.70582	1730.77206	810.478391	371.653477	180.193512	96.8353773	63.2519282	44.4424177	25.4969211	21.1858925	16.8748639]./CF;

%95th Percentile
MaxFlux95_AE9Omni= [91558289.4	41140003.2	21777989.6	8588339.04	2658679.36	821612.364	310389.614	138246.272	62029.605	28853.2118	13085.9523	6466.13666	3657.18038	2383.31583	1430.15343	1220.62378	1011.13351]./CF;
CenFlux95_AE9Omni= [82438977.2	40152604.7	21041509.5	8128343.18	2453176.17	740869.067	281208.686	124800.997	54691.1672	23963.8502	9361.59773	4105.59847	2160.44605	1471.82903	872.166313	729.849019	587.531725]./CF;
ApoFlux95_AE9Omni= [18390202.3	5481304.18	1997029.96	647185.244	186317.874	53039.0122	18715.0165	7344.80766	3270.54056	1486.48666	743.047183	432.730168	297.843794	218.480884	134.831044	114.986834	95.1426237]./CF;

%Mean
MaxFluxMean_AE9Omni= [29719534	10294550.9	5145356.75	1972260.08	591967.488	180958.825	67634.9545	29785.0655	13309.3245	6179.96504	2801.18165	1385.17797	791.346898	521.590481	314.831095	269.160117	223.703666]./CF;
CenFluxMean_AE9Omni= [20975174.9	10008126	5076315.94	1896363.99	556804.147	164682.291	61532.9153	27007.8912	11780.2207	5153.02796	2022.54965	889.781998	467.824031	316.607283	186.893584	156.369112	125.84464]./CF;
ApoFluxMean_AE9Omni= [5682405.63	1525791.21	504437.068	156112.947	43709.6229	12295.5287	4282.62169	1678.31053	759.102006	346.00969	171.227699	97.4892416	66.3249877	48.2615504	29.4447683	25.0441188	20.6434693]./CF;

%Calculate flux value
% 25th Percetile
CenFlux25Omni = interp1(x,CenFlux25_AE9Omni,r);
MaxFlux25Omni = interp1(x,MaxFlux25_AE9Omni,r);
ApoFlux25Omni = interp1(x,ApoFlux25_AE9Omni,r);

%50th Percetile
CenFlux50Omni = interp1(x,CenFlux50_AE9Omni,r);
MaxFlux50Omni = interp1(x,MaxFlux50_AE9Omni,r);
ApoFlux50Omni = interp1(x,ApoFlux50_AE9Omni,r);

%75th Percetile
CenFlux75Omni = interp1(x,CenFlux75_AE9Omni,r);
MaxFlux75Omni = interp1(x,MaxFlux75_AE9Omni,r);
ApoFlux75Omni = interp1(x,ApoFlux75_AE9Omni,r);

%95th Percetile
CenFlux95Omni = interp1(x,CenFlux95_AE9Omni,r);
MaxFlux95Omni = interp1(x,MaxFlux95_AE9Omni,r);
ApoFlux95Omni = interp1(x,ApoFlux95_AE9Omni,r);

%Average
CenFluxMeanOmni = interp1(x,CenFluxMean_AE9Omni,r);
MaxFluxMeanOmni = interp1(x,MaxFluxMean_AE9Omni,r);
ApoFluxMeanOmni = interp1(x,ApoFluxMean_AE9Omni,r);




end


