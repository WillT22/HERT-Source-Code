function [MaxFlux25Omni,MinFlux25Omni,ApoFlux25Omni,MaxFlux50Omni,MinFlux50Omni,ApoFlux50Omni,MaxFlux75Omni,MinFlux75Omni,ApoFlux75Omni,MaxFlux95Omni,MinFlux95Omni,ApoFlux95Omni,MaxFluxMeanOmni,MinFluxMeanOmni,ApoFluxMeanOmni] = FluxRateOmni(r)
%FluxRate Summary of this function goes here
%   Inputs:
%     X- energy at flux
%     Outputs:
%     Flux50:Expected average flux at energy X
%     Flux75:Expected 75th percentile flux at energy X
%     Flux95:Expected 95th percentile flux at energy X

x = [0.5,0.75,1:0.5:7];

% Conversion Factor to /sr from Omni
CF = 4*pi;

%% Omnidirectional
%Values from AE9
%25th Percentile
MaxFlux25_AE9Omni= [8111977.47	1573149.04	570832.259	145269.2	29218.1317	6903.77267	2063.2227	742.263646	266.679309	109.348528	49.3155003	24.100156	13.0634938	7.06531883	3.28402254]./CF;
MinFlux25_AE9Omni= [119336.699	8225.82878	2706.75572	1120.70882	441.600604	172.554702	70.32912	28.4551312	12.1755939	5.72643673	3.11791022	1.9170406	1.23915444	0.804382395	0.466120109]./CF;
ApoFlux25_AE9Omni= [1329043.75	252562.735	54297.3124	12851.4779	2957.10699	756.164159	234.966797	90.8712388	47.1895398	22.0176569	10.0054712	4.57520635	2.71797272	1.77125091	0.887628471]./CF;

%50th Percentile
MaxFlux50_AE9Omni= [19893090.7	5097010.91	2097907.42	624573.54	154628.957	40909.9701	13350.1255	5049.95988	2026.19355	857.250336	375.862092	170.394972	85.9078482	50.603342	26.4107426]./CF;
MinFlux50_AE9Omni= [831159.082	75215.3761	24933.8101	10050.2331	4041.88421	1583.18182	688.277682	298.083931	133.99576	67.5783659	38.6913855	24.7525623	16.1691211	10.7614667	6.38568385]./CF;
ApoFlux50_AE9Omni= [3560799.93	801096.429	209891.941	55885.6218	13987.2765	3724.35075	1213.94357	471.837509	231.23887	106.863573	50.3731833	25.199803	15.7532854	10.6828812	5.7595443]./CF;

%75th Percentile
MaxFlux75_AE9Omni= [40883035.7	12922351.4	6049174.9	2077265.16	590046.286	166606.194	58924.2062	23885.7512	10125.2594	4511.16945	1978.85991	914.656732	473.531305	291.004977	163.06214]./CF;
MinFlux75_AE9Omni= [3846756.52	434054.751	144124.757	56744.2736	23170.7812	9094.12836	4159.27749	1901.0335	888.89288	473.988988	282.537685	186.68609	123.24219	84.0403946	51.0708825]./CF;
ApoFlux75_AE9Omni= [7745852.95	1990748.82	609664.598	178140.465	47641.6501	13094.5472	4431.71299	1730.04776	810.168494	371.511692	180.126638	96.8010134	63.2327947	44.4305583	25.490612]./CF;

%95th Percentile
MaxFlux95_AE9Omni= [91558289.4	41140003.2	21777989.6	8588339.04	2658679.36	821612.364	310389.614	138246.272	62029.605	28853.2118	13085.9523	6466.13666	3657.18038	2383.31583	1430.15343]./CF;
MinFlux95_AE9Omni= [21177994.6	3071319.81	1017724.32	389492.874	161615.176	63576.8611	30758.273	14929.137	7296.80974	4142.65193	2584.25855	1771.65054	1185.81287	834.132484	522.86174]./CF;
ApoFlux95_AE9Omni= [18385065.5	5479128.5	1996072.17	646850.622	186225.504	53014.2791	18706.5562	7341.4144	3269.1016	1485.84951	742.754015	432.57455	297.754015	218.423896	134.799101]./CF;

%Mean
MaxFluxMean_AE9Omni= [29719534	10294550.9	5145356.75	1972260.08	591967.488	180958.825	67634.9545	29785.0655	13309.3245	6179.96504	2801.18165	1385.17797	791.346898	521.590481	314.831095]./CF;
MinFluxMean_AE9Omni= [4590162.32	658963.392	218115.238	83384.8825	34597.477	13610.6814	6596.79335	3214.55575	1577.85667	902.813172	567.118013	391.331833	262.762424	186.186808	117.636928]./CF;
ApoFluxMean_AE9Omni= [5680714.31	1525170.91	504193.498	156031.64	43687.6255	12289.7226	4280.69018	1677.55596	758.781532	345.866514	171.16136	97.4543412	66.3050963	48.2490523	29.4378704]./CF;

%Calculate flux value
% 25th Percetile
MinFlux25Omni = interp1(x,MinFlux25_AE9Omni,r);
MaxFlux25Omni = interp1(x,MaxFlux25_AE9Omni,r);
ApoFlux25Omni = interp1(x,ApoFlux25_AE9Omni,r);

%50th Percetile
MinFlux50Omni = interp1(x,MinFlux50_AE9Omni,r);
MaxFlux50Omni = interp1(x,MaxFlux50_AE9Omni,r);
ApoFlux50Omni = interp1(x,ApoFlux50_AE9Omni,r);

%75th Percetile
MinFlux75Omni = interp1(x,MinFlux75_AE9Omni,r);
MaxFlux75Omni = interp1(x,MaxFlux75_AE9Omni,r);
ApoFlux75Omni = interp1(x,ApoFlux75_AE9Omni,r);

%95th Percetile
MinFlux95Omni = interp1(x,MinFlux95_AE9Omni,r);
MaxFlux95Omni = interp1(x,MaxFlux95_AE9Omni,r);
ApoFlux95Omni = interp1(x,ApoFlux95_AE9Omni,r);

%Average
MinFluxMeanOmni = interp1(x,MinFluxMean_AE9Omni,r);
MaxFluxMeanOmni = interp1(x,MaxFluxMean_AE9Omni,r);
ApoFluxMeanOmni = interp1(x,ApoFluxMean_AE9Omni,r);


%% User created interpol
% if r>= 0.5 && r<= 0.75
%     MinFlux50Omni = -5264851.86*r + 4021596.26;%y = -5,264,851.86x + 4,021,596.26
%     
%     ApoFlux50Omni = -18630863.76*r+ 15952251.27;
%     
%     
%     MaxFlux50Omni = -70332511.16*r + 58073941.28;
%     MaxFlux75Omni = -130533391.20*r + 111367423.10;%y = -130,533,391.20x + 111,367,423.10
%     MaxFlux95Omni = -234642452*r + 218019366;
%     
% elseif r>= 0.75 && r<= 1.0
%     MinFlux50Omni = -221590.29*r + 239150.09; %y = -221,590.29x + 239,150.09
%     
%     ApoFlux50Omni = -5477135.85*r+ 6086955.34;
%     
%     MaxFlux50Omni = -12459253.88*r + 14668998.32;
%     MaxFlux75Omni = -28791595.52*r + 35061076.34;%y = -28,791,595.52x + 35,061,076.34
%     MaxFlux95Omni = -78888027.2*r + 101203547;
%     
% elseif r>= 1.0 && r<= 1.5
%     MinFlux50Omni = -23672.54*r + 41232.33;%y = -23,672.54x + 41,232.33
%     
%     ApoFlux50Omni = -845354.54*r+ 1455174.03;
%     
%     MaxFlux50Omni = -3104268.57*r + 5314013.01;
%     MaxFlux75Omni = -8225808.24*r + 14495289.06;%y = -8,225,808.24x + 14,495,289.06
%     MaxFlux95Omni = -26999680.40*r + 49315200.60;
%     
% elseif r>= 1.5 && r<= 2.0
%     MinFlux50Omni = -7444.22*r + 16889.85;%y = -7,444.22x + 16,889.85
%     
%     ApoFlux50Omni = -270345.05*r+ 592656.79;
%     
%     MaxFlux50Omni = -992812.41*r + 2146828.77;
%     MaxFlux75Omni = -3091070.44*r+ 6796182.35;%y = -3,091,070.44x + 6,793,182.35
%     MaxFlux95Omni = -12185602*r+27094083;
%     
% elseif r>= 2.0 && r<= 2.5
%     MinFlux50Omni = -2656.43*r + 7314.26; %y = -2,656.43x + 7,314.26
%     
%     ApoFlux50Omni = -74030.92*r+ 200031.53;
%     
%     
%     MaxFlux50Omni = -236816.57*r + 634837.08;
%     MaxFlux75Omni = -875112.76*r + 2361266.99;%y = -875,112.76x + 2,361,266.99
%     MaxFlux95Omni = -3757378*r + 10237635;
%     
% elseif r>= 2.5 && r<= 3.0
%     MinFlux50Omni = -810.78*r + 2700.15;%y = -810.78x + 2,700.15
%     
%     ApoFlux50Omni = -20697.8*r+ 66698.73;
%     
%     MaxFlux50Omni = -57669.67*r + 186969.84;
%     MaxFlux75Omni = -224890.53*r + 735711.43;%y = -224,890.53x + 735,711.43
%     MaxFlux95Omni = -1049424.2*r + 3467750.50;
%     
% elseif r>= 3.0 && r<= 3.5
%     MinFlux50Omni = -307.70*r + 1181.91;%y = -304.70x + 1,181.91
%     
%     ApoFlux50Omni = -6149.2*r+ 23052.95;
%     
%     MaxFlux50Omni = -17434.76*r + 66265.12;
%     MaxFlux75Omni = -72402.04*r + 278245.97;%y = -72,402.04x + 278,245.97
%     MaxFlux95Omni = -354126.40*r + 1381857.10;
%     
% elseif r>= 3.5 && r<= 4.0
%     MinFlux50Omni = -125.08*r +553.23;%y = -125.08x + 553.23
%     
%     ApoFlux50Omni = -1870.65*r+ 8078.02;
%     
%     MaxFlux50Omni = -6315.43*r + 27347.45;
%     MaxFlux75Omni = -28900.74*r + 125991.42;%y = -28,900.74x + 125,991.42
%     MaxFlux95Omni = -157786.20*r + 694666.40;
%     
% elseif r>= 4.0 && r<= 4.5
%     MinFlux50Omni = -50.83*r + 256.23;%y = -50.83x + 256.23
%     
%     ApoFlux50Omni = -701.91*r+ 3403.04;
%     
%     MaxFlux50Omni = -2413.47*r + 11739.61;
%     MaxFlux75Omni = -11572.09*r + 56676.81;%y = -11,572.09x + 56,676.81
%     MaxFlux95Omni = -68310.70*r + 336764.40;
%     
% elseif r>= 4.5 && r<= 5.0
%     MinFlux50Omni = -23.13*r + 131.59; %y = -23.13x + 131.59
%     
%     ApoFlux50Omni = -266.24*r+ 1442.54;
%     
%     MaxFlux50Omni = -982.49*r + 5300.21;
%     MaxFlux75Omni = -5118.29*r + 27634.71;%y = -5,118.29x + 27,634.71
%     MaxFlux95Omni = -32040.82*r + 173549.94;
%     
% elseif r>= 5.0 && r<= 5.5
%     MinFlux50Omni = -11.74*r+74.64;%y = -11.74x + 74.64
%     
%     ApoFlux50Omni = -109.93*r+ 660.95;
%     
%     MaxFlux50Omni = -419.40*r + 2484.75;
%     MaxFlux75Omni = -2170.11*r + 12893.83;%y = -2,170.11x + 12,893.83
%     MaxFlux95Omni = -13367.48*r + 80183.26;
%     
% elseif r>= 5.5 && r<= 6.0
%     MinFlux50Omni = -6.77*r + 47.29;%y = -6.77x + 47.29
%     
%     ApoFlux50Omni = -47.9*r+ 319.8;
%     
%     MaxFlux50Omni = -173.75*r + 1133.68;
%     MaxFlux75Omni = -921.19*r + 6024.73;%y = -921.19x + 6,024.73
%     MaxFlux95Omni = -5867.98*r + 38935.98;
%     
% elseif r>= 6.0 && r<= 6.5
%     MinFlux50Omni = -4.07*r + 31.13;%y = -4.07x + 31.13
%     
%     ApoFlux50Omni = -24.12*r+ 177.14;
%     
%     MaxFlux50Omni = -74.75*r + 539.66;
%     MaxFlux75Omni = -385.44*r + 2810.25;%y = -385.44x + 2,810.25
%     MaxFlux95Omni = -2642.26*r + 19581.67;
%     
% elseif r>= 6.5 && r<= 7.0
%     MinFlux50Omni = -3.64*r + 28.31;%y = -3.64x + 28.31
%     
%     ApoFlux50Omni = -19.71*r+ 148.48;
%     
%     MaxFlux50Omni = -51.35*r + 387.59;
%     MaxFlux75Omni = -272.16*r + 2073.91;%y = -272.16x + 2,073.91
%     MaxFlux95Omni = -1929.60*r + 14949.40;
%     
% end



end

