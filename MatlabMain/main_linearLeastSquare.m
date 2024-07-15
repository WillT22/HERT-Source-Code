%% Routine: Use linear least square algorithms (Selesnick, 4.1) to determine flux and compare that with AE9/AP9 flux spectra
% by LYK 
% History:
% Initiate (Sep 9, 2020)
% 10/9/2020 : adapted for extreme flux input
% 10/13/2020: Clean up
% Input: 
% input file yield from routine_Part1_input file (should be in the
% same InverseProblem file
% Action: 
% 1) for normal spectra (comment out jE, n)
%    for other spectra (update the n, jE to reflect the correct spectra)
% 2) Update the save file to reflect the correct spectra 
% and make sure the save section is comment out when experimenting with the optimal
% del, sigmaM
% 3) to find optimal del, sigmaM: first try with a large del and reasonable
% sigmaM and slowly decrease it. 

%housekeeping
clear all;
clc;
close all;

cd '/Users/lengying/Documents/GEANT4_2017/InverseProb'
%% Step 0: reference Sr-90 spectrum
load('/Users/lengying/Documents/Research/projects/CIRBE/Calibration:Test/strontium/Vaugh_Sr90EnergySpectra');
load('/Users/lengying/Documents/Research/projects/CIRBE/Calibration:Test/strontium/Euramet_Sr90EnergySpectra');
eSr = [0.005,0.041,0.062,0.093,0.124,0.155,0.207,0.269,0.31,0.341,0.362,0.424,0.496,0.651,0.816,1.002,1.178,1.395,1.56,1.684,1.88,2.247];
jSr = [0,0.096,0.225,0.508,0.667,0.879,0.992,0.888,0.671,0.513,0.388,0.254,0.188,0.217,0.225,0.196,0.15,0.096,0.05,0.029,0.004,0.008];

%% Step 1: Read Sr results (count rate, flux, geometric factor)
h5File = '/Users/lengying/Documents/Research/projects/CIRBE/Calibration:Test/strontium/Sr90_REPTile2_IntegratedCPT.h5'
h5disp(h5File);
cr = h5read(h5File,'/countRate');
dt = h5readatt(h5File,'/countRate','Approximate total time (s)');
cr(cr<0) = 0; % remove count rates with negative values

%% Step 2: Read geometric factor 
h5file = '/Users/lengying/Documents/Research/projects/CIRBE/Response Function/GF_incidentE_15000k_EMZ_T0_1MeV_Tg0_671MeV.h5';
GF = h5read(h5file,'/geometricFactor'); Emin = h5readatt(h5file,'/geometricFactor/','Emin'); Emax = h5readatt(h5file,'/geometricFactor/','Emax'); dE = h5readatt(h5file,'/geometricFactor/','dE');
Tg = h5readatt(h5file,'/geometricFactor/','Tg'); T = h5readatt(h5file,'/geometricFactor/','T'); T1s = h5readatt(h5file,'/geometricFactor/','T1s'); T2s = h5readatt(h5file,'/geometricFactor/','T2s'); T3s = h5readatt(h5file,'/geometricFactor/','T3s'); T4s = h5readatt(h5file,'/geometricFactor/','T4s');
dE = 0.01;
energy = [Emin:dE:Emax]; ener = energy(1:end-1)+dE/2;

g = GF(:,61:110); cr = cr(61:110);
delx = dE;indE = ener<3; % only take GF from certain energy
%% Step 3: Setup
% Cd0d0 is the poisson error of counting rate (hence it is not just a
% simple sqrt(counts))
d0 = cr;
numr = length(cr); numx = sum(indE);
g = g(indE,:); 
% x = (ener(indE));
x = log(ener(indE));
Cm0 = zeros(numx,numx);
Cd0 = zeros(numr,numr);
j0 = zeros(numx,1);
rerr = sqrt(cr*dt+1)/dt;
for nr = 1:1:numr
   Cd0d0(nr,nr) = rerr(nr)^2; 
end

%%
fig1 = figure(1); fig1.Position = [1424 353 557 267];
plot(Sr90_energyEuramet/1000,Sr90_countsEuramet,'k-','LineWidth',2); hold on;
ylabel('Sr-90 spectrum'); xlabel('Energy (MeV)');
set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on'); set(gca,'FontSize',12);

fig2 = figure(2); fig2.Position = [1424 353 557 267];
errorbar([1:50],d0,sqrt(diag(Cd0d0)),'k*'); hold on;
ylabel('Count rates'); xlabel('Electron channels');ylim([0 2.5]);
set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on'); set(gca,'FontSize',12);
%% Step 4: Establish priori model, Cm
p0 = j0; sigmaM = 500; del = 0.5;
%sigmaM = 500; del = 0.6;
for numX = 1:1:length(x)
%     Cp0p0(:,numX) = sigmaM^2*exp(-0.5*(x-x(numX)).^2./del^2);
    Cp0p0(:,numX) = sigmaM^2*exp(-0.5*(x-x(numX)).^2./del^2);
end

for numCh = 1:1:length(cr)
    Gn(numCh,:) = (g(:,numCh).*delx);
end

%% Step 5: Use equations
Sn = (Cd0d0+Gn*Cp0p0*Gn');
p = p0+Cp0p0*Gn'*(Sn^(-1))*(d0-Gn*p0);
Cpp = Cp0p0 - Cp0p0*Gn'*(Sn)^(-1)*Gn*Cp0p0;
dp = Gn*p;
p_chi2 = sum((d0-dp).^2./rerr.^2);
%% Flux method: Bowtie
load('/Users/lengying/Documents/Research/projects/CIRBE/Response Function/Bowtie/Bowtie_GF_incidentE_15000k_EMZ_T0_1MeV_Tg0_671MeV_Claudepierre');
% Ebt = EigE(1:50); jbt = cr./gE(1:50);
gE = gEoptAvg(1:50);
Ebt = Eopt(1:50); jbt = cr./gE;
% Step 2b: Poisson distributed errors for Bowtie results
jbt_err = rerr./(cr).*jbt;
% jbt_err = sqrt(jbt);
% jbt_err = 1/gE(1:50).*rerr';
%% plot results
figure();Ek = ener(indE);
subplot(2,1,2);
% plot(eSr,400*jSr,'r-','LineWidth',2); hold on;
% plot(Sr90_energy,300*Sr90_counts,'k-','LineWidth',2); hold on;
plot(Ek,p,'r.','LineWidth',2,'MarkerSize',10); hold on;
% plot(Ebt,jbt,'b.','LineWidth',2,'MarkerSize',12); 
errorbar(Ebt,jbt,jbt_err','b.','LineWidth',2,'MarkerSize',12); hold on;
% plot(Sr90_energyEuramet/1000,5500*Sr90_countsEuramet,':','Color','k','LineWidth',2.5); 
plot(Sr90_energyEuramet/1000,5000*Sr90_countsEuramet,':','Color',[125,125,125]./255,'LineWidth',2.5); 
plot(Ek,p+sqrt(diag(Cpp)),'r:','LineWidth',1); 
plot(Ek,p-sqrt(diag(Cpp)),'r:','LineWidth',1); 
% plot(Ebt,jbt+sqrt(diag(jbt_err)),'b:','LineWidth',1); hold on;
% plot(Ebt,jbt-sqrt(diag(jbt_err)),'b:','LineWidth',1); hold on;
set(gca,'YScale','log');
leg=legend('Linear least square','Bowtie','Sr-90 spectrum x 5000'); leg.Box = 'off'; %'Linear least square - boresight'
leg.Location = 'southwest';
xlim([0,2.5]);
ylim([0.01 1000]); set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on'); 
xlabel('Energy (MeV)'); ylabel('j, #/(cm^2 sr s MeV)');
set(gca,'FontSize',12);
% title(['Boresight energy response function \sigma=',num2str(sigmaM),' \Delta=',num2str(del)]);
subplot(2,1,1);
errorbar([61:110],d0,sqrt(diag(Cd0d0)),'k*'); hold on;
plot([61:110],dp,'r*','LineWidth',1.5); hold on;
leg = legend('observed count rates','estimated count rates'); leg.Box = 'off';
leg.Location = 'northeast';
set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on'); 
xlabel('Electron channels'); ylabel('count rate'); ylim([-0.1 2]);
set(gca,'FontSize',12);xlim([60 110]);


% figure();
% plot(Ek,p,'b.-','LineWidth',2,'MarkerSize',10); hold on;
% plot(Ek,p+sqrt(diag(Cpp)),'c--','LineWidth',1); hold on;
% plot(Ek,p-sqrt(diag(Cpp)),'c--','LineWidth',1); hold on;
% set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on'); 
% xlabel('Energy (MeV)'); ylabel('j, #/cm^2 sr s MeV');
% set(gca,'FontSize',12);
