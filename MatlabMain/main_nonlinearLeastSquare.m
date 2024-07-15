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
% close all;

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
cr(cr<=0) = 0; % Change for nonlinear: remove count rates with negative values

%% Step 2: Read geometric factor 
h5file = '/Users/lengying/Documents/Research/projects/CIRBE/Response Function/GF_Boresight_incidentE_5000k_EMZ_T0_1MeV_Tg0_671MeV.h5';
GF = h5read(h5file,'/geometricFactor'); Emin = h5readatt(h5file,'/geometricFactor/','Emin'); Emax = h5readatt(h5file,'/geometricFactor/','Emax'); dE = h5readatt(h5file,'/geometricFactor/','dE');
Tg = h5readatt(h5file,'/geometricFactor/','Tg'); T = h5readatt(h5file,'/geometricFactor/','T'); T1s = h5readatt(h5file,'/geometricFactor/','T1s'); T2s = h5readatt(h5file,'/geometricFactor/','T2s'); T3s = h5readatt(h5file,'/geometricFactor/','T3s'); T4s = h5readatt(h5file,'/geometricFactor/','T4s');
dE =0.01;energy = [Emin:dE:Emax]; ener = energy(1:end-1); %ener = energy(1:end-1)+dE/2;

g = GF(:,61:110); cr = cr(61:110);
indE = ener<=3; % only take GF from certain energy

%% Step 3: Setup
% Cd0d0 is the poisson error of counting rate (hence it is not just a
% simple sqrt(counts))

% setup for background (energy, geometric factor)
g = g(indE,:); x = log(ener');  x2 = log(energy); 
delx = x2(2:end)-x2(1:end-1); delx = delx(indE);
x = x(indE); numx = sum(indE);
 
% setup for priori model
j0 = zeros(numx,1)+30;  
Cm0 = zeros(numx,numx);

% setup for priori count rates/observations
C = 10^-5;
d0 = log(cr+C); numr = length(cr);
Cd0 = zeros(numr,numr);
rerr = sqrt(cr*dt+1)/dt; dsig = rerr./(cr+C);
for nr = 1:1:numr
   Cd(nr,nr) = (rerr(nr)/(cr(nr)+C)).^2; 
end

%% Step 4: Establish priori model, Cm
m0 = log(j0); mn=m0; sigmaM = 10; del = 0.5;%sigmaM = 10; del = 0.5;
for numX = 1:1:length(x)
    Cm(:,numX) = sigmaM^2*exp(-0.5*(x-x(numX)).^2./del^2);
end

for numCh = 1:1:length(cr)
    dn(numCh,1) = log(sum(g(:,numCh).*exp(mn).*dE')+C);
    rn(numCh,1) = exp(dn(numCh,1))-C;
    Gn(numCh,:) = 1/(rn(numCh)+C).*g(:,numCh).*exp(mn).*dE';
end

%% Step 5: Use equations
Sn = Cd + Gn*Cm*Gn';
B = Cm*Gn'*Sn^-1;
mn1 = m0 + Cm*Gn'*Sn^(-1)*(d0-dn+Gn*(mn-m0));
Cmn1 = Cm - Cm*Gn'*Sn^(-1)*Gn*Cm;
count = 1;
while max(abs(mn1-mn))>0.1
    mn = []; Gn = []; dn = []; Sn = []; B = [];
     mn = mn1; mn1 = []; 
    for numCh = 1:1:length(cr)
        dn(numCh,1) = log(sum(g(:,numCh).*exp(mn).*dE')+C);
        rn(numCh,1) = exp(dn(numCh,1))-C;
        Gn(numCh,:) = 1/(rn(numCh)+C).*g(:,numCh).*exp(mn).*dE';
    end
    Sn = Cd + Gn*Cm*Gn';
    B = Cm*Gn'*Sn^-1;
    mn1 = m0 + B*(d0-dn+Gn*(mn-m0));
    Cmm =  Cm - Cm*Gn'*(Sn^-1)*Gn*Cm;
    jn = exp(mn1);
    jsig = jn.*sqrt(diag(Cmm));
    figure(10);
    plot(ener(indE),j0,'k-','LineWidth',2); hold on;
    plot(ener(indE),jn,'b-','LineWidth',2); hold on;
    plot(ener(indE),jn+jsig,'c--','LineWidth',2); hold on;
    plot(ener(indE),jn-jsig,'c--','LineWidth',2); hold on;
    set(gca,'YScale','log'); xlabel('Energy (MeV)'); ylabel('j');
    set(gca,'FontSize',12);
    maxMn1 = max(mn1)
    maxM0 = max(m0)
    maxdif = max(abs(mn1-mn))
    count
%     pause;
    count = count + 1;
end
jn = exp(mn1);
%% Step 7: Compute Cm
Gn = []; 
dn = [];
Sn = [];
B = [];
for numCh = 1:1:length(cr)
    dn(numCh,1) = log(sum(g(:,numCh).*exp(mn1).*dE)+C);
    rn(numCh,1) = exp(dn(numCh,1))-C;
    Gn(numCh,:) = 1/(rn(numCh)+C).*g(:,numCh).*exp(mn1).*dE;
end
Sn = Cd + Gn*Cm*Gn';
B = Cm*Gn'*Sn^-1;
Cmm =  Cm - Cm*Gn'*(Sn^-1)*Gn*Cm;
jn = exp(mn1); jsig = jn.*sqrt(diag(Cmm));
%% Goodness of test
for df = 1:1000
Chi2(df,:) = chi2inv(1-[0.995 0.99 0.975 0.95 0.9 0.1 0.05 0.025 0.01 0.005],df);
end
disp(Chi2);
df = 1:1000;
% dfcomp = (length(dp)-1)*(length(p)-1);
v = length(dn)%-length(p);
indDf = (df==v);
chi2Exp = Chi2(indDf,9)
p_chi2 = (dn-d0)'*(Gn*Cm*Gn'+Cd)^(-1)*(dn-d0)

%% Flux method: Bowtie
load('/Users/lengying/Documents/GEANT4_2017/Analysis/Bowtie_GF_incidentE_5000k_EMZ_T0_1MeV_Tg0_671MeV');
Ebt = EigE(1:50); jbt = cr./gE(1:50);

% Step 2b: Poisson distributed errors for Bowtie results
jbt_err = sqrt(rerr)./(cr).*jbt;

%%
figure();Ek = ener(indE);
plot(Ek,exp(mn1),'b.-','LineWidth',2,'MarkerSize',10); hold on;
plot(Ek,exp(mn1+sqrt(diag(Cmm))),'b--','LineWidth',1); hold on;
plot(Ek,exp(mn1-sqrt(diag(Cmm))),'b--','LineWidth',1); hold on;

Ek = ener(indE); jn = exp(mn1); jsig = jn.*sqrt(diag(Cmm));
plot(Ek,jn,'r.-','LineWidth',2,'MarkerSize',10); hold on;
plot(Ek,jn+jsig,'m--','LineWidth',1); hold on;
plot(Ek,jn-jsig,'m--','LineWidth',1); hold on;
set(gca,'YScale','log');
%% plot results
figure();Ek = ener(indE);
jsig = jn.*sqrt(diag(Cmm));
% plot(eSr,400*jSr,'r-','LineWidth',2); hold on;
subplot(2,1,2);
% plot(Sr90_energyEuramet/1000,6000*Sr90_countsEuramet,'k-','LineWidth',2); hold on;
% plot(Ebt,jbt,'bx','LineWidth',2); hold on;
plot(Ek,jn,'m.','LineWidth',2,'MarkerSize',10); hold on;
errorbar(Ebt,jbt,jbt_err,'b.','LineWidth',2,'MarkerSize',12); 
plot(Sr90_energyEuramet/1000,5500*Sr90_countsEuramet,':','Color',[125,125,125]./255,'LineWidth',2.5); 
plot(Ek,jn-jsig,'m:','LineWidth',1); hold on;
plot(Ek,jn+jsig,'m:','LineWidth',1); hold on;
set(gca,'YScale','log');
leg=legend('Nonlinear least square','Bowtie','Sr-90 spectrum x 5000'); leg.Box = 'off'; 
leg.Location = 'southwest';
% leg=legend('Linear least square','Bowtie','Sr-90 spectrum','Nonlinear least square'); leg.Box = 'off'; %'Linear least square - boresight'
xlim([0 2.5]); set(gca,'YScale','log');
ylim([10^-2 1000]); set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on'); 
xlabel('Energy (MeV)'); ylabel('j, #/cm^2 sr s MeV');
set(gca,'FontSize',12);  
% title(['Boresight: \sigma=',num2str(sigmaM),' \Delta=',num2str(del),'  if \chi2 = ',num2str(p_chi2,'%.01f'),...
%     '<',num2str(chi2Exp,'%.01f'), '  then the fit is acceptable.']);

subplot(2,1,1);
errorbar([61:110],exp(d0)-C,rerr,'k*'); hold on;
plot([61:110],exp(dn)-C,'r*'); hold on;
leg = legend('observed count rates','estimated count rates'); leg.Box = 'off';
set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on'); 
xlabel('Electron channels'); ylabel('count rate'); ylim([-0.1 2]);
set(gca,'FontSize',12); xlim([60 110]);
