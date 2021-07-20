%% SpikePCA

% addpath('General Functions')
clear
clc
pathD = 'C:\Users\tomas\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\LegTracking\Dark\';
pathL = 'C:\Users\tomas\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\WT TB\';

% pathL = '';
params = GetParams();
% get the angular spikes matrices and PCs in light and darkness
[SPKVRD, SPKVRRandD, SPKVRDistD, SPKVFD,pcaD, ptypesD] = GetSpikePCA(pathD, params);
[SPKVRL, SPKVRRandL, SPKVRDistL, SPKVFL,pcaL, ptypesL] = GetSpikePCA(pathL, params);
%%
nd = 1; % trial type dark
nl = 5; %trial type light

% Calculate PC1 grand means 
gmD = pcaD.PC1{nd}*pcaD.NSP{nd}/sum(pcaD.NSP{nd});
semD = sqrt(((bsxfun(@minus, pcaD.PC1{nd}, gmD).*...
    bsxfun(@minus, pcaD.PC1{nd}, gmD)*pcaD.NSP{nd}))/sum(pcaD.NSP{nd}));
gmL = pcaL.PC1{nl}*pcaL.NSP{nl}/sum(pcaL.NSP{nl});
semL = sqrt(((bsxfun(@minus, pcaL.PC1{nl}, gmL).*...
    bsxfun(@minus, pcaL.PC1{nl}, gmL)*pcaL.NSP{nl}))/sum(pcaL.NSP{nl}));

% Calculate explained variance grand means
gmEL = pcaL.EXP{nl}*pcaL.NSP{nl}/sum(pcaL.NSP{nl});
gsdEL = sqrt(((bsxfun(@minus, pcaL.EXP{nl}, gmEL).*...
    bsxfun(@minus, pcaL.EXP{nl}, gmEL)*pcaL.NSP{nl}))/sum(pcaL.NSP{nl}));
gmED = pcaD.EXP{nd}*pcaD.NSP{nd}/sum(pcaD.NSP{nd});
gsdED = sqrt(((bsxfun(@minus, pcaD.EXP{nd}, gmED).*...
    bsxfun(@minus, pcaD.EXP{nd}, gmED)*pcaD.NSP{nd}))/sum(pcaD.NSP{nd}));

figure,
% plot PC1 shape
subplot(1,2,1)
t = (-20:20)/60;
hold on
plot(t,gmD, 'k', 'linewidth',3)
plot(t,gmD-semD, 'k', 'linewidth',1)
plot(t,gmD+semD, 'k', 'linewidth',1)
plot(t,gmL, 'r', 'linewidth',3)
plot(t,gmL-semL, 'r', 'linewidth',1)
plot(t,gmL+semL, 'r', 'linewidth',1)
xlabel('Time (s)')
ylabel('PC1')
axis([-0.3 0.3 -0.1 0.5])
subplot(1,2,2)
% plot Explained Variance
scatter(1.5*ones(size(pcaL.EXP{nl},2),1), pcaL.EXP{nl}(1,:), 'r')
hold on
errorbar(1.5, gmEL(1), gsdEL(1), 'color', 'r', ...
    'linewidth', 3,'marker', 'o', 'markersize', 10, 'markerfacecolor','r')
scatter(1*ones(size(pcaD.EXP{nd},2),1), pcaD.EXP{nd}(1,:), 'k')
hold on
errorbar(1, gmED(1), gsdED(1), 'color', 'k', ...
    'linewidth', 3,'marker', 'o', 'markersize', 10, 'markerfacecolor','k')
axis([0.5 2 0 100])
ylabel('Var. Explained PC1')
xt={'Dark' ; 'Light'} ; 
set(gca,'xtick',[1 1.5]); 
set(gca,'xticklabel',xt);

dd = pcaD.EXP{nd}(1,:);
ll = pcaL.EXP{nl}(1,:);
% Stat test of explained variance across conditions
p = mwwtest(dd(pcaD.NSP{nd}>20),ll(pcaL.NSP{nl}>20),0);
p.p(2)

%%
% Calculate grand mean amplitude distributions
gmDistL = pcaL.Dist{nl}'*pcaL.NSP{nl}/sum(pcaL.NSP{nl});
gsdDistL = sqrt(((bsxfun(@minus, pcaL.Dist{nl}', gmDistL).*...
    bsxfun(@minus, pcaL.Dist{nl}', gmDistL)*pcaL.NSP{nl}))/sum(pcaL.NSP{nl}))/sqrt(length(pcaL.NSP{nl}));

gmDistD = pcaD.Dist{nd}'*pcaD.NSP{nd}/sum(pcaD.NSP{nd});
gsdDistD = sqrt(((bsxfun(@minus, pcaD.Dist{nd}', gmDistD).*...
    bsxfun(@minus, pcaD.Dist{nd}', gmDistD)*pcaD.NSP{nd}))/sum(pcaD.NSP{nd}))/sqrt(length(pcaD.NSP{nd}));

Vmdar = sum(bsxfun(@times, pcaD.Dist{nd}, pcaD.VrCents),2);
Vmlig = sum(bsxfun(@times, pcaL.Dist{nl}, pcaD.VrCents),2);

mwwtest(Vmdar,Vmlig,1);

figure,
plot(pcaD.VrCents, (gmDistL), 'r', 'linewidth', 3)
hold on
plot(pcaD.VrCents, (gmDistL) + gsdDistL, 'r', 'linewidth', 1)
hold on
plot(pcaD.VrCents, (gmDistL) - gsdDistL, 'r', 'linewidth', 1)
hold on
plot(pcaD.VrCents, (gmDistD), 'k', 'linewidth', 3)
hold on
plot(pcaD.VrCents, (gmDistD) + gsdDistD, 'k', 'linewidth', 1)
hold on
plot(pcaD.VrCents, (gmDistD) - gsdDistD, 'k', 'linewidth', 1)
xlabel('Peak velocity (º/s)')
ylabel('Distribution')
axis([110 1200 0 0.15])
