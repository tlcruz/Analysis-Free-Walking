%% SpikePCA
addpath('Saccade PCA')
clear
clc
pathC1 = '\Split T4T5 Cnt\';
pathC2 = '\EmptySplit 10xKir\';
pathS = '\SplitT4T5 Kir\';
% get the angular spikes matrices and PCs in silenced and control
% conditions
params = GetParams();
params.spkTempPath = 'SpikeTemplateL.mat';
params.cutoff = 0.15;
[SPKVRC1, SPKVRRandC1, SPKVRDistC1, SPKVFC1, pcaC1, ptypesC1] = ...
    GetSpikePCA(pathC1, params);
[SPKVRC2, SPKVRRandC2, SPKVRDistC2, SPKVFC2, pcaC2, ptypesC2] = ...
    GetSpikePCA(pathC2, params);
[SPKVRS, SPKVRRandS, SPKVRDistS, SPKVFS, pcaS, ptypesS] = ...
    GetSpikePCA(pathS, params);
%%
% Calculate PC1 grand means 
n = 1;
gmC1 = pcaC1.PC1{n}*pcaC1.NSP{n}/sum(pcaC1.NSP{n});
semC1 = sqrt(((bsxfun(@minus, pcaC1.PC1{n}, gmC1).*...
    bsxfun(@minus, pcaC1.PC1{n}, gmC1)*pcaC1.NSP{n}))/sum(pcaC1.NSP{n}));
gmC2 = pcaC2.PC1{n}*pcaC2.NSP{n}/sum(pcaC2.NSP{n});
semC2 = sqrt(((bsxfun(@minus, pcaC2.PC1{n}, gmC2).*...
    bsxfun(@minus, pcaC2.PC1{n}, gmC2)*pcaC2.NSP{n}))/sum(pcaC2.NSP{n}));
gmS = pcaS.PC1{n}*pcaS.NSP{n}/sum(pcaS.NSP{n});
semS = sqrt(((bsxfun(@minus, pcaS.PC1{n}, gmS).*...
    bsxfun(@minus, pcaS.PC1{n}, gmS)*pcaS.NSP{n}))/sum(pcaS.NSP{n}));

% Calculate explained variance grand means
gmEC1 = pcaC1.EXP{n}*pcaC1.NSP{n}/sum(pcaC1.NSP{n});
gsdEC1 = sqrt(((bsxfun(@minus, pcaC1.EXP{n}, gmEC1).*...
    bsxfun(@minus, pcaC1.EXP{n}, gmEC1)*pcaC1.NSP{n}))/sum(pcaC1.NSP{n}));
gmEC2 = pcaC2.EXP{n}*pcaC2.NSP{n}/sum(pcaC2.NSP{n});
gsdEC2 = sqrt(((bsxfun(@minus, pcaC2.EXP{n}, gmEC2).*...
    bsxfun(@minus, pcaC2.EXP{n}, gmEC2)*pcaC2.NSP{n}))/sum(pcaC2.NSP{n}));
gmES = pcaS.EXP{n}*pcaS.NSP{n}/sum(pcaS.NSP{n});
gsdES = sqrt(((bsxfun(@minus, pcaS.EXP{n}, gmES).*...
    bsxfun(@minus, pcaS.EXP{n}, gmES)*pcaS.NSP{n}))/sum(pcaS.NSP{n}));

figure,
% plot PC1 shape
subplot(1,2,1)
t = (-20:20)/60;
hold on
plot(t,gmC1, 'color', [0 0 0.7], 'linewidth',3)
plot(t,gmC1-semC1, 'color', [0 0 0.7], 'linewidth',1)
plot(t,gmC1+semC1, 'color', [0 0 0.7], 'linewidth',1)
plot(t,gmC2, 'color', [0.3 0.3 1], 'linewidth',3)
plot(t,gmC2-semC2, 'color', [0.3 0.3 1], 'linewidth',1)
plot(t,gmC2+semC2, 'color', [0.3 0.3 1], 'linewidth',1)
plot(t,gmS, 'color', [1 0.5 0], 'linewidth',3)
plot(t,gmS-semS, 'color', [1 0.5 0], 'linewidth',1)
plot(t,gmS+semS, 'color', [1 0.5 0], 'linewidth',1)
xlabel('Time (s)')
ylabel('PC1')
title(['T4T5 ' ptypesS{n}])
axis([-0.3 0.3 -0.1 0.5])

% plot Explained Variance
subplot(1,2,2)
hold on
errorbar(1, gmEC1(1), gsdEC1(1), 'color', [0 0 0.7], ...
    'linewidth', 3,'marker', 'o', 'markersize', 10, 'markerfacecolor',[0 0 0.7])
errorbar(1.5, gmES(1), gsdES(1), 'color', [1 0.5 0], ...
    'linewidth', 3,'marker', 'o', 'markersize', 10, 'markerfacecolor',[1 0.5 0])
errorbar(2, gmEC2(1), gsdEC2(1), 'color', [0.3 0.3 1], ...
    'linewidth', 3,'marker', 'o', 'markersize', 10, 'markerfacecolor',[0.3 0.3 1])
axis([0.5 2.5 0 100])
ylabel('Var. Explained PC1')
xt={'Control 1' ; 'Silenced'; 'Control 2'} ; 
set(gca,'xtick',[1 1.5 2]); 
set(gca,'xticklabel',xt);

% Stat test of explained variance across conditions
c1e = pcaC1.EXP{n}(1,:);
c2e = pcaC2.EXP{n}(1,:);
se = pcaS.EXP{n}(1,:);
pS1 = mwwtest(se(pcaS.NSP{n}>20),c1e(pcaC1.NSP{n}>20),0);
pS2 = mwwtest(se(pcaS.NSP{n}>20),c2e(pcaC2.NSP{n}>20),0);
p12 = mwwtest(c2e(pcaC2.NSP{n}>20),c1e(pcaC1.NSP{n}>20),0);

title(['p(S-C1):' num2str(0.0001*floor(10000*pS1.p(2))) ...
    ' p(S-C2):' num2str(0.0001*floor(10000*pS2.p(2))) ...
    ' p(C1-C2):' num2str(0.0001*floor(10000*p12.p(2))) ])