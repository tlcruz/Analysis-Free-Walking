%% SpikePCA
clear
clc
pathD = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Dark\Dark WTTB\';
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Windmill\NG-RG\';
pathL = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\WT TB\';
params = GetParams();
params.spkTempPath = 'SpikeTemplateD.mat';
params.cutoff = 0.15;
[SPKVRD, SPKVRRandD, SPKVRDistD, SPKVFD,pcaD, ptypesD] = GetSpikePCA(pathD, params);
params.spkTempPath = 'SpikeTemplateL.mat';
params.cutoff = 0.15;
[SPKVRL, SPKVRRandL, SPKVRDistL, SPKVFL,pcaL, ptypesL] = GetSpikePCA(pathL, params);
%%
nd = 1;
nl = 1;
gmD = pcaD.PC1{nd}*pcaD.NSP{nd}/sum(pcaD.NSP{nd});
semD = sqrt(((bsxfun(@minus, pcaD.PC1{nd}, gmD).*...
    bsxfun(@minus, pcaD.PC1{nd}, gmD)*pcaD.NSP{nd}))/sum(pcaD.NSP{nd}));
gmL = pcaL.PC1{nl}*pcaL.NSP{nl}/sum(pcaL.NSP{nl});
semL = sqrt(((bsxfun(@minus, pcaL.PC1{nl}, gmD).*...
    bsxfun(@minus, pcaL.PC1{nl}, gmD)*pcaL.NSP{nl}))/sum(pcaL.NSP{nl}));

gmEL = pcaL.EXP{nl}*pcaL.NSP{nl}/sum(pcaL.NSP{nl});
gsdEL = sqrt(((bsxfun(@minus, pcaL.EXP{nl}, gmEL).*...
    bsxfun(@minus, pcaL.EXP{nl}, gmEL)*pcaL.NSP{nl}))/sum(pcaL.NSP{nl}));
gmED = pcaD.EXP{nd}*pcaD.NSP{nd}/sum(pcaD.NSP{nd});
gsdED = sqrt(((bsxfun(@minus, pcaD.EXP{nd}, gmED).*...
    bsxfun(@minus, pcaD.EXP{nd}, gmED)*pcaD.NSP{nd}))/sum(pcaD.NSP{nd}));
figure,
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

p = mwwtest(dd(pcaD.NSP{nd}>20),ll(pcaL.NSP{nl}>20),0);
p.p(2)

%%
gmDistL = pcaL.Dist{nl}'*pcaL.NSP{nl}/sum(pcaL.NSP{nl});
gsdDistL = sqrt(((bsxfun(@minus, pcaL.Dist{nl}', gmDistL).*...
    bsxfun(@minus, pcaL.Dist{nl}', gmDistL)*pcaL.NSP{nl}))/sum(pcaL.NSP{nl}))/sqrt(length(pcaL.NSP{nl}));

gmDistD = pcaD.Dist{nd}'*pcaD.NSP{nd}/sum(pcaD.NSP{nd});
gsdDistD = sqrt(((bsxfun(@minus, pcaD.Dist{nd}', gmDistD).*...
    bsxfun(@minus, pcaD.Dist{nd}', gmDistD)*pcaD.NSP{nd}))/sum(pcaD.NSP{nd}))/sqrt(length(pcaD.NSP{nd}));

Vmdar = sum(bsxfun(@times, pcaD.Dist{nd}, pcaD.VrCents),2);
Vmlig = sum(bsxfun(@times, pcaL.Dist{nl}, pcaD.VrCents),2);

mwwtest(Vmdar,Vmlig,1)

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


%%
figure,
vr = [];
vf = [];
for n = 1 : size(SPKVRL,1)
    for j = 1 : size(SPKVRL{n,1}, 2)
        vr = horzcat(vr, SPKVRL{n,1}(:,j));
        vf = horzcat(vf, SPKVFL{n,1}(:,j));
    end
end

[B, I] = sort(vr(20,:));
cmap = jet(size(vr,2));
v = randi(length(I),length(I));

subplot(1,2,1)
hold on
for n = 1 : length(I)
    plot((-20:20)/60,vr(:,I(v(n))), 'color', cmap(v(n),:))
end
xlabel('Time (s)')
ylabel('Vr (º/s)')
axis([-0.3 0.3 -1500 1500])

subplot(1,2,2)
hold on
for n = 1 : length(I)
    plot((-20:20)/60,vf(:,I(v(n))), 'color', cmap(v(n),:))
end
xlabel('Time (s)')
ylabel('Vf (mm/s)')
axis([-0.3 0.3 -20 50])





% %%
% figure,
% hold on
% templ = PC1*NSP/sum(PC1*NSP);
% v = [];
% for n = 1 : 72
%     t = abs(SPKVR{nn}(:,n));
%     c = (t'*templ)/sqrt(t'*t);
%     v = vertcat(v, c);
%     
%     if c < 0.15
%         plot(t)
%     end
% end
% plot(v)

% axis([1 41 -10 1000])
% 

% save('SpikeTemplateD.mat', 'templ')




