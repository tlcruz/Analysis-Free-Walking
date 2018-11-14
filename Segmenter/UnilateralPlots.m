path1 = 'C:\Users\Tomas\Dropbox (Sensorimotor)\ChiappeLabNew\data\TOMÁS\Free Walking VR\Data\Unilateral Balanced\WTTB\';
path2 = 'C:\Users\tomas\Dropbox (Sensorimotor)\ChiappeLabNew\data\TOMÁS\Free Walking VR\Data\Unilateral Stimulus\WT TB\';
% [ot1] = UnilPSSVals(path1, 2);
% [ot2] = UnilPSSVals(path2, 1);
% [ot1] = UnilPSSValsVr(path1, 2);
% [ot2] = UnilPSSValsVr(path2, 1);
[ot1] = UnilPSSValsVrForwBouts(path1, 2);
[ot2] = UnilPSSValsVrForwBouts(path2, 1);
ot1s = ot1;
ot2s = ot2;
%%
clear
clc
pathExampleNG = 'C:\Users\tomas\Dropbox (Sensorimotor)\ChiappeLabNew\data\TOMÁS\Free Walking VR\Data\Unilateral Stimulus\WT TB\Fly 1886\DataLowRes.mat';
dtNG = load(pathExampleNG);
dtNG = dtNG.Flies;
partsNG = [7];

pathExampleRG = 'C:\Users\tomas\Dropbox (Sensorimotor)\ChiappeLabNew\data\TOMÁS\Free Walking VR\Data\Unilateral Stimulus\WT TB\Fly 1886\DataLowRes.mat';
dtRG = load(pathExampleRG);
dtRG = dtRG.Flies;
partsRG = [8];
%

tiNG = dtNG.Protocol(partsNG(1));
tfNG = dtNG.Protocol(partsNG(end)+1);
tNG = tiNG:tfNG;
tiRG = dtRG.Protocol(partsRG(1));
tfRG = dtRG.Protocol(partsRG(end)+1);
tRG = tiRG:tfRG;

VrNG = [];
VfNG = [];
VsNG = [];
VtNG = [];
XNG = [];
YNG = [];
actstNG = [];
for n = 1 : length(partsNG)
    VrNG = vertcat(VrNG, dtNG.Data{partsNG(n)}.Vr);
    VfNG = vertcat(VfNG, dtNG.Data{partsNG(n)}.Vf);
    VsNG = vertcat(VsNG, dtNG.Data{partsNG(n)}.Vs);
    VtNG = vertcat(VtNG, dtNG.Data{partsNG(n)}.Vt);
    XNG = vertcat(XNG, dtNG.Data{partsNG(n)}.X);
    YNG = vertcat(YNG, dtNG.Data{partsNG(n)}.Y);
    actstNG = vertcat(actstNG, dtNG.Data{partsNG(n)}.actState);
end

VrRG = [];
VfRG = [];
VsRG = [];
VtRG = [];
XRG = [];
YRG = [];
actstRG = [];
for n = 1 : length(partsRG)
    VrRG = vertcat(VrRG, dtRG.Data{partsRG(n)}.Vr);
    VfRG = vertcat(VfRG, dtRG.Data{partsRG(n)}.Vf);
    VsRG = vertcat(VsRG, dtRG.Data{partsRG(n)}.Vs);
    VtRG = vertcat(VtRG, dtRG.Data{partsRG(n)}.Vt);
    XRG = vertcat(XRG, dtRG.Data{partsRG(n)}.X);
    YRG = vertcat(YRG, dtRG.Data{partsRG(n)}.Y);
    actstRG = vertcat(actstRG, dtRG.Data{partsRG(n)}.actState);
end

clear i nSeqs path
[params] = GetParams();
[locs, ~, cmhSong, thr, dx] = SortSpikes(VrNG, actstNG, params);
[locsF, pksIF, pksEF, ~] = CullSpikesBasedOnVf(VrNG, locs, cmhSong, thr, VfNG, params);
[locsFNG, pksIFNG, pksEFNG, scores] = TemplateCompSpike(VrNG,locsF,pksIF, pksEF, cmhSong, thr, params);
[FBoutsNG, FBoutsNGAll] = GetFbouts(actstNG, VfNG, locsFNG, pksIFNG,pksEFNG, params);

[params] = GetParams();
[locs, ~, cmhSong, thr] = SortSpikes(VrRG, actstRG, params);
[locsF, pksIF,pksEF, ~] = CullSpikesBasedOnVf(VrRG, locs, cmhSong, thr, VfRG, params);
[locsFRG, pksIFRG, pksEFRG, ~] = TemplateCompSpike(VrRG,locsF,pksIF, pksEF, cmhSong, thr, params);
[FBoutsRG, FBoutsRGAll] = GetFbouts(actstRG, VfRG, locsFRG, pksIFRG, pksEFRG, params);

%%
figure,
subplot(2,3,1)
ttNG = tNG((60*60):(75*60))'-tNG(1);
plot(XNG(ttNG), YNG(ttNG), 'k', 'linewidth', 2)
axis square
axis([-40 40 -40 40])
axis off

subplot(2,3,4)
ttRG = tRG((70*60):(85*60))'-tRG(1);
plot(XRG(ttRG), YRG(ttRG), 'k', 'linewidth', 2)
axis square
axis([-40 40 -40 40])
axis off

subplot(2,3,[2 3])
plot(ttNG/60, VrNG(ttNG), 'k', 'linewidth', 2)
hold on
for i = 1 : length(FBoutsNG)
    plot((FBoutsNG{i})/60, VrNG((FBoutsNG{i})),'r', 'linewidth', 2)
    hold on
end
for i = 1 : length(locsFNG)
    dt = ((locsFNG(i)-pksIFNG(i)):(locsFNG(i)+pksEFNG(i)));
    plot(dt/60, VrNG(dt),'b', 'linewidth', 2)
    hold on
end
plot(ttNG/60, 10*VfNG(ttNG)-1500, 'k', 'linewidth', 1)
for i = 1 : length(FBoutsNG)
    plot((FBoutsNG{i})/60, 10*VfNG((FBoutsNG{i}))-1500,'r', 'linewidth', 1)
    hold on
end
for i = 1 : length(locsFNG)
    dt = ((locsFNG(i)-pksIFNG(i)):(locsFNG(i)+pksEFNG(i)));
    plot(dt/60, 10*VfNG(dt)-1500,'b', 'linewidth', 1)
    hold on
end
axis([ttNG(1)/60 ttNG(end)/60 -1600 1000])
box off

subplot(2,3,[5 6])
plot(ttRG/60, VrRG(ttRG), 'k', 'linewidth', 2)
hold on
for i = 1 : length(FBoutsRGAll)
    plot((FBoutsRGAll{i})/60, VrRG((FBoutsRGAll{i})),'r', 'linewidth', 2)
    hold on
end
for i = 1 : length(locsFRG)
    dt = ((locsFRG(i)-pksIFRG(i)):(locsFRG(i)+pksEFRG(i)));
    plot(dt/60, VrRG(dt),'b', 'linewidth', 2)
    hold on
end
plot(ttRG/60, 10*VfRG(ttRG)-1500, 'k', 'linewidth', 1)
for i = 1 : length(FBoutsRGAll)
    plot((FBoutsRGAll{i})/60, 10*VfRG((FBoutsRGAll{i}))-1500,'r', 'linewidth', 1)
    hold on
end
for i = 1 : length(locsFRG)
    dt = ((locsFRG(i)-pksIFRG(i)):(locsFRG(i)+pksEFRG(i)));
    plot(dt/60, 10*VfRG(dt)-1500,'b', 'linewidth', 1)
    hold on
end
axis([ttRG(1)/60 ttRG(end)/60 -1600 1000])
box off

%%
ot1 = ot1s;
ot2 = ot2s;

ot2.PSSV = cell(6, size(ot2.PSSF,2), size(ot2.PSSF,3));
ot2.NPV = cell(6, size(ot2.PSSF,2), size(ot2.PSSF,3));
ot1.PSSV = cell(6, size(ot1.PSSF,2), size(ot1.PSSF,3));
ot1.NPV = cell(6, size(ot1.PSSF,2), size(ot1.PSSF,3));

if length(ot1.pTypes) == 12
    vectOr = [5 6 1 2 3 4];
    ot1.PSSV = cell(6, size(ot1.PSSF,2), size(ot1.PSSF,3));
    ot1.NPV = cell(6, size(ot1.PSSF,2), size(ot1.PSSF,3));
    for i = 1 : 6
        for j = 1 : size(ot1.PSSF,2)
            for k = 1 : size(ot1.PSSF,3)
                 ot1.PSSV{vectOr(i),j,k} = vertcat(ot1.PSSF{2*i-1,j,k}, ot1.PSSF{2*i,j,k});
                 ot1.NPV{vectOr(i),j,k} = vertcat(ot1.NPF{2*i-1,j,k}, ot1.NPF{2*i,j,k});
            end
        end
    end
end
if length(ot2.pTypes) == 8
    ot2.PSSV = cell(6, size(ot2.PSSF,2), size(ot2.PSSF,3));
    ot2.NPV = cell(6, size(ot2.PSSF,2), size(ot2.PSSF,3));
    vectOr = [5 6 1 2 3 4];
    vectIr = [1 2 5 6 7 8];
    for i = 1:6
        for j = 1 : size(ot2.PSSF,2)
            for k = 1 : size(ot2.PSSF,3)
                ot2.PSSV{vectOr(i),j,k} = ot2.PSSF{vectIr(i),j,k};
                ot2.NPV{vectOr(i),j,k} = ot2.NPF{vectIr(i),j,k};
            end
        end
    end
end

PSS = cell(6, 2, 2);
NP = cell(6, 2, 2);
for i = 1 : 6
    for j = 1 : size(ot2.PSSV,2)
        for k = 1 : size(ot2.PSSV,3)
            PSS{i,j,k} = vertcat(ot1.PSSV{i,j,k}, ot2.PSSV{i,j,k});%vertcat(ot2.PSSV{i,j,k});%, ot2.PSSV{i,j,k});%
            NP{i,j,k} = vertcat(ot1.NPV{i,j,k}, ot2.NPV{i,j,k});%vertcat(ot2.NPV{i,j,k});%, ot2.NPV{i,j,k});%
        end
    end
end

%

figure,
hold on
plot([0 4], [0.5 0.5], '--g', 'linewidth', 1.5)
vDiff = [];
evDiff = [];
for i = 1 : 6
    if mod(i,2) == 1
        pssFTB = vertcat(PSS{i, 1, 1}, PSS{i, 2, 2});
        npFTB = vertcat(NP{i, 1, 1}, NP{i, 2, 2});
        npFTB = npFTB(~isnan(pssFTB));
        pssFTB = pssFTB(~isnan(pssFTB));
        gmFTB = pssFTB'*npFTB/sum(npFTB);
        semFTB = sqrt(((pssFTB- gmFTB).*(pssFTB- gmFTB))'*...
            npFTB/sum(npFTB))/sqrt(50);
        pssBTF = vertcat(PSS{i, 2, 1}, PSS{i, 1, 2});
        npBTF = vertcat(NP{i, 2, 1}, NP{i, 1, 2});
        npBTF = npBTF(~isnan(pssBTF));
        pssBTF = pssBTF(~isnan(pssBTF));
        gmBTF = pssBTF'*npBTF/sum(npBTF);
        semBTF = sqrt(((pssBTF- gmBTF).*(pssBTF- gmBTF))'*...
            npBTF/sum(npBTF))/sqrt(50);
        
        vDiff = vertcat(vDiff, (gmFTB-gmBTF)/(gmFTB+gmBTF));
        evDiff = vertcat(evDiff, sqrt(((2*gmBTF/((gmFTB+gmBTF)^2))*semFTB)^2 +((-2*gmFTB/((gmFTB+gmBTF)^2))*semBTF)^2));
        errorbar(ceil(i/2), gmFTB, semFTB, 'o', 'color', [1 0 0], 'linewidth', 3)
        errorbar(ceil(i/2), gmBTF, semBTF, 'o', 'color', [0 0 1], 'linewidth', 3)
    else
        pssFTB = vertcat(PSS{i, 2, 1}, PSS{i, 1, 2});
        npFTB = vertcat(NP{i, 2, 1}, NP{i, 1, 2});
        npFTB = npFTB(~isnan(pssFTB));
        pssFTB = pssFTB(~isnan(pssFTB));
        gmFTB = pssFTB'*npFTB/sum(npFTB);
        semFTB = sqrt(((pssFTB- gmFTB).*(pssFTB- gmFTB))'*...
            npFTB/sum(npFTB))/sqrt(50);
        pssBTF = vertcat(PSS{i, 1, 1}, PSS{i, 2, 2});
        npBTF = vertcat(NP{i, 1, 1}, NP{i, 2, 2});
        npBTF = npBTF(~isnan(pssBTF));
        pssBTF = pssBTF(~isnan(pssBTF));
        gmBTF = pssBTF'*npBTF/sum(npBTF);
        semBTF = sqrt(((pssBTF- gmBTF).*(pssBTF- gmBTF))'*...
            npBTF/sum(npBTF))/sqrt(50);
        
        vDiff = vertcat(vDiff, (gmFTB-gmBTF)/(gmFTB+gmBTF));
        evDiff = vertcat(evDiff, sqrt(((2*gmBTF/((gmFTB+gmBTF)^2))*semFTB)^2 +((-2*gmFTB/((gmFTB+gmBTF)^2))*semBTF)^2));
        errorbar(ceil(i/2), gmFTB, semFTB, 'o', 'color', [0.5 0 0], 'linewidth', 3)
        errorbar(ceil(i/2), gmBTF, semBTF, 'o', 'color', [0 0 0.5], 'linewidth', 3)
    end
end
% axis([0.5 3.5 0 120])
% axis([0 4 0.3 0.8])
xlabel('Dot Size (º)')
% ylabel('PSS|VisFeedback')
ylabel('<Vr(º/s)>|VisFeedback - FBouts')


figure,

%%
ot1 = ot1s;
ot2 = ot2s;

ot2.SPKPROMV = cell(6, size(ot2.PSSF,2), size(ot2.PSSF,3));
ot2.NSPKV = cell(6, size(ot2.PSSF,2), size(ot2.PSSF,3));
ot1.SPKPROMV = cell(6, size(ot1.PSSF,2), size(ot1.PSSF,3));
ot1.NSPKV = cell(6, size(ot1.PSSF,2), size(ot1.PSSF,3));

if length(ot1.pTypes) == 12
    vectOr = [5 6 1 2 3 4];
    ot1.SPKPROMV = cell(6, size(ot1.SPKSZ,2), size(ot1.SPKSZ,3));
    ot1.NSPKV = cell(6, size(ot1.NSPK,2), size(ot1.NSPK,3));
    for i = 1 : 6
        for j = 1 : size(ot1.SPKSZ,2)
            for k = 1 : size(ot1.SPKSZ,3)
%                  ot1.SPKPROMV{vectOr(i),j,k} = vertcat(ot1.SPKSZ{2*i-1,j,k});%,...
% %                      ot1.SPKSZ{2*i,j,k});
%                  ot1.NSPKV{vectOr(i),j,k} = vertcat(ot1.NSPK{2*i-1,j,k});%,...
% %                      ot1.NSPK{2*i,j,k});
                 ot1.SPKPROMV{vectOr(i),j,k} = vertcat(ot1.SPKSZ{2*i-1,j,k},...
                     ot1.SPKSZ{2*i,j,k});
                 ot1.NSPKV{vectOr(i),j,k} = vertcat(ot1.NSPK{2*i-1,j,k},...
                     ot1.NSPK{2*i,j,k});
            end
        end
    end
end
if length(ot2.pTypes) == 8
    ot2.SPKPROMV = cell(6, size(ot2.SPKSZ,2), size(ot2.SPKSZ,3));
    ot2.NSPKV = cell(6, size(ot2.NSPK,2), size(ot2.NSPK,3));
    vectOr = [5 6 1 2 3 4];
    vectIr = [1 2 5 6 7 8];
    for i = 1:6
        for j = 1 : size(ot2.SPKSZ,2)
            for k = 1 : size(ot2.SPKSZ,3)
                ot2.SPKPROMV{vectOr(i),j,k} = ot2.SPKSZ{vectIr(i),j,k};
                ot2.NSPKV{vectOr(i),j,k} = ot2.NSPK{vectIr(i),j,k};
            end
        end
    end
end

SSPK = cell(6, 2, 2);
NSPK = cell(6, 2, 2);
for i = 1 : 6
    for j = 1 : size(ot2.SPKPROMV,2)
        for k = 1 : size(ot2.SPKPROMV,3)
            SSPK{i,j,k} = vertcat(ot1.SPKPROMV{i,j,k}, ot2.SPKPROMV{i,j,k});%vertcat(ot2.PSSV{i,j,k});%, ot2.PSSV{i,j,k});%
            NSPK{i,j,k} = vertcat(ot1.NSPKV{i,j,k}, ot2.NSPKV{i,j,k});%vertcat(ot2.NPV{i,j,k});%, ot2.NPV{i,j,k});%
        end
    end
end

% SSPK = NSPK;
% NSPK = NP;
%
figure,
hold on
% plot([0 4], [0.5 0.5], '--g', 'linewidth', 1.5)
for i = 1 : 6
    if mod(i,2) == 1
        pssFTB = vertcat(SSPK{i, 1, 1}, SSPK{i, 2, 2});
        npFTB = vertcat(NSPK{i, 1, 1}, NSPK{i, 2, 2});
        indsnanFTB = find(isnan(pssFTB) == 1 | npFTB < 40);
        pssBTF = vertcat(SSPK{i, 2, 1}, SSPK{i, 1, 2});
        npBTF = vertcat(NSPK{i, 2, 1}, NSPK{i, 1, 2});
        indsnanBTF = find(isnan(pssBTF) == 1 | npBTF < 40);
        indsnan = union(indsnanFTB, indsnanBTF);
        npFTB(indsnan) = [];
        pssFTB(indsnan) = [];
        npBTF(indsnan) = [];
        pssBTF(indsnan) = [];
        
%         npFTB = npFTB(~isnan(pssFTB));
%         pssFTB = pssFTB(~isnan(pssFTB));
        gmFTB = pssFTB'*npFTB/sum(npFTB);
        semFTB = sqrt(((pssFTB- gmFTB).*(pssFTB- gmFTB))'*...
            npFTB/sum(npFTB))/sqrt(50);

%         npBTF = npBTF(~isnan(pssBTF));
%         pssBTF = pssBTF(~isnan(pssBTF));
        gmBTF = pssBTF'*npBTF/sum(npBTF);
        semBTF = sqrt(((pssBTF- gmBTF).*(pssBTF- gmBTF))'*...
            npBTF/sum(npBTF))/sqrt(50);
%         scatter((0.1+ceil(i/2))*ones(length(pssFTB),1), pssFTB, 50, [1 0 0])
%         scatter((-0.1+ceil(i/2))*ones(length(pssBTF),1), pssBTF, 50, [0 0 1])
%         plot([(0.1+ceil(i/2)) (-0.1+ceil(i/2))],[pssFTB pssBTF]', 'k')
        errorbar(ceil(i/2), gmFTB, semFTB, 'o', 'color', [1 0 0], 'linewidth', 3)
        errorbar(ceil(i/2), gmBTF, semBTF, 'o', 'color', [0 0 1], 'linewidth', 3)
%         errorbar(ceil(i/2), mean(pssFTB), std(pssFTB)/sqrt(53), 'o', 'color', [1 0 0], 'linewidth', 3)
%         errorbar(ceil(i/2), mean(pssBTF), std(pssBTF)/sqrt(53), 'o', 'color', [0 0 1], 'linewidth', 3)
    else
        pssFTB = vertcat(SSPK{i, 2, 1}, SSPK{i, 1, 2});
        npFTB = vertcat(NSPK{i, 2, 1}, NSPK{i, 1, 2});
        indsnanFTB = find(isnan(pssFTB) == 1 | npFTB < 40);
        pssBTF = vertcat(SSPK{i, 1, 1}, SSPK{i, 2, 2});
        npBTF = vertcat(NSPK{i, 1, 1}, NSPK{i, 2, 2});
        indsnanBTF = find(isnan(pssBTF) == 1 | npBTF < 40);
        indsnan = union(indsnanFTB, indsnanBTF);
        npFTB(indsnan) = [];
        pssFTB(indsnan) = [];
        npBTF(indsnan) = [];
        pssBTF(indsnan) = [];
        
        
        
        
%         pssFTB = vertcat(SSPK{i, 2, 1}, SSPK{i, 1, 2});
%         npFTB = vertcat(NSPK{i, 2, 1}, NSPK{i, 1, 2});
%         npFTB = vertcat(NSPK{i, 2, 1}, NSPK{i, 1, 2});
%         pssFTB = pssFTB./npFTB;
%         npFTB = npFTB(~isnan(pssFTB));
%         pssFTB = pssFTB(~isnan(pssFTB));
        gmFTB = pssFTB'*npFTB/sum(npFTB);
        semFTB = sqrt(((pssFTB- gmFTB).*(pssFTB- gmFTB))'*...
            npFTB/sum(npFTB))/sqrt(50);
%         pssBTF = vertcat(SSPK{i, 1, 1}, SSPK{i, 2, 2});
%         npBTF = vertcat(NSPK{i, 1, 1}, NSPK{i, 2, 2});
%         npBTF = vertcat(NSPK{i, 1, 1}, NSPK{i, 2, 2});
%         pssBTF = pssBTF./npBTF;
        npBTF = npBTF(~isnan(pssBTF));
        pssBTF = pssBTF(~isnan(pssBTF));
        gmBTF = pssBTF'*npBTF/sum(npBTF);
        semBTF = sqrt(((pssBTF- gmBTF).*(pssBTF- gmBTF))'*...
            npBTF/sum(npBTF))/sqrt(50);
        
        errorbar(ceil(i/2), gmFTB, semFTB, 'o', 'color', [0.5 0 0], 'linewidth', 3)
        errorbar(ceil(i/2), gmBTF, semBTF, 'o', 'color', [0 0 0.5], 'linewidth', 3)

%         errorbar(ceil(i/2), mean(pssFTB), std(pssFTB)/sqrt(53), 'o', 'color', [.5 0 0], 'linewidth', 3)
%         errorbar(ceil(i/2), mean(pssBTF), std(pssBTF)/sqrt(53), 'o', 'color', [0 0 .5], 'linewidth', 3)
    end
end
% axis([0.5 3.5 0 120])
axis([0 4 200 500])
% axis([0 4 0 20])
% axis([0 4 10 90])
xlabel('Dot Size (º)')
% ylabel('Size Spike Turns (º/s)')
ylabel('Spike Prominence (a.u.)')
% ylabel('# Spike Turns')
% ylabel('<Vr(º/s)>|VisFeedback - FBouts')


%%
STRR = zeros(3,1);
SEMSTR = zeros(3,1);
NSTR = zeros(3,1);
if length(ot1.pTypes) == 12
    vectIr = [5 9 1];
    for i = 1 : 3 %size(ot1.STR,1)
        na = [];
        stra = [];
        for n = 1 : size(ot1.STR,2)
            if~isempty(ot1.STR{vectIr(i),n}) && ~isempty(ot1.STR{vectIr(i)+1,n})
                na = vertcat(na, sum(ot1.N{vectIr(i),n}));
                na = vertcat(na, sum(ot1.N{vectIr(i)+1,n}));
                stra = vertcat(stra, ot1.STR{vectIr(i),n}'*ot1.N{vectIr(i),n}/sum(ot1.N{vectIr(i),n}));
                stra = vertcat(stra, ot1.STR{vectIr(i)+1,n}'*ot1.N{vectIr(i)+1,n}/sum(ot1.N{vectIr(i)+1,n}));
            end
        end
        str = stra'*na/sum(na);
        semStr = sqrt(((stra-str).*(stra-str))'*na/sum(na))/sqrt(length(na));
        STRR(i) = str;
        SEMSTR(i) = semStr;
        NSTR(i) = sum(na);
    end
end
if length(ot2.pTypes) == 8
    vectIr = [5 7 1];
    for i = 1 : 3
        na = [];
        stra = [];
        for n = 1 : size(ot2.STR,2)
            if~isempty(ot2.STR{vectIr(i),n})
                na = vertcat(na, sum(ot2.N{vectIr(i),n}));
                stra = vertcat(stra, ot2.STR{vectIr(i),n}'*ot2.N{vectIr(i),n}/sum(ot2.N{vectIr(i),n}));
            end
        end
        str = stra'*na/sum(na);
        semStr = sqrt(((stra-str).*(stra-str))'*na/sum(na))/sqrt(length(na));
        
        STRR(i) = (STRR(i)*NSTR(i) + str*sum(na))/(NSTR(i) + sum(na));
        SEMSTR(i) = (SEMSTR(i)*NSTR(i) + semStr*sum(na))/(NSTR(i) + sum(na));
        NSTR(i) = NSTR(i) + sum(na);
        
    end
end

%%



figure,

subplot(1,2,1)
hold on
for j = 1 : size(ot1.STR,1)
    na = [];
    stra = [];
    for n = 1 : size(ot1.STR,2)
        if~isempty(ot1.STR{j,n})
            na = vertcat(na, sum(ot1.N{j,n}));
            stra = vertcat(stra, ot1.STR{j,n}'*ot1.N{j,n}/sum(ot1.N{j,n}));
        end
    end
    str = stra'*na/sum(na);
    scatter(ceil(j/2), str, 100, 'k')
end
subplot(1,2,2)
hold on
for j = 1 : size(ot2.STR,1)
    na = [];
    stra = [];
    for n = 1 : size(ot2.STR,2)
        if~isempty(ot2.STR{j,n})
            na = vertcat(na, sum(ot2.N{j,n}));
            stra = vertcat(stra, ot2.STR{j,n}'*ot2.N{j,n}/sum(ot2.N{j,n}));
        end
    end
    str = stra'*na/sum(na);
    scatter(ceil(j/2), str, 100, 'k')
end






























%%
path = 'C:\Users\Tomas\Dropbox (Sensorimotor)\ChiappeLabNew\data\TOMÁS\Free Walking VR\Data\Unilateral Stimulus\R39VT05\LexAPOKir\';
% path = 'C:\Users\Tomas\Dropbox (Sensorimotor)\ChiappeLabNew\data\TOMÁS\Free Walking VR\Data\Unilateral Stimulus\R39VT05\R39VT05Cnt\';
% path = 'C:\Users\Tomas\Dropbox (Sensorimotor)\ChiappeLabNew\data\TOMÁS\Free Walking VR\Data\Unilateral Stimulus\R39VT05\R39VT05Kir\';

% [ot] = UnilPSSVals(path, 1);
[ot] = UnilPSSValsVr(path, 1);
%%
ot.PSSV = cell(6, size(ot.PSSF,2), size(ot.PSSF,3));
ot.NPV = cell(6, size(ot.PSSF,2), size(ot.PSSF,3));
vectOr = [5 6 1 2 3 4];
vectIr = [1 2 5 6 7 8];
for i = 1:6
    for j = 1 : size(ot.PSSF,2)
        for k = 1 : size(ot.PSSF,3)
            ot.PSSV{vectOr(i),j,k} = ot.PSSF{vectIr(i),j,k};
            ot.NPV{vectOr(i),j,k} = ot.NPF{vectIr(i),j,k};
        end
    end
end
PSS = cell(6, 2, 2);
NP = cell(6, 2, 2);
for i = 1 : 6
    for j = 1 : size(ot.PSSV,2)
        for k = 1 : size(ot.PSSV,3)
            PSS{i,j,k} = vertcat(ot.PSSV{i,j,k});%vertcat(ot2.PSSV{i,j,k});%, ot2.PSSV{i,j,k});%
            NP{i,j,k} = vertcat(ot.NPV{i,j,k});%vertcat(ot2.NPV{i,j,k});%, ot2.NPV{i,j,k});%
        end
    end
end
figure,
hold on
plot([0 4], [0.5 0.5], '--g', 'linewidth', 1.5)
for i = 1 : 6
    if mod(i,2) == 1
        pssFTB = vertcat(PSS{i, 1, 1}, PSS{i, 2, 2});
        npFTB = vertcat(NP{i, 1, 1}, NP{i, 2, 2});
        npFTB = npFTB(~isnan(pssFTB));
        pssFTB = pssFTB(~isnan(pssFTB));
        gmFTB = pssFTB'*npFTB/sum(npFTB);
        semFTB = sqrt(((pssFTB- gmFTB).*(pssFTB- gmFTB))'*...
            npFTB/sum(npFTB))/sqrt(50);
        pssBTF = vertcat(PSS{i, 2, 1}, PSS{i, 1, 2});
        npBTF = vertcat(NP{i, 2, 1}, NP{i, 1, 2});
        npBTF = npBTF(~isnan(pssBTF));
        pssBTF = pssBTF(~isnan(pssBTF));
        gmBTF = pssBTF'*npBTF/sum(npBTF);
        semBTF = sqrt(((pssBTF- gmBTF).*(pssBTF- gmBTF))'*...
            npBTF/sum(npBTF))/sqrt(50);
        
        errorbar(ceil(i/2), gmFTB, semFTB, 'o', 'color', [1 0 0], 'linewidth', 3)
        errorbar(ceil(i/2), gmBTF, semBTF, 'o', 'color', [0 0 1], 'linewidth', 3)
    else
        pssFTB = vertcat(PSS{i, 2, 1}, PSS{i, 1, 2});
        npFTB = vertcat(NP{i, 2, 1}, NP{i, 1, 2});
        npFTB = npFTB(~isnan(pssFTB));
        pssFTB = pssFTB(~isnan(pssFTB));
        gmFTB = pssFTB'*npFTB/sum(npFTB);
        semFTB = sqrt(((pssFTB- gmFTB).*(pssFTB- gmFTB))'*...
            npFTB/sum(npFTB))/sqrt(50);
        pssBTF = vertcat(PSS{i, 1, 1}, PSS{i, 2, 2});
        npBTF = vertcat(NP{i, 1, 1}, NP{i, 2, 2});
        npBTF = npBTF(~isnan(pssBTF));
        pssBTF = pssBTF(~isnan(pssBTF));
        gmBTF = pssBTF'*npBTF/sum(npBTF);
        semBTF = sqrt(((pssBTF- gmBTF).*(pssBTF- gmBTF))'*...
            npBTF/sum(npBTF))/sqrt(50);
        
        errorbar(ceil(i/2), gmFTB, semFTB, 'o', 'color', [0.7 0 0], 'linewidth', 3)
        errorbar(ceil(i/2), gmBTF, semBTF, 'o', 'color', [0 0 0.7], 'linewidth', 3)
    end
end
% axis([0.5 3.5 0 120])
axis([0.5 3.5 50 250])
% axis([0 4 0.3 0.8])
xlabel('Dot Size (º)')
% ylabel('PSS|VisFeedback')
ylabel('<Vr(º/s)>|VisFeedback - FBouts')
%%
























figure,

% subplot(1,2,1)
hold on
plot([0 4], [0.5 0.5], '--g', 'linewidth', 1.5)
for i = 1 : length(pSeqNG)
% + NG
    gmRp = PS{pSeqNG(i),stR}'*NP{pSeqNG(i),stR}/sum(NP{pSeqNG(i),stR});
    np = sum(NP{pSeqNG(i),stR});
    semRp = sqrt(((PS{pSeqNG(i),stR}- gmRp).*(PS{pSeqNG(i),1}- gmRp))'*...
        NP{pSeqNG(i),stR}/sum(NP{pSeqNG(i),stR}))/sqrt(length(NP{pSeqNG(i),stR}));
    gmRm = PS{mSeqNG(i),stR}'*NP{mSeqNG(i),stR}/sum(NP{mSeqNG(i),stR});
    nm = sum(NP{mSeqNG(i),stR});
    semRm = sqrt(((PS{mSeqNG(i),stR}- gmRm).*(PS{mSeqNG(i),stR}- gmRm))'*...
        NP{mSeqNG(i),stR}/sum(NP{mSeqNG(i),stR}))/sqrt(length(NP{mSeqNG(i),stR}));
    gmR = (gmRp*np+gmRm*nm)/(np+nm);
    semR = (semRp*np+semRm*nm)/(np+nm);
    errorbar(i, gmR, semR, 'o', 'color', [0 0 1], 'linewidth', 3)
    gmLp = PS{pSeqNG(i),stL}'*NP{pSeqNG(i),stL}/sum(NP{pSeqNG(i),stL});
    np = sum(NP{pSeqNG(i),stL});
    semLp = sqrt(((PS{pSeqNG(i),stL}- gmLp).*(PS{pSeqNG(i),stL}- gmLp))'*...
        NP{pSeqNG(i),stL}/sum(NP{pSeqNG(i),stL}))/sqrt(length(NP{pSeqNG(i),stL}));
    gmLm = PS{mSeqNG(i),stL}'*NP{mSeqNG(i),stL}/sum(NP{mSeqNG(i),stL});
    nm = sum(NP{mSeqNG(i),stL});
    semLm = sqrt(((PS{mSeqNG(i),stL}- gmLm).*(PS{mSeqNG(i),2}- gmLm))'*...
        NP{mSeqNG(i),stL}/sum(NP{mSeqNG(i),stL}))/sqrt(length(NP{mSeqNG(i),stL}));
    gmL = (gmLp*np+gmLm*nm)/(np+nm);
    semL = (semLp*np+semLm*nm)/(np+nm);
    errorbar(i, gmL, semL, 'o', 'color', [1 0 0], 'linewidth', 3)
end
ylabel('Probability turn stimulus side')
% axis([0.5 3.5 0 1])
axis([0.5 3.5 0 250])
%%
% subplot(1,2,2)
figure,
hold on
plot([0 4], [0.5 0.5], '--g', 'linewidth', 1.5)
for i = 1 : length(pSeqNG)
% + NG
    gmR = PSS{pSeqNG(i),stR}'*NP{pSeqNG(i),stR}/sum(NP{pSeqNG(i),stR});
    semR = sqrt(((PSS{pSeqNG(i),stR}- gmR).*(PSS{pSeqNG(i),stR}- gmR))'*NP{pSeqNG(i),stR}/sum(NP{pSeqNG(i),stR}))/sqrt(length(NP{pSeqNG(i),stR}));
    errorbar(i, gmR, semR, 'o', 'color', [0.5 0 0], 'linewidth', 3)
    gmL = PSS{pSeqNG(i),stL}'*NP{pSeqNG(i),stL}/sum(NP{pSeqNG(i),stL});
    semL = sqrt(((PSS{pSeqNG(i),stL}- gmL).*(PSS{pSeqNG(i),stL}- gmL))'*NP{pSeqNG(i),stL}/sum(NP{pSeqNG(i),stL}))/sqrt(length(NP{pSeqNG(i),stL}));
    errorbar(i, gmL, semL, 'o', 'color', [1 0 0], 'linewidth', 3)
% - NG   
    gmR = PSS{mSeqNG(i),stR}'*NP{mSeqNG(i),stR}/sum(NP{mSeqNG(i),stR});
    semR = sqrt(((PSS{mSeqNG(i),1}- gmR).*(PSS{mSeqNG(i),stR}- gmR))'*NP{mSeqNG(i),stR}/sum(NP{mSeqNG(i),stR}))/sqrt(length(NP{mSeqNG(i),stR}));
    errorbar(i, gmR, semR, 'o', 'color', [0 0 0.5], 'linewidth', 3)
    gmL = PSS{mSeqNG(i),stL}'*NP{mSeqNG(i),stL}/sum(NP{mSeqNG(i),stL});
    semL = sqrt(((PSS{mSeqNG(i),stL}- gmL).*(PSS{mSeqNG(i),stL}- gmL))'*NP{mSeqNG(i),stL}/sum(NP{mSeqNG(i),stL}))/sqrt(length(NP{mSeqNG(i),stL}));
    errorbar(i, gmL, semL, 'o', 'color', [0 0 1], 'linewidth', 3)
% + RG
    gmR = PSS{pSeqRG(i),stR}'*NP{pSeqRG(i),stR}/sum(NP{pSeqRG(i),stR});
    semR = sqrt(((PSS{pSeqRG(i),stR}- gmR).*(PSS{pSeqRG(i),stR}- gmR))'*NP{pSeqRG(i),stR}/sum(NP{pSeqRG(i),stR}))/sqrt(length(NP{pSeqRG(i),stR}));
    errorbar(i, gmR, semR, 'o', 'color', [0.5 0 0.5], 'linewidth', 3)
    gmL = PSS{pSeqRG(i),stL}'*NP{pSeqRG(i),stL}/sum(NP{pSeqRG(i),stL});
    semL = sqrt(((PSS{pSeqRG(i),stL}- gmL).*(PSS{pSeqRG(i),stL}- gmL))'*NP{pSeqRG(i),stL}/sum(NP{pSeqRG(i),stL}))/sqrt(length(NP{pSeqRG(i),stL}));
    errorbar(i, gmL, semL, 'o', 'color', [1 0 1], 'linewidth', 3)   
% - RG
    gmR = PSS{mSeqRG(i),stR}'*NP{mSeqRG(i),stR}/sum(NP{mSeqRG(i),stR});
    semR = sqrt(((PSS{mSeqRG(i),stR}- gmR).*(PSS{mSeqRG(i),stR}- gmR))'*NP{mSeqRG(i),stR}/sum(NP{mSeqRG(i),stR}))/sqrt(length(NP{mSeqRG(i),stR}));
    errorbar(i, gmR, semR, 'o', 'color', [0 0.5 0.5], 'linewidth', 3)
    gmL = PSS{mSeqRG(i),stL}'*NP{mSeqRG(i),stL}/sum(NP{mSeqRG(i),stL});
    semL = sqrt(((PSS{mSeqRG(i),stL}- gmL).*(PSS{mSeqRG(i),stL}- gmL))'*NP{mSeqRG(i),stL}/sum(NP{mSeqRG(i),stL}))/sqrt(length(NP{mSeqRG(i),stL}));
    errorbar(i, gmL, semL, 'o', 'color', [0 1 1], 'linewidth', 3)     
end
ylabel('Probability turn same side')
% axis([0.5 3.5 0.3 0.9])


%%

figure,
% pSeqNG = [5 9 1];
% pSeqRG = [7 11 3];
% mSeqNG = [6 10 2];
% mSeqRG = [8 12 4];
% subplot(1,2,1)
hold on
plot([0 4], [0.5 0.5], '--g', 'linewidth', 1.5)
for i = 1 : length(pSeqNG)
% + NG 
    gmR = PSS1{pSeqNG(i),stR}'*NP{pSeqNG(i),stR}/sum(NP{pSeqNG(i),stR});
    semR = sqrt(((PSS1{pSeqNG(i),stR}- gmR).*(PSS1{pSeqNG(i),stR}- gmR))'*NP{pSeqNG(i),stR}/sum(NP{pSeqNG(i),stR}))/sqrt(length(NP{pSeqNG(i),stR}));
    errorbar(i, gmR, semR, 'o', 'color', [0.5 0 0], 'linewidth', 3)
    gmL = PSS1{pSeqNG(i),stL}'*NP{pSeqNG(i),stL}/sum(NP{pSeqNG(i),stL});
    semL = sqrt(((PSS1{pSeqNG(i),stL}- gmL).*(PSS1{pSeqNG(i),stL}- gmL))'*NP{pSeqNG(i),stL}/sum(NP{pSeqNG(i),stL}))/sqrt(length(NP{pSeqNG(i),stL}));
    errorbar(i, gmL, semL, 'o', 'color', [1 0 0], 'linewidth', 3)
% - NG   
    gmR = PSS1{mSeqNG(i),stR}'*NP{mSeqNG(i),stR}/sum(NP{mSeqNG(i),stR});
    semR = sqrt(((PSS1{mSeqNG(i),stR}- gmR).*(PSS1{mSeqNG(i),stR}- gmR))'*NP{mSeqNG(i),stR}/sum(NP{mSeqNG(i),stR}))/sqrt(length(NP{mSeqNG(i),stR}));
    errorbar(i+0.2, gmR, semR, 'o', 'color', [0 0 0.5], 'linewidth', 3)
    gmL = PSS1{mSeqNG(i),stL}'*NP{mSeqNG(i),stL}/sum(NP{mSeqNG(i),stL});
    semL = sqrt(((PSS1{mSeqNG(i),stL}- gmL).*(PSS1{mSeqNG(i),stL}- gmL))'*NP{mSeqNG(i),stL}/sum(NP{mSeqNG(i),stL}))/sqrt(length(NP{mSeqNG(i),stL}));
    errorbar(i+0.2, gmL, semL, 'o', 'color', [0 0 1], 'linewidth', 3)
% + RG
    gmR = PS{pSeqRG(i),stR}'*NP{pSeqRG(i),stR}/sum(NP{pSeqRG(i),stR});
    semR = sqrt(((PS{pSeqRG(i),stR}- gmR).*(PS{pSeqRG(i),stR}- gmR))'*NP{pSeqRG(i),stR}/sum(NP{pSeqRG(i),stR}))/sqrt(length(NP{pSeqRG(i),stR}));
%     errorbar(i, gmR, semR, 'o', 'color', [0.5 0 0.5], 'linewidth', 3)
    gmL = PS{pSeqRG(i),stL}'*NP{pSeqRG(i),stL}/sum(NP{pSeqRG(i),stL});
    semL = sqrt(((PS{pSeqRG(i),stL}- gmL).*(PS{pSeqRG(i),stL}- gmL))'*NP{pSeqRG(i),stL}/sum(NP{pSeqRG(i),stL}))/sqrt(length(NP{pSeqRG(i),stL}));
%     errorbar(i, gmL, semL, 'o', 'color', [1 0 1], 'linewidth', 3)   
% - RG
    gmR = PS{mSeqRG(i),stR}'*NP{mSeqRG(i),stR}/sum(NP{mSeqRG(i),stR});
    semR = sqrt(((PS{mSeqRG(i),stR}- gmR).*(PS{mSeqRG(i),stR}- gmR))'*NP{mSeqRG(i),stR}/sum(NP{mSeqRG(i),stR}))/sqrt(length(NP{mSeqRG(i),stR}));
%     errorbar(i+0.2, gmR, semR, 'o', 'color', [0 0.5 0.5], 'linewidth', 3)
    gmL = PS{mSeqRG(i),stL}'*NP{mSeqRG(i),stL}/sum(NP{mSeqRG(i),stL});
    semL = sqrt(((PS{mSeqRG(i),stL}- gmL).*(PS{mSeqRG(i),stL}- gmL))'*NP{mSeqRG(i),stL}/sum(NP{mSeqRG(i),stL}))/sqrt(length(NP{mSeqRG(i),stL}));
%     errorbar(i+0.2, gmL, semL, 'o', 'color', [0 1 1], 'linewidth', 3)     
end
ylabel('Probability turn stimulus side')

% axis([0.5 3.5 0.3 0.9])


%%

figure,
hold on
plot([0 4], [0.25 0.25], '--g', 'linewidth', 1.5)
for i = 1 : length(pSeqNG)
% R NG 
    gmR = PSS1{pSeqNG(i),stR}'*NP{pSeqNG(i),stR}/sum(NP{pSeqNG(i),stR});
    semR = sqrt(((PSS1{pSeqNG(i),stR}- gmR).*(PSS1{pSeqNG(i),stR}- gmR))'*NP{pSeqNG(i),stR}/sum(NP{pSeqNG(i),stR}))/sqrt(length(NP{pSeqNG(i),stR}));
    errorbar(i, gmR, semR, 'o', 'color', [0.5 0 0], 'linewidth', 3)
    gmL = PSS1{pSeqNG(i),stL}'*NP{pSeqNG(i),stL}/sum(NP{pSeqNG(i),stL});
    semL = sqrt(((PSS1{pSeqNG(i),stL}- gmL).*(PSS1{pSeqNG(i),stL}- gmL))'*NP{pSeqNG(i),stL}/sum(NP{pSeqNG(i),stL}))/sqrt(length(NP{pSeqNG(i),stL}));
    errorbar(i, gmL, semL, 'o', 'color', [1 0 0], 'linewidth', 3)
% L NG   
    gmR = PSS2{pSeqNG(i),stR}'*NP{pSeqNG(i),stR}/sum(NP{pSeqNG(i),stR});
    semR = sqrt(((PSS2{pSeqNG(i),stR}- gmR).*(PSS2{pSeqNG(i),stR}- gmR))'*NP{pSeqNG(i),stR}/sum(NP{mSeqNG(i),stR}))/sqrt(length(NP{mSeqNG(i),stR}));
    errorbar(i+0.2, gmR, semR, 'o', 'color', [0 0 0.5], 'linewidth', 3)
    gmL = PSS2{pSeqNG(i),stL}'*NP{pSeqNG(i),stL}/sum(NP{pSeqNG(i),stL});
    semL = sqrt(((PSS2{pSeqNG(i),stL}- gmL).*(PSS2{pSeqNG(i),stL}- gmL))'*NP{pSeqNG(i),stL}/sum(NP{mSeqNG(i),stL}))/sqrt(length(NP{mSeqNG(i),stL}));
    errorbar(i+0.2, gmL, semL, 'o', 'color', [0 0 1], 'linewidth', 3)
% R RG
    gmR = PSS1{pSeqRG(i),stR}'*NP{pSeqRG(i),stR}/sum(NP{pSeqRG(i),stR});
    semR = sqrt(((PSS1{pSeqRG(i),stR}- gmR).*(PSS1{pSeqRG(i),stR}- gmR))'*NP{pSeqRG(i),stR}/sum(NP{pSeqRG(i),stR}))/sqrt(length(NP{pSeqRG(i),stR}));
    errorbar(i, gmR, semR, 'o', 'color', [0.5 0 0.5], 'linewidth', 3)
    gmL = PSS1{pSeqRG(i),stL}'*NP{pSeqRG(i),stL}/sum(NP{pSeqRG(i),stL});
    semL = sqrt(((PSS1{pSeqRG(i),stL}- gmL).*(PSS1{pSeqRG(i),stL}- gmL))'*NP{pSeqRG(i),stL}/sum(NP{pSeqRG(i),stL}))/sqrt(length(NP{pSeqRG(i),stL}));
    errorbar(i, gmL, semL, 'o', 'color', [1 0 1], 'linewidth', 3)   
% L RG
    gmR = PSS2{pSeqRG(i),stR}'*NP{pSeqRG(i),stR}/sum(NP{pSeqRG(i),stR});
    semR = sqrt(((PSS2{pSeqRG(i),stR}- gmR).*(PSS2{pSeqRG(i),stR}- gmR))'*NP{pSeqRG(i),stR}/sum(NP{mSeqRG(i),stR}))/sqrt(length(NP{mSeqRG(i),stR}));
    errorbar(i+0.2, gmR, semR, 'o', 'color', [0 0.5 0.5], 'linewidth', 3)
    gmL = PSS2{pSeqRG(i),stL}'*NP{pSeqRG(i),stL}/sum(NP{pSeqRG(i),stL});
    semL = sqrt(((PSS2{pSeqRG(i),stL}- gmL).*(PSS2{pSeqRG(i),stL}- gmL))'*NP{pSeqRG(i),stL}/sum(NP{mSeqRG(i),stL}))/sqrt(length(NP{mSeqRG(i),stL}));
    errorbar(i+0.2, gmL, semL, 'o', 'color', [0 1 1], 'linewidth', 3)     
end
ylabel('Probability turn stimulus side')

% axis([0.5 3.5 0 0.5])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
path1 = 'C:\Users\Tomas\Dropbox (Sensorimotor)\ChiappeLabNew\data\TOMÁS\Free Walking VR\Data\Unilateral Balanced\WTTB\';
path2 = 'C:\Users\tomas\Dropbox (Sensorimotor)\ChiappeLabNew\data\TOMÁS\Free Walking VR\Data\Unilateral Stimulus\WT TB\';
[params] = GetParams();
[STRALL1, DSTALL1, NSTRALL1, PSSALL1, NPSSALL1, pTypes1] = GetStrAndVisInf(path1, params);
[STRALL2, DSTALL2, NSTRALL2, PSSALL2, NPSSALL2, pTypes2] = GetStrAndVisInf(path2, params);
%%
params.lthr = 350;
params.pTypes = pTypes1;
[MSTR1, SEMSTR1, MF1, NF1] = GetGMSEM(STRALL1, DSTALL1, NSTRALL1, params);
params.pTypes = pTypes2;
[MSTR2, SEMSTR2, MF2, NF2] = GetGMSEM(STRALL2, DSTALL2, NSTRALL2, params);
%%
STRV = cell(3,1);
NSTRV = cell(3,1);
vectOr = [5 6 9 10 1 2];
for i = 1 : 3
    for j = 1 : size(MF1,2)
        STRV{i} = vertcat(STRV{i},MF1(vectOr(i),j));
        STRV{i} = vertcat(STRV{i},MF1(vectOr(i+1),j));
        NSTRV{i} = vertcat(NSTRV{i},NF1(vectOr(i),j));
        NSTRV{i} = vertcat(NSTRV{i},NF1(vectOr(i+1),j));
    end
end
vectOr = [5 7 1];
for i = 1 : 3
    for j = 1 : size(MF2,2)
        STRV{i} = vertcat(STRV{i},MF2(vectOr(i),j));
        NSTRV{i} = vertcat(NSTRV{i},NF2(vectOr(i),j));
    end
end
%%
MSTR = zeros(3,1);
SEMSTR = zeros(3,1);
thr = 600;
figure,
% for n = 1 : 16
%     thr = n *100;
%     subplot(4,4,n)
    for i = 1 : 3
        auxSt = STRV{i};
        auxNSt = NSTRV{i};
        auxNSt(isnan(auxSt)) = [];
        auxSt(isnan(auxSt)) = [];
        auxSt(auxNSt<thr) = [];
        auxNSt(auxNSt<thr) = [];
        MSTR(i) = auxSt'*auxNSt/sum(auxNSt);
        SEMSTR(i) = sqrt(((auxSt-MSTR(i)).*(auxSt-MSTR(i)))'*auxNSt/sum(auxNSt))/sqrt(length(auxSt));
    end
    
    hold on
    errorbar([1 5 10], [36.5 47.5 52.5], [2 2 2], 'or')
    errorbar([1 5 10],MSTR,SEMSTR, 'ok')
    axis([0 11 20 55])
%     title(num2str(thr))
% end
%%
cmp = hot(10);
% vecSt = [5 6 9 10 1 2];
vecSt = [6 7 1];
figure,
hold on
for n = 1 : length(vecSt)
    errorbar(n, MSTR1(vecSt(n)), SEMSTR1(vecSt(n)), 'color', cmp(n+2,:), 'linewidth', 3, 'marker', 'o', 'markersize', 10, 'markerfacecolor',cmp(n+2,:));
end
axis([0 7 30 55])
xlabel('Dot Size')
ylabel('Straightness (a.u.)')