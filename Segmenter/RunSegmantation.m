%% Load Data Just One Fly
% clear
clc
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Dark\Dark WTTB\Fly 1659\DataLowRes.mat';
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Windmill\NG-RG\Fly 1666\DataLowRes.mat';
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\VTR39\VTFLPKir\Fly 2623\DataLowRes.mat';
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\VTR39\VTR39FLPKir\Fly 2650\DataLowRes.mat';
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Unilateral Stimulus\R39VT05\R39VT05Cnt\Fly 2284\DataLowRes.mat';
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Unilateral Stimulus\R39VT05\R39VT05Kir\Fly 2303\DataLowRes.mat';
path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\WT TB\Fly 1906\DataLowRes.mat';

dt = load(path);
dt = dt.Flies;
nSeqs = length(dt.Data);
partsCtr = find(strcmp([dt.Seq], '10NG'))';
Vf = [];
Vs = [];
Vr = [];
acSt = [];
flp = [];
xx = [];
yy = [];
for i = partsCtr
    xx = vertcat(xx, dt.Data{i}.X);
    yy = vertcat(yy, dt.Data{i}.Y);
    flp = vertcat(flp, dt.Data{i}.flp);
    Vf = vertcat(Vf, dt.Data{i}.Vf);
    Vs = vertcat(Vs, dt.Data{i}.Vs);
    Vr = vertcat(Vr, dt.Data{i}.Vr);
    acSt = vertcat(acSt,dt.Data{i}.actState);
end
clear i nSeqs path
[params] = GetParams();
[locs, ~, cmhSong, thr, dx] = SortSpikes(Vr, acSt, params);
[locsF, pksIF, pksEF, pV] = CullSpikesBasedOnVf(Vr, locs, cmhSong, thr, Vf, params);
[locsF2, pksIF2, pksEF2, scores] = TemplateCompSpike(Vr,locsF,pksIF, pksEF, cmhSong, thr, params);
[FWBoutStr, FWBout] = GetFbouts(acSt, Vf, locsF2, pksIF2,pksEF2, params);
%
vr1 = Vr;
vr1(acSt==1) = nan;
vf1 = Vf;
vf1(acSt==1) = nan;
FVect = ones(length(Vr),1);
for i = 1 : length(FWBoutStr)
    if (length(FWBoutStr{i}) > 10)
        FVect(FWBoutStr{i}) = 0;
    end
end
vr3 = Vr;
vr3(FVect==1) = nan;
vf3 = Vf;
vf3(FVect==1) = nan;

vect = ones(length(Vr),1);
for i = 1 : length(locsF)
    vect((locsF(i)-pksIF(i)):(locsF(i)+pksEF(i))) = 0;
end
vr2 = Vr;
vr2(vect==1) = nan;
vf2 = Vf;
vf2(vect==1) = nan;

FVect = ones(length(Vr),1);
for i = 1 : length(FWBout)
    FVect(FWBout{i}) = 0;
end
vr4 = Vr;
vr4(FVect==1) = nan;
vf4 = Vf;
vf4(FVect==1) = nan;



%
VISpk = [];
VISpkH = [];
VISpkVf = [];
tw = 2;
w = 6;
for i = 1 : length(locsF)
    if (locsF(i)-pksIF(i)) > w+3
        auxVr = Vr((locsF(i)-pksIF(i) - w - tw):(locsF(i)-pksIF(i) - w + tw));
        auxVf = Vf((locsF(i)-pksIF(i) - w - tw):(locsF(i)-pksIF(i) - w + tw));
        
        if abs(mean(auxVf))>6
            if sign(mean(auxVr)) == sign(Vr(locsF(i)))
                VISpk = vertcat(VISpk, 1);
            else
                VISpk = vertcat(VISpk, 0);
            end
        end
    end
end
% disp(num2str(sum(VISpk)/length(VISpk)))

% disp(['% Spikes: ' num2str(100*sum(abs(vect-1))/sum(abs(acSt)))])
% disp(['% ForwBouts: ' num2str(100*sum(abs(FVect-1))/sum(abs(acSt)))])

%
% vr2 = Vr;
% vr2(cmhSong<thr) = nan;
% vf2 = Vf;
% vf2(cmhSong<thr) = nan;
%
% spksVect = zeros(length(Vr),1);
% for i = 1 : length(locsF)
%     spksVect((locsF(i)-pksIF(i)):(locsF(i)+pksEF(i))) = 1;
% end
% vr3 = Vr;
% vr3(spksVect==1) = nan;
% vf3 = Vf;
% vf3(spksVect==1) = nan;

figure,
plot([0 length(Vr)], [0 0], '--k')
hold on
plot((1:length(Vr))/60, Vr, 'b', 'linewidth', 1.5)
hold on
plot((1:length(Vr))/60,vr1,'k','linewidth', 1.5)
hold on
plot((1:length(Vr))/60,vr2,'r', 'linewidth', 1.5)
hold on
plot((1:length(Vr))/60,vr4,'color', [0 0.5 0], 'linewidth', 1.5)
hold on
plot((1:length(Vr))/60,vr3,'color', [0 1 0], 'linewidth', 1.5)
hold on

plot([0 length(Vr)], -1200+[0 0], '--k')
hold on
plot((1:length(Vf))/60, -1200+10*Vf, 'b', 'linewidth', 1.5)
hold on
plot((1:length(Vr))/60, -1200+10*vf1,'k','linewidth', 1.5)
hold on
plot((1:length(Vf))/60, -1200+10*vf2, 'r', 'linewidth', 1.5)
hold on
plot((1:length(Vr))/60, -1200+10*vf4,'color', [0 0.5 0], 'linewidth', 1.5)
hold on
plot((1:length(Vr))/60, -1200+10*vf3,'color', [0 1 0], 'linewidth', 1.5)
hold on
plot((1:length(cmhSong))/60,1300-cmhSong, 'color', [0.2 0.2 0], 'linewidth', 1.5)
plot((1:length(Vr))/60, -1300+0.5*flp)
% scatter(locsF/60,Vr(locsF),100, 'g', 'filled');
% scatter(dx.locsF/60,1300-cmhSong(dx.locsF),100, 'm', 'linewidth', 2);

hold on
plot([0 length(Vr)],-[thr thr]+1300, 'color', [0 0.7 0])
hold on
locsF2 = locsF2(locsF2>47*60);
locsF2 = locsF2(locsF2<67*60);
cmap = jet(length(locsF2));
% scatter(locsF/60,Vr(locsF),100, 'g', 'filled');
% scatter(dx.locsV/60,Vr(dx.locsV),100, 'g', 'linewidth', 2);
% scatter(dx.locsF/60,Vr(dx.locsF),100, 'm', 'linewidth', 2);
% scatter(locs/60,Vr(locs),100, 'k', 'filled');
% scatter(locs/60,Vr(locs),100, 'c', 'filled');
scatter(locsF2/60,Vr(locsF2),100, cmap, 'filled');
% cmap = jet(101);
% scores(scores<0.2) = 0.2;
% for i = 1 : length(scores)
%     sc = (scores(i)-min(scores))/(max(scores)-min(scores));
%     scatter(locsF(i)/60,Vr(locsF(i)),100, cmap(floor(100*sc+1),:), 'filled');
% end
axis([47 67 -1400 1400])
% axis([16.7 20.5 -1400 1400])
%%

ti = 47*60;
tf = 67*60;

figure,
hold on
plot(xx(ti:tf), yy(ti:tf), 'k', 'linewidth', 1.5)
scatter(xx(locsF2),yy(locsF2),100, cmap, 'filled');
% scatter(xx(ti),yy(ti),100,'r','filled')
% for r = 1 : 2
%     ti2 = ti + (r-1)*60;
%     tf2 = ti + r*60;
%     plot(xx(ti2:tf2), yy(ti2:tf2), 'color', cmap(r,:), 'linewidth', 1.5)
% end
axis([-40 40 -40 40])
axis square
xlabel('X (mm)')
ylabel('Y (mm)')


%%
BtsL  = [];
BtsN  = [];


path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\WT TB\';
flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowRes.mat'], 'Flies');
seq = dt.Flies.Seq;
pTypes = unique(seq);

for n = 1 : length(flies)
    disp(num2str(n))
    dt = load([path flies(n).name '\DataLowRes.mat'], 'Flies');
    seq = dt.Flies.Seq;
    dt = dt.Flies;
    for k = 1 : length(dt.Data)
        for j = 1 : length(pTypes)
            Vf = [];
            Vs = [];
            Vr = [];
            acSt = [];
            flp = [];
            switch seq{k}
                case pTypes{j}
                    flp = vertcat(flp, dt.Data{k}.flp);
                    Vf = vertcat(Vf, dt.Data{k}.Vf);
                    Vs = vertcat(Vs, dt.Data{k}.Vs);
                    Vr = vertcat(Vr, dt.Data{k}.Vr);
                    acSt = vertcat(acSt,dt.Data{k}.actState);
                    BtsN = vertcat(BtsN, length(dt.Data{k}.Bouts));
                    for i = 1 : length(dt.Data{k}.Bouts)
                        BtsL = vertcat(BtsL, length(dt.Data{k}.Bouts{i}));
                    end
            end
%             [params] = GetParams();
%             [locs, ~, cmhSong, thr, dx] = SortSpikes(Vr, acSt, params);
%             [locsF, pksIF, pksEF, pV] = CullSpikesBasedOnVf(Vr, locs, cmhSong, thr, Vf, params);
%             [locsF2, pksIF2, pksEF2, scores] = TemplateCompSpike(Vr,locsF,pksIF, pksEF, cmhSong, thr, params);
%             [FWBoutStr, FWBout] = GetFbouts(acSt, Vf, locsF2, pksIF2,pksEF2, params);
            
        end
        
    end
    
end


path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Dark\Dark WTTB\';
flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowRes.mat'], 'Flies');
seq = dt.Flies.Seq;
pTypes = unique(seq);

for n = 1 : length(flies)
    disp(num2str(n))
    dt = load([path flies(n).name '\DataLowRes.mat'], 'Flies');
    seq = dt.Flies.Seq;
    dt = dt.Flies;
    for k = 1 : length(dt.Data)
        for j = 1 : length(pTypes)
            Vf = [];
            Vs = [];
            Vr = [];
            acSt = [];
            flp = [];
            switch seq{k}
                case pTypes{j}
                    flp = vertcat(flp, dt.Data{k}.flp);
                    Vf = vertcat(Vf, dt.Data{k}.Vf);
                    Vs = vertcat(Vs, dt.Data{k}.Vs);
                    Vr = vertcat(Vr, dt.Data{k}.Vr);
                    acSt = vertcat(acSt,dt.Data{k}.actState);
                    BtsN = vertcat(BtsN, length(dt.Data{k}.Bouts));
                    for i = 1 : length(dt.Data{k}.Bouts)
                        BtsL = vertcat(BtsL, length(dt.Data{k}.Bouts{i}));
                    end
            end
%             [params] = GetParams();
%             [locs, ~, cmhSong, thr, dx] = SortSpikes(Vr, acSt, params);
%             [locsF, pksIF, pksEF, pV] = CullSpikesBasedOnVf(Vr, locs, cmhSong, thr, Vf, params);
%             [locsF2, pksIF2, pksEF2, scores] = TemplateCompSpike(Vr,locsF,pksIF, pksEF, cmhSong, thr, params);
%             [FWBoutStr, FWBout] = GetFbouts(acSt, Vf, locsF2, pksIF2,pksEF2, params);
            
        end
        
    end
    
end

%%
figure,

lcents = 0:5:10000;
distBoutSize = hist(BtsL,lcents);
plot(lcents, distBoutSize)
save('DistBoutLength.mat', 'lcents', 'distBoutSize') 


















%%
figure,
hcents = -5:2:400;
hs = hist(cmhSong(acSt==0),hcents);
hm = hist(cmhSong,hcents);
plot(hcents, hs/sum(hs), 'k', 'linewidth', 2)
hold on
plot(hcents, hm/sum(hm), 'b', 'linewidth', 2)
hold on
plot([thr thr], [0 0.15], '--g')
axis([0 200 0 0.15])



%% Load Data Just One Fly
clear
clc
% path = 'C:\Users\Tomas\Dropbox (Sensorimotor)\ChiappeLabNew\data\TOMÁS\Free Walking VR\Data\Dark\Dark WTTB\Fly 1641\DataLowRes.mat';
path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Windmill\NG-RG\Fly 1664\DataLowRes.mat';
dt = load(path);
dt = dt.Flies;
nSeqs = length(dt.Data);
Vf = [];
Vs = [];
Vr = [];
acSt = [];
for i = 1 %: nSeqs
    for j = 1 : length(dt.Data{i}.Bouts)
        Vf = vertcat(Vf, dt.Data{i}.Vf(dt.Data{i}.Bouts{j}));
        Vs = vertcat(Vs, dt.Data{i}.Vs(dt.Data{i}.Bouts{j}));
        Vr = vertcat(Vr, dt.Data{i}.Vr(dt.Data{i}.Bouts{j}));
        acSt = vertcat(acSt,dt.Data{i}.actState(dt.Data{i}.Bouts{j}));
    end
end
clear i nSeqs path
%%
[params] = GetParams();
[locs, ~, cmhSong, thr, dx] = SortSpikes(Vr, acSt, params);
[locsF, pksIF, pksEF, ~] = CullSpikesBasedOnVf(Vr, locs, cmhSong, thr, Vf, params);
[locsF, pksIF, pksEF, scores] = TemplateCompSpike(Vr,locsF,pksIF, pksEF, cmhSong, thr, params);
[FWBoutStr, FWBout] = GetFbouts(acSt, Vf, locsF, pksIF,pksEF, params);
% [locs, p, cmhSong, thr] = SortSpikes(Vr, acSt, params);
% [locsF, pksIF,pksEF, pV] = CullSpikesBasedOnVf(Vr, locs, cmhSong, thr, Vf, params);
% [FWBout] = GetFbouts(acSt, Vf, locsF, pksIF,pksEF, params);
% [locsF, pksIF, pksEF, score] = TemplateCompSpike(Vr,locsF,pksIF, pksEF,params);
% pksIF = locsF(score>params.cutoff);
% pksEF = locsF(score>params.cutoff);
% locsF = locsF(score>params.cutoff);

%%
vffb = [];
vrrb = [];
cmhb = [];
for j = 1 : length(FWBout)
    vffb = vertcat(vffb, Vf(FWBout{j}));
    vrrb = vertcat(vrrb, Vr(FWBout{j}));
    cmhb = vertcat(cmhb, cmhSong(FWBout{j})');
end

[mapSumT, ~, ~] = SpeedMapOccupancy(Vr, Vf, ones(size(Vr)));
[mapSumF, ~, ~] = SpeedMapOccupancy(vrrb, vffb, ones(size(Vr)));

%%
figure,
cents = 1:1:200;
plot(cents, hist(cmhb,cents)/length(cmhb), 'r')
hold on
plot(cents, hist(cmhSong,cents)/length(cmhSong), 'k')
figure,
mapToPlot=(mapSumT/sum(sum(mapSumT))-mapSumF/sum(sum(mapSumF)));
mapToPlot(isnan(mapToPlot)) = 0;
imagesc(-1:10:1000, -10:1:40, mapToPlot)
axis([0 1000 -10 40])
caxis(0.003*[-1 1])
colormap redblue
set(gca,'Ydir','normal')

%%

% pks = findpeaks(abs(Vr));
%
% figure,
% subplot(1,2,1)
% hold on
% scatter(pV.vfstd, pV.peakVal, 100, 'k')
% scatter(pV1.vfstd, pV1.peakVal, 100, 'r')
% subplot(1,2,2)
% plot(0:20:900, log(hist(pks, 0:20:900)), 'k')
% hold on
% plot(0:20:900, log(hist(pV1.peakVal, 0:20:900)), 'r')

figure,
imagesc(-1200:20:1200, -10:1:40, mapSum)
caxis([1 70])
colormap jet
set(gca,'Ydir','normal')
%%






%% Load Data All Flies
clear
clc
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Dark\Dark WTTB\';
path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\WT TB\';
flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowRes.mat'], 'Flies');
seq = dt.Flies.Seq;
pTypes = unique(seq);

VrAll = cell(length(flies),1);
VfAll = cell(length(flies),1);
CmAll = cell(length(flies),1);
VrForw = cell(length(flies),1);
VfForw = cell(length(flies),1);
CmForw = cell(length(flies),1);
VrSpikes = cell(length(flies),1);
VfSpikes = cell(length(flies),1);
CmSpikes = cell(length(flies),1);
aaa = [];
bbb = [];
pSpks = cell(length(flies),length(pTypes));
pFwd = cell(length(flies),length(pTypes));
pT = cell(length(flies),length(pTypes));

VISP = cell(length(flies),length(pTypes));
NSP = cell(length(flies),length(pTypes));
for n = 1 : length(flies)
    disp(num2str(n))
    dt = load([path flies(n).name '\DataLowRes.mat'], 'Flies');
    seq = dt.Flies.Seq;
    dt = dt.Flies.Data;
    for k = 1 : length(dt)
        for j = 1 : length(pTypes)
                        vra = [];
                        vfa = [];
                        cma = [];
                        vrf = [];
                        vff = [];
                        cmf = [];
                        vrs = [];
                        vfs = [];
                        cms = [];
            ps = [];
            pf = [];
            pt = [];
            vispk = [];
            switch seq{k}
                case pTypes{j}
                    Vf = dt{k}.Vf;
                    Vs = dt{k}.Vs;
                    Vr = dt{k}.Vr;
                    acSt = dt{k}.actState;
                    
                    [params] = GetParams();
                    
                    [locs, ~, cmhSong, thr, dx] = SortSpikes(Vr, acSt, params);
                    [locsF, pksIF, pksEF, ~] = CullSpikesBasedOnVf(Vr, locs, cmhSong, thr, Vf, params);
                    [locsF, pksIF, pksEF, scores] = TemplateCompSpike(Vr,locsF,pksIF, pksEF, cmhSong, thr, params);
                    [FWBoutStr, FWBout] = GetFbouts(acSt, Vf, locsF, pksIF,pksEF, params);
                    
                    FVect = zeros(length(Vr),1);
                    for i = 1 : length(FWBout)
                        FVect(FWBout{i}) = 1;
                    end
                    
                    
                                        vect = zeros(length(Vr),1);
                                        tw = 2;
                                        w = 6;
                                        for i = 1 : length(locsF)
                                            if (locsF(i)-pksIF(i)) > w+3
                                                auxVr = Vr((locsF(i)-pksIF(i) - w - tw):(locsF(i)-pksIF(i) - w + tw));
                                                auxVf = Vf((locsF(i)-pksIF(i) - w - tw):(locsF(i)-pksIF(i) - w + tw));
                                                if abs(auxVr) > 200%mean(auxVf)>6
                                                    if sign(mean(auxVr)) == sign(Vr(locsF(i)))
                                                        vispk = vertcat(vispk, 1);
                                                    else
                                                        vispk = vertcat(vispk, 0);
                                                    end
                                                end
                                            end
                                            vect((locsF(i)-pksIF(i)):(locsF(i)+pksEF(i))) = 1;
                                        end
                    
                    npss = [];
                    pss = [];
                    for i = 1 : length(FWBout)
%                         if length(FWBout{i}) > params.thr/2
                            vrb = Vr(FWBout{i});
                            vfb = Vf(FWBout{i});
                            [BoutVr,~,BoutP,~] = GetProbVec(vrb,vfb, params.delta,3,150,inf);
                            npss = vertcat(npss,length(BoutP));
                            pss = vertcat(pss,sum(BoutP)/(length(BoutP)+1));
                            %                             pss = vertcat(pss,mean(abs(BoutVr)));
%                         end
                    end
                    
                    
                    
                    ps = vertcat(ps, sum(abs(vect))/sum(abs(acSt)));
                    pf = vertcat(pf, sum(abs(FVect))/sum(abs(acSt)));
                    pt = vertcat(pt, sum(abs(acSt)));
%                     ps = vertcat(ps, pss'*npss/sum(npss));
%                     pf = vertcat(pf, sum(abs(FVect))/sum(abs(acSt)));
%                     pt = vertcat(pt, sum(npss));%
                    
                    for btf = 1 : length(FWBout)
                        vrf = vertcat(vrf, Vr(FWBout{btf}));
                        vff = vertcat(vff, Vf(FWBout{btf}));
                        cmf = vertcat(cmf, cmhSong(FWBout{btf})');
                    end
                    for bt = 1 : length(dt{k}.Bouts)
                        vra = vertcat(vra, Vr(dt{k}.Bouts{bt}));
                        vfa = vertcat(vfa, Vf(dt{k}.Bouts{bt}));
                        cma = vertcat(cma, cmhSong(dt{k}.Bouts{bt})');
                    end
                    for st = 1 : length(locsF)
                        vrs = vertcat(vrs, Vr((locsF(st)-pksIF(st)):(locsF(st)-pksEF(st))));
                        vfs = vertcat(vfs, Vf((locsF(st)-pksIF(st)):(locsF(st)-pksEF(st))));
                        cms = vertcat(cms, cmhSong((locsF(st)-pksIF(st)):(locsF(st)-pksEF(st)))');
                    end
            end
            
            %             pSpks{n,j} = vertcat(pSpks{n,j}, ps);
            %             pFwd{n,j} = vertcat(pFwd{n,j}, pf);
            %             pT{n,j} = vertcat(pT{n,j}, pt);
            VISP{n,j} = vertcat(VISP{n,j}, ps);
%             pFwd{n,j} = vertcat(pFwd{n,j}, pf);
            NSP{n,j} = vertcat(NSP{n,j}, pt);
            %             VISP{n,j} = vertcat(VISP{n,j}, sum(vispk)/length(vispk));
            %             NSP{n,j} = vertcat(NSP{n,j}, length(vispk));
        end
                VrAll{n} = vertcat(VrAll{n}, vra);
                VfAll{n} = vertcat(VfAll{n}, vfa);
                CmAll{n} = vertcat(CmAll{n}, cma);
                VrForw{n} = vertcat(VrForw{n}, vrf);
                VfForw{n} = vertcat(VfForw{n}, vff);
                CmForw{n} = vertcat(CmForw{n}, cmf);
                VrSpikes{n} = vertcat(VrSpikes{n}, vrs);
                VfSpikes{n} = vertcat(VfSpikes{n}, vfs);
                CmSpikes{n} = vertcat(CmSpikes{n}, cms);
    end
end
VISP2 = VISP;
NSP2 = NSP;
%%
VISP = VISP2;
NSP = NSP2;
figure,
hold on
inds = [5 6 7 8 9 10 1 2 3 4];
vs = cell(length(pTypes), 1);
vf = cell(length(pTypes), 1);
for j = 1 : length(pTypes)-2
    vaux = [];
    vT = [];
    a = inds(j);
    for n = 1 : length(flies)
        NSP{n,a}(isnan(VISP{n,a})) = 1;
        VISP{n,a}(isnan(VISP{n,a})) = 0;
%         NSP{n,a}(isempty(VISP{n,a})) = 1;
%         VISP{n,a}(isempty(VISP{n,a})) = 0;
        vaux = vertcat(vaux, VISP{n,a}'*NSP{n,a}/sum(NSP{n,a}));
        vT = vertcat(vT, sum(NSP{n,a}));
    end
    vT(isnan(vaux)) = [];
    vaux(isnan(vaux)) = [];
    vaux = vaux(1:length(vT));
    vs{j} = vaux;
    if mod(j,2) == 1
        errorbar(j, vaux'*vT/sum(vT), sqrt(sum((vaux-vaux'*vT/sum(vT)).*(vaux-vaux'*vT/sum(vT)))/length(flies))/sqrt(length(flies)), 'or')
    else
        errorbar(j-1, vaux'*vT/sum(vT), sqrt(sum((vaux-vaux'*vT/sum(vT)).*(vaux-vaux'*vT/sum(vT)))/length(flies))/sqrt(length(flies)), 'ob')
    end
end
axis([0 11 0 1])

%%
figure,
hold on
inds = [5 6 7 8 9 10 1 2 3 4];
vs = cell(length(pTypes), 1);
vf = cell(length(pTypes), 1);
for j = 1 : length(pTypes)-2
    vauxS = [];
    vauxF = [];
    vT = [];
    a = inds(j);
    for n = 1 : length(flies)
        pSpks{n,a}(isnan(pSpks{n,a})) = 0;
        pFwd{n,a}(isnan(pFwd{n,a})) = 0;
        vauxS = vertcat(vauxS, pSpks{n,a}'*pT{n,a}/sum(pT{n,a}));
        vauxF = vertcat(vauxF, pFwd{n,a}'*pT{n,a}/sum(pT{n,a}));
        vT = vertcat(vT, sum(pT{n,a}));
    end
    vs{j} = vauxS;
    vf{j} = vauxF;
    if mod(j,2) == 0
        errorbar(j, vauxS'*vT/sum(vT), sqrt(sum((vauxS-vauxS'*vT/sum(vT)).*(vauxS-vauxS'*vT/sum(vT)))/length(flies))/sqrt(length(flies)), 'or')
        errorbar(j, vauxF'*vT/sum(vT), sqrt(sum((vauxF-vauxF'*vT/sum(vT)).*(vauxF-vauxF'*vT/sum(vT)))/length(flies))/sqrt(length(flies)), 'ob')
    end
end
axis([0 11 0 1])


%%
% figure,
% for n = 1 : length(flies)
%     subplot(4, ceil(length(flies)/4), n)
%     [mapSum, ~, ~] = SpeedMapOccupancy(Vrs{n}, Vfs{n}, ones(size(Vrs{n})));
%     imagesc(-1200:20:1200, -10:1:40, mapSum./sum(sum(mapSum)))
%     caxis([0 0.005])
%     colormap jet
%     set(gca,'Ydir','normal')
%     disp(num2str(n))
% end
%%
vrra = [];
vffa = [];
vrrf = [];
vfff = [];
vrrs = [];
vffs = [];
cmma = [];
cmmf = [];
cmms = [];
for n = 1 : length(flies)
    vrra = vertcat(vrra, VrAll{n});
    vffa = vertcat(vffa, VfAll{n});
    cmma = vertcat(cmma, CmAll{n});
    vrrf = vertcat(vrrf, VrForw{n});
    vfff = vertcat(vfff, VfForw{n});
    cmmf = vertcat(cmmf, CmForw{n});
    vrrs = vertcat(vrrs, VrSpikes{n});
    vffs = vertcat(vffs, VfSpikes{n});
    cmms = vertcat(cmms, CmSpikes{n});
end
vrcent = linspace(-5,800,50);
vfcent = linspace(-5,35,40);
[mapSumAll, ~, ~] = SpeedMapOccupancy(abs(vrra), vffa, ones(size(vrra)), vrcent, vfcent);
[mapSumForw, ~, ~] = SpeedMapOccupancy(abs(vrrf), vfff, ones(size(vrrf)), vrcent, vfcent);
[mapSumSpikes, ~, ~] = SpeedMapOccupancy(abs(vrrs), vffs, ones(size(vrrs)), vrcent, vfcent);
% peido com molho
%%
figure,
subplot(2,2,1)
imagesc(-1:10:1000, -10:1:40, (mapSumSpikes)/sum(sum(mapSumSpikes)))
axis([0 1000 -10 40])
caxis([0 0.005])
colormap jet
set(gca,'Ydir','normal')
subplot(2,2,2)
imagesc(-1:10:1000, -10:1:40, (mapSumForw)/sum(sum(mapSumAll)))
axis([0 1000 -10 40])
caxis([0 0.005])
colormap jet
set(gca,'Ydir','normal')
subplot(2,2,3)
mapToPlot=(mapSumAll/sum(sum(mapSumAll))-mapSumSpikes/sum(sum(mapSumSpikes)));
mapToPlot(isnan(mapToPlot)) = 0;
imagesc(-1:10:1000, -10:1:40, mapToPlot)
axis([0 1000 -10 40])
caxis(0.003*[-1 1])
colormap redblue
set(gca,'Ydir','normal')
subplot(2,2,4)
mapToPlot=(mapSumAll/sum(sum(mapSumAll))-mapSumForw/sum(sum(mapSumForw)));
mapToPlot(isnan(mapToPlot)) = 0;
imagesc(-1:10:1000, -10:1:40, mapToPlot)
axis([0 1000 -10 40])
caxis(0.003*[-1 1])
colormap redblue
set(gca,'Ydir','normal')

%%
figure,
mapToPlot=(mapSumAll/sum(sum(mapSumAll))-mapSumForw/sum(sum(mapSumForw)));
% mapToPlot(isnan(mapToPlot)) = 0;
imagesc(-1:10:1000, -10:1:40, mapToPlot)
axis([0 800 -10 40])
caxis(0.004*[-1 1])
colormap redblue
set(gca,'Ydir','normal')


%%
% vrcent = linspace(-5,800,50);
% vfcent = linspace(-5,35,40);
% [mapSumAll, ~, ~] = SpeedMapOccupancy(abs(vrra), vffa, ones(size(vrra)), vrcent, vfcent);
%%
figure,
imagesc(vrcent, vfcent, ((mapSumAll)/sum(sum(mapSumAll))))
axis([0 800 -5 35])
caxis(0.0035*[0 1])
% caxis([-inf -4.5])
colormap hot
set(gca,'Ydir','normal')

%%
figure,
mapToPlot=mapSumForw/sum(sum(mapSumForw))-mapSumSpikes/sum(sum(mapSumSpikes));
% mapToPlot(isnan(mapToPlot)) = 0;
[X,Y] = meshgrid(vrcent,vfcent);
mapToCont = mapToPlot;
mapToCont(mapToCont > 0.0001) = 1;
mapToCont(mapToCont < -0.00005) = -1;
imagesc(vrcent, vfcent, mapToPlot)
hold on
contour(X, Y,mapToCont,1, 'k', 'linewidth',2)
axis([0 800 -5 35])
caxis(0.004*[-1 1])
colormap redblue
set(gca,'Ydir','normal')
