% load optomotor data and protocol types
addpath('Saccade Segmenter')
path = '\Optomotor 1D 10D\';
flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowRes.mat'], 'Flies');
seq = dt.Flies.Seq;
pTypes = unique(seq);
VR = cell(size(pTypes));
VF = cell(size(pTypes));
params = GetParams();
for n = 1: length(flies)
    pathF = [path flies(n).name];
    pathi = [pathF '\DataLowHighRes.mat'];
    if exist(pathi, 'file') == 2
        dt = load(pathi);
        dt = dt.Flies;
        seq = dt.Seq;
        pTypes = unique(seq);
        dt = dt.Data;
        for k = 2 : length(dt)-1
            for l = 1 : length(pTypes)
                switch seq{k}
                    case pTypes{l}
                        vrb = dt{k}.Vr;
                        vsb = dt{k}.Vs;
                        vfb = dt{k}.Vf;
                        actst = dt{k}.actState;
                        % segment saccades and forward runs
                        [locs, ~, cmhSong, thr] = SortSpikes(vrb, actst, params);
                        [locsF, pksIF,pksEF] = CullSpikesBasedOnVf(vrb, locs, cmhSong, thr, vfb, params);
                        [locsF, pksIF, pksEF] = TemplateCompSpike(vrb,locsF,pksIF, pksEF, cmhSong, thr, params);
                        [~,FBouts] = GetFbouts(actst, vfb, locsF, pksIF,pksEF, params);
                        
                        bouts = FBouts;
                        auxvr = zeros(size(vrb));
                        auxvf = zeros(size(vrb));
                        for i = 1 : length(bouts)
                            auxvr(bouts{i}) = 1;
                        end
%                         for i = 1 : length(locsF)
% %                             auxvr((locsF(i)-pksIF(i)):(locsF(i)+pksEF(i))) = 1;%vrb(bouts{i});
%                             auxvf((locsF(i)-pksIF(i)):(locsF(i)+pksEF(i))) = 1;%vfb(bouts{i});
%                         end
                        VR{l} = horzcat(VR{l}, auxvr);
                        VF{l} = horzcat(VF{l}, auxvf);
                end
            end
        end
    end
    disp(['Fly ' num2str(n) ' Done'])
end
% plot angular velocity under different conditions
figure,
for i = 1 : 6
subplot(2,3,i)
hold on
plot((1:601)/60,nanmean(VR{i},2), 'k', 'linewidth', 3)
plot((1:601)/60,nanmean(VR{i},2)+nanstd(VF{i},1,2)/sqrt(30), 'k', 'linewidth', 1)
plot((1:601)/60,nanmean(VR{i},2)-nanstd(VF{i},1,2)/sqrt(30), 'k', 'linewidth', 1)
axis([0 10 -5 50])
title(pTypes{i})
end






















