%% Requires Saccade segmenter
function [SPKVR, SPKVRRand, SPKVRDist,SPKVF, pcav, pTypes, dtout] = GetSpikePCA(path, params)
addpath('Saccade Segmenter')
% get the list of flies to analyse
flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowRes.mat'], 'Flies');
seq = dt.Flies.Seq;
pTypes = unique(seq);

% initialize variables for each fly and trial type 
SPKVR = cell(length(flies),length(pTypes));
SPKVF = cell(length(flies),length(pTypes));
SPKVS = cell(length(flies),length(pTypes));
SPKVRRand = cell(length(flies),length(pTypes));
SPKVFRand = cell(length(flies),length(pTypes));
SPKVSRand = cell(length(flies),length(pTypes));
SPKVRDist = cell(length(flies),length(pTypes));

% define spike window
wi = 20;
dm = 7;
for n = 1 : length(flies)
    % load data for a specific fly
    dt = load([path flies(n).name '\DataLowRes.mat'], 'Flies');
    seq = dt.Flies.Seq;
    dt = dt.Flies.Data;
    % iterate across protocol segments
    for k = 1 : length(dt)
        % iterate across trial types
        for j = 1 : length(pTypes)
            spkvr = [];
            spkvf = [];
            spkvs = [];
            spkvrrand = [];
            spkvfrand = [];
            spkvsrand = [];
            spkdist = [];
            switch seq{k}
                % when trial type matches a protocol segment
                case pTypes{j}
                    % Load walking parameters for this protocol segment
                    Vf = dt{k}.Vf;
                    Vs = dt{k}.Vs;
                    Vr = dt{k}.Vr;
                    acSt = dt{k}.actState;
                    % segment angular spikes
                    [locs, ~, cmhSong, thr] = SortSpikes(Vr, acSt, params);
                    [locsF, pksIF,pksEF, ~] = CullSpikesBasedOnVf(Vr, locs, cmhSong, thr, Vf, params);
                    [locsF, ~, ~] = TemplateCompSpike(Vr,locsF,pksIF, pksEF, cmhSong, thr, params);
                    locsRand = randi([wi+dm+1 length(Vr)-(wi+dm+1)], length(locsF),1);
                    % isolate a window around each spike
                    if(length(locsF) > 3)
                        for st = 1 : length(locsF)
                            if st == 1
                                if locsF(st)>wi+dm
                                    if (locsF(st+1)-locsF(st) > wi+dm)
                                        spkvr = horzcat(spkvr, Vr((locsF(st)-wi):(locsF(st)+wi)));
                                        spkvf = horzcat(spkvf, Vf((locsF(st)-wi):(locsF(st)+wi)));
                                        spkvs = horzcat(spkvs, Vs((locsF(st)-wi):(locsF(st)+wi)));
                                        spkvrrand = horzcat(spkvrrand,  Vr((locsRand(st)-wi):(locsRand(st)+wi)));
                                        spkvfrand = horzcat(spkvfrand,  Vf((locsRand(st)-wi):(locsRand(st)+wi)));
                                        spkvsrand = horzcat(spkvsrand,  Vs((locsRand(st)-wi):(locsRand(st)+wi)));
                                        spkdist = horzcat(spkdist, Vr(locsF(st)));
                                    end
                                end
                            elseif st == length(locsF)
                                if locsF(st)+wi+dm < length(Vr)
                                    if (locsF(st)-locsF(st-1) > wi+dm)
                                        spkvr = horzcat(spkvr, Vr((locsF(st)-wi):(locsF(st)+wi)));
                                        spkvf = horzcat(spkvf, Vf((locsF(st)-wi):(locsF(st)+wi)));
                                        spkvs = horzcat(spkvs, Vs((locsF(st)-wi):(locsF(st)+wi)));
                                        spkvrrand = horzcat(spkvrrand,  Vr((locsRand(st)-wi):(locsRand(st)+wi)));
                                        spkvfrand = horzcat(spkvfrand,  Vf((locsRand(st)-wi):(locsRand(st)+wi)));
                                        spkvsrand = horzcat(spkvsrand,  Vf((locsRand(st)-wi):(locsRand(st)+wi)));
                                        spkdist = horzcat(spkdist, Vr(locsF(st)));
                                    end
                                end
                            else
                                if (locsF(st)-locsF(st-1) > wi+dm && locsF(st+1)-locsF(st) > wi+dm)
                                    spkvr = horzcat(spkvr, Vr((locsF(st)-wi):(locsF(st)+wi)));
                                    spkvf = horzcat(spkvf, Vf((locsF(st)-wi):(locsF(st)+wi)));
                                    spkvs = horzcat(spkvs, Vs((locsF(st)-wi):(locsF(st)+wi)));
                                    spkvrrand = horzcat(spkvrrand,  Vr((locsRand(st)-wi):(locsRand(st)+wi)));
                                    spkvfrand = horzcat(spkvfrand,  Vf((locsRand(st)-wi):(locsRand(st)+wi)));
                                    spkvsrand = horzcat(spkvsrand,  Vs((locsRand(st)-wi):(locsRand(st)+wi)));
                                    spkdist = horzcat(spkdist, Vr(locsF(st)));
                                end
                            end
                        end
                    end
            end
            % store spikes in matrices
            SPKVR{n,j} = horzcat(SPKVR{n,j}, spkvr);
            SPKVF{n,j} = horzcat(SPKVF{n,j}, spkvf);
            SPKVS{n,j} = horzcat(SPKVS{n,j}, spkvs);
            SPKVRRand{n,j} = horzcat(SPKVRRand{n,j}, spkvrrand);
            SPKVFRand{n,j} = horzcat(SPKVFRand{n,j}, spkvfrand);
            SPKVSRand{n,j} = horzcat(SPKVSRand{n,j}, spkvsrand);
            SPKVRDist{n,j} = horzcat(SPKVRDist{n,j}, spkdist);
        end
    end
end
dtout.SPKVR = SPKVR;
dtout.SPKVF = SPKVF;
dtout.SPKVS = SPKVS;
dtout.SPKVRRand = SPKVRRand;
dtout.SPKVFRand = SPKVFRand;
dtout.SPKVSRand = SPKVSRand;
dtout.SPKVRDist = SPKVRDist;

% run PCA in the spike matrix
pcav.EXP = cell(length(pTypes),1);
pcav.NSP = cell(length(pTypes),1);
pcav.PC1 = cell(length(pTypes),1);
pcav.SCO = cell(length(flies),length(pTypes));
pcav.Dist = cell(length(pTypes),1);
pcav.VrCents = -10:50:1500;
for j = 1 : length(pTypes)
    for nn = 1 : length(flies)
        pcav.NSP{j} = vertcat(pcav.NSP{j}, size(SPKVR{nn,j},2));
        [coeff,score,~,~,explained,~] = pca(SPKVR{nn,j}');
        if size(coeff,2) > 0
            pcav.PC1{j} = horzcat(pcav.PC1{j}, coeff(:,1));
            pcav.SCO{nn,j} = horzcat(pcav.SCO{nn,j}, score(:,1));
            expl = zeros(6,1);
            expl(1:min(6,length(explained))) = explained(1:min(6,length(explained)));
            pcav.EXP{j} = horzcat(pcav.EXP{j}, expl);
        end
        pcav.Dist{j} = vertcat(pcav.Dist{j}, hist(abs(SPKVRDist{nn,j}), ...
            pcav.VrCents)/sum(hist(abs(SPKVRDist{nn,j}), pcav.VrCents)));
    end
end

end