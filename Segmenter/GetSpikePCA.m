function [SPKVR, SPKVRRand, SPKVRDist,SPKVF, pcav, pTypes] = GetSpikePCA(path, params)

flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowRes.mat'], 'Flies');
seq = dt.Flies.Seq;
pTypes = unique(seq);

SPKVR = cell(length(flies),length(pTypes));
SPKVF = cell(length(flies),length(pTypes));
SPKVRRand = cell(length(flies),length(pTypes));
SPKVFRand = cell(length(flies),length(pTypes));
SPKVRDist = cell(length(flies),length(pTypes));
wi = 20;
dm = 7;
for n = 1 : length(flies)
    dt = load([path flies(n).name '\DataLowRes.mat'], 'Flies');
    seq = dt.Flies.Seq;
    dt = dt.Flies.Data;
    for k = 1 : length(dt)
        for j = 1 : length(pTypes)
            spkvr = [];
            spkvf = [];
            spkvrrand = [];
            spkvfrand = [];
            spkdist = [];
            switch seq{k}
                case pTypes{j}%'10NG'%
                    Vf = dt{k}.Vf;
%                     Vs = dt{k}.Vs;
                    Vr = dt{k}.Vr;
                    acSt = dt{k}.actState;
%                     [params] = GetParams();
                    [locs, ~, cmhSong, thr] = SortSpikes(Vr, acSt, params);
                    [locsF, pksIF,pksEF, ~] = CullSpikesBasedOnVf(Vr, locs, cmhSong, thr, Vf, params);
                    [locsF, ~, ~] = TemplateCompSpike(Vr,locsF,pksIF, pksEF, cmhSong, thr, params);
%                     [FWBout] = GetFbouts(acSt, Vf, locsF, pksIF,pksEF, params);
                    locsRand = randi([wi+dm+1 length(Vr)-(wi+dm+1)], length(locsF),1);
                    if(length(locsF) > 3)
                        for st = 1 : length(locsF)
                            if st == 1
                                if locsF(st)>wi+dm
                                    if (locsF(st+1)-locsF(st) > wi+dm)
                                        spkvr = horzcat(spkvr, Vr((locsF(st)-wi):(locsF(st)+wi)));
                                        spkvf = horzcat(spkvf, Vf((locsF(st)-wi):(locsF(st)+wi)));
                                        spkvrrand = horzcat(spkvrrand,  Vr((locsRand(st)-wi):(locsRand(st)+wi)));
                                        spkvfrand = horzcat(spkvfrand,  Vf((locsRand(st)-wi):(locsRand(st)+wi)));
                                        spkdist = horzcat(spkdist, Vr(locsF(st)));
                                    end
                                end
                            elseif st == length(locsF)
                                if locsF(st)+wi+dm < length(Vr)
                                    if (locsF(st)-locsF(st-1) > wi+dm)
                                        spkvr = horzcat(spkvr, Vr((locsF(st)-wi):(locsF(st)+wi)));
                                        spkvf = horzcat(spkvf, Vf((locsF(st)-wi):(locsF(st)+wi)));
                                        spkvrrand = horzcat(spkvrrand,  Vr((locsRand(st)-wi):(locsRand(st)+wi)));
                                        spkvfrand = horzcat(spkvfrand,  Vf((locsRand(st)-wi):(locsRand(st)+wi)));
                                        spkdist = horzcat(spkdist, Vr(locsF(st)));
                                    end
                                end
                            else
                                if (locsF(st)-locsF(st-1) > wi+dm && locsF(st+1)-locsF(st) > wi+dm)
                                    spkvr = horzcat(spkvr, Vr((locsF(st)-wi):(locsF(st)+wi)));
                                    spkvf = horzcat(spkvf, Vf((locsF(st)-wi):(locsF(st)+wi)));
                                    spkvrrand = horzcat(spkvrrand,  Vr((locsRand(st)-wi):(locsRand(st)+wi)));
                                    spkvfrand = horzcat(spkvfrand,  Vf((locsRand(st)-wi):(locsRand(st)+wi)));
                                    spkdist = horzcat(spkdist, Vr(locsF(st)));
                                end
                            end
                            
                        end
                    end
            end
            SPKVR{n,j} = horzcat(SPKVR{n,j}, spkvr);
            SPKVF{n,j} = horzcat(SPKVF{n,j}, spkvf);
            SPKVRRand{n,j} = horzcat(SPKVRRand{n,j}, spkvrrand);
            SPKVFRand{n,j} = horzcat(SPKVFRand{n,j}, spkvfrand);
            SPKVRDist{n,j} = horzcat(SPKVRDist{n,j}, spkdist);
        end
    end
end

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
        pcav.PC1{j} = horzcat(pcav.PC1{j}, coeff(:,1));
        pcav.SCO{nn,j} = horzcat(pcav.SCO{nn,j}, score(:,1));
        expl = zeros(6,1);
        expl(1:min(6,length(explained))) = explained(1:min(6,length(explained)));
        pcav.EXP{j} = horzcat(pcav.EXP{j}, expl);
        pcav.Dist{j} = vertcat(pcav.Dist{j}, hist(abs(SPKVRDist{nn,j}), ...
            pcav.VrCents)/sum(hist(abs(SPKVRDist{nn,j}), pcav.VrCents)));
    end
end

end