%% function to isolate leg movments during individual forward runs
function [LegMovForwSeg, pTypes] = GetLegMovForwSeg(path, params)
addpath('Saccade Segmenter')
flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowRes.mat'], 'Flies');
seq = dt.Flies.Seq;
pTypes = unique(seq);
LegMovForwSeg = cell(length(pTypes),length(flies));
for n = 1 : length(flies)-1
    pathF = [path flies(n).name];
    pathi = [pathF '\DataLowHighRes.mat'];
    if exist(pathi, 'file') == 2
        dt = load(pathi);
        dt = dt.Flies;
        seq = dt.Seq;
        pTypes = unique(seq);
        dt = dt.Data;
        for k = 1 : length(dt)
            for l = 1 : length(pTypes)
                switch seq{k}
                    case pTypes{l}
                        vrb = dt{k}.Vr;
                        vsb = dt{k}.Vs;
                        vfb = dt{k}.Vf;
                        actst = dt{k}.actState;
                        % Get saccades and forward runs
                        [locs, ~, cmhSong, thr] = SortSpikes(vrb, actst, params);
                        [locsF, pksIF,pksEF] = CullSpikesBasedOnVf(vrb, locs, cmhSong, thr, vfb, params);
                        [locsF, pksIF, pksEF] = TemplateCompSpike(vrb,locsF,pksIF, pksEF, cmhSong, thr, params);
                        [FBouts] = GetFbouts(actst, vfb, locsF, pksIF,pksEF, params);
                        
                        cSac = length(LegMovForwSeg{l,n});
                        for i = 1 : length(FBouts)
                            if (FBouts{i}(1)) > 0 && (FBouts{i}(end)) < length(dt{k}.FramesC2)
                                frms = (dt{k}.FramesC2(FBouts{i}(1)):dt{k}.FramesC2(FBouts{i}(end)))-dt{k}.FramesC2(1);
                                if min(frms) > 0 && max(frms) < length(dt{k}.LeftFrontLegX)
                                    LegMovForwSeg{l,n}{i+cSac}.VrLR = vrb(FBouts{i});
                                    LegMovForwSeg{l,n}{i+cSac}.VfLR = vfb(FBouts{i});
                                    LegMovForwSeg{l,n}{i+cSac}.VsLR = vsb(FBouts{i});
                                    LegMovForwSeg{l,n}{i+cSac}.FLX = dt{k}.LeftFrontLegX(frms);
                                    LegMovForwSeg{l,n}{i+cSac}.FLY = dt{k}.LeftFrontLegY(frms);
                                    LegMovForwSeg{l,n}{i+cSac}.FLErr = dt{k}.LeftFrontLegErr(frms);
                                    LegMovForwSeg{l,n}{i+cSac}.FRX = dt{k}.RightFrontLegX(frms);
                                    LegMovForwSeg{l,n}{i+cSac}.FRY = dt{k}.RightFrontLegY(frms);
                                    LegMovForwSeg{l,n}{i+cSac}.FRErr = dt{k}.RightFrontLegErr(frms);
                                    LegMovForwSeg{l,n}{i+cSac}.MLX = dt{k}.LeftMiddleLegX(frms);
                                    LegMovForwSeg{l,n}{i+cSac}.MLY = dt{k}.LeftMiddleLegY(frms);
                                    LegMovForwSeg{l,n}{i+cSac}.MLErr = dt{k}.LeftMiddleLegErr(frms);
                                    LegMovForwSeg{l,n}{i+cSac}.MRX = dt{k}.RightMiddleLegX(frms);
                                    LegMovForwSeg{l,n}{i+cSac}.MRY = dt{k}.RightMiddleLegY(frms);
                                    LegMovForwSeg{l,n}{i+cSac}.MRErr = dt{k}.RightMiddleLegErr(frms);
                                    LegMovForwSeg{l,n}{i+cSac}.HLX = dt{k}.LeftHindLegX(frms);
                                    LegMovForwSeg{l,n}{i+cSac}.HLY = dt{k}.LeftHindLegY(frms);
                                    LegMovForwSeg{l,n}{i+cSac}.HLErr = dt{k}.LeftHindLegErr(frms);
                                    LegMovForwSeg{l,n}{i+cSac}.HRX = dt{k}.RightHindLegX(frms);
                                    LegMovForwSeg{l,n}{i+cSac}.HRY = dt{k}.RightHindLegY(frms);
                                    LegMovForwSeg{l,n}{i+cSac}.HRErr = dt{k}.RightHindLegErr(frms);
                                end
                            end
                        end
                end
            end
        end
    end
    disp(['Fly ' num2str(n) ' Done'])
end
end