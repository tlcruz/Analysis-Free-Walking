function [locs, p, cmhSong, thr, dx] = SortSpikes(Vr, acSt, params)
% get the frequency power signal from Vr in a predefined frequency band 
[cmhSong] = GetMainFreqs(Vr, params);
cmhSong = double(cmhSong);
sig4Test = cmhSong;
% get noise level by looking at power signal during no locomotion
sig4Test2 = sig4Test;
sig4Test2(acSt==1) = [];
thr = mean(sig4Test2) + params.nstds*std(sig4Test2);
% cut power signal at noise level
sig4Test(sig4Test<thr) = thr;

% if ther is too much locomotion use pre-defined value
if (sum(acSt)/length(acSt) > 0.75) || thr > 75
    thr = 45;
end

% get peaks in speed and frequency spaces
[~,locsV,~,pV] = findpeaks(abs(Vr),'MinPeakDistance',params.MinPeakDist,...
    'MinPeakProminence',params.MinPeakPromV);
[~,locsF,~,~] = findpeaks(sig4Test,'MinPeakDistance',params.MinPeakDist,...
    'MinPeakProminence',params.MinPeakPromF);

dx.locsV = locsV;
dx.pV = pV;
dx.locsF = locsF;

% compare frequency and speed peaks and keep only when they match
p = [];
locs = [];
for i = 1 : length(locsF)
    auxL = 10000;
    locA = 0;
    pA = 0;
    for j = 1 : length(locsV)
        if(abs(locsV(j) - locsF(i)) < auxL)
            auxL = abs(locsV(j) - locsF(i));
            locA = locsV(j);
            pA = pV(j);
        end
    end
    if auxL < params.thrDist
        locs = vertcat(locs, locA);
        p = vertcat(p, pA);
    end
end
locsV = locsV(pV > 350);
for i = 1 : length(locsV)
    if(isempty(find(locs == locsV(i), 1)))
       locs =  vertcat(locs,locsV(i));
       p = vertcat(p, pV(i));
    end
end
locs = unique(locs);

for i = 1 : length(locs)
    if acSt(locs(i)) == 0
        locs(i) = nan;
        p(i) = nan;
    end
end
p(isnan(locs)) = [];
locs(isnan(locs)) = [];
end