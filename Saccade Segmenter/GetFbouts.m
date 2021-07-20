function [FWBoutStr, FWBout] = GetFbouts(actState, Vf, locsF, pksIF,pksEF, params)
% get vector with spike timmings
TSSub = zeros(size(actState));
for i = 1 : length(locsF)
    TSSub((floor(locsF(i)-pksIF(i))):(ceil(locsF(i)+pksEF(i)))) = 1;
end
% get vector with locomotion but no angular spikes
TSSub2 = zeros(size(Vf));
TSSub2(Vf > params.vft) = 1;
FWST = mod(actState+TSSub,2);
FWST = mod(FWST+TSSub2,3);
FWST(FWST ~= 2) = 0;
FWST(FWST == 2) = 1;
% transfor vector to bouts with a bout size threshold
[FWBoutStr] = TimeSeriesToBout(FWST, 20);
[FWBout] = TimeSeriesToBout(FWST, 0);
end