function [indsStance, indsSwing] = GetLegMovStepPars(FLLY)
peakProminence = 40;
maxStanceTime = 15;
minStanceTime = 2;
maxSwingTime = 15;
minSwingTime = 2;
maxStepLength = 130;
minStepLength = 0;
[~, indFLSwing] = findpeaks(FLLY, 'MinPeakProminence', peakProminence);
[~, indFLStance] = findpeaks(-FLLY, 'MinPeakProminence', peakProminence);

indsStance = [];
indsSwing = [];
% indsDoubleSwing = [];
% indsDoubleStance = [];
for j = 1 : length(indFLStance)
    indStance = indFLStance(j);
    indPrevSwing = find(indFLSwing<indStance);
    if ~isempty(indPrevSwing)
        indPrevSwing = indFLSwing(indPrevSwing(end));
        swingTime = indStance - indPrevSwing;
        stepLength = -(FLLY(indStance) - FLLY(indPrevSwing));
        indNextSwing = find(indFLSwing>indFLStance(j), 1);
        if  ~isempty(indNextSwing)
            indNextSwing = indFLSwing(indNextSwing);
            stanceTime = indNextSwing - indStance;
            
            if stanceTime < maxStanceTime && stanceTime > minStanceTime  && ...
                    swingTime < maxSwingTime && swingTime > minSwingTime && ...
                    stepLength < maxStepLength && stepLength > minStepLength
            
                indsStance = vertcat(indsStance, indStance);
                indsSwing = vertcat(indsSwing, indPrevSwing);
                
%                 indNextStance = find(indFLStance>indStance);
%                 if ~isempty(indNextStance)
%                     indNextStance = indFLStance(indNextStance(1));
%                     swingTime2 = indNextStance - indNextSwing;
%                     if swingTime2 < maxSwingTime && swingTime2 > minSwingTime
%                         indsDoubleStance = vertcat(indsDoubleStance, [indStance indNextStance]);
%                         indsDoubleSwing = vertcat(indsDoubleSwing, [indPrevSwing indNextSwing]);
%                     end
%                 end
                
            end
        end
    end
end
end