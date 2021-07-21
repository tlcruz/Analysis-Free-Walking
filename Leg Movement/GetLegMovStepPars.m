function [indsStance, indsSwing] = GetLegMovStepPars(FLLY)
% criteria for high quality tracking steps
peakProminence = 40;
maxStanceTime = 15;
minStanceTime = 2;
maxSwingTime = 15;
minSwingTime = 2;
maxStepLength = 130;
minStepLength = 0;
% find the timestamps of transitions from stance to swing
[~, indFLSwing] = findpeaks(FLLY, 'MinPeakProminence', peakProminence);
[~, indFLStance] = findpeaks(-FLLY, 'MinPeakProminence', peakProminence);

indsStance = [];
indsSwing = [];
for j = 1 : length(indFLStance)
    % get the previous swing
    indStance = indFLStance(j);
    indPrevSwing = find(indFLSwing<indStance);
    if ~isempty(indPrevSwing)
        indPrevSwing = indFLSwing(indPrevSwing(end));
        swingTime = indStance - indPrevSwing;
        stepLength = -(FLLY(indStance) - FLLY(indPrevSwing));
        % get the next swing
        indNextSwing = find(indFLSwing>indFLStance(j), 1);
        if  ~isempty(indNextSwing)
            indNextSwing = indFLSwing(indNextSwing);
            stanceTime = indNextSwing - indStance;
            % check if it fills the criteria
            if stanceTime < maxStanceTime && stanceTime > minStanceTime  && ...
                    swingTime < maxSwingTime && swingTime > minSwingTime && ...
                    stepLength < maxStepLength && stepLength > minStepLength
            
                indsStance = vertcat(indsStance, indStance);
                indsSwing = vertcat(indsSwing, indPrevSwing);

            end
        end
    end
end
end