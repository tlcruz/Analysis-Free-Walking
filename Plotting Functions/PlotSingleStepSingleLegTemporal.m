clc
addpath('Leg Movement')
% load the low and high resolution data relative to the forward segments
path = '\';
params = GetParams();
[LegMovForwSeg] = GetLegMovForwSeg(path, params);

thr = 0.02;
VrCents = -20:1:20;

% Individual Leg Step Spatial Parameters
LeftLegSwingTime = [];
NLeftLegSwingTime = [];
RightLegSwingTime = [];
NRightLegSwingTime = [];
LeftLegStanceTime = [];
NLeftLegStanceTime = [];
RightLegStanceTime = [];
NRightLegStanceTime = [];

% Iterate across trial types
for j = 1 : 2 : size(LegMovForwSeg, 1)
    % matrices to put all flies
    leftLegSwingTime = [];
    nleftLegSwingTime = [];
    rightLegSwingTime = [];
    nrightLegSwingTime = [];
    leftLegStanceTime = [];
    nleftLegStanceTime = [];
    rightLegStanceTime = [];
    nrightLegStanceTime = [];
    % iterate across flies
    for n = 1 : size(LegMovForwSeg, 2)
        [lD] = GetCleanLegData(LegMovForwSeg{j,n}, thr);
        % vectors to put all step data
        stepLeftLegSwingTime = [];
        stepLeftLegStanceTime = [];
        stepRightLegSwingTime = [];
        stepRightLegStanceTime = [];
        stepVr = [];
        for i = 1 : length(lD.VR)
            if ~isempty(lD.FLLY{i}) && ~isempty(lD.FRLY{i})
                % get the step parameters for specific legs
                stepPars1 = GetStepParameters(lD.FLLY{i}, lD.FRLY{i}, lD.FLLX{i}, lD.FRLX{i}, lD.VR{i}, lD.VF{i});
                for k = 1 : length(stepPars1.MVRL1)
                    stepLeftLegSwingTime = vertcat(stepLeftLegSwingTime, ...
                        stepPars1.L1SwingTime(k));
                    stepLeftLegStanceTime = vertcat(stepLeftLegStanceTime, ...
                        stepPars1.L1StanceTime(k));
                    stepRightLegSwingTime = vertcat(stepRightLegSwingTime, ...
                        stepPars1.L2SwingTime(k));
                    stepRightLegStanceTime = vertcat(stepRightLegStanceTime, ...
                        stepPars1.L2StanceTime(k));
                    stepVr = vertcat(stepVr, (stepPars1.MVRL1(k)+stepPars1.MVRL2(k))/2);
                end
                % get the step parameters for specific legs
                stepPars2 = GetStepParameters(lD.FRLY{i}, lD.FLLY{i}, lD.FRLX{i}, lD.FLLX{i}, lD.VR{i}, lD.VF{i});
                for k = 1 : length(stepPars2.MVRL1)
                    stepLeftLegSwingTime = vertcat(stepLeftLegSwingTime, ...
                        stepPars2.L2SwingTime(k));
                    stepLeftLegStanceTime = vertcat(stepLeftLegStanceTime, ...
                        stepPars2.L2StanceTime(k));
                    stepRightLegSwingTime = vertcat(stepRightLegSwingTime, ...
                        stepPars2.L1SwingTime(k));
                    stepRightLegStanceTime = vertcat(stepRightLegStanceTime, ...
                        stepPars2.L1StanceTime(k));
                    stepVr = vertcat(stepVr, (stepPars2.MVRL1(k)+stepPars2.MVRL2(k))/2);
                end
            end
        end
        
        % calculate distributions
        lSwingTimeAux = nan;
        nlSwingTimeAux = 0;
        lStanceTimeAux = nan;
        nlStanceTimeAux = 0;
        rSwingTimeAux = nan;
        nrSwingTimeAux = 0;
        rStanceTimeAux = nan;
        nrStanceTimeAux = 0;
        
        for i = 1 : length(VrCents)-1
            inds = find(stepVr > VrCents(i) & stepVr < VrCents(i+1));
            if~isempty(inds)
                lSwingTimeAux = vertcat(lSwingTimeAux, mean(stepLeftLegSwingTime(inds)));
                nlSwingTimeAux = vertcat(nlSwingTimeAux, length(inds));
                lStanceTimeAux = vertcat(lStanceTimeAux, mean(stepLeftLegStanceTime(inds)));
                nlStanceTimeAux = vertcat(nlStanceTimeAux, length(inds));
                rSwingTimeAux = vertcat(rSwingTimeAux, mean(stepRightLegSwingTime(inds)));
                nrSwingTimeAux = vertcat(nrSwingTimeAux, length(inds));
                rStanceTimeAux = vertcat(rStanceTimeAux, mean(stepRightLegStanceTime(inds)));
                nrStanceTimeAux = vertcat(nrStanceTimeAux, length(inds));
            else
                lSwingTimeAux = vertcat(lSwingTimeAux, nan);
                nlSwingTimeAux = vertcat(nlSwingTimeAux, 0);
                lStanceTimeAux = vertcat(lStanceTimeAux, nan);
                nlStanceTimeAux = vertcat(nlStanceTimeAux, 0);
                rSwingTimeAux = vertcat(rSwingTimeAux, nan);
                nrSwingTimeAux = vertcat(nrSwingTimeAux, 0);
                rStanceTimeAux = vertcat(rStanceTimeAux, nan);
                nrStanceTimeAux = vertcat(nrStanceTimeAux, 0);
            end
        end
        leftLegSwingTime = horzcat(leftLegSwingTime, lSwingTimeAux);
        nleftLegSwingTime = horzcat(nleftLegSwingTime, nlSwingTimeAux);
        rightLegSwingTime = horzcat(rightLegSwingTime, rSwingTimeAux);
        nrightLegSwingTime = horzcat(nrightLegSwingTime, nrSwingTimeAux);
        leftLegStanceTime = horzcat(leftLegStanceTime, lStanceTimeAux);
        nleftLegStanceTime = horzcat(nleftLegStanceTime, nlStanceTimeAux);
        rightLegStanceTime = horzcat(rightLegStanceTime, rStanceTimeAux);
        nrightLegStanceTime = horzcat(nrightLegStanceTime, nrStanceTimeAux);
    end
    
    LeftLegSwingTime = horzcat(LeftLegSwingTime, leftLegSwingTime);
    NLeftLegSwingTime = horzcat(NLeftLegSwingTime, nleftLegSwingTime);
    RightLegSwingTime = horzcat(RightLegSwingTime, rightLegSwingTime);
    NRightLegSwingTime = horzcat(NRightLegSwingTime, nrightLegSwingTime);
    LeftLegStanceTime = horzcat(LeftLegStanceTime, leftLegStanceTime);
    NLeftLegStanceTime = horzcat(NLeftLegStanceTime, nleftLegStanceTime);
    RightLegStanceTime = horzcat(RightLegStanceTime, rightLegStanceTime);
    NRightLegStanceTime = horzcat(NRightLegStanceTime, nrightLegStanceTime);
    
end
% get the grand mean of the step parameters
[gmLeftLegSwingTime, semLeftLegSwingTime] = GetGMSEMFromMatrix(LeftLegSwingTime, NLeftLegSwingTime,4);
[gmRightLegSwingTime, semRightLegSwingTime] = GetGMSEMFromMatrix(RightLegSwingTime, NRightLegSwingTime,4);
[gmLeftLegStanceTime, semLeftLegStanceTime] = GetGMSEMFromMatrix(LeftLegStanceTime, NLeftLegStanceTime,4);
[gmRightLegStanceTime, semRightLegStanceTime] = GetGMSEMFromMatrix(RightLegStanceTime, NRightLegStanceTime,4);

% plot the swing time for the pair of legs
figure,
hold on
% plot([VrCents(1) VrCents(end)], [0 0], '--g')
plot([0 0], [0 100], '--m')
plot(VrCents, 1000*gmLeftLegSwingTime, 'k', 'linewidth', 2)
plot(VrCents, 1000*gmLeftLegSwingTime+1000*semLeftLegSwingTime, 'k', 'linewidth', 1)
plot(VrCents, 1000*gmLeftLegSwingTime-1000*semLeftLegSwingTime, 'k', 'linewidth', 1)
axis([-15 15 0 100])
xlabel('Angular Deflection (º)')
ylabel('Left Leg Swing Time (ms)')
title('All conditions')

figure,
hold on
% plot([VrCents(1) VrCents(end)], [0 0], '--g')
plot([0 0], [0 100], '--m')
plot(VrCents, 1000*gmRightLegSwingTime, 'k', 'linewidth', 2)
plot(VrCents, 1000*gmRightLegSwingTime+1000*semRightLegSwingTime, 'k', 'linewidth', 1)
plot(VrCents, 1000*gmRightLegSwingTime-1000*semRightLegSwingTime, 'k', 'linewidth', 1)
axis([-15 15 0 100])
xlabel('Angular Deflection (º)')
ylabel('Right Leg Swing Time (ms)')
title('All conditions')

% plot the stance time for the pair of legs
figure,
hold on
% plot([VrCents(1) VrCents(end)], [0 0], '--g')
plot([0 0], [0 100], '--m')
plot(VrCents, 1000*gmLeftLegStanceTime, 'k', 'linewidth', 2)
plot(VrCents, 1000*gmLeftLegStanceTime+1000*semLeftLegStanceTime, 'k', 'linewidth', 1)
plot(VrCents, 1000*gmLeftLegStanceTime-1000*semLeftLegStanceTime, 'k', 'linewidth', 1)
axis([-15 15 0 100])
xlabel('Angular Deflection (º)')
ylabel('Left Leg Stance Time (ms)')
title('All conditions')

figure,
hold on
% plot([VrCents(1) VrCents(end)], [0 0], '--g')
plot([0 0], [0 100], '--m')
plot(VrCents, 1000*gmRightLegStanceTime, 'k', 'linewidth', 2)
plot(VrCents, 1000*gmRightLegStanceTime+1000*semRightLegStanceTime, 'k', 'linewidth', 1)
plot(VrCents, 1000*gmRightLegStanceTime-1000*semRightLegStanceTime, 'k', 'linewidth', 1)
axis([-15 15 0 100])
xlabel('Angular Deflection (º)')
ylabel('Right Leg Stance Time (ms)')
title('All conditions')