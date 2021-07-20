clc
addpath('Leg Movement')
% load the low and high resolution data relative to the forward segments
path = '\';
params = GetParams();
[LegMovForwSeg] = GetLegMovForwSeg(path, params);

thr = 0.02;
VrCents = -20:1:20;

% Individual Leg Step Spatial Parameters
LeftLegLatMov = [];
NLeftLegLatMov = [];
RightLegLatMov = [];
NRightLegLatMov = [];
LeftLegStepLen = [];
NLeftLegStepLen = [];
RightLegStepLen = [];
NRightLegStepLen = [];
vvr = [];
ll = [];

% Iterate across trial types
for j = 1 : 2 : size(LegMovForwSeg, 1)
    % matrices to put all flies
    leftLegLatMov = [];
    nleftLegLatMov = [];
    rightLegLatMov = [];
    nrightLegLatMov = [];
    leftLegStepLen = [];
    nleftLegStepLen = [];
    rightLegStepLen = [];
    nrightLegStepLen = [];
    % iterate across flies
    for n = 1 : size(LegMovForwSeg, 2)
        [lD] = GetCleanLegData(LegMovForwSeg{j,n}, thr);
        % vectors to put all step data
        stepLeftLegLatMov = [];
        stepLeftLegStepLen = [];
        stepRightLegLatMov = [];
        stepRightLegStepLen = [];
        stepVr = [];
        for i = 1 : length(lD.VR)
            if ~isempty(lD.HLLY{i}) && ~isempty(lD.HRLY{i})
                % get the step parameters for specific legs
                stepPars1 = GetStepParameters(lD.HLLY{i}, lD.HRLY{i}, lD.HLLX{i}, lD.HRLX{i}, lD.VR{i}, lD.VF{i});
                for k = 1 : length(stepPars1.MVRL1)
                    stepLeftLegLatMov = vertcat(stepLeftLegLatMov, ...
                        stepPars1.L1SwingLateral(k));
                    stepLeftLegStepLen = vertcat(stepLeftLegStepLen, ...
                        stepPars1.L1SwingLength(k));
                    stepRightLegLatMov = vertcat(stepRightLegLatMov, ...
                        stepPars1.L2SwingLateral(k));
                    stepRightLegStepLen = vertcat(stepRightLegStepLen, ...
                        stepPars1.L2SwingLength(k));
                    stepVr = vertcat(stepVr, (stepPars1.MVRL1(k)+stepPars1.MVRL2(k))/2);
                end
                % get the step parameters for specific legs
                stepPars2 = GetStepParameters(lD.HRLY{i}, lD.HLLY{i}, lD.HRLX{i}, lD.HLLX{i}, lD.VR{i}, lD.VF{i});
                for k = 1 : length(stepPars2.MVRL1)
                    stepLeftLegLatMov = vertcat(stepLeftLegLatMov, ...
                        stepPars2.L2SwingLateral(k));
                    stepLeftLegStepLen = vertcat(stepLeftLegStepLen, ...
                        stepPars2.L2SwingLength(k));
                    stepRightLegLatMov = vertcat(stepRightLegLatMov, ...
                        stepPars2.L1SwingLateral(k));
                    stepRightLegStepLen = vertcat(stepRightLegStepLen, ...
                        stepPars2.L1SwingLength(k));
                    stepVr = vertcat(stepVr, (stepPars2.MVRL1(k)+stepPars2.MVRL2(k))/2);
                end
            end
        end
        ll = vertcat(ll, stepRightLegLatMov);
        vvr = vertcat(vvr, stepVr);
        
        % calculate distributions
        lLatMovAux = nan;
        nlLatMovAux = 0;
        lSLenAux = nan;
        nlSLenAux = 0;
        rLatMovAux = nan;
        nrLatMovAux = 0;
        rSLenAux = nan;
        nrSLenAux = 0;
        
        for i = 1 : length(VrCents)-1
            inds = find(stepVr > VrCents(i) & stepVr < VrCents(i+1));
            if~isempty(inds)
                lLatMovAux = vertcat(lLatMovAux, mean(stepLeftLegLatMov(inds)));
                nlLatMovAux = vertcat(nlLatMovAux, length(inds));
                lSLenAux = vertcat(lSLenAux, mean(stepLeftLegStepLen(inds)));
                nlSLenAux = vertcat(nlSLenAux, length(inds));
                rLatMovAux = vertcat(rLatMovAux, mean(stepRightLegLatMov(inds)));
                nrLatMovAux = vertcat(nrLatMovAux, length(inds));
                rSLenAux = vertcat(rSLenAux, mean(stepRightLegStepLen(inds)));
                nrSLenAux = vertcat(nrSLenAux, length(inds));
            else
                lLatMovAux = vertcat(lLatMovAux, nan);
                nlLatMovAux = vertcat(nlLatMovAux, 0);
                lSLenAux = vertcat(lSLenAux, nan);
                nlSLenAux = vertcat(nlSLenAux, 0);
                rLatMovAux = vertcat(rLatMovAux, nan);
                nrLatMovAux = vertcat(nrLatMovAux, 0);
                rSLenAux = vertcat(rSLenAux, nan);
                nrSLenAux = vertcat(nrSLenAux, 0);
            end
        end
        leftLegLatMov = horzcat(leftLegLatMov, lLatMovAux);
        nleftLegLatMov = horzcat(nleftLegLatMov, nlLatMovAux);
        rightLegLatMov = horzcat(rightLegLatMov, rLatMovAux);
        nrightLegLatMov = horzcat(nrightLegLatMov, nrLatMovAux);
        leftLegStepLen = horzcat(leftLegStepLen, lSLenAux);
        nleftLegStepLen = horzcat(nleftLegStepLen, nlSLenAux);
        rightLegStepLen = horzcat(rightLegStepLen, rSLenAux);
        nrightLegStepLen = horzcat(nrightLegStepLen, nrSLenAux);
    end
    
    LeftLegLatMov = horzcat(LeftLegLatMov, leftLegLatMov);
    NLeftLegLatMov = horzcat(NLeftLegLatMov, nleftLegLatMov);
    RightLegLatMov = horzcat(RightLegLatMov, rightLegLatMov);
    NRightLegLatMov = horzcat(NRightLegLatMov, nrightLegLatMov);
    LeftLegStepLen = horzcat(LeftLegStepLen, leftLegStepLen);
    NLeftLegStepLen = horzcat(NLeftLegStepLen, nleftLegStepLen);
    RightLegStepLen = horzcat(RightLegStepLen, rightLegStepLen);
    NRightLegStepLen = horzcat(NRightLegStepLen, nrightLegStepLen);
    
end
% get the grand mean of the step parameters
[gmLeftLegLatMov, semLeftLegLatMov] = GetGMSEMFromMatrix(LeftLegLatMov, NLeftLegLatMov,4);
[gmRightLegLatMov, semRightLegLatMov] = GetGMSEMFromMatrix(RightLegLatMov, NRightLegLatMov,4);
[gmLeftLegStepLen, semLeftLegStepLen] = GetGMSEMFromMatrix(LeftLegStepLen, NLeftLegStepLen,4);
[gmRightLegStepLen, semRightLegStepLenv] = GetGMSEMFromMatrix(RightLegStepLen, NRightLegStepLen,4);

%% Lateral Movement
% plot the movement in the X direction (lateral) for the pair of legs
figure,
hold on
plot([VrCents(1) VrCents(end)], [0 0], '--g')
plot([0 0], [-0.3 0.3], '--m')
plot(VrCents, gmLeftLegLatMov, 'k', 'linewidth', 2)
plot(VrCents, gmLeftLegLatMov+semLeftLegLatMov, 'k', 'linewidth', 1)
plot(VrCents, gmLeftLegLatMov-semLeftLegLatMov, 'k', 'linewidth', 1)
axis([-15 15 -0.2 0.2])
xlabel('Angular Deflection (º)')
ylabel('Left Leg Lateral Movement (BL)')
title('All conditions')

figure,
hold on
plot([VrCents(1) VrCents(end)], [0 0], '--g')
plot([0 0], [-0.3 0.3], '--m')
plot(VrCents, gmRightLegLatMov, 'k', 'linewidth', 2)
plot(VrCents, gmRightLegLatMov+semRightLegLatMov, 'k', 'linewidth', 1)
plot(VrCents, gmRightLegLatMov-semRightLegLatMov, 'k', 'linewidth', 1)
axis([-15 15 -0.2 0.2])
xlabel('Angular Deflection (º)')
ylabel('Right Leg Lateral Movement (BL)')
title('All conditions')

%% Step Length
% plot the movement in the Y direction (length) for the pair of legs
figure,
hold on
% plot([VrCents(1) VrCents(end)], [0 0], '--g')
plot([0 0], [0 1], '--m')
plot(VrCents, gmLeftLegStepLen, 'k', 'linewidth', 2)
plot(VrCents, gmLeftLegStepLen+semLeftLegStepLen, 'k', 'linewidth', 1)
plot(VrCents, gmLeftLegStepLen-semLeftLegStepLen, 'k', 'linewidth', 1)
axis([-15 15 0 1])
xlabel('Angular Deflection (º)')
ylabel('Left Leg Step Length (BL)')
title('All conditions')

figure,
hold on
% plot([VrCents(1) VrCents(end)], [0 0], '--g')
plot([0 0], [0 1], '--m')
plot(VrCents, gmRightLegStepLen, 'k', 'linewidth', 2)
plot(VrCents, gmRightLegStepLen+semRightLegStepLenv, 'k', 'linewidth', 1)
plot(VrCents, gmRightLegStepLen-semRightLegStepLenv, 'k', 'linewidth', 1)
axis([-15 15 0 1])
xlabel('Angular Deflection (º)')
ylabel('Right Leg Step Length (BL)')
title('All conditions')