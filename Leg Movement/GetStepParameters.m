function stepPars = GetStepParameters(Leg1Y, Leg2Y, Leg1X, Leg2X, VR, VF, Time)
maxStepSize = 130;
minSwingSize = 2;
maxSwingSize = 15;
minStanceSize = 2;
maxStanceSize = 14;
peakProminence = 40;

[~, indSwingL1] = findpeaks(Leg1Y, 'MinPeakProminence', peakProminence);
[~, indStanceL1] = findpeaks(-Leg1Y, 'MinPeakProminence', peakProminence);
[~, indSwingL2] = findpeaks(Leg2Y, 'MinPeakProminence', peakProminence);
[~, indStanceL2] = findpeaks(-Leg2Y, 'MinPeakProminence', peakProminence);

stepPars.StepDurationL1 = [];
stepPars.StepDurationL2 = [];
stepPars.MVRL1 = [];
stepPars.MVFL1 = [];
stepPars.MVRL2 = [];
stepPars.MVFL2 = [];
stepPars.L1SwingTime = [];
stepPars.L2SwingTime = [];
stepPars.L1StanceTime = [];
stepPars.L2StanceTime = [];
stepPars.L1SwingLateral = [];
stepPars.L2SwingLateral = [];
stepPars.L1SwingLength = [];
stepPars.L2SwingLength = [];
stepPars.L1StanceLateral = [];
stepPars.L2StanceLateral = [];
stepPars.L1StanceLength = [];
stepPars.L2StanceLength = [];
stepPars.StepLag = [];
stepPars.StepTime = [];

nst1 = length(indStanceL1)-2;
% Extract swing-stance-swing frames for both legs
for j = 1 : length(indStanceL1)
    % Get frame Leg1 stance
    indL1Stance = indStanceL1(j);
    % Get frame Leg2 stance
    indL2Stance = find(indStanceL2>=indL1Stance-3, 1);
    
    if ~isempty(indL2Stance)
        indL2Stance = indStanceL2(indL2Stance);
        indL1PrevSwing = find(indSwingL1<indL1Stance);
        % Get frame L1 PREVIOUS swing
        if ~isempty(indL1PrevSwing)
            indL1PrevSwing = indSwingL1(indL1PrevSwing(end));
            % Get frame L2 previous swing
            indL2PrevSwing = find(indSwingL2<indL2Stance);
            if ~isempty(indL2PrevSwing)
                indL2PrevSwing = indSwingL2(indL2PrevSwing(end));
                
                swingTimeL1 = indL1Stance - indL1PrevSwing;
                swingTimeL2 = indL2Stance - indL2PrevSwing;
                
                stepLag = indL2Stance - indL1Stance;
                
                swingLengthL1 = -(Leg1Y(indL1Stance) - Leg1Y(indL1PrevSwing));
                swingLengthL2 = -(Leg2Y(indL2Stance) - Leg2Y(indL2PrevSwing));
                swingLateralL1 = (Leg1X(indL1Stance) - Leg1X(indL1PrevSwing));
                swingLateralL2 = (Leg2X(indL2Stance) - Leg2X(indL2PrevSwing));
                
                % Get frame ipsi next swing
                indL1NextSwing = find(indSwingL1>indL1Stance, 1);
                if  ~isempty(indL1NextSwing)
                    indL1NextSwing = indSwingL1(indL1NextSwing);
                    % Get frame contra next swing
                    indL2NextSwing = find(indSwingL2>indL2Stance, 1);
                    if  ~isempty(indL2NextSwing)
                        indL2NextSwing = indSwingL2(indL2NextSwing);
                        
                        stanceTimeL1 = indL1NextSwing - indL1Stance;
                        stanceTimeL2 = indL2NextSwing - indL2Stance;

                        stanceLengthL1 = -(Leg1Y(indL1NextSwing) - Leg1Y(indL1Stance));
                        stanceLengthL2 = -(Leg2Y(indL2NextSwing) - Leg2Y(indL2Stance));
                        stanceLateralL1 = (Leg1X(indL1NextSwing) - Leg1X(indL1Stance));
                        stanceLateralL2 = (Leg2X(indL2NextSwing) - Leg2X(indL2Stance));
                        
                        % If all these frames are within non error margins
                        if stanceTimeL1 < maxStanceSize && stanceTimeL1 > minStanceSize  && ...
                                stanceTimeL2 < maxStanceSize && stanceTimeL2 > minStanceSize  && ...
                                swingTimeL1 < maxSwingSize && swingTimeL1 > minSwingSize  && ...
                                swingTimeL2 < maxSwingSize && swingTimeL2 > minSwingSize && ...
                                swingLengthL1 < maxStepSize && swingLengthL2 < maxStepSize
                            
                            stepPars.StepLag = vertcat(stepPars.StepLag, stepLag/120);
                            stepPars.StepDurationL1 = vertcat(stepPars.StepDurationL1, (indL1NextSwing-indL1PrevSwing)/120); % Vr step L1
                            stepPars.StepDurationL2 = vertcat(stepPars.StepDurationL2, (indL2NextSwing-indL2PrevSwing)/120); % Vr step L1
                            stepPars.MVRL1 = vertcat(stepPars.MVRL1, mean(VR(indL1PrevSwing:indL1NextSwing))*stepPars.StepDurationL1(end)); % Vr step L1
                            stepPars.MVFL1 = vertcat(stepPars.MVFL1, mean(VF(indL1PrevSwing:indL1NextSwing))*stepPars.StepDurationL1(end)); % Vf step L1
                            stepPars.MVRL2 = vertcat(stepPars.MVRL2, mean(VR(indL2PrevSwing:indL2NextSwing))*stepPars.StepDurationL2(end)); % Vr step L2
                            stepPars.MVFL2 = vertcat(stepPars.MVFL2, mean(VF(indL2PrevSwing:indL2NextSwing))*stepPars.StepDurationL2(end)); % Vf step L2
                            stepPars.StepTime = vertcat(stepPars.StepTime, Time(indL1Stance));
                            
                            stepPars.L1SwingTime = vertcat(stepPars.L1SwingTime, swingTimeL1/120); % Duration Swing L1
                            stepPars.L2SwingTime = vertcat(stepPars.L2SwingTime, swingTimeL2/120); % Duration Swing L2
                            stepPars.L1StanceTime = vertcat(stepPars.L1StanceTime, stanceTimeL1/120); % Duration Stance L1
                            stepPars.L2StanceTime = vertcat(stepPars.L2StanceTime, stanceTimeL2/120); % Duration Stance L2
                            
                            stepPars.L1SwingLateral = vertcat(stepPars.L1SwingLateral, swingLateralL1/180); % Lateral Mov L1 (swing)
                            stepPars.L1SwingLength = vertcat(stepPars.L1SwingLength, swingLengthL1/180); % Length Mov L1 (swing)
                            stepPars.L2SwingLateral = vertcat(stepPars.L2SwingLateral, swingLateralL2/180); % Lateral Mov L2 (swing)
                            stepPars.L2SwingLength = vertcat(stepPars.L2SwingLength, swingLengthL2/180); % Length Mov L2 (swing)
                            
                            stepPars.L1StanceLateral = vertcat(stepPars.L1StanceLateral, stanceLateralL1/180); % Lateral Mov L1 (swing)
                            stepPars.L1StanceLength = vertcat(stepPars.L1StanceLength, stanceLengthL1/180); % Length Mov L1 (swing)
                            stepPars.L2StanceLateral = vertcat(stepPars.L2StanceLateral, stanceLateralL2/180); % Lateral Mov L2 (swing)
                            stepPars.L2StanceLength = vertcat(stepPars.L2StanceLength, stanceLengthL2/180); % Length Mov L2 (swing)
                        end
                    end
                end
            end
        end
    end
end
stepPars.PStepQ = length(stepPars.StepLag)/nst1;
