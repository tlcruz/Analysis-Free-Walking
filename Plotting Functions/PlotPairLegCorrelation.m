clc
% load the low and high resolution data relative to the forward segments
path = '\';
params = GetParams();
[LegMovForwSeg] = GetForwSegLegData(path, params);
%
thr = 0.02;
vrcnt = [0 0.5 2 4 8 1000];

dXCents = -0.3:0.02:0.3;
cmap = summer(length(vrcnt)-1);
figure,
hold on
for vrr = 1 : length(vrcnt)-1
% Individual Leg Step Spatial Parameters
LatL1p1 = [];
NLatL1p1 = [];
LatL1p2 = [];
NLatL1p2 = [];
LatL2p1 = [];
NLatL2p1 = [];
LatL2p2 = [];
NLatL2p2 = [];
FPopL1 = [];
FPopL2 = [];
% iterate across trial types
for j = 1 : 2 : size(LegMovForwSeg, 1)
    % matrices to put all flies
    latl1p1 = [];
    nlatl1p1 = [];
    latl2p1 = [];
    nlatl2p1 = [];
    latl1p2 = [];
    nlatl1p2 = [];
    latl2p2 = [];
    nlatl2p2 = [];
    % iterate across flies
    for n = 1 : size(LegMovForwSeg, 2)
        [lD] = GetCleanLegData(LegMovForwSeg{j,n}, thr);
        % vectors to put all step data
        steplatl1 = [];
        steplatl2 = [];
        stepVr = [];
        for i = 1 : length(lD.VR)
            if ~isempty(lD.FLLY{i}) && ~isempty(lD.MRLY{i})
                % get parameters between a particular pair of legs
                stepPars1 = GetStepParameters(lD.FLLY{i}, lD.FRLY{i}, lD.FLLX{i}, lD.FRLX{i}, lD.VR{i}, lD.VF{i});
                % store the movements in the X direction
                for k = 1 : length(stepPars1.MVRL1)
                    steplatl1 = vertcat(steplatl1, ...
                        stepPars1.L1SwingLateral(k));
                    steplatl2 = vertcat(steplatl2, ...
                        stepPars1.L2SwingLateral(k));
                    stepVr = vertcat(stepVr, (stepPars1.MVRL1(k)+stepPars1.MVRL2(k))/2);
                end
                % get parameters between a particular pair of legs
                stepPars2 = GetStepParameters(lD.FRLY{i}, lD.FLLY{i}, lD.FRLX{i}, lD.FLLX{i}, lD.VR{i}, lD.VF{i});
                % store the movements in the X direction
                for k = 1 : length(stepPars2.MVRL1)
                    steplatl1 = vertcat(steplatl1, ...
                        stepPars2.L1SwingLateral(k));
                    steplatl2 = vertcat(steplatl2, ...
                        stepPars2.L2SwingLateral(k));
                    stepVr = vertcat(stepVr, (stepPars2.MVRL1(k)+stepPars2.MVRL2(k))/2);
                end
            end
        end
        
        % calculate distributions
        latl1Auxp1 = nan;
        nlatl1Auxp1 = 0;
        latl1Auxp2 = nan;
        nlatl1Auxp2 = 0;
        latl2Auxp1 = nan;
        nlatl2Auxp1 = 0;
        latl2Auxp2 = nan;
        nlatl2Auxp2 = 0;

        % subtract mean position leg
        steplatl1 = steplatl1 - nanmean(steplatl1);
        steplatl2 = steplatl2 - nanmean(steplatl2);
        
        % select for a specific bin of angular deflection
        inds1 = find(abs(stepVr) > vrcnt(vrr) & abs(stepVr) < vrcnt(vrr+1));
        steplatl1 = steplatl1(inds1);
        steplatl2 = steplatl2(inds1);
        
        % iterate across leg movement bins
        for i = 1 : length(dXCents)-1
            inds = find(steplatl1 > dXCents(i) & steplatl1 < dXCents(i+1));
            if~isempty(inds)
                auxl1 = steplatl1(inds);
                auxl2 = steplatl2(inds);
                latl1Auxp1 = vertcat(latl1Auxp1, median(auxl1));
                nlatl1Auxp1 = vertcat(nlatl1Auxp1, length(inds));
                latl2Auxp1 = vertcat(latl2Auxp1, median(auxl2));
                nlatl2Auxp1 = vertcat(nlatl2Auxp1, length(inds));
                
            else
                latl1Auxp1 = vertcat(latl1Auxp1, nan);
                nlatl1Auxp1 = vertcat(nlatl1Auxp1, 0);
                latl2Auxp1 = vertcat(latl2Auxp1, nan);
                nlatl2Auxp1 = vertcat(nlatl2Auxp1, 0);
            end
        end
        latl1p1 = horzcat(latl1p1, latl1Auxp1);
        nlatl1p1 = horzcat(nlatl1p1, nlatl1Auxp1);
        latl2p1 = horzcat(latl2p1, latl2Auxp1);
        nlatl2p1 = horzcat(nlatl2p1, nlatl2Auxp1);
    end
    LatL1p1 = horzcat(LatL1p1, latl1p1);
    NLatL1p1 = horzcat(NLatL1p1, nlatl1p1);
    LatL2p1 = horzcat(LatL2p1, latl2p1);
    NLatL2p1 = horzcat(NLatL2p1, nlatl2p1);
end

% calcuate grand mean of leg movement
[gmLatL1p1, semLatL1p1] = GetGMSEMFromMatrix(LatL1p1, NLatL1p1,4);
[gmLatL2p1, semLatL2p1] = GetGMSEMFromMatrix(LatL2p1, NLatL2p1,4);
gmLatL1 = (gmLatL1p1(2:end));
semLatL1 =(semLatL1p1(2:end));
gmLatL2 = (gmLatL2p1(2:end));
semLatL2 =(semLatL2p1(2:end));

% plot leg movement
hold on
plot([-0.1 0.1], [0 0], '--g')
plot([0 0], [-0.2 0.2], '--m')
plot(dXCents(2:end), gmLatL2, 'color', cmap(vrr,:), 'linewidth', 3)
plot(dXCents(2:end), gmLatL2+semLatL2, 'color', cmap(vrr,:), 'linewidth', 1)
plot(dXCents(2:end), gmLatL2-semLatL2, 'color', cmap(vrr,:), 'linewidth', 1)
end
axis([-0.1 0.1 -0.12 0.12])