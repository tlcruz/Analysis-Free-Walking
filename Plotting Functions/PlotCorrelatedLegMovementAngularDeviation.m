%% Plot correlation between angular deviation and correlated leg placement

clear
clc
path = 'D:\Dropbox\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\LegTracking\Dark\';
% load leg movements during forward runs 
params = GetParams();
[LegMovForwSegD, pTypesD] = GetLegMovForwSeg(path, params);
thr = 0.02;
vrcnt = [-5000 5000];

% Individual Leg Step Spatial Parameters
LLp1 = cell(size(LegMovForwSegD, 1),size(LegMovForwSegD, 2));
VRLp1 = cell(size(LegMovForwSegD, 1),size(LegMovForwSegD, 2));
StepTs = cell(size(LegMovForwSegD, 1),size(LegMovForwSegD, 2));
for j = 1 : 1 : size(LegMovForwSegD, 1)
    llp = [];
    nllp = [];
    for n = 1 : size(LegMovForwSegD, 2)
        n2 = n;
        [lD] = GetCleanLegData(LegMovForwSegD{j,n}, thr);
        % vectors to put all step data
        steplatl1 = [];
        steplatl2 = [];
        stepVr = [];
        steplat = [];
        stepT = [];
        for i = 1 : length(lD.VR)
            if ~isempty(lD.FLLY{i}) && ~isempty(lD.FRLY{i})
                stepPars1 = GetStepParameters(lD.FLLY{i}, lD.FRLY{i}, lD.FLLX{i}, lD.FRLX{i}, lD.VR{i}, lD.VF{i}, lD.Time{i});
                for k = 1 : length(stepPars1.MVRL1)
                    if  sign(stepPars1.L1SwingLateral(k)) > 0
                        LLp1{n2} = vertcat(LLp1{n2}, stepPars1.L1SwingLateral(k)*stepPars1.L2SwingLateral(k));
                        StepTs{n2} = vertcat(StepTs{n2}, stepPars1.StepTime(k));
                        VRLp1{n2} = vertcat(VRLp1{n2}, -(stepPars1.MVRL1(k)+stepPars1.MVRL2(k))/2);
                    end
                end
                stepPars2 = GetStepParameters(lD.FRLY{i}, lD.FLLY{i}, lD.FRLX{i}, lD.FLLX{i}, lD.VR{i}, lD.VF{i},lD.Time{i});
                for k = 1 : length(stepPars2.MVRL1)
                    if  sign(stepPars2.L1SwingLateral(k)) > 0
                        LLp1{n2} = vertcat(LLp1{n2}, stepPars2.L1SwingLateral(k)*stepPars2.L2SwingLateral(k));
                        VRLp1{n2} = vertcat(VRLp1{n2}, -(stepPars2.MVRL1(k)+stepPars2.MVRL2(k))/2);
                        StepTs{n2} = vertcat(StepTs{n2}, stepPars2.StepTime(k));
                    end
                end
            end
        end
    end
end

% plot angular deviations vs correlated leg placement
figure,
hold on
step = 0.003;
auxCents = -0.03:step:0.06;
aux = 1;
cmap = hot(5);
kk = [];
nn = [];
for n = 1 : size(LLp1,1)
    ssp = VRLp1{n};
    probX = zeros(length(auxCents)-1, 1);
    nnX = zeros(length(auxCents)-1, 1);
    for k = 1 : length(auxCents)-1
        inds = find(LLp1{n} > auxCents(k) & LLp1{n} <= auxCents(k+1));
        if ~isempty(inds)
            probX(k) = nanmean(ssp(inds));
            nnX(k) = length(inds);
        else
            probX(k) = 0;
            nnX(k) = 0;
        end
    end
    kk = horzcat(kk, probX);
    nn = horzcat(nn, nnX);
end
gm = [];
sem = [];
for k = 1 : length(auxCents)-1
    if sum(nn(k,:)) > 50
        gm(k) = kk(k,:)*nn(k,:)'/sum(nn(k,:));
        sem(k) = sqrt(((kk(k,:)-gm(k)).*(kk(k,:)-gm(k)))*nn(k,:)'/sum(nn(k,:)))/sqrt(length(find(nn(k,:)~=0)));
    else
        gm(k) = nan;
        sem(k) = nan;
    end
end
hold on
plot([-0.05 0.05], [0 0], '--g')
plot([0 0], [-10 5], '--m')
plot(auxCents(1:(end-1)), (gm), 'color', cmap(aux,:), 'linewidth', 3)
plot(auxCents(1:(end-1)), (gm+sem), 'color', cmap(aux,:), 'linewidth', 1)
plot(auxCents(1:(end-1)), (gm-sem), 'color', cmap(aux,:), 'linewidth', 1)
title(pTypes{j})
xlabel('<Lxi*Lxi+1>')
ylabel('Signed Corrected Vr (º/s)')
axis([-0.03 0.04 -10 10])
aux = aux + 1;
