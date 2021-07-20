clear
clc
% load the low and high resolution data relative to the forward segments
addpath('Leg Movement')
params = GetParams();
path = '\';
[LegMovForwSeg] = GetLegMovForwSeg(path, params);
%%
clc
% binning windows for the angular displacements
DevCents = [-1000 -3.75 -2.5 -1.25 0 1.25 2.5 3.75 1000];
thr = 0.02;
smth = 7;
% initialize variables to store data based on the angular displacement bin
MatAll = cell(length(DevCents)-1,15);
FLLX = cell(length(DevCents)-1,1);
FLLY = cell(length(DevCents)-1,1);  
FRLX = cell(length(DevCents)-1,1);
FRLY = cell(length(DevCents)-1,1);
MLLX = cell(length(DevCents)-1,1);
MLLY = cell(length(DevCents)-1,1);
MRLX = cell(length(DevCents)-1,1);
MRLY = cell(length(DevCents)-1,1);
HLLX = cell(length(DevCents)-1,1);
HLLY = cell(length(DevCents)-1,1);
HRLX = cell(length(DevCents)-1,1);
HRLY = cell(length(DevCents)-1,1);
SS = cell(length(DevCents)-1,1);
VR = cell(length(DevCents)-1,1);
VF = cell(length(DevCents)-1,1);
VS = cell(length(DevCents)-1,1);
NB = cell(length(DevCents)-1,1);
% iterate across angular displacement bins 
for jk = 1 : length(DevCents)-1
    minAngDev = DevCents(jk);
    maxAngDev = DevCents(jk+1);
    % iterate across the NG trial types
    for l = 1 : 2 : size(LegMovForwSeg,1)
        % iterate across flies
        for n = 1 : size(LegMovForwSeg,2)
            % iterate across forward segments
            for i = 1 : length(LegMovForwSeg{l,n})
                if ~isempty(LegMovForwSeg{l,n}{i})
                    % load low resolution data for this forward segment
                    vr = LegMovForwSeg{l,n}{i}.VrLR;
                    vf = LegMovForwSeg{l,n}{i}.VfLR;
                    vs = LegMovForwSeg{l,n}{i}.VsLR;
                    vt = sqrt(vf.*vf + vs.*vs);
                    angDev = mean(vr./vt);
                    % if the angular displacements belong to the bin
                    if angDev > minAngDev && angDev < maxAngDev
                        % store data appropriately
                        VR{jk} = vertcat(VR{jk}, vr);
                        VF{jk} = vertcat(VF{jk}, vf);
                        VS{jk} = vertcat(VS{jk}, vs);
                        NB{jk} = vertcat(NB{jk}, length(vs));
                        
                        % threshold high resolution data based on the
                        % tracking reliability
                        fflx = LegMovForwSeg{l,n}{i}.FLX;
                        ffly = LegMovForwSeg{l,n}{i}.FLY;
                        fflerr = LegMovForwSeg{l,n}{i}.FLErr;
                        fflx(fflerr > thr) = nan;
                        ffly(fflerr > thr) = nan;
                        FLLX{jk} = vertcat(FLLX{jk}, fflx);
                        FLLY{jk} = vertcat(FLLY{jk}, ffly);

                        ffrx = LegMovForwSeg{l,n}{i}.FRX;
                        ffry = LegMovForwSeg{l,n}{i}.FRY;
                        ffrerr = LegMovForwSeg{l,n}{i}.FRErr;
                        ffrx(ffrerr > thr) = nan;
                        ffry(ffrerr > thr) = nan;
                        FRLX{jk} = vertcat(FRLX{jk}, ffrx);
                        FRLY{jk} = vertcat(FRLY{jk}, ffry);

                        mmlx = LegMovForwSeg{l,n}{i}.MLX;
                        mmly = LegMovForwSeg{l,n}{i}.MLY;
                        mmlerr = LegMovForwSeg{l,n}{i}.MLErr;
                        mmlx(mmlerr > thr) = nan;
                        mmly(mmlerr > thr) = nan;
                        MLLX{jk} = vertcat(MLLX{jk}, mmlx);
                        MLLY{jk} = vertcat(MLLY{jk}, mmly);
    
                        mmrx = LegMovForwSeg{l,n}{i}.MRX;
                        mmry = LegMovForwSeg{l,n}{i}.MRY;
                        mmrerr = LegMovForwSeg{l,n}{i}.MRErr;
                        mmrx(mmrerr > thr) = nan;
                        mmry(mmrerr > thr) = nan;
                        MRLX{jk} = vertcat(MRLX{jk}, mmrx);
                        MRLY{jk} = vertcat(MRLY{jk}, mmry);

                        hhlx = LegMovForwSeg{l,n}{i}.HLX;
                        hhly = LegMovForwSeg{l,n}{i}.HLY;
                        hhlerr = LegMovForwSeg{l,n}{i}.HLErr;
                        hhlx(hhlerr > thr) = nan;
                        hhly(hhlerr > thr) = nan;
                        HLLX{jk} = vertcat(HLLX{jk}, hhlx);
                        HLLY{jk} = vertcat(HLLY{jk}, hhly);
                        
                        hhrx = LegMovForwSeg{l,n}{i}.HRX;
                        hhry = LegMovForwSeg{l,n}{i}.HRY;
                        hhrerr = LegMovForwSeg{l,n}{i}.HRErr;
                        hhrx(hhrerr > thr) = nan;
                        hhry(hhrerr > thr) = nan;
                        HRLX{jk} = vertcat(HRLX{jk}, hhrx);
                        HRLY{jk} = vertcat(HRLY{jk}, hhry);
                        
                        % get the swing stance pattern for each leg
                        swingstance = [];
                        naninds = [];
                        naninds = unique(vertcat(naninds, find(isnan(ffly))));
                        naninds = unique(vertcat(naninds, find(isnan(ffry))));
                        naninds = unique(vertcat(naninds, find(isnan(mmly))));
                        naninds = unique(vertcat(naninds, find(isnan(mmry))));
                        naninds = unique(vertcat(naninds, find(isnan(hhly))));
                        naninds = unique(vertcat(naninds, find(isnan(hhry))));

                        vleg = diff(smooth(ffly, smth/length(ffly), 'lowess'));
                        ss = ones(length(vleg)+1,1);
                        ss(vleg>-2) = 0;
                        ss(naninds) = [];
                        swingstance = horzcat(swingstance,ss);
                        
                        vleg = diff(smooth(ffry, smth/length(ffry), 'lowess'));
                        ss =ones(length(vleg)+1,1);
                        ss(vleg>-2) = 0;
                        ss(naninds) = [];
                        swingstance = horzcat(swingstance,ss);
                        
                        vleg = diff(smooth(mmly, smth/length(mmly), 'lowess'));
                        ss =ones(length(vleg)+1,1);
                        ss(vleg>-2) = 0;
                        ss(naninds) = [];
                        swingstance = horzcat(swingstance,ss);
                        
                        vleg = diff(smooth(mmry, smth/length(mmry), 'lowess'));
                        ss =ones(length(vleg)+1,1);
                        ss(vleg>-2) = 0;
                        ss(naninds) = [];
                        swingstance = horzcat(swingstance,ss);
                        
                        vleg = diff(smooth(hhly, smth/length(hhly), 'lowess'));
                        ss =ones(length(vleg)+1,1);
                        ss(vleg>-2) = 0;
                        ss(naninds) = [];
                        swingstance = horzcat(swingstance,ss);
                        
                        vleg = diff(smooth(hhry, smth/length(hhry), 'lowess'));
                        ss =ones(length(vleg)+1,1);
                        ss(vleg>-2) = 0;
                        ss(naninds) = [];
                        swingstance = horzcat(swingstance,ss);
                        
                        SS{jk} = vertcat(SS{jk}, swingstance);
                        
                        % calculate the phase difference between every pair of
                        legs
                        MatAll{jk,1} = vertcat(MatAll{jk,1}, ...
                            GetPhaseSignal(LegMovForwSeg{l,n}{i}.FLY, LegMovForwSeg{l,n}{i}.FRY));
                        MatAll{jk,2} = vertcat(MatAll{jk,2}, ...
                            GetPhaseSignal(LegMovForwSeg{l,n}{i}.FLY, LegMovForwSeg{l,n}{i}.MLY));
                        MatAll{jk,3} = vertcat(MatAll{jk,3}, ...
                            GetPhaseSignal(LegMovForwSeg{l,n}{i}.FLY, LegMovForwSeg{l,n}{i}.MRY));
                        MatAll{jk,4} = vertcat(MatAll{jk,4}, ...
                            GetPhaseSignal(LegMovForwSeg{l,n}{i}.FLY, LegMovForwSeg{l,n}{i}.HLY));
                        MatAll{jk,5} = vertcat(MatAll{jk,5}, ...
                            GetPhaseSignal(LegMovForwSeg{l,n}{i}.FLY, LegMovForwSeg{l,n}{i}.HRY));
                        MatAll{jk,6} = vertcat(MatAll{jk,6}, ...
                            GetPhaseSignal(LegMovForwSeg{l,n}{i}.FRY, LegMovForwSeg{l,n}{i}.MLY));
                        MatAll{jk,7} = vertcat(MatAll{jk,7}, ...
                            GetPhaseSignal(LegMovForwSeg{l,n}{i}.FRY, LegMovForwSeg{l,n}{i}.MRY));
                        MatAll{jk,8} = vertcat(MatAll{jk,8}, ...
                            GetPhaseSignal(LegMovForwSeg{l,n}{i}.FRY, LegMovForwSeg{l,n}{i}.HLY));
                        MatAll{jk,9} = vertcat(MatAll{jk,9}, ...
                            GetPhaseSignal(LegMovForwSeg{l,n}{i}.FRY, LegMovForwSeg{l,n}{i}.HRY));
                        MatAll{jk,10} = vertcat(MatAll{jk,10}, ...
                            GetPhaseSignal(LegMovForwSeg{l,n}{i}.MLY, LegMovForwSeg{l,n}{i}.MRY));
                        MatAll{jk,11} = vertcat(MatAll{jk,11}, ...
                            GetPhaseSignal(LegMovForwSeg{l,n}{i}.MLY, LegMovForwSeg{l,n}{i}.HLY));
                        MatAll{jk,12} = vertcat(MatAll{jk,12}, ...
                            GetPhaseSignal(LegMovForwSeg{l,n}{i}.MLY, LegMovForwSeg{l,n}{i}.HRY));
                        MatAll{jk,13} = vertcat(MatAll{jk,13}, ...
                            GetPhaseSignal(LegMovForwSeg{l,n}{i}.MRY, LegMovForwSeg{l,n}{i}.HLY));
                        MatAll{jk,14} = vertcat(MatAll{jk,14}, ...
                            GetPhaseSignal(LegMovForwSeg{l,n}{i}.MRY, LegMovForwSeg{l,n}{i}.HRY));
                        MatAll{jk,15} = vertcat(MatAll{jk,14}, ...
                            GetPhaseSignal(LegMovForwSeg{l,n}{i}.HLY, LegMovForwSeg{l,n}{i}.HRY));

                    end
                end
            end
        end
    end
end

% calculate number of legs in air
figure,
hold on
for jk = 1 : length(DevCents)-1
    bar(jk, nanmean(nansum(SS{jk},2)))
    errorbar(jk,  nanmean(nansum(SS{jk},2)), nanstd(nansum(SS{jk},2)))
end
axis([0 9 0 3.5])

%%
% Plot centered trajectories
xT = cell(1,1);
yT = cell(1,1);
indd = cell(1,1);
a = 1;
cmap1 = jet(length(VR));
% iterate through the angular deviation bins
for i = 2 : length(NB)-1
   nbs = NB{i};
   nbScs = [1; cumsum(nbs)];
   % iterate through the forward segmtnes
   for j = 1 : length(nbs)
       vr = VR{i}(nbScs(j):nbScs(j+1))/60;
       vf = VF{i}(nbScs(j):nbScs(j+1))/60;
       vs = VS{i}(nbScs(j):nbScs(j+1))/60;
       x = [];
       y = [];
       th = [];
       if length(vr) < 180
       % get centered X and Y
       for ij = 1 : length(vr)
           if isempty(x)
               th = vertcat(th, pi*vr(ij)/180);
               x = vertcat(x, vf(ij)*sin(0)+vs(ij)*cos(0));
               y = vertcat(y, vf(ij)*cos(0)-vs(ij)*sin(0));
           else
               th = vertcat(th, th(end) + pi*vr(ij)/180);
               x = vertcat(x, x(end) + vf(ij)*sin(th(ij-1))+vs(ij)*cos(th(ij-1)));
               y = vertcat(y, y(end) + vf(ij)*cos(th(ij-1))-vs(ij)*sin(th(ij-1)));
           end
       end
       end
       xT{a} = x;
       yT{a} = y;
       indd{a} = i;
       a = a + 1;
   end
end
% plot a color coded sample of the forward segments
inds = randperm(length(xT));
inds = inds(1:500);
figure,
hold on
for i = 1 : 500
    plot(yT{inds(i)}, xT{inds(i)}, 'color', cmap1(indd{inds(i)},:))
end
axis([-2 30 -25 25])

%%
% calculate leg phases
% pair of leg phase comparisons
pT = {'FL-FR', 'FL-ML', 'FL-MR', 'FL-HL', 'FL-HR', 'FR-ML', 'FR-MR', 'FR-HL', 'FR-HR' ...
    , 'ML-MR', 'ML-HL', 'ML-HR', 'MR-HL', 'MR-HR', 'HL-HR'};
% vector signaling legs that belong to same triangle
TF = [0 0 1 1 0 1 0 0 1 0 0 1 1 0 0];
Cent = 9;
figure,
hold on
for i = 1 : 15
%     subplot(3,5,i)
    [tout,rout] = rose(pi*MatAll{Cent,i}/180+pi/6, (-pi:pi/12:pi));
    tout = tout(1:end-1);
    rout = rout(1:end-1);
    inds = find(tout==0);
    tout(inds) = [];
    rout(inds) = [];
    tout = tout(1:2:end);
    rout = rout(1:2:end)/sum(rout);
    if TF(i) > 0
        polar(tout,smooth(rout,3)','b')
    else
        polar(tout,smooth(rout,3)','r')
    end
    title(pT{i})
end
axis square
axis([-0.15 0.15 -0.15 0.15])
%%
% plot landing distributions
% define landing position bins
centsX = 0 : 3 : 400;
centsY = 0 : 3 : 400;
aB = 1;
cmap = jet(length(DevCents)-1);
cmap2 = flipud(gray(50));

% parameters to define swing stance
maxStepSize = 130;
minSwingSize = 2;
maxSwingSize = 15;

figure,
for i = 1 : length(DevCents)-1
    % look for trajectories local minima and maxima to detect swing stance
    % transitions
    [~, indSwingIFL] = findpeaks(FLLY{i},'MinPeakProminence',40);
    [~, indSwingFFL] = findpeaks(-FLLY{i},'MinPeakProminence',40);
    stepLenghtFL = [];
    landingXFL = [];
    landingYFL = [];
    liftoffXFL = [];
    liftoffYFL = [];
    % iterate for all the local maxima
    for j = 1 : length(indSwingIFL)
        indF = find(indSwingFFL>indSwingIFL(j),1);
        if ~isempty(indF)
            swingTime = indSwingFFL(indF)-indSwingIFL(j);
            stepSize = FLLY{i}(indSwingFFL(indF)) - FLLY{i}(indSwingIFL(j));
            % if swing time and step size are within the parameters 
            if swingTime < maxSwingSize && swingTime > minSwingSize  && stepSize > - maxStepSize
                stepLenghtFL = vertcat(stepLenghtFL, stepSize);
                landingXFL = vertcat(landingXFL, FLLX{i}(indSwingFFL(indF)));
                landingYFL = vertcat(landingYFL, 400-FLLY{i}(indSwingFFL(indF)));
                liftoffXFL = vertcat(liftoffXFL, FLLX{i}(indSwingIFL(j)));
                liftoffYFL = vertcat(liftoffYFL, 400-FLLY{i}(indSwingIFL(j)));
            end
        end
    end
    stepLenghtFL = -stepLenghtFL;
    
    % look for trajectories local minima and maxima to detect swing stance
    % transitions
    [~, indSwingIFR] = findpeaks(FRLY{i},'MinPeakProminence',40);
    [~, indSwingFFR] = findpeaks(-FRLY{i},'MinPeakProminence',40);
    stepLenghtFR = [];
    landingXFR = [];
    landingYFR = [];
    liftoffXFR = [];
    liftoffYFR = [];
    % iterate for all the local maxima
    for j = 1 : length(indSwingIFR)
        indF = find(indSwingFFR>indSwingIFR(j),1);
        if ~isempty(indF)
            swingTime = indSwingFFR(indF)-indSwingIFR(j);
            stepSize = FRLY{i}(indSwingFFR(indF)) - FRLY{i}(indSwingIFR(j));
            % if swing time and step size are within the parameters 
            if swingTime < maxSwingSize && swingTime > minSwingSize  && stepSize > - maxStepSize
                stepLenghtFR = vertcat(stepLenghtFR, stepSize);
                landingXFR = vertcat(landingXFR, FRLX{i}(indSwingFFR(indF)));
                landingYFR = vertcat(landingYFR, 400-FRLY{i}(indSwingFFR(indF)));
                liftoffXFR = vertcat(liftoffXFR, FRLX{i}(indSwingIFR(j)));
                liftoffYFR = vertcat(liftoffYFR, 400-FRLY{i}(indSwingIFR(j)));
            end
        end
    end
    stepLenghtFR = -stepLenghtFR;
    
    % look for trajectories local minima and maxima to detect swing stance
    % transitions
	[~, indSwingIML] = findpeaks(MLLY{i},'MinPeakProminence',40);
    [~, indSwingFML] = findpeaks(-MLLY{i},'MinPeakProminence',40);
    stepLenghtML = [];
    landingXML = [];
    landingYML = [];
    liftoffXML = [];
    liftoffYML = [];
    % iterate for all the local maxima
    for j = 1 : length(indSwingIML)
        indF = find(indSwingFML>indSwingIML(j),1);
        if ~isempty(indF)
            swingTime = indSwingFML(indF)-indSwingIML(j);
            stepSize = MLLY{i}(indSwingFML(indF)) - MLLY{i}(indSwingIML(j));
            % if swing time and step size are within the parameters
            if swingTime < maxSwingSize && swingTime > minSwingSize  && stepSize > - maxStepSize
                stepLenghtML = vertcat(stepLenghtML, stepSize);
                landingXML = vertcat(landingXML, MLLX{i}(indSwingFML(indF)));
                landingYML = vertcat(landingYML, 400-MLLY{i}(indSwingFML(indF)));
                liftoffXML = vertcat(liftoffXML, MLLX{i}(indSwingIML(j)));
                liftoffYML = vertcat(liftoffYML, 400-MLLY{i}(indSwingIML(j)));
            end
        end
    end
    stepLenghtML = -stepLenghtML;
    
    % look for trajectories local minima and maxima to detect swing stance
    % transitions
	[~, indSwingIMR] = findpeaks(MRLY{i},'MinPeakProminence',40);
    [~, indSwingFMR] = findpeaks(-MRLY{i},'MinPeakProminence',40);
    stepLenghtMR = [];
    landingXMR = [];
    landingYMR = [];
    liftoffXMR = [];
    liftoffYMR = [];
    % iterate for all the local maxima
    for j = 1 : length(indSwingIMR)
        indF = find(indSwingFMR>indSwingIMR(j),1);
        if ~isempty(indF)
            swingTime = indSwingFMR(indF)-indSwingIMR(j);
            stepSize = MRLY{i}(indSwingFMR(indF)) - MRLY{i}(indSwingIMR(j));
            % if swing time and step size are within the parameters
            if swingTime < maxSwingSize && swingTime > minSwingSize  && stepSize > - maxStepSize
                stepLenghtMR = vertcat(stepLenghtMR, stepSize);
                landingXMR = vertcat(landingXMR, MRLX{i}(indSwingFMR(indF)));
                landingYMR = vertcat(landingYMR, 400-MRLY{i}(indSwingFMR(indF)));
                liftoffXMR = vertcat(liftoffXMR, MRLX{i}(indSwingIMR(j)));
                liftoffYMR = vertcat(liftoffYMR, 400-MRLY{i}(indSwingIMR(j)));
            end
        end
    end
    stepLenghtMR = -stepLenghtMR;    
    
    % look for trajectories local minima and maxima to detect swing stance
    % transitions
    [~, indSwingIHL] = findpeaks(HLLY{i},'MinPeakProminence',40);
    [~, indSwingFHL] = findpeaks(-HLLY{i},'MinPeakProminence',40);
    stepLenghtHL = [];
    landingXHL = [];
    landingYHL = [];
    liftoffXHL = [];
    liftoffYHL = [];
    % iterate for all the local maxima
    for j = 1 : length(indSwingIHL)
        indF = find(indSwingFHL>indSwingIHL(j),1);
        if ~isempty(indF)
            swingTime = indSwingFHL(indF)-indSwingIHL(j);
            stepSize = HLLY{i}(indSwingFHL(indF)) - HLLY{i}(indSwingIHL(j));
            % if swing time and step size are within the parameters
            if swingTime < maxSwingSize && swingTime > minSwingSize  && stepSize > - maxStepSize
                stepLenghtHL = vertcat(stepLenghtHL, stepSize);
                landingXHL = vertcat(landingXHL, HLLX{i}(indSwingFHL(indF)));
                landingYHL = vertcat(landingYHL, 400-HLLY{i}(indSwingFHL(indF)));
                liftoffXHL = vertcat(liftoffXHL, HLLX{i}(indSwingIHL(j)));
                liftoffYHL = vertcat(liftoffYHL, 400-HLLY{i}(indSwingIHL(j)));
            end
        end
    end
    stepLenghtHL = -stepLenghtHL;  
    
    % look for trajectories local minima and maxima to detect swing stance
    % transitions
    [~, indSwingIHR] = findpeaks(HRLY{i},'MinPeakProminence',40);
    [~, indSwingFHR] = findpeaks(-HRLY{i},'MinPeakProminence',40);
    stepLenghtHR = [];
    landingXHR = [];
    landingYHR = [];
    liftoffXHR = [];
    liftoffYHR = [];
    % iterate for all the local maxima
    for j = 1 : length(indSwingIHR)
        indF = find(indSwingFHR>indSwingIHR(j),1);
        if ~isempty(indF)
            swingTime = indSwingFHR(indF)-indSwingIHR(j);
            stepSize = HRLY{i}(indSwingFHR(indF)) - HRLY{i}(indSwingIHR(j));
            % if swing time and step size are within the parameters
            if swingTime < maxSwingSize && swingTime > minSwingSize  && stepSize > - maxStepSize
                stepLenghtHR = vertcat(stepLenghtHR, stepSize);
                landingXHR = vertcat(landingXHR, HRLX{i}(indSwingFHR(indF)));
                landingYHR = vertcat(landingYHR, 400-HRLY{i}(indSwingFHR(indF)));
                liftoffXHR = vertcat(liftoffXHR, HRLX{i}(indSwingIHR(j)));
                liftoffYHR = vertcat(liftoffYHR, 400-HRLY{i}(indSwingIHR(j)));
            end
        end
    end
    stepLenghtHR = -stepLenghtHR;
    
    % Get landing distributions for the X direction
    landingXFL(isnan(landingXFL)) = [];
    landingXFR(isnan(landingXFR)) = [];
    landingXML(isnan(landingXML)) = [];
    landingXMR(isnan(landingXMR)) = [];
    landingXHL(isnan(landingXHL)) = [];
    landingXHR(isnan(landingXHR)) = [];
    hhflx = smooth(hist(landingXFL, centsX)/length(landingXFL), 4/length(centsX), 'lowess');
    hhfrx = smooth(hist(landingXFR, centsX)/length(landingXFR), 4/length(centsX), 'lowess');
    hhmlx = smooth(hist(landingXML, centsX)/length(landingXML), 4/length(centsX), 'lowess');
    hhmrx = smooth(hist(landingXMR, centsX)/length(landingXMR), 4/length(centsX), 'lowess');
    hhhlx = smooth(hist(landingXHL, centsX)/length(landingXHL), 4/length(centsX), 'lowess');
    hhhrx = smooth(hist(landingXHR, centsX)/length(landingXHR), 4/length(centsX), 'lowess');
    
    % Get landing distributions for the Y direction
    landingYFL(isnan(landingYFL)) = [];
    landingYFR(isnan(landingYFR)) = [];
    landingYML(isnan(landingYML)) = [];
    landingYMR(isnan(landingYMR)) = [];
    landingYHL(isnan(landingYHL)) = [];
    landingYHR(isnan(landingYHR)) = [];
    hhfly = smooth(hist(landingYFL, centsY)/length(landingYFL), 4/length(centsY), 'lowess');
    hhfry = smooth(hist(landingYFR, centsY)/length(landingYFR), 4/length(centsY), 'lowess');
    hhmly = smooth(hist(landingYML, centsY)/length(landingYML), 4/length(centsY), 'lowess');
    hhmry = smooth(hist(landingYMR, centsY)/length(landingYMR), 4/length(centsY), 'lowess');
    hhhly = smooth(hist(landingYHL, centsY)/length(landingYHL), 4/length(centsY), 'lowess');
    hhhry = smooth(hist(landingYHR, centsY)/length(landingYHR), 4/length(centsY), 'lowess');

%     % calculate 2D histogram for leg landing position
    mapLanding = Occupancy(...
        vertcat(landingXFL, landingXFR, landingXML, landingXMR, landingXHL, landingXHR),...
        vertcat(landingYFL, landingYFR, landingYML, landingYMR, landingYHL, landingYHR),...
        centsX, centsY);
    mapLanding = flipud(mapLanding'./sum(sum(mapLanding)));   

%  calculate 2D histogram for leg liftoff position
%     mapLanding = Occupancy(...
%         vertcat(liftoffXFL, liftoffXFR, liftoffXML, liftoffXMR, liftoffXHL, liftoffXHR),...
%         vertcat(liftoffYFL, liftoffYFR, liftoffYML, liftoffYMR, liftoffYHL, liftoffYHR),...
%         centsX, centsY);
%     mapLanding = flipud(mapLanding'./sum(sum(mapLanding)));

    % plot 2D histogram
    ax = subplot(2, ceil((length(DevCents)-1)/2), i);
    hold on
    cmap3 = horzcat(cmap2(:,1).*cmap(i,1), cmap2(:,2).*cmap(i,2), cmap2(:,3).*cmap(i,3));
    cmap3(1,:) = [1 1 1];
    caxis([0 0.002])
    [c,h]=contourf(centsX,centsY,imgaussfilt(mapLanding,1), [0.0001:0.0005:0.002]);
    set(h,'LineColor','none')
    title(['[' num2str(DevCents(i)) ' ' num2str(DevCents(i+1)) ']'])
    colormap(ax, cmap3)
    title(['[' num2str(DevCents(i)) ' ' num2str(DevCents(i+1)) ']'])
end

