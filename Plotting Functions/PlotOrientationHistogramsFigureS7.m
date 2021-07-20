clear
clc
% Load data
[params] = GetParams();
path = '\';
flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowRes.mat'], 'Flies');
seq = dt.Flies.Seq;
pTypes = unique(seq);

% initialize variables to store orientation and angular speed across flies
% and trial types
Ang = cell(length(flies), length(pTypes));
vr = cell(length(flies), length(pTypes));
NAng = cell(length(flies), length(pTypes));
for n = 1 : length(flies)
    % load data for a specific fly
    dt = load([path flies(n).name '\DataLowRes.mat'], 'Flies');
    seq = dt.Flies.Seq;
    dt = dt.Flies.Data;
    % iterate across protocol segments
    for k = 1 : length(dt)
        % iterate across trial types
        for j = 1 : length(pTypes)
            switch seq{k}
                % when trial type matches a protocol segment
                case pTypes{j}
                    % Load walking parameters for this protocol segment
                    vrb = dt{k}.Vr;
                    actst = dt{k}.actState;
                    ang = dt{k}.Angle;
                    actst(wd < params.mDistWall) = 0;
                    
                    FWB = dt{k}.Bouts;
                    % for each walking bout extract angular speed and orientation 
                    for i = 1 : length(FWB)
                        Ang{n,j} = vertcat(Ang{n,j}, ang(FWB{i}));
                        vr{n,j} = vertcat(vr{n,j}, vrb(FWB{i}));
                        NAng{n,j} = horzcat(NAng{n,j}, length(FWB{i}));
                    end

                    disp(['Fly:' num2str(n) '  Protocol:' pTypes{j}])
            end
        end
        
    end
end

%%
% define orientation window
% angCents = 0.01:1:90.09;
angCents = 0.01:1:360.09;
thrN = 000;
figure,
sampD = [];
% iterate across NG trials
for j = 1 :2: length(pTypes)
    subplot(1,5,(j+1)/2)
    hold on
    HHA = [];
    NA = [];
    for n = 1 : length(flies)
        if sum(NAng{n,j}) > thrN
            angF = mod(Ang{n,j}+10, 360); %185
%             angF = mod(Ang{n,j}+25, 90); %185
            sampD = vertcat(sampD,  angF);
            hA = smooth(hist(angF, angCents)/sum(hist(angF, angCents)))';
            HHA = vertcat(HHA, hA);
            NA = horzcat(NA, sum(NAng{n,j}));
        end
    end
    % calculate weighted mean across flies, using walking time as weight
    gmA = HHA'*NA'/sum(NA);
    semA = zeros(1,length(gmA));
    for n = 1 : size(HHA,1)
        semA = semA + (HHA(n,:)-gmA').*(HHA(n,:)-gmA')*NA(n);
    end
    semA = sqrt(semA/sum(NA));
    
    plot(angCents, gmA, 'k', 'linewidth', 3)
    plot(angCents, gmA+semA', 'k', 'linewidth', 1)
    plot(angCents, gmA-semA', 'k', 'linewidth', 1)
    axis([0 floor(angCents(end)) 0 0.006])
end