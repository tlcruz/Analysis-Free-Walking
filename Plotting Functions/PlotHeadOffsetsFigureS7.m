clear
clc

path = '\';

% get the list of flies to analyse
flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowAndHighRes.mat']);
pTypes = dt.Fly.Seq;
pTypes = unique(pTypes);

% initialize variables for each fly and trial type 
HASALLL = cell(length(pTypes), length(flies));
VRALLL = cell(length(pTypes), length(flies));
HASALLR = cell(length(pTypes), length(flies));
VRALLR = cell(length(pTypes), length(flies));

for n = 1 : length(flies)
    % load data for a specific fly
    if exist([path flies(n).name '\DataLowAndHighRes.mat'], 'file') == 2
    dt = load([path flies(n).name '\DataLowAndHighRes.mat']);
    seq = dt.Fly.Seq;
    dt = dt.Fly.Data;
    btthr = 150;
    % iterate across protocol segments
    for k = 1 : length(dt)
        % iterate across trial types
        for l = 1 : length(pTypes)
            HASL = [];
            VRSL = [];
            HASR = [];
            VRSR = [];
            switch seq{k}
                % when trial type matches a protocol segment
                case pTypes{l}
                    % Load walking parameters for this protocol segment
                    difBt = dt{k}.EndBout-dt{k}.StartBout; % bout duration
                    bsha = nanmean(dt{k}.HangDS); % average head angle though trial (baseline)
                    for i = 1 : length(difBt)
                        if difBt(i) > 200
                            if (dt{k}.EndBout(i)>20) && (dt{k}.EndBout(i)+20<length(dt{k}.HangDS))
                                if mean(dt{k}.Vr((dt{k}.EndBout(i)-20):(dt{k}.EndBout(i)))) >= 0 &&  mean(dt{k}.Vr((dt{k}.EndBout(i)-20):(dt{k}.EndBout(i)))) <= inf
                                    % save head angle and angular speed
                                    % around the bout termination for Vr < 0
                                    HASL = horzcat(HASL, dt{k}.HangDS((dt{k}.EndBout(i)-20):(dt{k}.EndBout(i)+20))-bsha);
                                    VRSL = horzcat(VRSL, dt{k}.Vr((dt{k}.EndBout(i)-20):(dt{k}.EndBout(i)+20)));
                                    
                                elseif mean(dt{k}.Vr((dt{k}.EndBout(i)-20):(dt{k}.EndBout(i)))) <= 0 && mean(dt{k}.Vr((dt{k}.EndBout(i)-20):(dt{k}.EndBout(i)))) >= -inf
                                    % save head angle and angular speed
                                    % around the bout termination for Vr > 0
                                    HASR = horzcat(HASR, dt{k}.HangDS((dt{k}.EndBout(i)-20):(dt{k}.EndBout(i)+20))-bsha);
                                    VRSR = horzcat(VRSR, dt{k}.Vr((dt{k}.EndBout(i)-20):(dt{k}.EndBout(i)+20)));
                                end
                            end
                        end
                    end
                    HASALLL{l,n} = horzcat(HASALLL{l,n}, HASL);
                    VRALLL{l,n} = horzcat(VRALLL{l,n}, VRSL);
                    HASALLR{l,n} = horzcat(HASALLR{l,n}, HASR);
                    VRALLR{l,n} = horzcat(VRALLR{l,n}, VRSR);
            end
        end
        
    end
    end
end

% Plot average head angle around walking termination divided by angular
% direction at bout termination
figure,
for i = 1 : length(pTypes)
    subplot(2,ceil(length(pTypes)/2), i)
    HAL = [];
    HAR = [];
    VRL = [];
    VRR = [];
    for k = 1 : size(HASALLL,2)
        HAL = horzcat(HAL, HASALLL{i,k});
        HAR = horzcat(HAR, HASALLR{i,k});
        VRL = horzcat(VRL, VRALLL{i,k});
        VRR = horzcat(VRR, VRALLR{i,k});
    end
    hold on
    if (~isempty(HAL)) && (~isempty(HAR))
        plot((-20:20)/60,nanmean(HAL,2),'color', [0 0 1],'linewidth', 3)
        plot((-20:20)/60,nanmean(HAL,2)+nanstd(HAL,1,2)/sqrt(length(path)),'color', [0 0 1],'linewidth', 1)
        plot((-20:20)/60,nanmean(HAL,2)-nanstd(HAL,1,2)/sqrt(length(path)),'color', [0 0 1],'linewidth', 1)
        plot((-20:20)/60,nanmean(HAR,2),'color', [1 0 0],'linewidth', 3)
        plot((-20:20)/60,nanmean(HAR,2)+nanstd(HAR,1,2)/sqrt(length(path)),'color', [1 0 0],'linewidth', 1)
        plot((-20:20)/60,nanmean(HAR,2)-nanstd(HAR,1,2)/sqrt(length(path)),'color', [1 0 0],'linewidth', 1)
        hold on
        plot((-20:20)/60,0.01*nanmean(VRL,2)-15,'color', [0 0 1],'linewidth', 2)
        plot((-20:20)/60,0.01*nanmean(VRL,2)+0.01*nanstd(VRL,1,2)/sqrt(length(path))-15,'color', [0 0 1],'linewidth', 1)
        plot((-20:20)/60,0.01*nanmean(VRL,2)-0.01*nanstd(VRL,1,2)/sqrt(length(path))-15,'color', [0 0 1],'linewidth', 1)
        plot((-20:20)/60,0.01*nanmean(VRR,2)-15,'color', [1 0 0],'linewidth', 2)
        plot((-20:20)/60,0.01*nanmean(VRR,2)+0.01*nanstd(VRR,1,2)/sqrt(length(path))-15,'color', [1 0 0],'linewidth', 2)
        plot((-20:20)/60,0.01*nanmean(VRR,2)-0.01*nanstd(VRR,1,2)/sqrt(length(path))-15,'color', [1 0 0],'linewidth', 2)
        hal1 = nanmean(HAL,2);
        hal2 = nanmean(HAR,2);
        title([pTypes{i} '   ' num2str(mean(hal1(25:40))-mean(hal2(25:40)))])
        axis([-20/60 20/60 -25 10])
    end
end