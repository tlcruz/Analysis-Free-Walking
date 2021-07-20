%% Plot gaze deviations
path = '\';

% get paramenters and protocol types
params = GetParams();
flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowAndHighRes.mat']);
pTypes = dt.Fly.Seq;
pTypes = unique(pTypes);

mbody = cell(length(flies), length(pTypes));
nbody = cell(length(flies), length(pTypes));
mgaze = cell(length(flies), length(pTypes));
ngaze = cell(length(flies), length(pTypes));

for n = 1 : length(flies)
    % if there is both high and low resolution data
    if exist([path flies(n).name '\DataLowAndHighRes.mat'], 'file') == 2
        dt = load([path flies(n).name '\DataLowAndHighRes.mat']);
        seq = dt.Fly.Seq;
        dt = dt.Fly.Data;
        for k = 1 : length(dt)-1
            for l = 1 : length(pTypes)
                switch seq{k}
                    case pTypes{l}
                        vrb = dt{k}.Vr;
                        vfb = dt{k}.Vf;
                        vsb = dt{k}.Vs;
                        actst = dt{k}.actState;
                        hA = dt{k}.HangDS;
                        % segment angular spikes and forward segments
                        [locs, ~, cmhSong, thr] = SortSpikes(vrb, actst, params);
                        [locsF, pksIF,pksEF] = CullSpikesBasedOnVf(vrb, locs, cmhSong, thr, vfb, params);
                        [locsF, pksIF, pksEF] = TemplateCompSpike(vrb,locsF,pksIF, pksEF, cmhSong, thr, params);
                        [FBouts] = GetFbouts(actst, vfb, locsF, pksIF,pksEF, params);
                        
                        avg = [];
                        avb = [];
                        nn = [];
                        % iterate across forward runs
                        for j = 1 : length(FBouts)
                            if length(hA) > FBouts{j}(end)
                                vr = vrb(FBouts{j});
                                bb = cumsum(vr)/60; % body angle during forward segment
                                hh = hA(FBouts{j}); % head angle during forward segment
                                inds = find(abs(diff(hh))>5);
                                if ~isempty(inds)
                                    if (inds(1) > 2)
                                        inds = vertcat(inds(1)-2,inds(1)-1, inds, inds(1) + 1);
                                    end
                                end
                                hh(inds) = nan;
                                gz = bb + hh; % gaze angle during forward segment
                                
                                % compute gaze, body and head deviations
%                                 if(max(abs(bb)) < 15 && length(find(isnan(gz)==0)) > 20 && max(abs(hh))<10)
                                if(max(abs(bb)) < 30 && length(find(isnan(gz)==0)) > 20  && max(abs(hh))<10)
                                    avg = vertcat(avg, nanstd(gz));
                                    avb = vertcat(avb, nanstd(bb));
                                    nn = vertcat(nn, length(bb));
                                end
                            end
                        end

                        if ~isempty(avg)
                            mbody{n,l} = vertcat(mbody{n,l}, nanmean(avb));
                            nbody{n,l} = vertcat(nbody{n,l}, sum(nn));
                            mgaze{n,l} = vertcat(mgaze{n,l}, nanmean(avg));
                            ngaze{n,l} = vertcat(ngaze{n,l}, sum(nn));
                        end
                end
            end
        end
    end
    disp([num2str(n) ' out of ' num2str(length(flies))])
end

% plot body and gaze deviations
cfs = jet(length(flies));
figure,
for j = 1 : ceil(length(pTypes)/2)
    subplot(1, ceil(length(pTypes)/2), j)
    hold on
    l=2*j-1;
    ab = [];
    ag = [];
    for n = 1: length(flies)
        inds = find(isnan(mgaze{n,l}) == 0);
        if ~isempty(inds)
            if nbody{n,l}(inds) > 0
                avb = mbody{n,l}(inds)'*nbody{n,l}(inds)/sum(nbody{n,l}(inds));
                semvb = sqrt(((mbody{n,l}(inds)-avb).*(mbody{n,l}(inds)-avb))'*nbody{n,l}(inds)/sum(nbody{n,l}(inds)))/sqrt(length(flies));
                avg = mgaze{n,l}(inds)'*ngaze{n,l}(inds)/sum(ngaze{n,l}(inds));
                semvg = sqrt(((mgaze{n,l}(inds)-avg).*(mgaze{n,l}(inds)-avg))'*ngaze{n,l}(inds)/sum(ngaze{n,l}(inds)))/sqrt(length(flies));
                ab = vertcat(ab, avb);
                ag = vertcat(ag, avg);
                errorbar([1 2], [avb avg], [semvb semvg], '-o', 'color', cfs(n,:))
            end
        end
    end
    % do pairwise comparison
    [p] = signrank(ab, ag);
    title([pTypes{l} ' ' num2str(p)])
    axis([0 3 0 10])
    ylabel('Average angle change')
end

%% Plot examples of head, body and gaze traces for one fly
avg = [];
avb = [];
for j = 1: length(FBouts)
    vr = vrb(FBouts{j});
    bb = cumsum(vr)/60;
    hh = hA(FBouts{j});
    gz = bb + hh;
    
    if(max(abs(bb)) < 10 && length(find(isnan(gz)==0)) > 10 && max(abs(hh))<10 ) 
        figure,
        hold on
        plot([0 length(bb)/60], [0 0], '--', 'color', [0.6 0.6 0.6])
        plot((1:length(hh))/60, hh, 'color', [0 0.8 0])
        plot((1:length(bb))/60, bb, 'color', [0 0 0.8])
        plot((1:length(gz))/60, gz, 'color', [0.8 0 0]);
        title([num2str(abs(nanmean(bb))-abs(nanmean(gz)))])
        axis([0 length(bb)/60 -7 7])
        avg = vertcat(avg, abs(nanmean(gz)));
        avb = vertcat(avb, abs(nanmean(bb)));
    end 
end



