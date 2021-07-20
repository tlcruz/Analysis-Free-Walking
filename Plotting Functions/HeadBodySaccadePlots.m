%% plot head, gaze and body movements during saccades
clear
clc
path = '\';
% get paramenters and protocol types
params = GetParams();
flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowAndHighRes.mat'], 'Fly');
seq = dt.Fly.Seq;
pTypes = unique(seq);

SPKHA = cell(length(flies),length(pTypes));
SPKVR = cell(length(flies),length(pTypes));
SPKVF = cell(length(flies),length(pTypes));
SPKVS = cell(length(flies),length(pTypes));

% define window around saccade
wi = 20;
dm = 7;
vrmin = 100;
vrmax = 200;
errThr = 0.4;
for n = 1 : length(flies)
    dt = load([path flies(n).name '\DataLowAndHighRes.mat'], 'Fly');
    seq = dt.Fly.Seq;
    dt = dt.Fly.Data;
    for k = 1 : length(dt)
        for j = 1 : length(pTypes)
            spkha = [];
            spkvr = [];
            spkvf = [];
            spkvs = [];
            spkvrrand = [];
            spkvfrand = [];
            spkvsrand = [];
            spkdist = [];
            switch seq{k}
                case pTypes{j}
                    Vf = dt{k}.Vf;
                    Vs = dt{k}.Vs;
                    Vr = dt{k}.Vr;
                    Hang = (dt{k}.HangDS-nanmedian(dt{k}.HangDS));
                    HangErr = dt{k}.ErrTrkHRDS;
                    acSt = dt{k}.actState;
                    % segment the saccades
                    [locs, ~, cmhSong, thr] = SortSpikes(Vr, acSt, params);
                    [locsF, pksIF,pksEF, ~] = CullSpikesBasedOnVf(Vr, locs, cmhSong, thr, Vf, params);
                    [locsF, ~, ~] = TemplateCompSpike(Vr,locsF,pksIF, pksEF, cmhSong, thr, params);
                    locsRand = randi([wi+dm+1 length(Vr)-(wi+dm+1)], length(locsF),1);
                    % crop window around the saccades
                    if(length(locsF) > 3)
                        for st = 1 : length(locsF)
                            if locsF(st)-wi < length(HangErr)
                                if st == 1
                                    if locsF(st)>wi+dm && abs(Vr(locsF(st))) > vrmin && abs(Vr(locsF(st))) < vrmax
                                        if (locsF(st+1)-locsF(st) > wi+dm)
                                            if isempty(find(HangErr((locsF(st)-wi):(locsF(st)+wi)) > errThr))
                                                spkha = horzcat(spkha, Hang((locsF(st)-wi):(locsF(st)+wi)));
                                                spkvr = horzcat(spkvr, Vr((locsF(st)-wi):(locsF(st)+wi)));
                                                spkvf = horzcat(spkvf, Vf((locsF(st)-wi):(locsF(st)+wi)));
                                                spkvs = horzcat(spkvs, Vs((locsF(st)-wi):(locsF(st)+wi)));
                                            end
                                        end
                                    end
                                elseif st == length(locsF)
                                    if locsF(st)+wi+dm < length(Vr) && abs(Vr(locsF(st))) > vrmin && abs(Vr(locsF(st))) < vrmax
                                        if (locsF(st)-locsF(st-1) > wi+dm)
                                            if isempty(find(HangErr((locsF(st)-wi):(locsF(st)+wi)) > errThr))
                                                spkha = horzcat(spkha, Hang((locsF(st)-wi):(locsF(st)+wi)));
                                                spkvr = horzcat(spkvr, Vr((locsF(st)-wi):(locsF(st)+wi)));
                                                spkvf = horzcat(spkvf, Vf((locsF(st)-wi):(locsF(st)+wi)));
                                                spkvs = horzcat(spkvs, Vs((locsF(st)-wi):(locsF(st)+wi)));
                                            end
                                        end
                                    end
                                else
                                    if (locsF(st)-locsF(st-1) > wi+dm && locsF(st+1)-locsF(st) > wi+dm) && abs(Vr(locsF(st))) > vrmin && abs(Vr(locsF(st))) < vrmax
                                        if length(find(HangErr((locsF(st)-wi):(locsF(st)+wi)) > errThr)) == 0
                                            spkha = horzcat(spkha, Hang((locsF(st)-wi):(locsF(st)+wi)));
                                            spkvr = horzcat(spkvr, Vr((locsF(st)-wi):(locsF(st)+wi)));
                                            spkvf = horzcat(spkvf, Vf((locsF(st)-wi):(locsF(st)+wi)));
                                            spkvs = horzcat(spkvs, Vs((locsF(st)-wi):(locsF(st)+wi)));
                                        end
                                    end
                                end
                            end
                        end
                    end
            end
            SPKHA{n,j} = horzcat(SPKHA{n,j}, spkha);
            SPKVR{n,j} = horzcat(SPKVR{n,j}, spkvr);
            SPKVF{n,j} = horzcat(SPKVF{n,j}, spkvf);
            SPKVS{n,j} = horzcat(SPKVS{n,j}, spkvs);
        end
    end
    disp(['Fly ' num2str(n) ' of ' num2str(length(flies)) ' Done'])
end
% combined plot of head, body ang gaze position and velocity
figure,
aux = 1;
npl = 4;
for cond = [5 7 9 1 3]%
    haa = [];
    hvv = [];
    aar = [];
    vvr = [];
    gzaa = [];
    gzvv = [];
    vvf = [];
    haRVfH = [];
    haRVfL = [];
    hvRVfH = [];
    hvRVfL = [];
    vrRVfH = [];
    vrRVfL = [];
    arRVfH = [];
    arRVfL = [];
    gzaRVfH = [];
    gzaRVfL = [];
    gzvRVfH = [];
    gzvRVfL = [];
    vfRVfH = [];
    vfRVfL = [];
    haLVfH = [];
    haLVfL = [];
    hvLVfH = [];
    hvLVfL = [];
    vrLVfH = [];
    vrLVfL = [];
    arLVfH = [];
    arLVfL = [];
    gzaLVfH = [];
    gzaLVfL = [];
    gzvLVfH = [];
    gzvLVfL = [];
    vfLVfH = [];
    vfLVfL = [];
    for n = 1 : size(SPKHA,1)
        for ns = 1 : size(SPKHA{n,cond},2)
            ha = (SPKHA{n,cond}(:,ns)-SPKHA{n,cond}(1,ns));
            hv = vertcat(0,60*diff(SPKHA{n,cond}(:,ns)));
            vr = SPKVR{n,cond}(:,ns);
            ar = (cumsum(SPKVR{n,cond}(:,ns))-SPKVR{n,cond}(1,ns))/60;
            vf = SPKVF{n,cond}(:,ns);
            gza = ar + ha;
            gzv = vr + hv;
            
            if vr(20) > 0 && nanmean(vf(1:15)) < 5
                haa = horzcat(haa, ha);
                hvv = horzcat(hvv, hv);
                aar = horzcat(aar, ar);
                vvr = horzcat(vvr, vr);
                gzaa = horzcat(gzaa, gza);
                gzvv = horzcat(gzvv, gzv);
                vvf = horzcat(vvf, vf);
                haRVfL = horzcat(haRVfL, ha);
                hvRVfL = horzcat(hvRVfL, hv);
                vrRVfL = horzcat(vrRVfL, vr);
                arRVfL = horzcat(arRVfL, ar);
                gzaRVfL = horzcat(gzaRVfL, gza);
                gzvRVfL = horzcat(gzvRVfL, gzv);
                vfRVfL = horzcat(vfRVfL, vf);
            elseif vr(20) > 0 && nanmean(vf(1:15)) > 5
                haa = horzcat(haa, ha);
                hvv = horzcat(hvv, hv);
                aar = horzcat(aar, ar);
                vvr = horzcat(vvr, vr);
                gzaa = horzcat(gzaa, gza);
                gzvv = horzcat(gzvv, gzv);
                vvf = horzcat(vvf, vf);
                haRVfH = horzcat(haRVfH, ha);
                hvRVfH = horzcat(hvRVfH, hv);
                vrRVfH = horzcat(vrRVfH, vr);
                arRVfH = horzcat(arRVfH, ar);
                gzaRVfH = horzcat(gzaRVfH, gza);
                gzvRVfH = horzcat(gzvRVfH, gzv);
                vfRVfH = horzcat(vfRVfH, vf);
            elseif vr(20) < 0 && nanmean(vf(1:15)) < 5
                haa = horzcat(haa, -ha);
                hvv = horzcat(hvv, -hv);
                aar = horzcat(aar, -ar);
                vvr = horzcat(vvr, -vr);
                gzaa = horzcat(gzaa, -gza);
                gzvv = horzcat(gzvv, -gzv);
                vvf = horzcat(vvf, vf);
                haLVfL = horzcat(haLVfL, ha);
                hvLVfL = horzcat(hvLVfL, hv);
                vrLVfL = horzcat(vrLVfL, vr);
                arLVfL = horzcat(arLVfL, ar);
                gzaLVfL = horzcat(gzaLVfL, gza);
                gzvLVfL = horzcat(gzvLVfL, gzv);
                vfLVfL = horzcat(vfLVfL, vf);
            elseif vr(20) < 0 && nanmean(vf(1:15)) > 5
                haa = horzcat(haa, -ha);
                hvv = horzcat(hvv, -hv);
                aar = horzcat(aar, -ar);
                vvr = horzcat(vvr, -vr);
                gzaa = horzcat(gzaa, -gza);
                gzvv = horzcat(gzvv, -gzv);
                vvf = horzcat(vvf, vf);
                haLVfH = horzcat(haLVfH, ha);
                hvLVfH = horzcat(hvLVfH, hv);
                vrLVfH = horzcat(vrLVfH, vr);
                arLVfH = horzcat(arLVfH, ar);
                gzaLVfH = horzcat(gzaLVfH, gza);
                gzvLVfH = horzcat(gzvLVfH, gzv);
                vfLVfH = horzcat(vfLVfH, vf);
            end
        end
    end

    t = 1000*(-20:20)/60;
    subplot(2,npl,aux)
    hold on
    plot(t,nanmean(haa,2), 'm', 'linewidth', 2)
    plot(t,nanmean(haa,2)+nanstd(haa,1,2), 'm', 'linewidth', 0.5)
    plot(t,nanmean(haa,2)-nanstd(haa,1,2), 'm', 'linewidth', 0.5)
    plot(t,nanmean(aar,2), 'r', 'linewidth', 2)
    plot(t,nanmean(aar,2)+nanstd(aar,1,2), 'r', 'linewidth', 0.5)
    plot(t,nanmean(aar,2)-nanstd(aar,1,2), 'r', 'linewidth', 0.5)
    plot(t,nanmean(gzaa,2), 'c', 'linewidth', 2)
    plot(t,nanmean(gzaa,2)+nanstd(gzaa,1,2), 'c', 'linewidth', 0.5)
    plot(t,nanmean(gzaa,2)-nanstd(gzaa,1,2), 'c', 'linewidth', 0.5)
    axis([-330 330 -20 100])
    xlabel('Time (ms)')
    ylabel('Angle (º)')
    title([pTypes{cond}])
    subplot(2,npl, aux + npl)
    hold on
    plot(t,nanmean(hvv,2), 'm', 'linewidth', 2)
    plot(t,nanmean(hvv,2)+nanstd(hvv,1,2), 'm', 'linewidth', 0.5)
    plot(t,nanmean(hvv,2)-nanstd(hvv,1,2), 'm', 'linewidth', 0.5)
    plot(t,nanmean(vvr,2), 'r', 'linewidth', 2)
    plot(t,nanmean(vvr,2)+nanstd(vvr,1,2), 'r', 'linewidth', 0.5)
    plot(t,nanmean(vvr,2)-nanstd(vvr,1,2), 'r', 'linewidth', 0.5)
    plot(t,nanmean(gzvv,2), 'c', 'linewidth', 2)
    plot(t,nanmean(gzvv,2)+nanstd(gzvv,1,2), 'c', 'linewidth', 0.5)
    plot(t,nanmean(gzvv,2)-nanstd(gzvv,1,2), 'c', 'linewidth', 0.5)
    axis([-330 330 -150 750])
    xlabel('Time (ms)')
    ylabel('Speed (º/s)')
    aux = aux + 1;
end
