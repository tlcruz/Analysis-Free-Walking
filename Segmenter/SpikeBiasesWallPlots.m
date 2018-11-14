clear
clc
path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\WT TB\';
params = GetParams();

flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowRes.mat'], 'Flies');
seq = dt.Flies.Seq;
pTypes = unique(seq);

SPKVAL = cell(length(flies),length(pTypes));
SPKWDIST = cell(length(flies),length(pTypes));
SPKPROB = cell(length(flies),length(pTypes));
AngleP = cell(length(flies),length(pTypes)); 
AngleA = cell(length(flies),length(pTypes));

wDCents = 0:2:50;
wdbin = 2;
for n = 1 : length(flies)
    dt = load([path flies(n).name '\DataLowRes.mat'], 'Flies');
    seq = dt.Flies.Seq;
    dt = dt.Flies.Data;
    for k = 1 : length(dt)
        for j = 1 %: length(pTypes)
            spkVal = [];
            wDist = [];
            ap = [];
            aa = [];
            distr = [];
            switch seq{k}
                case pTypes{j}
                    Vf = dt{k}.Vf;
                    Vr = dt{k}.Vr;
                    X = dt{k}.X;
                    Y = dt{k}.Y;
                    WD = dt{k}.WallDist;
                    acSt = dt{k}.actState;
                    Angle = dt{k}.Angle;
                    [locs, ~, cmhSong, thr] = SortSpikes(Vr, acSt, params);
                    [locsF, pksIF,pksEF, ~] = CullSpikesBasedOnVf(Vr, locs, cmhSong, thr, Vf, params);
                    [locsF, pksIF, pksEF] = TemplateCompSpike(Vr,locsF,pksIF, pksEF, cmhSong, thr, params);
                    spkSt = zeros(length(acSt),1);
                    for i = 1 : length(locsF)
                        if X(locsF(i)) < 0 && Y(locsF(i)) > 0
                            spkVal = vertcat(spkVal, Vr(locsF(i)));
                            wDist = vertcat(wDist, WD(locsF(i)));
                            ap = vertcat(ap, Angle(locsF(i)-pksIF(i)));
                            aa = vertcat(aa, Angle(locsF(i)+pksIF(i)));
                            spkSt(locsF(i)) = 1;
                        end
                    end
                    WDAct = WD(acSt==1);
                    WDSpk = WD(spkSt==1);
                    distr = [];
                    
                    for x = 1 : length(wDCents)
                        indsa = find(WDAct >= wDCents(1)+(x-1)*wdbin & WDAct < x*wdbin);
                        indsb = find(WDSpk >= wDCents(1)+(x-1)*wdbin & WDSpk < x*wdbin);
                        distr = vertcat(distr, length(indsb)/(length(indsa)+1));
                    end 
            end
            SPKPROB{n,j} = horzcat(SPKPROB{n,j}, distr);
            SPKVAL{n,j} = vertcat(SPKVAL{n,j}, spkVal);
            SPKWDIST{n,j} = vertcat(SPKWDIST{n,j}, wDist);
            AngleP{n,j} = vertcat(AngleP{n,j}, ap);
            AngleA{n,j} = vertcat(AngleA{n,j}, aa);
        end
    end
end

%%
distSpkP = [];
for k = 1 : size(SPKPROB,2)
    if mod(k, 2) == 1
        for n = 1 : size(SPKPROB,2)
            distSpkP = horzcat(distSpkP,SPKPROB{n,k});
        end
    end
end

plot(wDCents, mean(distSpkP,2))

%%
figure,
wDCents = 0:1:50;

for k = 1 : size(SPKVAL,2)
    spkAmp = [];
    spkWD = [];
    spkA = [];
    for n = 1 : size(SPKVAL,1)
        spkAmp = vertcat(spkAmp, abs(SPKVAL{n,k}));
        spkWD = vertcat(spkWD, abs(SPKWDIST{n,k}));
        wdf = abs(SPKWDIST{n,k});
        amf = abs(SPKVAL{n,k});
        sp = zeros(size(wDCents));
        for x = 1 : length(wDCents)
            if x < length(wDCents)
                wdbin = wDCents(x+1)-wDCents(x);
            else
                wdbin = -wDCents(x-1)+wDCents(x);
            end
            tmpa=find(wdf >= wDCents(1)+(x-1)*wdbin);
            tmpb=find(wdf < wDCents(1)+x*wdbin);
            tmpx=tmpa(ismembc(tmpa,tmpb));
            if(length(tmpx)>00)
                sp(x) = mean(amf(tmpx));
            else
                sp(x) = nan;
            end
        end
        spkA = vertcat(spkA, sp);
    end

    subplot(2,4,k)
    hold on
    errorbar(wDCents, nanmean(spkA), nanstd(spkA,1,1)/sqrt(size(SPKVAL,1)),'ok')
    title(pTypes{k})
    axis([0 45 0 1000])
    xlabel('Dist. Wall (mm)')
    ylabel('Spike Turn Size (º/s)')
    disp(num2str(k))
end


%%
dwauxThrmin = 0;
dwauxThrmax = 40;

figure,
wDCents = -180:25:180;
for k = 1 : size(SPKVAL,2)
    if mod(k,2) == 1
    ap = [];
    spkAmp = [];
    spkA = [];
    for n = 1 : size(SPKVAL,1)
        dwaux = SPKWDIST{n,k};
        apaux = mod(AngleA{n,k}+180,360)-180;
        spaux = sign(SPKVAL{n,k});
        inds = find(dwaux<dwauxThrmax & dwaux>dwauxThrmin);

        apaux = apaux(inds);
        spaux = spaux(inds);
        
        ap = vertcat(ap, apaux);
        spkAmp = vertcat(spkAmp, spaux);
        
        sp = zeros(size(wDCents));
        for x = 1 : length(wDCents)
            if x < length(wDCents)
                wdbin = wDCents(x+1)-wDCents(x);
            else
                wdbin = -wDCents(x-1)+wDCents(x);
            end
            tmpa=find(apaux >= wDCents(1)+(x-1)*wdbin);
            tmpb=find(apaux < wDCents(1)+x*wdbin);
            tmpx=tmpa(ismembc(tmpa,tmpb));
            if(length(tmpx)>0)
                sp(x) = mean(spaux(tmpx));
            else
                sp(x) = nan;
            end
        end
        spkA = vertcat(spkA, sp);
    end

    subplot(2,5,k)
    hold on
    errorbar(wDCents, nanmean(spkA), nanstd(spkA,1,1)/sqrt(size(SPKVAL,1)),'ok')
    title(pTypes{k})
    axis([-180 180 -1 1])
%     xlabel('Dist. Wall (mm)')
%     ylabel('Spike Turn Size (º/s)')
    disp(num2str(k))
    end
end





            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
