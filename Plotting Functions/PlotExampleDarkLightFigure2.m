%% Raw Examples - Requires saccade segmenter
clear
clc
% Load the data in dark conditions
pathDark = '\Data\Dark\Dark WTTB\Fly 1643\DataLowRes.mat';
dtDark = load(pathDark);
dtDark = dtDark.Flies;
partsDark = [3 4];

VrDark = [];
VfDark = [];
VsDark = [];
VtDark = [];
XDark = [];
YDark = [];
actstDark = [];
for n = 1 : length(partsDark)
    VrDark = vertcat(VrDark, dtDark.Data{partsDark(n)}.Vr);
    VfDark = vertcat(VfDark, dtDark.Data{partsDark(n)}.Vf);
    VsDark = vertcat(VsDark, dtDark.Data{partsDark(n)}.Vs);
    VtDark = vertcat(VtDark, dtDark.Data{partsDark(n)}.Vt);
    XDark = vertcat(XDark, dtDark.Data{partsDark(n)}.X);
    YDark = vertcat(YDark, dtDark.Data{partsDark(n)}.Y);
    actstDark = vertcat(actstDark, dtDark.Data{partsDark(n)}.actState);
end



% Load the data in light conditions
pathLight = '\Data\Bilateral Stimulus Size\WT TB\Fly 1902\DataLowRes.mat';
dtLight = load(pathLight);
dtLight = dtLight.Flies;
partsLight = find(strcmp([dtLight.Seq], '10NG'))';

VrLight = [];
VfLight = [];
VsLight = [];
VtLight = [];
XLight = [];
YLight = [];
actstLight = [];
for n = 1 : length(partsLight)
    VrLight = vertcat(VrLight, dtLight.Data{partsLight(n)}.Vr);
    VfLight = vertcat(VfLight, dtLight.Data{partsLight(n)}.Vf);
    VsLight = vertcat(VsLight, dtLight.Data{partsLight(n)}.Vs);
    VtLight = vertcat(VtLight, dtLight.Data{partsLight(n)}.Vt);
    XLight = vertcat(XLight, dtLight.Data{partsLight(n)}.X);
    YLight = vertcat(YLight, dtLight.Data{partsLight(n)}.Y);
    actstLight = vertcat(actstLight, dtLight.Data{partsLight(n)}.actState);
end

% define time windows to plot
tiDark = dtDark.Protocol(partsDark(1));
tfDark = dtDark.Protocol(partsDark(end)+1);
tDark = tiDark:tfDark;

tiLight = dtLight.Protocol(partsLight(1));
tfLight = dtLight.Protocol(partsLight(end)+1);
tLight = tiLight:tfLight;

% Get the angular spikes and the forward segments in light conditions
[params] = GetParams();
[locs, ~, cmhSong, thr, dx] = SortSpikes(VrLight, actstLight, params);
[locsF, pksIF, pksEF, ~] = CullSpikesBasedOnVf(VrLight, locs, cmhSong, thr, VfLight, params);
[locsF, pksIF, pksEF, scores] = TemplateCompSpike(VrLight,locsF,pksIF, pksEF, cmhSong, thr, params);
[FBoutsLight, FBoutsLightAll] = GetFbouts(actstLight, VfLight, locsF, pksIF,pksEF, params);

% Get the angular spikes and the forward segments in dark conditions
[params] = GetParams();
[locs, ~, cmhSong, thr] = SortSpikes(VrDark, actstDark, params);
[locsF, pksIF,pksEF, ~] = CullSpikesBasedOnVf(VrDark, locs, cmhSong, thr, VfDark, params);
[locsF, pksIF, pksEF, ~] = TemplateCompSpike(VrDark,locsF,pksIF, pksEF, cmhSong, thr, params);
[FBoutsDark, FBoutsDarkAll] = GetFbouts(actstDark, VfDark, locsF, pksIF,pksEF, params);


addpath('cbrewer')
ncolor = 50;
bss = -10;
[cmap]=cbrewer('div', 'PiYG', ncolor, 'PiYG');

figure,
% XY plot Light
subplot(2,3,1)
plot(XLight, YLight,'k', 'linewidth', 1.5)
axis square
axis([-40 40 -40 40])
set(gca,'xtick',[-40 -20 0 20 40]); 
set(gca,'ytick',[-40 -20 0 20 40]);
xlabel('X (mm)')
ylabel('Y (mm)')

% label the forward segments inside the plotting interval
subplot(2,3,2)
for k = 1 : length(FBoutsLightAll)
    hold on
    plot(XLight(FBoutsLightAll{k}), YLight(FBoutsLightAll{k}), 'b', 'linewidth', 1.5)
end
axis square
axis([-40 40 -40 40])
set(gca,'xtick',[-40 -20 0 20 40]); 
set(gca,'ytick',[-40 -20 0 20 40]);
xlabel('X (mm)')
ylabel('Y (mm)')

% calculate the straightness and color code it
subplot(2,3,3)
for k = 1 : length(FBoutsLight)
    vrr = VrLight(FBoutsLight{k});
    vff = VfLight(FBoutsLight{k});
    vss = VsLight(FBoutsLight{k});
    xaux = XLight(FBoutsLight{k});
    yaux = YLight(FBoutsLight{k});
    vtaux = VtLight(FBoutsLight{k});
    nv = length(xaux) - params.windStr;
    dm = [];
    ss = [];
    for l = 1 : nv
        pi = [xaux(l),yaux(l),0];
        pf = [xaux(l+params.windStr),yaux(l+params.windStr),0];
        dis = sum(vtaux((l):(l+params.windStr)))/60;
        pt = [xaux(l+floor(params.windStr/2)),yaux(l+floor(params.windStr/2)),0];
        dr = point_to_line(pt,pi,pf);
        dm = vertcat(dm, dr);
        ss = vertcat(ss, dis);
    end
    if size(dm,1)>1
        strt = sum(ss)/sum(dm);
    hold on
    plot(XLight(FBoutsLight{k}), YLight(FBoutsLight{k}), 'color',cmap(min(ncolor,floor(max(1,1*strt+bss))),:),...
        'linewidth', 1.5)
    end
end
axis square
axis([-40 40 -40 40])
set(gca,'xtick',[-40 -20 0 20 40]); 
set(gca,'ytick',[-40 -20 0 20 40]);
xlabel('X (mm)')
ylabel('Y (mm)')
colormap(cmap);
c = colorbar('Ticks',[0, 0.25, 0.5, 0.75, 1], ...
    'TickLabels',{'10', '20', '30', '40', '50'});
c.Label.String = 'Straightness (a.u.)';

% XY plot Dark
subplot(2,3,4)
plot(XDark, YDark,'k', 'linewidth', 1.5)
axis square
axis([-40 40 -40 40])
set(gca,'xtick',[-40 -20 0 20 40]); 
set(gca,'ytick',[-40 -20 0 20 40]);
xlabel('X (mm)')
ylabel('Y (mm)')

% label the forward segments inside the plotting interval
subplot(2,3,5)
for k = 1 : length(FBoutsDarkAll)
    hold on
    plot(XDark(FBoutsDarkAll{k}), YDark(FBoutsDarkAll{k}), 'b', 'linewidth', 1.5)
end
axis square
axis([-40 40 -40 40])
set(gca,'xtick',[-40 -20 0 20 40]); 
set(gca,'ytick',[-40 -20 0 20 40]);
xlabel('X (mm)')
ylabel('Y (mm)')

% calculate the straightness and color code it
subplot(2,3,6)
for k = 1 : length(FBoutsDark)
    vrr = VrDark(FBoutsDark{k});
    vff = VfDark(FBoutsDark{k});
    vss = VsDark(FBoutsDark{k});
    xaux = XDark(FBoutsDark{k});
    yaux = YDark(FBoutsDark{k});
    vtaux = VtDark(FBoutsDark{k});
    nv = length(xaux) - params.windStr;
    dm = [];
    ss = [];
    for l = 1 : nv
        pi = [xaux(l),yaux(l),0];
        pf = [xaux(l+params.windStr),yaux(l+params.windStr),0];
        dis = sum(vtaux((l):(l+params.windStr)))/60;
        pt = [xaux(l+floor(params.windStr/2)),yaux(l+floor(params.windStr/2)),0];
        dr = point_to_line(pt,pi,pf);
        dm = vertcat(dm, dr);
        ss = vertcat(ss, dis);
    end
    if size(dm,1)>1
        strt = sum(ss)/sum(dm);
        hold on
        plot(XDark(FBoutsDark{k}), YDark(FBoutsDark{k}), 'color',cmap(min(ncolor,floor(max(1,1*strt+bss))),:),...
            'linewidth', 1.5)
    end
end
axis square
axis([-40 40 -40 40])
set(gca,'xtick',[-40 -20 0 20 40]); 
set(gca,'ytick',[-40 -20 0 20 40]);
xlabel('X (mm)')
ylabel('Y (mm)')
colormap(cmap);
c = colorbar('Ticks',[0, 0.25, 0.5, 0.75, 1], ...
    'TickLabels',{'10', '20', '30', '40', '50'});
c.Label.String = 'Straightness (a.u.)';