%% Raw Examples
clear
clc
% Load the data in 1Deg conditions
path1Deg = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\WT TB\Fly 1904\DataLowRes.mat';
dt1Deg = load(path1Deg);
dt1Deg = dt1Deg.Flies;
parts1Deg = find(strcmp([dt1Deg.Seq], '1NG'))';

Vr1Deg = [];
Vf1Deg = [];
Vs1Deg = [];
Vt1Deg = [];
X1Deg = [];
Y1Deg = [];
actst1Deg = [];
for n = 1 : length(parts1Deg)
    Vr1Deg = vertcat(Vr1Deg, dt1Deg.Data{parts1Deg(n)}.Vr);
    Vf1Deg = vertcat(Vf1Deg, dt1Deg.Data{parts1Deg(n)}.Vf);
    Vs1Deg = vertcat(Vs1Deg, dt1Deg.Data{parts1Deg(n)}.Vs);
    Vt1Deg = vertcat(Vt1Deg, dt1Deg.Data{parts1Deg(n)}.Vt);
    X1Deg = vertcat(X1Deg, dt1Deg.Data{parts1Deg(n)}.X);
    Y1Deg = vertcat(Y1Deg, dt1Deg.Data{parts1Deg(n)}.Y);
    actst1Deg = vertcat(actst1Deg, dt1Deg.Data{parts1Deg(n)}.actState);
end

% Load the data in 10Deg conditions
path10Deg = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\WT TB\Fly 2116\DataLowRes.mat';
dt10Deg = load(path10Deg);
dt10Deg = dt10Deg.Flies;
parts10Deg = find(strcmp([dt10Deg.Seq], '10NG'))';

Vr10Deg = [];
Vf10Deg = [];
Vs10Deg = [];
Vt10Deg = [];
X10Deg = [];
Y10Deg = [];
actst10Deg = [];
for n = 1 : length(parts10Deg)
    Vr10Deg = vertcat(Vr10Deg, dt10Deg.Data{parts10Deg(n)}.Vr);
    Vf10Deg = vertcat(Vf10Deg, dt10Deg.Data{parts10Deg(n)}.Vf);
    Vs10Deg = vertcat(Vs10Deg, dt10Deg.Data{parts10Deg(n)}.Vs);
    Vt10Deg = vertcat(Vt10Deg, dt10Deg.Data{parts10Deg(n)}.Vt);
    X10Deg = vertcat(X10Deg, dt10Deg.Data{parts10Deg(n)}.X);
    Y10Deg = vertcat(Y10Deg, dt10Deg.Data{parts10Deg(n)}.Y);
    actst10Deg = vertcat(actst10Deg, dt10Deg.Data{parts10Deg(n)}.actState);
end

% define time windows to plot
ti1Deg = dt1Deg.Protocol(parts1Deg(1));
tf1Deg = dt1Deg.Protocol(parts1Deg(end)+1);
t1Deg = ti1Deg:tf1Deg;

ti10Deg = dt10Deg.Protocol(parts10Deg(1));
tf10Deg = dt10Deg.Protocol(parts10Deg(end)+1);
t10Deg = ti10Deg:tf10Deg;

% Get the angular spikes and the forward segments in 10Deg conditions
[params] = GetParams();
[locs, ~, cmhSong, thr] = SortSpikes(Vr10Deg, actst10Deg, params);
[locsF, pksIF, pksEF, ~] = CullSpikesBasedOnVf(Vr10Deg, locs, cmhSong, thr, Vf10Deg, params);
[locsF, pksIF, pksEF, ~] = TemplateCompSpike(Vr10Deg, locsF, pksIF, pksEF, cmhSong, thr, params);
[~,FBouts10Deg] = GetFbouts(actst10Deg, Vf10Deg, locsF, pksIF, pksEF, params);

% Get the angular spikes and the forward segments in 1Deg conditions
[params] = GetParams();
[locs, ~, cmhSong, thr] = SortSpikes(Vr1Deg, actst1Deg, params);
[locsF, pksIF,pksEF, ~] = CullSpikesBasedOnVf(Vr1Deg, locs, cmhSong, thr, Vf1Deg, params);
[locsF, pksIF, pksEF, ~] = TemplateCompSpike(Vr1Deg,locsF,pksIF, pksEF, cmhSong, thr, params);
[~,FBouts1Deg] = GetFbouts(actst1Deg, Vf1Deg, locsF, pksIF,pksEF, params);


addpath('cbrewer')
ncolor = 50;
bss = -10;
[cmap]=cbrewer('div', 'PiYG', ncolor, 'PiYG');

figure,
% XY plot 10Deg
subplot(2,3,1)
plot(X10Deg, Y10Deg,'k', 'linewidth', 1.5)
axis square
axis([-40 40 -40 40])
set(gca,'xtick',[-40 -20 0 20 40]); 
set(gca,'ytick',[-40 -20 0 20 40]);
xlabel('X (mm)')
ylabel('Y (mm)')

% label the forward segments inside the plotting interval
subplot(2,3,2)
for k = 1 : length(FBouts10Deg)
    hold on
    plot(X10Deg(FBouts10Deg{k}), Y10Deg(FBouts10Deg{k}), 'b', 'linewidth', 1.5)
end
axis square
axis([-40 40 -40 40])
set(gca,'xtick',[-40 -20 0 20 40]); 
set(gca,'ytick',[-40 -20 0 20 40]);
xlabel('X (mm)')
ylabel('Y (mm)')

% calculate the straightness and color code it
subplot(2,3,3)
for k = 1 : length(FBouts10Deg)
    vrr = Vr10Deg(FBouts10Deg{k});
    vff = Vf10Deg(FBouts10Deg{k});
    vss = Vs10Deg(FBouts10Deg{k});
    xaux = X10Deg(FBouts10Deg{k});
    yaux = Y10Deg(FBouts10Deg{k});
    vtaux = Vt10Deg(FBouts10Deg{k});
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
    plot(X10Deg(FBouts10Deg{k}), Y10Deg(FBouts10Deg{k}), 'color',cmap(min(ncolor,floor(max(1,1*strt+bss))),:),...
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

% XY plot 1Deg
subplot(2,3,4)
plot(X1Deg, Y1Deg,'k', 'linewidth', 1.5)
axis square
axis([-40 40 -40 40])
set(gca,'xtick',[-40 -20 0 20 40]); 
set(gca,'ytick',[-40 -20 0 20 40]);
xlabel('X (mm)')
ylabel('Y (mm)')

% label the forward segments inside the plotting interval
subplot(2,3,5)
for k = 1 : length(FBouts1Deg)
    hold on
    plot(X1Deg(FBouts1Deg{k}), Y1Deg(FBouts1Deg{k}), 'b', 'linewidth', 1.5)
end
axis square
axis([-40 40 -40 40])
set(gca,'xtick',[-40 -20 0 20 40]); 
set(gca,'ytick',[-40 -20 0 20 40]);
xlabel('X (mm)')
ylabel('Y (mm)')

% calculate the straightness and color code it
subplot(2,3,6)
for k = 1 : length(FBouts1Deg)
    vrr = Vr1Deg(FBouts1Deg{k});
    vff = Vf1Deg(FBouts1Deg{k});
    vss = Vs1Deg(FBouts1Deg{k});
    xaux = X1Deg(FBouts1Deg{k});
    yaux = Y1Deg(FBouts1Deg{k});
    vtaux = Vt1Deg(FBouts1Deg{k});
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
        plot(X1Deg(FBouts1Deg{k}), Y1Deg(FBouts1Deg{k}), 'color',cmap(min(ncolor,floor(max(1,1*strt+bss))),:),...
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

