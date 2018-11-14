fc = 2639;
pathCtr = ['C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\VTR39\VTR39Kir\Fly ' num2str(fc) '\DataLowRes.mat'];
dtCtr = load(pathCtr);
dtCtr = dtCtr.Flies;
partsCtr = find(strcmp([dtCtr.Seq], '10NG'))';

fk = 2659;
pathKir = ['C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\VTR39\VTR39FLPKir\Fly ' num2str(fk) '\DataLowRes.mat'];
dtKir = load(pathKir);
dtKir = dtKir.Flies;
partsKir = find(strcmp([dtKir.Seq], '10NG'))';


tiCtr = dtCtr.Protocol(partsCtr(1));
tfCtr = dtCtr.Protocol(partsCtr(end)+1);
tCtr = tiCtr:tfCtr;

tiKir = dtKir.Protocol(partsKir(1));
tfKir = dtKir.Protocol(partsKir(end)+1);
tKir = tiKir:tfKir;

VrCtr = [];
VfCtr = [];
VsCtr = [];
VtCtr = [];
XCtr = [];
YCtr = [];
actstCtr = [];
for n = 1 : length(partsCtr)
    VrCtr = vertcat(VrCtr, dtCtr.Data{partsCtr(n)}.Vr);
    VfCtr = vertcat(VfCtr, dtCtr.Data{partsCtr(n)}.Vf);
    VsCtr = vertcat(VsCtr, dtCtr.Data{partsCtr(n)}.Vs);
    VtCtr = vertcat(VtCtr, dtCtr.Data{partsCtr(n)}.Vt);
    XCtr = vertcat(XCtr, dtCtr.Data{partsCtr(n)}.X);
    YCtr = vertcat(YCtr, dtCtr.Data{partsCtr(n)}.Y);
    actstCtr = vertcat(actstCtr, dtCtr.Data{partsCtr(n)}.actState);
end

VrKir = [];
VfKir = [];
VsKir = [];
VtKir = [];
XKir = [];
YKir = [];
actstKir = [];
for n = 1 : length(partsKir)
    VrKir = vertcat(VrKir, dtKir.Data{partsKir(n)}.Vr);
    VfKir = vertcat(VfKir, dtKir.Data{partsKir(n)}.Vf);
    VsKir = vertcat(VsKir, dtKir.Data{partsKir(n)}.Vs);
    VtKir = vertcat(VtKir, dtKir.Data{partsKir(n)}.Vt);
    XKir = vertcat(XKir, dtKir.Data{partsKir(n)}.X);
    YKir = vertcat(YKir, dtKir.Data{partsKir(n)}.Y);
    actstKir = vertcat(actstKir, dtKir.Data{partsKir(n)}.actState);
end

clear i nSeqs path
[params] = GetParams();
[locs, ~, cmhSong, thr] = SortSpikes(VrKir, actstKir, params);
[locsF, pksIF, pksEF, ~] = CullSpikesBasedOnVf(VrKir, locs, cmhSong, thr, VfKir, params);
[locsF, pksIF, pksEF, ~] = TemplateCompSpike(VrKir, locsF, pksIF, pksEF ,cmhSong, thr, params);
[~,FBoutsKir] = GetFbouts(actstKir, VfKir, locsF, pksIF, pksEF, params);

[params] = GetParams();
[locs, ~, cmhSong, thr] = SortSpikes(VrCtr, actstCtr, params);
[locsF, pksIF,pksEF, ~] = CullSpikesBasedOnVf(VrCtr, locs, cmhSong, thr, VfCtr, params);
[locsF, pksIF, pksEF, ~] = TemplateCompSpike(VrCtr,locsF,pksIF, pksEF, cmhSong, thr, params);
[~,FBoutsCtr] = GetFbouts(actstCtr, VfCtr, locsF, pksIF,pksEF, params);


addpath('cbrewer')
ncolor = 50;
bss = -10;
[cmap]=cbrewer('div', 'PiYG', ncolor, 'PiYG');

figure,
subplot(2,3,1)
plot(XKir, YKir,'k', 'linewidth', 1.5)
axis square
axis([-40 40 -40 40])
set(gca,'xtick',[-40 -20 0 20 40]); 
set(gca,'ytick',[-40 -20 0 20 40]);
xlabel('X (mm)')
ylabel('Y (mm)')
title(num2str(fk))
subplot(2,3,2)
for k = 1 : length(FBoutsKir)
    hold on
    plot(XKir(FBoutsKir{k}), YKir(FBoutsKir{k}), 'b', 'linewidth', 1.5)
end
axis square
axis([-40 40 -40 40])
set(gca,'xtick',[-40 -20 0 20 40]); 
set(gca,'ytick',[-40 -20 0 20 40]);
xlabel('X (mm)')
ylabel('Y (mm)')

subplot(2,3,3)
for k = 1 : length(FBoutsKir)
    vrr = VrKir(FBoutsKir{k});
    vff = VfKir(FBoutsKir{k});
    vss = VsKir(FBoutsKir{k});
    xaux = XKir(FBoutsKir{k});
    yaux = YKir(FBoutsKir{k});
    vtaux = VtKir(FBoutsKir{k});
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
    plot(XKir(FBoutsKir{k}), YKir(FBoutsKir{k}), 'color',cmap(min(ncolor,floor(max(1,1*strt+bss))),:),...
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

subplot(2,3,4)
plot(XCtr, YCtr,'k', 'linewidth', 1.5)
axis square
axis([-40 40 -40 40])
set(gca,'xtick',[-40 -20 0 20 40]); 
set(gca,'ytick',[-40 -20 0 20 40]);
xlabel('X (mm)')
ylabel('Y (mm)')
title(num2str(fc))
subplot(2,3,5)
for k = 1 : length(FBoutsCtr)
    hold on
    plot(XCtr(FBoutsCtr{k}), YCtr(FBoutsCtr{k}), 'b', 'linewidth', 1.5)
end
axis square
axis([-40 40 -40 40])
set(gca,'xtick',[-40 -20 0 20 40]); 
set(gca,'ytick',[-40 -20 0 20 40]);
xlabel('X (mm)')
ylabel('Y (mm)')

subplot(2,3,6)
for k = 1 : length(FBoutsCtr)
    vrr = VrCtr(FBoutsCtr{k});
    vff = VfCtr(FBoutsCtr{k});
    vss = VsCtr(FBoutsCtr{k});
    xaux = XCtr(FBoutsCtr{k});
    yaux = YCtr(FBoutsCtr{k});
    vtaux = VtCtr(FBoutsCtr{k});
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
        plot(XCtr(FBoutsCtr{k}), YCtr(FBoutsCtr{k}), 'color',cmap(min(ncolor,floor(max(1,1*strt+bss))),:),...
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
