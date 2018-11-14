function [] = ZMap(mMap1, mMap2, stdMap1)

CRange=[-5 5];
aux1 = abs(mMap1 - mMap2);

inds = find(aux1 ~= 0);
aux1(inds) = aux1(inds)./stdMap1(inds);
MaskMat=aux1;
MaskMat(MaskMat>CRange(2))=CRange(2);
MaskMat(MaskMat==0)=1000;
vfStart=-5;
vfbin=1;
vfEnd=50;
vrStart=-800;
vrbin=10;
vrEnd=800;
vrx=vrStart:vrbin:vrEnd;
vfx=vfStart:vfbin:vfEnd;

figure,
hold on
h=imagesc(MaskMat);
colormap redblue
set(gca,'Ytick',1:10:length(vfx))
set(gca,'YtickLabel',vfx(1):10*vfbin:vfx(end))
set(gca,'Xtick',1:10:length(vrx))
set(gca,'XtickLabel',vrx(1):10*vrbin:vrx(end))
set(gca,'box','off')
colorbar
xlabel('Rot Speed [deg/s]')
ylabel('Forw Speed [mm/s]')
colordata = colormap;
colordata = [colordata; 1 1 1];
colormap(colordata);
axis square
caxis([CRange(1) CRange(2)+0.1])

hold on

aux1 = abs(mMap1 - mMap2);
MaskMat=aux1;
MaskMat(MaskMat>CRange(2))=CRange(2);
MaskMat(MaskMat==0)=1000;
aux1 = MaskMat;
inds = find(aux1 ~= 0);
aux1(inds) = aux1(inds)./stdMap1(inds)-1.5;
aux1(abs(aux1) > 5 & abs(aux1) < 500) = 5;
aux1(abs(aux1) > 4 & abs(aux1) < 5) = 4;
aux1(abs(aux1) > 3 & abs(aux1) < 4) = 3;
aux1(abs(aux1) > 2 & abs(aux1) < 3) = 2;
aux1(abs(aux1) > 1 & abs(aux1) < 2) = 1;
aux1(find(aux1 < 1)) = 10;

aux1(aux1>500) = 0;

contour(1:1:length(vrx), 1:1:length(vfx), aux1, 1, 'k', 'linewidth', 2);

%%
end