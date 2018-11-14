%% StrPlotsFig3
clear
clc
[params] = GetParams();
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\WT TB\';
path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Density\WTTB 5 Deg\';
[STRALL1, DSTALL1, NSTRALL1, PSSALL1, NPSSALL1, pTypes1] = GetStrAndVisInf(path, params);
params.pTypes = pTypes1;
%%
params.lthr = 350;
[MSTR1, SEMSTR1, MF1, NF1] = GetGMSEM(STRALL1, DSTALL1, NSTRALL1, params);
params.vthr = 0;
[VI1, NVI1, MVI1, SEMVI1] = PSStoVisInf(PSSALL1, NPSSALL1, params);
% vecSt = [5 6 7 8 9 10 1 2]; %size
% vecVi = [3 4 5 1]; %size
% cmp = hot(10);
vecSt = [5 6 7 8 3 4]; %density
vecVi = [3 4 2]; %density
cmp = winter(6);
%%
X = [];
Y = [];
figure,
hold on
for n = 1 : floor(size(MF1,1)/2)-1
    st = MF1(vecSt(2*n-1),:);
    nst = NF1(vecSt(2*n-1),:);
    vi = VI1(vecVi(n),:);
    nvi = NVI1(vecVi(n),:);
    inds = find(vi~=0 & nvi > params.vthr & nst > params.lthr);
    st = st(inds);
    vi = vi(inds);
    
    scatter(vi, st, 100, cmp(n+2,:), 'filled')
    X = horzcat(X, vi);
    Y = horzcat(Y, st);
end
axis([-0.05 0.4 10 70])
[R,P] = corrcoef(X,Y);
disp(['C = ' num2str(R(1,2)) '  P = ' num2str(P(1,2))])
xlabel('Visual Influence (a.u.)')
ylabel('Straightness (a.u.)')
title(['C = ' num2str(R(1,2)) '  P = ' num2str(P(1,2))])
%%
figure,
hold on
for n = 1 : length(MVI1)-1
    errorbar(n, MVI1(vecVi(n)), SEMVI1(vecVi(n)), 'color', cmp(n+2,:), 'linewidth', 3, 'marker', 'o', 'markersize', 10, 'markerfacecolor',cmp(n+2,:));
end
% axis([0 5 -0.05 0.35])
% set(gca,'xtick',[1 2 3 4]); 
% set(gca,'xticklabel',{'1º','2.5º','5º','10º'});
% xlabel('Dot Size')
axis([0 4 -0.05 0.35])
set(gca,'xtick',[1 2 3]); 
set(gca,'xticklabel',{'5ºLD','5ºMD','5ºHD'});
xlabel('Dot Density')
ylabel('Visual Influence (a.u.)')
%%
figure,
hold on
for n = 1 : length(MSTR1)/2 -1
    errorbar(n, MSTR1(vecSt(2*n-1)), SEMSTR1(vecSt(2*n-1)), 'color', cmp(n+2,:), 'linewidth', 3, 'marker', 'o', 'markersize', 10, 'markerfacecolor',cmp(n+2,:));
end
axis([0 5 20 55])
set(gca,'xtick',[1 2 3 4]); 
set(gca,'xticklabel',{'1º','2.5º','5º','10º'});
xlabel('Dot Size')
% axis([0 4 20 55])
% set(gca,'xtick',[1 2 3]); 
% set(gca,'xticklabel',{'5ºLD','5ºMD','5ºHD'});
% xlabel('Dot Density')
ylabel('Straightness (a.u.)')
%% SIZE
% clc
disp('Straightness:')
a1 = MF1(vecSt(1),:);
a2 = MF1(vecSt(3),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º - 2.5º: ' num2str(p.p(2))])
a1 = MF1(vecSt(1),:);
a2 = MF1(vecSt(5),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º - 5º: ' num2str(p.p(2))])
a1 = MF1(vecSt(1),:);
a2 = MF1(vecSt(7),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º - 10º: ' num2str(p.p(2))])
a1 = MF1(vecSt(3),:);
a2 = MF1(vecSt(5),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['2.5º - 5º: ' num2str(p.p(2))])
a1 = MF1(vecSt(3),:);
a2 = MF1(vecSt(7),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['2.5º - 10º: ' num2str(p.p(2))])
a1 = MF1(vecSt(5),:);
a2 = MF1(vecSt(7),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['5º - 10º: ' num2str(p.p(2))])
disp(' ')

disp('Visual Influence:')
a1 = VI1(vecVi(1),:);
a2 = VI1(vecVi(2),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º - 2.5º: ' num2str(p.p(2))])
a1 = VI1(vecVi(1),:);
a2 = VI1(vecVi(3),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º - 5º: ' num2str(p.p(2))])
a1 = VI1(vecVi(1),:);
a2 = VI1(vecVi(4),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º - 10º: ' num2str(p.p(2))])
a1 = VI1(vecVi(2),:);
a2 = VI1(vecVi(3),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['2.5º - 5º: ' num2str(p.p(2))])
a1 = VI1(vecVi(2),:);
a2 = VI1(vecVi(4),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['2.5º - 10º: ' num2str(p.p(2))])
a1 = VI1(vecVi(3),:);
a2 = VI1(vecVi(4),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['5º - 10º: ' num2str(p.p(2))])
disp(' ')


%% Density
vecSt = [1 2 5 6 7 8 3 4]; %density
vecVi = [1 3 4 2]; %density
clc
disp('Straightness:')
a1 = MF1(vecSt(1),NF1(vecSt(1),:)>params.lthr);
a2 = MF1(vecSt(3),NF1(vecSt(3),:)>params.lthr);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º - 5Lº: ' num2str(p.p(2))])
a1 = MF1(vecSt(1),NF1(vecSt(1),:)>params.lthr);
a2 = MF1(vecSt(5),NF1(vecSt(5),:)>params.lthr);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º - 5Mº: ' num2str(p.p(2))])
a1 = MF1(vecSt(1),NF1(vecSt(1),:)>params.lthr);
a2 = MF1(vecSt(7),NF1(vecSt(7),:)>params.lthr);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º - 5Hº: ' num2str(p.p(2))])
a1 = MF1(vecSt(3),NF1(vecSt(3),:)>params.lthr);
a2 = MF1(vecSt(5),NF1(vecSt(5),:)>params.lthr);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['5Lº - 5Mº: ' num2str(p.p(2))])
a1 = MF1(vecSt(3),NF1(vecSt(3),:)>params.lthr);
a2 = MF1(vecSt(7),NF1(vecSt(7),:)>params.lthr);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['5Lº - 5Hº: ' num2str(p.p(2))])
a1 = MF1(vecSt(5),NF1(vecSt(5),:)>params.lthr);
a2 = MF1(vecSt(7),NF1(vecSt(7),:)>params.lthr);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['5Mº - 5Hº: ' num2str(p.p(2))])
disp(' ')

disp('Visual Influence:')
a1 = VI1(vecVi(1),:);
a2 = VI1(vecVi(2),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º - 5Lº: ' num2str(p.p(2))])
a1 = VI1(vecVi(1),:);
a2 = VI1(vecVi(3),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º - 5Mº: ' num2str(p.p(2))])
a1 = VI1(vecVi(1),:);
a2 = VI1(vecVi(4),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º - 5Hº: ' num2str(p.p(2))])
a1 = VI1(vecVi(2),:);
a2 = VI1(vecVi(3),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['5Lº - 5Mº: ' num2str(p.p(2))])
a1 = VI1(vecVi(2),:);
a2 = VI1(vecVi(4),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['5Lº - 5Hº: ' num2str(p.p(2))])
a1 = VI1(vecVi(3),:);
a2 = VI1(vecVi(4),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['5Mº - 5Hº: ' num2str(p.p(2))])
disp(' ')
 


%%
nbts = [];
for k1 = 1 : size(NPSSALL1,1)
    for k2 = 1 : size(NPSSALL1,2)
        nbts = vertcat(nbts, length(NPSSALL1{k1,k2}));
    end
end
sum(nbts)