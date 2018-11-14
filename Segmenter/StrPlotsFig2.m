%% Straightness plots figure 1
clear
clc
[params] = GetParams();
path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Dark\Dark WTTB\';
[STRALL1, DSTALL1, NSTRALL1, PSSALL1, NPSSALL1, pTypes1, ANGDALL1, VTALL1] = GetStrAndVisInf(path, params);
path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\WT TB\';
[STRALL2, DSTALL2, NSTRALL2, PSSALL2, NPSSALL2, pTypes2, ANGDALL2, VTALL2] = GetStrAndVisInf(path, params);
params.pTypes = pTypes1;
%%
params.lthr = 350;
params.pTypes = pTypes1;
[MSTR1, SEMSTR1, MF1, NF1] = GetGMSEM(STRALL1, DSTALL1, NSTRALL1, params);
params.pTypes = pTypes2;
[MSTR2, SEMSTR2, MF2, NF2] = GetGMSEM(STRALL2, DSTALL2, NSTRALL2, params);

figure,
errorbar(1, MSTR1(1), SEMSTR1(1), 'o', 'color', [0 0 0.6], 'markersize', 10, 'markerfacecolor',[0 0 0.6]);
hold on
errorbar(2, MSTR2(1), SEMSTR2(1), 'o', 'color', [1 0 0], 'markersize', 10, 'markerfacecolor',[1 0 0]);
axis([0 3 15 55])
xt={'Dark' ; 'Light'} ; 
set(gca,'xtick',[1 2]); 
set(gca,'xticklabel',xt);
ylabel('Straightness (a.u.)')
[p] = mwwtest(MF1,MF2(1,:),0);
disp(['Dark - Light: pVal = ' num2str(p.p(2))])



%% Rotation Comparison
params.lthr = 350;
thrv = 0;
params.pTypes = pTypes1;
[MSTR1, SEMSTR1, MF1, NF1] = GetGMSEM(STRALL1, DSTALL1, NSTRALL1, params);
params.pTypes = pTypes2;
[MSTR2, SEMSTR2, MF2, NF2] = GetGMSEM(STRALL2, DSTALL2, NSTRALL2, params);
params.pTypes = pTypes1;
[MVR1, SEMVR1, MR1, NR1] = GetGMSEM(ANGDALL1, DSTALL1, NSTRALL1, params);
params.pTypes = pTypes2;
[MVR2, SEMVR2, MR2, NR2] = GetGMSEM(ANGDALL2, DSTALL2, NSTRALL2, params);


MF1 = MF1(NF1>thrv);
MR1 = MR1(NR1>thrv);
MF210 = MF2(1,:);
MF210 = MF210(NF2(1,:)>thrv);
MR210 = MR2(1,:);
MR210 = MR210(NR2(1,:)>thrv);

figure,
hold on
scatter(MF1(MF1~=0), MR1(MR1~=0), 100, 'k', 'filled')
hold on
scatter(MF210(MF210~=0), MR210(MR210~=0), 100, 'r', 'filled')
axis([10 75 0 7])
xlabel('Straightness(a.u.)')
ylabel('Angular Deviations (rad/mm)')
[R,P] = corrcoef(vertcat(MF1(MF1~=0)', MF210(MF210~=0)'), vertcat(MR1(MR1~=0)', MR210(MR210~=0)'));
disp(['Corr AngDev - Str: Val = ' num2str(R(1,2)) ', pVal = ' num2str(P(1,2))])

%% Translation Comparison
params.lthr = 350;
thrv = 0;
params.pTypes = pTypes1;
[MSTR1, SEMSTR1, MF1, NF1] = GetGMSEM(STRALL1, DSTALL1, NSTRALL1, params);
params.pTypes = pTypes2;
[MSTR2, SEMSTR2, MF2, NF2] = GetGMSEM(STRALL2, DSTALL2, NSTRALL2, params);
params.pTypes = pTypes1;
[MVT1, SEMVR1, MT1, NT1] = GetGMSEM(VTALL1, DSTALL1, NSTRALL1, params);
params.pTypes = pTypes2;
[MVT2, SEMVR2, MT2, NT2] = GetGMSEM(VTALL2, DSTALL2, NSTRALL2, params);

MF1 = MF1(NF1>thrv);
MT1 = MT1(NT1>thrv);
MF210 = MF2(1,:);
MF210 = MF210(NF2(1,:)>thrv);
MT210 = MT2(1,:);
MT210 = MT210(NT2(1,:)>thrv);



figure,
hold on
scatter(MF1(MF1~=0), MT1(MT1~=0), 100, 'k', 'filled')
hold on
scatter(MF210(MF210~=0), MT210(MT210~=0), 100, 'r', 'filled')
axis([10 70 0 35])
xlabel('Straightness(a.u.)')
ylabel('Translational Speed (mm/s)')
[R,P] = corrcoef(vertcat(MF1(MF1~=0)', MF210(MF210~=0)'), vertcat(MT1(MT1~=0)', MT210(MT210~=0)'));
disp(['Corr AngDev - Str: Val = ' num2str(R(1,2)) ', pVal = ' num2str(P(1,2))])

%%
nbts = [];
for k1 = 1 %: size(NSTRALL2,1)
    for k2 = 1 : size(NSTRALL2,2)
        nbts = vertcat(nbts, length(NSTRALL2{k1,k2}));
    end
end
sum(nbts)
%%
figure,
hold on
plot([0 5], [0.5 0.5], '--g', 'linewidth', 2)
inds = [5 6 7 8 9 10 1 2 3 4];
for k = 1 : size(PSSALL2,1)
    mps = [];
    nps = [];
    for n = 1 : size(PSSALL2,2)
        pss = PSSALL2{inds(k),n};
        npss = NPSSALL2{inds(k),n};
        if ~isempty(pss'*npss/sum(npss)) && ~isnan(pss'*npss/sum(npss))
            mps = vertcat(mps, pss'*npss/sum(npss));
            nps = vertcat(nps, sum(npss));
        end
    end
    gm = mps'*nps/sum(nps);
    sem = (((mps-gm).*(mps-gm))'*nps/sum(nps));
    
    errorbar(floor((k+1)/2), gm, sem, 'ok', 'markersize', 10)
    axis([0.5 4.5 0 1])
end
xlabel('Dot Size (º)')
ylabel('Probability Same Side Turning')

%%
params.pTypes = pTypes1;
params.vthr = 0;
[VI1, NVI1, MVI1, SEMVI1] = PSStoVisInf(PSSALL1, NPSSALL1, params);
params.pTypes = pTypes2;
[VI2, NVI2, MVI2, SEMVI2] = PSStoVisInf(PSSALL2, NPSSALL2, params);

%%
vecSt = [5 6 7 8 9 10 1 2]; %size
vecVi = [3 4 5 1]; %size
cmp = hot(10);
figure,
hold on
errorbar(0, MSTR1, SEMSTR1, 'k', 'linewidth', 3, 'marker', 'o', 'markersize', 10, 'markerfacecolor','k');

for n = 1 : length(MSTR2)/2 -1
    errorbar(n, MSTR2(vecSt(2*n-1)), SEMSTR2(vecSt(2*n-1)), 'color', cmp(n+2,:), 'linewidth', 3, 'marker', 'o', 'markersize', 10, 'markerfacecolor',cmp(n+2,:));
end
axis([-1 5 20 55])
set(gca,'xtick',[0 1 2 3 4]); 
set(gca,'xticklabel',{'Dark', '1º','2.5º','5º','10º'});
xlabel('Dot Size')
ylabel('Straightness (a.u.)')




%%
vecSt = [5 6 7 8 9 10 1 2]; %size
vecVi = [3 4 5 1]; %size
cmp = hot(10);
figure,
hold on

errorbar(0, MVI1, SEMVI1, 'k', 'linewidth', 3, 'marker', 'o', 'markersize', 10, 'markerfacecolor','k');

for n = 1 : length(MVI2)-1
    errorbar(n, MVI2(vecVi(n)), SEMVI2(vecVi(n)), 'color', cmp(n+2,:), 'linewidth', 3, 'marker', 'o', 'markersize', 10, 'markerfacecolor',cmp(n+2,:));
end
axis([-1 5 -0.05 0.35])
set(gca,'xtick',[0 1 2 3 4]); 
set(gca,'xticklabel',{'Dark', '1º','2.5º','5º','10º'});
xlabel('Dot Size')
ylabel('Visual Influence (a.u.)')


%%
X = [];
Y = [];
figure,
hold on
for n = 1 : floor(size(MF2,1)/2)-1
    st = MF2(vecSt(2*n-1),:);
    nst = NF2(vecSt(2*n-1),:);
    vi = VI2(vecVi(n),:);
    nvi = NVI2(vecVi(n),:);
    inds = find(vi~=0 & nvi > params.vthr & nst > params.lthr);
    st = st(inds);
    vi = vi(inds);
    
    scatter(vi, st, 100, cmp(n+2,:), 'filled')
    X = horzcat(X, vi);
    Y = horzcat(Y, st);
end
st = MF1;
nst = NF1;
vi = VI1;
nvi = NVI1;
inds = find(vi~=0 & nvi > params.vthr & nst > params.lthr);
st = st(inds);
vi = vi(inds);

scatter(vi, st, 100,'k', 'filled')
X = horzcat(X, vi);
Y = horzcat(Y, st);

axis([-0.05 0.4 10 70])
[R,P] = corrcoef(X,Y);
disp(['C = ' num2str(R(1,2)) '  P = ' num2str(P(1,2))])
xlabel('Visual Influence (a.u.)')
ylabel('Straightness (a.u.)')
title(['C = ' num2str(R(1,2)) '  P = ' num2str(P(1,2))])

%%
disp('Straightness:')
a1 = MF2(vecSt(1),:);
a2 = MF1;
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['Darkº - 1º: ' num2str(p.p(2))])
a1 = MF2(vecSt(3),:);
a2 = MF1;
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['Darkº - 2.5º: ' num2str(p.p(2))])
a1 = MF2(vecSt(5),:);
a2 = MF1;
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['Darkº - 5º: ' num2str(p.p(2))])
a1 = MF2(vecSt(7),:);
a2 = MF1;
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['Darkº - 10º: ' num2str(p.p(2))])

a1 = MF2(vecSt(1),:);
a2 = MF2(vecSt(3),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º - 2.5º: ' num2str(p.p(2))])
a1 = MF2(vecSt(1),:);
a2 = MF2(vecSt(5),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º - 5º: ' num2str(p.p(2))])
a1 = MF2(vecSt(1),:);
a2 = MF2(vecSt(7),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º - 10º: ' num2str(p.p(2))])
a1 = MF2(vecSt(3),:);
a2 = MF2(vecSt(5),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['2.5º - 5º: ' num2str(p.p(2))])
a1 = MF2(vecSt(3),:);
a2 = MF2(vecSt(7),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['2.5º - 10º: ' num2str(p.p(2))])
a1 = MF2(vecSt(5),:);
a2 = MF2(vecSt(7),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['5º - 10º: ' num2str(p.p(2))])
disp(' ')

disp('Visual Influence:')
a1 = VI2(vecVi(1),:);
a2 = VI1;
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['Darkº - 1º: ' num2str(p.p(2))])
a1 = VI2(vecVi(2),:);
a2 = VI1;
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['Darkº - 2.5º: ' num2str(p.p(2))])
a1 = VI2(vecVi(3),:);
a2 = VI1;
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['Darkº - 5º: ' num2str(p.p(2))])
a1 = VI2(vecVi(4),:);
a2 = VI1;
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['Darkº - 10º: ' num2str(p.p(2))])

a1 = VI2(vecVi(1),:);
a2 = VI2(vecVi(2),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º - 2.5º: ' num2str(p.p(2))])
a1 = VI2(vecVi(1),:);
a2 = VI2(vecVi(3),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º - 5º: ' num2str(p.p(2))])
a1 = VI2(vecVi(1),:);
a2 = VI2(vecVi(4),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º - 10º: ' num2str(p.p(2))])
a1 = VI2(vecVi(2),:);
a2 = VI2(vecVi(3),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['2.5º - 5º: ' num2str(p.p(2))])
a1 = VI2(vecVi(2),:);
a2 = VI2(vecVi(4),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['2.5º - 10º: ' num2str(p.p(2))])
a1 = VI2(vecVi(3),:);
a2 = VI2(vecVi(4),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['5º - 10º: ' num2str(p.p(2))])
disp(' ')



