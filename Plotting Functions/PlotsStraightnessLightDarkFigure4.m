%% Straightness plots light/dark figure 1
clear
clc
% load dark and light data
% calculare straightness and visual influence
[params] = GetParams();
path = '\';
[STRALL1, DSTALL1, NSTRALL1, PSSALL1, VRSSALL1, NPSSALL1, pTypes1, ANGDALL1, VTALL1]  = GetStrAndVisInf(path, params);
path = '\';
[STRALL2, DSTALL2, NSTRALL2, PSSALL2, VRSSALL2, NPSSALL2, pTypes2, ANGDALL2, VTALL2]  = GetStrAndVisInf(path, params);
params.pTypes = pTypes1;    

%%
% calculate straigness weighted grand mean an standard error
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


%%
% calculate straigness weighted grand mean an standard error
params.lthr = 350;
params.pTypes = pTypes1;
[MSTR1, SEMSTR1, MF1, NF1] = GetGMSEM(ANGDALL1, DSTALL1, NSTRALL1, params);
params.pTypes = pTypes2;
[MSTR2, SEMSTR2, MF2, NF2] = GetGMSEM(ANGDALL2, DSTALL2, NSTRALL2, params);

figure,
errorbar(1, MSTR1(1), SEMSTR1(1), 'o', 'color', [0 0 0.6], 'markersize', 10, 'markerfacecolor',[0 0 0.6]);
hold on
errorbar(2, MSTR2(1), SEMSTR2(1), 'o', 'color', [1 0 0], 'markersize', 10, 'markerfacecolor',[1 0 0]);
axis([0 3 1 5])
xt={'Dark' ; 'Light'} ; 
set(gca,'xtick',[1 2]); 
set(gca,'xticklabel',xt);
ylabel('Ang. Dev. (º/mm)')
[p] = mwwtest(MF1,MF2(1,:),0);
disp(['Dark - Light: pVal = ' num2str(p.p(2))])

%% Rotation Comparison
% calculate angular deflection weighted grand mean an standard error
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
