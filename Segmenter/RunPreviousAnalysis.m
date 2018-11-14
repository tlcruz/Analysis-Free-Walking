%% Light vs Dark
clear
clc
[params] = GetParams();
path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Dark\Dark WTTB\';
[STRALL1, DSTALL1, NSTRALL1, PSSALL1, NPSSALL1, pTypes1] = GetStrAndVisInf(path, params);
path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Windmill\NG-RG\';
[STRALL2, DSTALL2, NSTRALL2, PSSALL2, NPSSALL2, pTypes2] = GetStrAndVisInf(path, params);
params.pTypes = pTypes1;
params.lthr = 350;
[MSTR1, SEMSTR1, MF1, NF1] = GetGMSEM(STRALL1, DSTALL1, NSTRALL1, params);
params.pTypes = pTypes2;
[MSTR2, SEMSTR2, MF2, NF2] = GetGMSEM(STRALL2, DSTALL2, NSTRALL2, params);
figure,
subplot(1,2,1)
errorbar(1, MSTR1(1), SEMSTR1(1), 'o', 'color', [0 0 0.6], 'markersize', 10, 'markerfacecolor',[0 0 0.6]);
hold on
errorbar(2, MSTR2(1), SEMSTR2(1), 'o', 'color', [1 0 0], 'markersize', 10, 'markerfacecolor',[1 0 0]);
axis([0 5 15 55])

%% Random Dots WT
clear
clc
[params] = GetParams();
path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Density\WTTB 5 Deg\';
[STRALL1, DSTALL1, NSTRALL1, PSSALL1, NPSSALL1, pTypes1] = GetStrAndVisInf(path, params);
params.pTypes = pTypes1;
%%
params.lthr = 350;
[MSTR1, SEMSTR1, MF1, NF1] = GetGMSEM(STRALL1, DSTALL1, NSTRALL1, params);
params.vthr = 400;
[VI1, NVI1, MVI1, SEMVI1] = PSStoVisInf(PSSALL1, NPSSALL1, params);

figure,
subplot(1,2,1)
inds = [1 5 7 3];
errorbar([1 2 3 4], MSTR1(inds), SEMSTR1(inds), 'o', 'color', [0 0 0.6], 'markersize', 10, 'markerfacecolor',[0 0 0.6]);
axis([0 5 15 55])
subplot(1,2,2)
indsvi = [1 3 4 2];
errorbar([1 2 3 4], MVI1(indsvi), SEMVI1(indsvi), 'o', 'color', [0 0 0.6], 'markersize', 10, 'markerfacecolor',[0 0 0.6]);
axis([0 5 -0.1 0.35])
%%
% Random Dots T4T5
clear
clc
ph = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\';
[params] = GetParams();
path = [ph 'Bilateral Stimulus Size\VTR39\VTR39FLPKir\'];
[STRALL1, DSTALL1, NSTRALL1, PSSALL1, NPSSALL1, pTypes1] = GetStrAndVisInf(path, params);

path = [ph 'Bilateral Stimulus Size\VTR39\VTFLPKir\'];
[STRALL2, DSTALL2, NSTRALL2, PSSALL2, NPSSALL2, pTypes2] = GetStrAndVisInf(path, params);

path = [ph 'Bilateral Stimulus Size\VTR39\VTR39Kir\'];
[STRALL3, DSTALL3, NSTRALL3, PSSALL3, NPSSALL3, pTypes3] = GetStrAndVisInf(path, params);

%%
params.pTypes = pTypes1;
params.lthr = 350;
[MSTR1, SEMSTR1, MF1, NF1] = GetGMSEM(STRALL1, DSTALL1, NSTRALL1, params);
[MSTR2, SEMSTR2, MF2, NF2] = GetGMSEM(STRALL2, DSTALL2, NSTRALL2, params);
[MSTR3, SEMSTR3, MF3, NF3] = GetGMSEM(STRALL3, DSTALL3, NSTRALL3, params);

params.vthr = 400;
[VI1, NVI1, MVI1, SEMVI1] = PSStoVisInf(PSSALL1, NPSSALL1, params);
[VI2, NVI2, MVI2, SEMVI2] = PSStoVisInf(PSSALL2, NPSSALL2, params);
[VI3, NVI3, MVI3, SEMVI3] = PSStoVisInf(PSSALL3, NPSSALL3, params);

%
th = 100;
inds = [3 5 7 1];
figure,
subplot(1,2,1)
errorbar([1 2 3 4], MSTR1(inds), SEMSTR1(inds), 'o', 'color', [1 0 0], 'markersize', 10, 'markerfacecolor',[1 0 0]);
hold on
errorbar([1 2 3 4], MSTR2(inds), SEMSTR2(inds), 'o', 'color', [0 0 0.6], 'markersize', 10, 'markerfacecolor',[0 0 0.6]);
hold on
errorbar([1 2 3 4], MSTR3(inds), SEMSTR3(inds), 'o', 'color', [0 0 1], 'markersize', 10, 'markerfacecolor',[0 0 1]);

% for i = 1 : 4
%    aa1 = MF1(inds(i),find(NF1(inds(i),:) > th));
%    scatter(i*ones(length(aa1),1)'+0.1, aa1, 0.05*NF1(find(NF1(inds(i),:) > th)) ,[1 0 0])
%    aa2 = MF2(inds(i),find(NF2(inds(i),:) > th));
%    scatter(i*ones(length(aa2),1)'-0.1, aa2, 0.05*(1+NF2(find(NF2(inds(i),:) > th))) ,[0 0 0.6])
%    aa3 = MF3(inds(i),find(NF3(inds(i),:) > th));
%    scatter(i*ones(length(aa3),1)', aa3, 0.05*NF3(find(NF3(inds(i),:) > th)) ,[0 0 1])
% end
axis([0 5 20 55])
set(gca,'xtick',[1 2 3 4]);
set(gca,'xticklabel',{'1º','2.5º','5º','10º'});
xlabel('Dot Size')
ylabel('Straightness (a.u.)')
subplot(1,2,2)
indsvi = [2 3 4 1];
errorbar([1 2 3 4], MVI1(indsvi), SEMVI1(indsvi), 'o', 'color', [1 0 0], 'markersize', 10, 'markerfacecolor',[1 0 0]);
hold on
errorbar([1 2 3 4], MVI2(indsvi), SEMVI2(indsvi), 'o', 'color', [0 0 0.6], 'markersize', 10, 'markerfacecolor',[0 0 0.6]);
hold on
errorbar([1 2 3 4], MVI3(indsvi), SEMVI3(indsvi), 'o', 'color', [0 0 1], 'markersize', 10, 'markerfacecolor',[0 0 1]);
% for i = 1 : 4
%    aa1 = VI1(indsvi(i),find(NVI1(indsvi(i),:) > th));
%    scatter(i*ones(length(aa1),1)'+0.1, aa1, 0.01*NVI1(find(NVI1(indsvi(i),:) > th)) ,[1 0 0])
%    aa2 = VI2(indsvi(i),find(NVI2(indsvi(i),:) > th));
%    scatter(i*ones(length(aa2),1)'-0.1, aa2, 0.01*(1+NVI2(find(NVI2(indsvi(i),:) > th))) ,[0 0 0.6])
%    aa3 = VI3(indsvi(i),find(NVI3(indsvi(i),:) > th));
%    scatter(i*ones(length(aa3),1)', aa3, 0.01*NVI1(find(NVI3(indsvi(i),:) > th)) ,[0 0 1])
% end
axis([0 5 -0.05 0.35])
set(gca,'xtick',[1 2 3 4]); 
set(gca,'xticklabel',{'1º','2.5º','5º','10º'});
xlabel('Dot Size')
ylabel('Visual Influence (a.u.)')

%%
disp(' ')
disp(' ')
disp(' ')
th = 100;

for j = 1:2: length(pTypes1)
    v1 = MF1(j,find(NF1(j,:) > th));
    v1(v1 == 0 | isnan(v1)) = [];
    v2 = MF3(j,find(NF3(j,:) > th));
    v2(v2 == 0 | isnan(v2)) = [];
    p = mwwtest(v1,v2,0);
    disp([pTypes1{j} ': ' num2str(p.p(1))])
end
disp(' ')
% clc
for j = 1:(length(pTypes1)/2)
    v1 = VI1(j,find(NVI1(j,:) > th));
    v1(v1 == 0 | isnan(v1)) = [];
    v2 = VI3(j,find(NVI3(j,:) > th));
    v2(v2 == 0 | isnan(v2)) = [];
    p = mwwtest(v1,v2,0);
    disp([pTypes1{2*j} ': ' num2str(p.p(1))])
end

disp(' ')
disp(' ')

for j = 1:2: length(pTypes1)
    v1 = MF1(j,find(NF1(j,:) > th));
    v1(v1 == 0 | isnan(v1)) = [];
    v2 = MF2(j,find(NF2(j,:) > th));
    v2(v2 == 0 | isnan(v2)) = [];
    p = mwwtest(v1,v2,0);
    disp([pTypes1{j} ': ' num2str(p.p(1))])
end
disp(' ')
% clc
for j = 1:(length(pTypes1)/2)
    v1 = VI1(j, find(NVI1(j,:) > th));
    v1(v1 == 0 | isnan(v1)) = [];
    v2 = VI2(j, find(NVI2(j,:) > th));
    v2(v2 == 0 | isnan(v2)) = [];
    p = mwwtest(v1,v2,0);
    disp([pTypes1{2*j} ': ' num2str(p.p(1))])
end























%% Light vs Dark

clear
clc
ph = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\';
[params] = GetParams();
path = [ph 'Head Movement Bilateral Stimulus size\SS02552\SS02552Kir\'];
[STRALL1, DSTALL1, NSTRALL1, PSSALL1, NPSSALL1, pTypes1] = GetStrAndVisInf(path, params);

path = [ph 'Head Movement Bilateral Stimulus size\SS02552\SS02552Cnt\'];
[STRALL2, DSTALL2, NSTRALL2, PSSALL2, NPSSALL2, pTypes2] = GetStrAndVisInf(path, params);

path = [ph 'Head Movement Bilateral Stimulus size\EmptySplitKir\'];
[STRALL3, DSTALL3, NSTRALL3, PSSALL3, NPSSALL3, pTypes3] = GetStrAndVisInf(path, params);

%
params.pTypes = pTypes1;
params.lthr = 350;
[MSTR1, SEMSTR1, MF1, NF1] = GetGMSEM(STRALL1, DSTALL1, NSTRALL1, params);
[MSTR2, SEMSTR2, MF2, NF2] = GetGMSEM(STRALL2, DSTALL2, NSTRALL2, params);
[MSTR3, SEMSTR3, MF3, NF3] = GetGMSEM(STRALL3, DSTALL3, NSTRALL3, params);

params.vthr = 400;
[VI1, NVI1, MVI1, SEMVI1] = PSStoVisInf(PSSALL1, NPSSALL1, params);
[VI2, NVI2, MVI2, SEMVI2] = PSStoVisInf(PSSALL2, NPSSALL2, params);
[VI3, NVI3, MVI3, SEMVI3] = PSStoVisInf(PSSALL3, NPSSALL3, params);

%
th = 100;
inds = [2 3 1];
figure,
subplot(1,2,1)
errorbar([1 2 3], MSTR1(inds), SEMSTR1(inds), 'o', 'color', [1 0 0], 'markersize', 10, 'markerfacecolor',[1 0 0]);
hold on
errorbar([1 2 3], MSTR2(inds), SEMSTR2(inds), 'o', 'color', [0 0 0.6], 'markersize', 10, 'markerfacecolor',[0 0 0.6]);
hold on
errorbar([1 2 3], MSTR3(inds), SEMSTR3(inds), 'o', 'color', [0 0 1], 'markersize', 10, 'markerfacecolor',[0 0 1]);

% for i = 1 : 4
%    aa1 = MF1(inds(i),find(NF1(inds(i),:) > th));
%    scatter(i*ones(length(aa1),1)'+0.1, aa1, 0.05*NF1(find(NF1(inds(i),:) > th)) ,[1 0 0])
%    aa2 = MF2(inds(i),find(NF2(inds(i),:) > th));
%    scatter(i*ones(length(aa2),1)'-0.1, aa2, 0.05*(1+NF2(find(NF2(inds(i),:) > th))) ,[0 0 0.6])
%    aa3 = MF3(inds(i),find(NF3(inds(i),:) > th));
%    scatter(i*ones(length(aa3),1)', aa3, 0.05*NF3(find(NF3(inds(i),:) > th)) ,[0 0 1])
% end
axis([0 5 20 55])
set(gca,'xtick',[1 2 3]);
set(gca,'xticklabel',{'1º','5º','10º'});
xlabel('Dot Size')
ylabel('Straightness (a.u.)')
subplot(1,2,2)
indsvi = [2 3 1];
errorbar([1 2 3], MVI1(indsvi), SEMVI1(indsvi), 'o', 'color', [1 0 0], 'markersize', 10, 'markerfacecolor',[1 0 0]);
hold on
errorbar([1 2 3], MVI2(indsvi), SEMVI2(indsvi), 'o', 'color', [0 0 0.6], 'markersize', 10, 'markerfacecolor',[0 0 0.6]);
hold on
errorbar([1 2 3], MVI3(indsvi), SEMVI3(indsvi), 'o', 'color', [0 0 1], 'markersize', 10, 'markerfacecolor',[0 0 1]);
% for i = 1 : 4
%    aa1 = VI1(indsvi(i),find(NVI1(indsvi(i),:) > th));
%    scatter(i*ones(length(aa1),1)'+0.1, aa1, 0.01*NVI1(find(NVI1(indsvi(i),:) > th)) ,[1 0 0])
%    aa2 = VI2(indsvi(i),find(NVI2(indsvi(i),:) > th));
%    scatter(i*ones(length(aa2),1)'-0.1, aa2, 0.01*(1+NVI2(find(NVI2(indsvi(i),:) > th))) ,[0 0 0.6])
%    aa3 = VI3(indsvi(i),find(NVI3(indsvi(i),:) > th));
%    scatter(i*ones(length(aa3),1)', aa3, 0.01*NVI1(find(NVI3(indsvi(i),:) > th)) ,[0 0 1])
% end
axis([0 5 -0.05 0.35])
set(gca,'xtick',[1 2 3]); 
set(gca,'xticklabel',{'1º','5º','10º'});
xlabel('Dot Size')
ylabel('Visual Influence (a.u.)')