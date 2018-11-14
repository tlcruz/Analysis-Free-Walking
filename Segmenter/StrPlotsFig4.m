%% StrPlotsFig4
clear
clc
[params] = GetParams();
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Dark\Dark SplitT4T5Kir\';
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\SplitT4T5 Kir\';
path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Size wTB\T4T5\wTBSplitT4T5OtdKIR\';
[STRALL1, DSTALL1, NSTRALL1, PSSALL1, NPSSALL1, pTypes1] = GetStrAndVisInf(path, params);

% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Dagoo I guessrk\Dark SplitT4T5\';
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\Split T4T5 Cnt\';
path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Size wTB\T4T5\wTBSplitT4T5KIR\';
[STRALL2, DSTALL2, NSTRALL2, PSSALL2, NPSSALL2, pTypes2] = GetStrAndVisInf(path, params);
NSTRALL22 = NSTRALL2;
STRALL22 = STRALL2;
NPSSALL22 = NPSSALL2;
PSSALL22 = PSSALL2;

% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Dark\Dark EmptyKir\';
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\EmptySplit 10xKir\';
path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Size wTB\T4T5\wTBEmptySplitOtdKIR\';
[STRALL3, DSTALL3, NSTRALL3, PSSALL3, NPSSALL3, pTypes3] = GetStrAndVisInf(path, params);

NSTRALL32 = NSTRALL3;
STRALL32 = STRALL3;
NPSSALL32 = NPSSALL3;
PSSALL32 = PSSALL3;

NSTRALL12 = NSTRALL1;
STRALL12 = STRALL1;
NPSSALL12 = NPSSALL1;
PSSALL12 = PSSALL1;

%%
STRALL1 = STRALL12;
NSTRALL1 = NSTRALL12;
PSSALL1 = PSSALL12;
NPSSALL1 = NPSSALL12;
STRALL2 = STRALL22;
NSTRALL2 = NSTRALL22;
PSSALL2 = PSSALL22;
NPSSALL2 = NPSSALL22;
STRALL3 = STRALL32;
NSTRALL3 = NSTRALL32;
PSSALL3 = PSSALL32;
NPSSALL3 = NPSSALL32;
% for j = 1:8
%     NSTRALL1{j,5} = zeros(length(NSTRALL1{j,5}),1);
%     NSTRALL1{j,12} = zeros(length(NSTRALL1{j,12}),1);
%     NSTRALL1{j,18} = zeros(length(NSTRALL1{j,18}),1);
%     NSTRALL1{j,23} = zeros(length(NSTRALL1{j,23}),1);
%     NSTRALL1{j,29} = zeros(length(NSTRALL1{j,29}),1);
%     NSTRALL1{j,30} = zeros(length(NSTRALL1{j,30}),1);
%     NSTRALL1{j,31} = zeros(length(NSTRALL1{j,31}),1);
%     NSTRALL1{j,34} = zeros(length(NSTRALL1{j,34}),1);
%     NSTRALL1{j,35} = zeros(length(NSTRALL1{j,35}),1);
%     NSTRALL1{j,38} = zeros(length(NSTRALL1{j,38}),1);
%     
%     NSTRALL2{j,8} = zeros(length(NSTRALL2{j,8}),1);
%     NSTRALL2{j,16} = zeros(length(NSTRALL2{j,16}),1);
% %     NSTRALL2{j,17} = zeros(length(NSTRALL2{j,17}),1);
%     
%     NSTRALL3{j,30} = zeros(length(NSTRALL3{j,30}),1);
%     NSTRALL3{j,12} = zeros(length(NSTRALL3{j,12}),1);
%     NSTRALL3{j,31} = zeros(length(NSTRALL3{j,31}),1);
%     
%     NPSSALL1{j,5} = zeros(length(NPSSALL1{j,5}),1);
%     NPSSALL1{j,12} = zeros(length(NPSSALL1{j,12}),1);
%     NPSSALL1{j,18} = zeros(length(NPSSALL1{j,18}),1);
%     NPSSALL1{j,23} = zeros(length(NPSSALL1{j,23}),1);
%     NPSSALL1{j,29} = zeros(length(NPSSALL1{j,29}),1);
%     NPSSALL1{j,30} = zeros(length(NPSSALL1{j,30}),1);
%     NPSSALL1{j,31} = zeros(length(NPSSALL1{j,31}),1);
%     
%     NPSSALL1{j,34} = zeros(length(NPSSALL1{j,34}),1);
%     NPSSALL1{j,35} = zeros(length(NPSSALL1{j,35}),1);
%     NPSSALL1{j,38} = zeros(length(NPSSALL1{j,38}),1);
%     
%     NPSSALL2{j,8} = zeros(length(NPSSALL2{j,8}),1);
%     NPSSALL2{j,16} = zeros(length(NPSSALL2{j,16}),1);
% %     NPSSALL2{j,17} = zeros(length(NPSSALL2{j,17}),1);
% %     NPSSALL2{j,28} = zeros(length(NPSSALL2{j,28}),1);
%     NPSSALL3{j,30} = zeros(length(NPSSALL3{j,30}),1);
%     NPSSALL3{j,31} = zeros(length(NPSSALL3{j,31}),1);
%     NPSSALL3{j,12} = zeros(length(NPSSALL3{j,12}),1);
% 
% end
params.pTypes = pTypes1;
params.lthr = 350;
[MSTR1, SEMSTR1, MF1, NF1] = GetGMSEM(STRALL1, DSTALL1, NSTRALL1, params);
[MSTR2, SEMSTR2, MF2, NF2] = GetGMSEM(STRALL2, DSTALL2, NSTRALL2, params);
[MSTR3, SEMSTR3, MF3, NF3] = GetGMSEM(STRALL3, DSTALL3, NSTRALL3, params);

params.vthr = 0;
[VI1, NVI1, MVI1, SEMVI1] = PSStoVisInf(PSSALL1, NPSSALL1, params);
[VI2, NVI2, MVI2, SEMVI2] = PSStoVisInf(PSSALL2, NPSSALL2, params);
[VI3, NVI3, MVI3, SEMVI3] = PSStoVisInf(PSSALL3, NPSSALL3, params);

%
inds = [3 5 7 1];
% inds = [5 7 9 1];
% inds = 1;
figure,
subplot(1,2,1)
errorbar(1:length(inds), MSTR1(inds), SEMSTR1(inds), 'o', 'color', [1 0 0], 'markersize', 10, 'markerfacecolor',[1 0 0]);
hold on
errorbar(1:length(inds), MSTR2(inds), SEMSTR2(inds), 'o', 'color', [0 0 0.6], 'markersize', 10, 'markerfacecolor',[0 0 0.6]);
hold on
errorbar(1:length(inds), MSTR3(inds), SEMSTR3(inds), 'o', 'color', [0 0 1], 'markersize', 10, 'markerfacecolor',[0 0 1]);
% for i = 1 : 4
%    aa1 = MF1(inds(i),find(NF1(inds(i),:) > th));
%    scatter((i+0.1)*ones(length(aa1),1)', aa1, 30 ,[1 0 0])
%    aa2 = MF2(inds(i),find(NF2(inds(i),:) > th));
%    scatter((i-0.1)*ones(length(aa2),1)', aa2, 30 ,[0 0 0.6])
%    aa3 = MF3(inds(i),find(NF3(inds(i),:) > th));
%    scatter(i*ones(length(aa3),1)', aa3, 30 ,[0 0 1])
% end
axis([0 5 20 55])
set(gca,'xtick',[1 2 3 4]);
set(gca,'xticklabel',{'1º','2.5º','5º','10º'});
xlabel('Dot Size')
ylabel('Straightness (a.u.)')
subplot(1,2,2)
indsvi = [2 3 4 1];
% indsvi = [3 4 5 1];
errorbar(1:4, MVI1(indsvi), SEMVI1(indsvi), 'o', 'color', [1 0 0], 'markersize', 10, 'markerfacecolor',[1 0 0]);
hold on
errorbar(1:4, MVI2(indsvi), SEMVI2(indsvi), 'o', 'color', [0 0 0.6], 'markersize', 10, 'markerfacecolor',[0 0 0.6]);
hold on
errorbar(1:4, MVI3(indsvi), SEMVI3(indsvi), 'o', 'color', [0 0 1], 'markersize', 10, 'markerfacecolor',[0 0 1]);
% for i = 1 : 4
%    aa1 = VI1(indsvi(i),find(NVI1(indsvi(i),:) > th));
%    scatter((i+0.1)*ones(length(aa1),1)', aa1, 30 ,[1 0 0])
%    aa2 = VI2(indsvi(i),find(NVI2(indsvi(i),:) > th));
%    scatter((i-0.1)*ones(length(aa2),1)', aa2, 30 ,[0 0 0.6])
%    aa3 = VI3(indsvi(i),find(NVI3(indsvi(i),:) > th));
%    scatter(i*ones(length(aa3),1)', aa3, 30 ,[0 0 1])
% end
axis([0 5 -0.05 0.35])
% axis([0 5 0 200])
set(gca,'xtick',[1 2 3 4]); 
set(gca,'xticklabel',{'1º','2.5º','5º','10º'});
xlabel('Dot Size')
ylabel('Visual Influence (a.u.)')

%
vecSt = [3 4 5 6 7 8 1 2];
% vecSt = 1;
vecVi = [2 3 4 1];
clc
th = 250;
thk = 250;
disp('Straightness: Kir - Cnt1')
a1 = MF1(vecSt(1),find(NF1(vecSt(1),:)> thk));
a2 = MF2(vecSt(1),find(NF2(vecSt(1),:)> th));
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º Kir - Cnt1: ' num2str(p.p(2))])
a1 = MF1(vecSt(3),find(NF1(vecSt(3),:)> thk));
a2 = MF2(vecSt(3),find(NF2(vecSt(3),:)> th));
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['2.5º Kir - Cnt1: ' num2str(p.p(2))])
a1 = MF1(vecSt(5),find(NF1(vecSt(5),:)> thk));
a2 = MF2(vecSt(5),find(NF2(vecSt(5),:)> th));
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['5º Kir - Cnt1: ' num2str(p.p(2))])
a1 = MF1(vecSt(7),find(NF1(vecSt(7),:)> thk));
a2 = MF2(vecSt(7),find(NF2(vecSt(7),:)> th));
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['10º Kir - Cnt1: ' num2str(p.p(2))])
disp(' ')
disp('Straightness: Kir - Cnt2')
a1 = MF1(vecSt(1),find(NF1(vecSt(1),:)> thk));
a2 = MF3(vecSt(1),find(NF3(vecSt(1),:)> th));
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º Kir - Cnt2: ' num2str(p.p(2))])
a1 = MF1(vecSt(3),find(NF1(vecSt(3),:)> thk));
a2 = MF3(vecSt(3),find(NF3(vecSt(3),:)> th));
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['2.5º Kir - Cnt2: ' num2str(p.p(2))])
a1 = MF1(vecSt(5),find(NF1(vecSt(5),:)> thk));
a2 = MF3(vecSt(5),find(NF3(vecSt(5),:)> th));
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['5º Kir - Cnt2: ' num2str(p.p(2))])
a1 = MF1(vecSt(7),find(NF1(vecSt(7),:)> thk));
a2 = MF3(vecSt(7),find(NF3(vecSt(7),:)> th));
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['10º Kir - Cnt2: ' num2str(p.p(2))])
disp(' ')
disp('Straightness: Cnt1 - Cnt2')
a1 = MF2(vecSt(1),find(NF2(vecSt(1),:)> th));
a2 = MF3(vecSt(1),find(NF3(vecSt(1),:)> th));
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º Cnt1 - Cnt2: ' num2str(p.p(2))])
a1 = MF2(vecSt(3),find(NF2(vecSt(3),:)> th));
a2 = MF3(vecSt(3),find(NF3(vecSt(3),:)> th));
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['2.5º Cnt1 - Cnt2: ' num2str(p.p(2))])
a1 = MF2(vecSt(5),find(NF2(vecSt(5),:)> th));
a2 = MF3(vecSt(5),find(NF3(vecSt(5),:)> th));
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['5º Cnt1 - Cnt2: ' num2str(p.p(2))])
a1 = MF2(vecSt(7),find(NF2(vecSt(7),:)> th));
a2 = MF3(vecSt(7),find(NF3(vecSt(7),:)> th));
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['10º Cnt1 - Cnt2: ' num2str(p.p(2))])
disp(' ')

tth = 350;
tthk = 350;

disp('Visual Influence: Kir - Cnt1')
a1 = VI1(vecVi(1),NVI1(vecVi(1),:)>tth);
a2 = VI2(vecVi(1),NVI2(vecVi(1),:)>tth);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º Kir - Cnt1: ' num2str(p.p(2))])
a1 = VI1(vecVi(2),NVI1(vecVi(2),:)>tth);
a2 = VI2(vecVi(2),NVI2(vecVi(2),:)>tth);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['2.5º Kir - Cnt1: ' num2str(p.p(2))])
a1 = VI1(vecVi(3),NVI1(vecVi(3),:)>tth);
a2 = VI2(vecVi(3),NVI2(vecVi(3),:)>tth);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['5º Kir - Cnt1: ' num2str(p.p(2))])
a1 = VI1(vecVi(4),NVI1(vecVi(4),:)>tth);
a2 = VI2(vecVi(4),NVI2(vecVi(4),:)>tth);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['10º Kir - Cnt1: ' num2str(p.p(2))])
disp(' ')
disp('Visual Influence: Kir - Cnt2')
a1 = VI1(vecVi(1),NVI1(vecVi(1),:)>tth);
a2 = VI3(vecVi(1),NVI3(vecVi(1),:)>tth);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º Kir - Cnt1: ' num2str(p.p(2))])
a1 = VI1(vecVi(2),NVI1(vecVi(2),:)>tth);
a2 = VI3(vecVi(2),NVI3(vecVi(2),:)>tth);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['2.5º Kir - Cnt1: ' num2str(p.p(2))])
a1 = VI1(vecVi(3),NVI1(vecVi(3),:)>tth);
a2 = VI3(vecVi(3),NVI3(vecVi(3),:)>tth);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['5º Kir - Cnt1: ' num2str(p.p(2))])
a1 = VI1(vecVi(4),NVI1(vecVi(4),:)>tth);
a2 = VI3(vecVi(4),NVI3(vecVi(4),:)>tth);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['10º Kir - Cnt1: ' num2str(p.p(2))])
disp(' ')
disp('Visual Influence: Cnt1 - Cnt2')
a1 = VI3(vecVi(1),NVI3(vecVi(1),:)>tth);
a2 = VI2(vecVi(1),NVI2(vecVi(1),:)>tth);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º Kir - Cnt1: ' num2str(p.p(2))])
a1 = VI3(vecVi(2),NVI3(vecVi(2),:)>tth);
a2 = VI2(vecVi(2),NVI2(vecVi(2),:)>tth);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['2.5º Kir - Cnt1: ' num2str(p.p(2))])
a1 = VI3(vecVi(3),NVI3(vecVi(3),:)>tth);
a2 = VI2(vecVi(3),NVI2(vecVi(3),:)>tth);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['5º Kir - Cnt1: ' num2str(p.p(2))])
a1 = VI3(vecVi(4),NVI3(vecVi(4),:)>tth);
a2 = VI2(vecVi(4),NVI2(vecVi(4),:)>tth);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['10º Kir - Cnt1: ' num2str(p.p(2))])

