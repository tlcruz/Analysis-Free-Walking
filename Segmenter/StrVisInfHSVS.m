%% StrPlotsFig4
clear
clc
[params] = GetParams();
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Dark\Dark SplitT4T5Kir\';
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\VTR39\VTR39FLPKir\';
path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Size wTB\HSVS\wTBR39FLPVT05KIR\';

[STRALL1, DSTALL1, NSTRALL1, PSSALL1, NPSSALL1, pTypes1] = GetStrAndVisInf(path, params);

% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Dark\Dark SplitT4T5\';
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\VTR39\VTR39Kir\';
path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Size wTB\HSVS\wTBR39VT05KIR\';
[STRALL2, DSTALL2, NSTRALL2, PSSALL2, NPSSALL2, pTypes2] = GetStrAndVisInf(path, params);

% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Dark\Dark EmptyKir\';
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\VTR39\VTFLPKir\';
path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Size wTB\HSVS\wTBFLPVT05KIR\';
[STRALL3, DSTALL3, NSTRALL3, PSSALL3, NPSSALL3, pTypes3] = GetStrAndVisInf(path, params);
NSTRALL12 = NSTRALL1;
STRALL12 = STRALL1;
NPSSALL12 = NPSSALL1;
PSSALL12 = PSSALL1;
NSTRALL22 = NSTRALL2;
STRALL22 = STRALL2;
NPSSALL22 = NPSSALL2;
PSSALL22 = PSSALL2;
NSTRALL32 = NSTRALL3;
STRALL32 = STRALL3;
NPSSALL32 = NPSSALL3;
PSSALL32 = PSSALL3;
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
%     NSTRALL1{j,2} = zeros(length(NSTRALL1{j,2}),1);
%     NSTRALL1{j,11} = zeros(length(NSTRALL1{j,11}),1);
% %     NSTRALL2{j,18} = zeros(length(NSTRALL2{j,18}),1);
% %     NSTRALL2{j,5} = zeros(length(NSTRALL2{j,5}),1);
%     NSTRALL2{j,25} = zeros(length(NSTRALL2{j,25}),1);
%     
% %     NSTRALL3{j,15} = zeros(length(NSTRALL3{j,15}),1);
% %     NSTRALL3{j,31} = zeros(length(NSTRALL3{j,31}),1);
% %     NSTRALL3{j,34} = zeros(length(NSTRALL3{j,33}),1);
% %     NSTRALL2{j,7} = zeros(length(NSTRALL2{j,7}),1);
% %     NSTRALL1{j,10} = zeros(length(NSTRALL1{j,10}),1);
%     
%     NPSSALL1{j,2} = zeros(length(NPSSALL1{j,2}),1);
%     NPSSALL1{j,11} = zeros(length(NPSSALL1{j,11}),1);
% %     NPSSALL1{j,27} = zeros(length(NPSSALL1{j,27}),1);
% %     NPSSALL2{j,18} = zeros(length(NPSSALL2{j,18}),1);
%     NPSSALL2{j,25} = zeros(length(NPSSALL2{j,25}),1);
%     
% %     NPSSALL3{j,15} = zeros(length(NPSSALL3{j,15}),1);
% %     NPSSALL3{j,31} = zeros(length(NPSSALL3{j,31}),1);
% %     NPSSALL2{j,20} = zeros(length(NPSSALL2{j,22}),1);
% %     NPSSALL1{j,10} = zeros(length(NPSSALL1{j,10}),1);
% end

params.pTypes = pTypes1;
params.lthr = 350;
[MSTR1, SEMSTR1, MF1, NF1] = GetGMSEM(STRALL1, DSTALL1, NSTRALL1, params);
[MSTR2, SEMSTR2, MF2, NF2] = GetGMSEM(STRALL2, DSTALL2, NSTRALL2, params);
[MSTR3, SEMSTR3, MF3, NF3] = GetGMSEM(STRALL3, DSTALL3, NSTRALL3, params);

params.vthr = 00;
[VI1, NVI1, MVI1, SEMVI1] = PSStoVisInf(PSSALL1, NPSSALL1, params);
[VI2, NVI2, MVI2, SEMVI2] = PSStoVisInf(PSSALL2, NPSSALL2, params);
[VI3, NVI3, MVI3, SEMVI3] = PSStoVisInf(PSSALL3, NPSSALL3, params);

th = 150;
%
inds = [3 5 7 1];
figure,
subplot(1,2,1)
errorbar(1:length(inds), MSTR1(inds), SEMSTR1(inds), 'o', 'color', [1 0 0], 'markersize', 10, 'markerfacecolor',[1 0 0]);
hold on
errorbar(1:length(inds), MSTR2(inds), SEMSTR2(inds), 'o', 'color', [0 0 0.6], 'markersize', 10, 'markerfacecolor',[0 0 0.6]);
hold on
errorbar(1:length(inds), MSTR3(inds), SEMSTR3(inds), 'o', 'color', [0 0 1], 'markersize', 10, 'markerfacecolor',[0 0 1]);

axis([0 5 20 55])
set(gca,'xtick',[1 2 3 4]);
set(gca,'xticklabel',{'1º','2.5º','5º','10º'});
xlabel('Dot Size')
ylabel('Straightness (a.u.)')
subplot(1,2,2)
indsvi = [2 3 4 1];
errorbar(1:4, MVI1(indsvi), SEMVI1(indsvi), 'o', 'color', [1 0 0], 'markersize', 10, 'markerfacecolor',[1 0 0]);
hold on
errorbar(1:4, MVI2(indsvi), SEMVI2(indsvi), 'o', 'color', [0 0 0.6], 'markersize', 10, 'markerfacecolor',[0 0 0.6]);
hold on
errorbar(1:4, MVI3(indsvi), SEMVI3(indsvi), 'o', 'color', [0 0 1], 'markersize', 10, 'markerfacecolor',[0 0 1]);

axis([0 5 -0.05 0.35])
% axis([0 5 0 200])
set(gca,'xtick',[1 2 3 4]); 
set(gca,'xticklabel',{'1º','2.5º','5º','10º'});
xlabel('Dot Size')
ylabel('Visual Influence (a.u.)')

%
vecSt = [3 4 5 6 7 8 1 2];
vecVi = [2 3 4 1];
clc
th = 350;
disp('Straightness: Kir - Cnt1')
a1 = MF1(vecSt(1),find(NF1(vecSt(1),:)> th));
a2 = MF2(vecSt(1),find(NF2(vecSt(1),:)> th));
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º Kir - Cnt1: ' num2str(p.p(2))])
a1 = MF1(vecSt(3),find(NF1(vecSt(3),:)> th));
a2 = MF2(vecSt(3),find(NF2(vecSt(3),:)> th));
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['2.5º Kir - Cnt1: ' num2str(p.p(2))])
a1 = MF1(vecSt(5),find(NF1(vecSt(5),:)> th));
a2 = MF2(vecSt(5),find(NF2(vecSt(5),:)> th));
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['5º Kir - Cnt1: ' num2str(p.p(2))])
a1 = MF1(vecSt(7),find(NF1(vecSt(7),:)> th));
a2 = MF2(vecSt(7),find(NF2(vecSt(7),:)> th));
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['10º Kir - Cnt1: ' num2str(p.p(2))])
disp(' ')
disp('Straightness: Kir - Cnt2')
a1 = MF1(vecSt(1),find(NF1(vecSt(1),:)> th));
a2 = MF3(vecSt(1),find(NF3(vecSt(1),:)> th));
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º Kir - Cnt2: ' num2str(p.p(2))])
a1 = MF1(vecSt(3),find(NF1(vecSt(3),:)> th));
a2 = MF3(vecSt(3),find(NF3(vecSt(3),:)> th));
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['2.5º Kir - Cnt2: ' num2str(p.p(2))])
a1 = MF1(vecSt(5),find(NF1(vecSt(5),:)> th));
a2 = MF3(vecSt(5),find(NF3(vecSt(5),:)> th));
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['5º Kir - Cnt2: ' num2str(p.p(2))])
a1 = MF1(vecSt(7),find(NF1(vecSt(7),:)> th));
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
disp('Visual Influence: Kir - Cnt1')
a1 = VI1(vecVi(1),:);
a2 = VI2(vecVi(1),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º Kir - Cnt1: ' num2str(p.p(2))])
a1 = VI1(vecVi(2),:);
a2 = VI2(vecVi(2),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['2.5º Kir - Cnt1: ' num2str(p.p(2))])
a1 = VI1(vecVi(3),:);
a2 = VI2(vecVi(3),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['5º Kir - Cnt1: ' num2str(p.p(2))])
a1 = VI1(vecVi(4),:);
a2 = VI2(vecVi(4),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['10º Kir - Cnt1: ' num2str(p.p(2))])
disp(' ')
disp('Visual Influence: Kir - Cnt2')
a1 = VI1(vecVi(1),:);
a2 = VI3(vecVi(1),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º Kir - Cnt1: ' num2str(p.p(2))])
a1 = VI1(vecVi(2),:);
a2 = VI3(vecVi(2),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['2.5º Kir - Cnt1: ' num2str(p.p(2))])
a1 = VI1(vecVi(3),:);
a2 = VI3(vecVi(3),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['5º Kir - Cnt1: ' num2str(p.p(2))])
a1 = VI1(vecVi(4),:);
a2 = VI3(vecVi(4),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['10º Kir - Cnt1: ' num2str(p.p(2))])
disp(' ')
disp('Visual Influence: Cnt1 - Cnt2')
a1 = VI3(vecVi(1),:);
a2 = VI2(vecVi(1),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['1º Kir - Cnt1: ' num2str(p.p(2))])
a1 = VI3(vecVi(2),:);
a2 = VI2(vecVi(2),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['2.5º Kir - Cnt1: ' num2str(p.p(2))])
a1 = VI3(vecVi(3),:);
a2 = VI2(vecVi(3),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['5º Kir - Cnt1: ' num2str(p.p(2))])
a1 = VI3(vecVi(4),:);
a2 = VI2(vecVi(4),:);
[p] = mwwtest(a1(a1~=0),a2(a2~=0),0);
disp(['10º Kir - Cnt1: ' num2str(p.p(2))])

%%
nbts = [];
for k1 = 1 : size(NPSSALL3,1)
    for k2 = 1 : size(NPSSALL3,2)
        nbts = vertcat(nbts, length(NPSSALL3{k1,k2}));
    end
end
sum(nbts)