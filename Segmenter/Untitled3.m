clear
clc
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Dark\Dark WTTB\';
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\WT TB\';
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\VTR39\VTFLPKir\';
path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Size wTB\HSVS\wTBFLPVT05KIR\';
[stfbcntHS1, pTypesHS] = GetSTandFB(path);
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\VTR39\VTR39Kir\';
path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Size wTB\HSVS\wTBR39VT05KIR\';
[stfbcntHS2, pTypesHS] = GetSTandFB(path);
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\VTR39\VTR39FLPKir\';
path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Size wTB\HSVS\wTBR39FLPVT05KIR\';
[stfbkirHS, pTypesHS] = GetSTandFB(path);

% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\EmptySplit 10xKir\';
path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Size wTB\T4T5\wTBEmptySplitOtdKIR\';
[stfbcntT4T51, pTypesT4T5] = GetSTandFB(path);
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\Split T4T5 Cnt\';
path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Size wTB\T4T5\wTBSplitT4T5KIR\';
[stfbcntT4T52, pTypesT4T5] = GetSTandFB(path);
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\SplitT4T5 Kir\';
path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Size wTB\T4T5\wTBSplitT4T5OtdKIR\';
[stfbkirT4T5, pTypesT4T5] = GetSTandFB(path);

% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\WT TB\';
% [stfb, pTypes] = GetSTandFB(path);

%% Plot %FB WT
figure,
hold on
inds = [5 6 7 8 9 10 1 2 3 4];
for j = 1 : length(pTypes)
    vaux = [];
    vT = [];
    a = inds(j);
    for n = 1 : size(stfb.pFB,1);
        stfb.pT{n,a}(isnan(stfb.pFB{n,a})) = 1;
        stfb.pFB{n,a}(isnan(stfb.pFB{n,a})) = 0;
        stfb.pT{n,a}(isempty(stfb.pFB{n,a})) = 1;
        stfb.pFB{n,a}(isempty(stfb.pFB{n,a})) = 0;
        vaux = vertcat(vaux, stfb.pFB{n,a}'*stfb.pT{n,a}/sum(stfb.pT{n,a}));
        vT = vertcat(vT, sum(stfb.pT{n,a}));
    end
    vT(isnan(vaux)) = [];
    vaux(isnan(vaux)) = [];
    if mod(j,2) == 1
        errorbar(j, vaux'*vT/sum(vT), sqrt(sum((vaux-vaux'*vT/sum(vT)).*(vaux-vaux'*vT/sum(vT)))/size(stfb.pFB,1))/sqrt(size(stfb.pFB,1)), 'or')
    else
        errorbar(j-1, vaux'*vT/sum(vT), sqrt(sum((vaux-vaux'*vT/sum(vT)).*(vaux-vaux'*vT/sum(vT)))/size(stfb.pFB,1))/sqrt(size(stfb.pFB,1)), 'ob')
    end
end
axis([0.5 8.5 0.1 0.5])
xlabel('Dot Size(º)')
ylabel('% Time FW')

%% Plot %Bias FB WT
figure,
hold on
inds = [5 6 7 8 9 10 1 2 3 4];
for j = 1 : length(pTypes)
    vaux = [];
    vT = [];
    a = inds(j);
    for n = 1 : size(stfb.biasFB,1);
        stfb.NbiasFB{n,a}(isnan(stfb.biasFB{n,a})) = 1;
        stfb.biasFB{n,a}(isnan(stfb.biasFB{n,a})) = 0;
        stfb.NbiasFB{n,a}(isempty(stfb.biasFB{n,a})) = 1;
        stfb.biasFB{n,a}(isempty(stfb.biasFB{n,a})) = 0;
        vaux = vertcat(vaux, max(stfb.biasFB{n,a}'*stfb.NbiasFB{n,a}/sum(stfb.NbiasFB{n,a})));
        vT = vertcat(vT, sum(stfb.pT{n,a}));
    end
    vT(isnan(vaux)) = [];
    vaux(isnan(vaux)) = [];
    
    if mod(j,2) == 1
        errorbar(j, vaux'*vT/sum(vT), ...
            sqrt(sum((vaux-vaux'*vT/sum(vT)).*(vaux-vaux'*vT/sum(vT)))...
            /size(stfb.biasFB,1))/sqrt(size(stfb.biasFB,1)), 'or')
    else
        errorbar(j-1, vaux'*vT/sum(vT), ...
            sqrt(sum((vaux-vaux'*vT/sum(vT)).*(vaux-vaux'*vT/sum(vT)))...
            /size(stfb.biasFB,1))/sqrt(size(stfb.biasFB,1)), 'ob')
    end
end
axis([0.5 8.5 0.5 1])
xlabel('Dot Size(º)')
ylabel('% Bias FW')

%% Plot Size Bias FB WT
figure,
hold on
inds = [5 6 7 8 9 10 1 2 3 4];
for j = 1 : length(pTypes)
    vaux = [];
    vT = [];
    a = inds(j);
    for n = 1 : size(stfb.valFB,1);
        stfb.NvalFB{n,a}(isnan(stfb.valFB{n,a})) = 1;
        stfb.valFB{n,a}(isnan(stfb.valFB{n,a})) = 0;
        stfb.NvalFB{n,a}(isempty(stfb.valFB{n,a})) = 1;
        stfb.valFB{n,a}(isempty(stfb.valFB{n,a})) = 0;
        vaux = vertcat(vaux, max(stfb.valFB{n,a}'*stfb.NvalFB{n,a}/sum(stfb.NvalFB{n,a})));
        vT = vertcat(vT, sum(stfb.NvalFB{n,a}));
    end
    vT(isnan(vaux)) = [];
    vaux(isnan(vaux)) = [];
    
    if mod(j,2) == 1
        errorbar(j, vaux'*vT/sum(vT), ...
            sqrt(sum((vaux-vaux'*vT/sum(vT)).*(vaux-vaux'*vT/sum(vT)))...
            /size(stfb.valFB,1))/sqrt(size(stfb.valFB,1)), 'or')
    else
        errorbar(j-1, vaux'*vT/sum(vT), ...
            sqrt(sum((vaux-vaux'*vT/sum(vT)).*(vaux-vaux'*vT/sum(vT)))...
            /size(stfb.valFB,1))/sqrt(size(stfb.valFB,1)), 'ob')
    end
end
% axis([0.5 8.5 0 180])
axis([0.5 8.5 0 250])
xlabel('Dot Size(º)')
ylabel('Drift (º/s)')

%% Plot %ST WT
figure,
hold on
inds = [5 6 7 8 9 10 1 2 3 4];
for j = 1 : length(pTypes)
    vaux = [];
    vT = [];
    a = inds(j);
    for n = 1 : size(stfb.pST,1);
        stfb.pT{n,a}(isnan(stfb.pST{n,a})) = 1;
        stfb.pST{n,a}(isnan(stfb.pST{n,a})) = 0;
        stfb.pT{n,a}(isempty(stfb.pST{n,a})) = 1;
        stfb.pST{n,a}(isempty(stfb.pST{n,a})) = 0;
        vaux = vertcat(vaux, stfb.pST{n,a}'*stfb.pT{n,a}/sum(stfb.pT{n,a}));
        vT = vertcat(vT, sum(stfb.pT{n,a}));
    end
    vT(isnan(vaux)) = [];
    vaux(isnan(vaux)) = [];
    if mod(j,2) == 1
        errorbar(j, vaux'*vT/sum(vT), sqrt(sum((vaux-vaux'*vT/sum(vT)).*(vaux-vaux'*vT/sum(vT)))/size(stfb.pFB,1))/sqrt(size(stfb.pFB,1)), 'or')
    else
        errorbar(j-1, vaux'*vT/sum(vT), sqrt(sum((vaux-vaux'*vT/sum(vT)).*(vaux-vaux'*vT/sum(vT)))/size(stfb.pFB,1))/sqrt(size(stfb.pFB,1)), 'ob')
    end
end
axis([0.5 8.5 0.3 0.7])

%% Plot %ST Bias
figure,
hold on
inds = [5 6 7 8 9 10 1 2 3 4];
for j = 1 : length(pTypes)
    vaux1 = [];
    vT1 = [];
    a = inds(j);
    for n = 1 : size(stfb.biasST1,1);
        stfb.NbiasST1{n,a}(isnan(stfb.biasST1{n,a})) = 1;
        stfb.biasST1{n,a}(isnan(stfb.biasST1{n,a})) = 0;
        stfb.NbiasST1{n,a}(isempty(stfb.biasST1{n,a})) = 1;
        stfb.biasST1{n,a}(isempty(stfb.biasST1{n,a})) = 0;
        vaux1 = vertcat(vaux1, stfb.biasST1{n,a}'*stfb.NbiasST1{n,a}/sum(stfb.NbiasST1{n,a}));
        vT1 = vertcat(vT1, sum(stfb.NbiasST1{n,a}));
    end
    vT1(isnan(vaux1)) = [];
    vaux1(isnan(vaux1)) = [];
    
    vaux20 = [];
    vT20 = [];
    a = inds(j);
    for n = 1 : size(stfb.biasST20,1);
        stfb.NbiasST20{n,a}(isnan(stfb.biasST20{n,a})) = 1;
        stfb.biasST20{n,a}(isnan(stfb.biasST20{n,a})) = 0;
        stfb.NbiasST20{n,a}(isempty(stfb.biasST20{n,a})) = 1;
        stfb.biasST20{n,a}(isempty(stfb.biasST20{n,a})) = 0;
        vaux20 = vertcat(vaux20, stfb.biasST20{n,a}'*stfb.NbiasST20{n,a}/sum(stfb.NbiasST20{n,a}));
        vT20 = vertcat(vT20, sum(stfb.NbiasST20{n,a}));
    end
    vT20(isnan(vaux20)) = [];
    vaux20(isnan(vaux20)) = [];
    
    vaux100 = [];
    vT100 = [];
    a = inds(j);
    for n = 1 : size(stfb.biasST100,1);
        stfb.NbiasST100{n,a}(isnan(stfb.biasST100{n,a})) = 1;
        stfb.biasST100{n,a}(isnan(stfb.biasST100{n,a})) = 0;
        stfb.NbiasST100{n,a}(isempty(stfb.biasST100{n,a})) = 1;
        stfb.biasST100{n,a}(isempty(stfb.biasST100{n,a})) = 0;
        vaux100 = vertcat(vaux100, stfb.biasST100{n,a}'*stfb.NbiasST100{n,a}/sum(stfb.NbiasST100{n,a}));
        vT100 = vertcat(vT100, sum(stfb.NbiasST100{n,a}));
    end
    vT100(isnan(vaux100)) = [];
    vaux100(isnan(vaux100)) = [];
    
    
    plot([0 9], [0.5 0.5], '--g')
    if mod(j,2) == 1
        errorbar(j, vaux1'*vT1/sum(vT1), ...
            sqrt(sum((vaux1-vaux1'*vT1/sum(vT1)).*(vaux1-vaux1'*vT1/sum(vT1)))/size(stfb.biasST1,1))...
            /sqrt(size(stfb.biasST1,1)), 'color', [1 0 0], 'marker', 'o')
        errorbar(j, vaux20'*vT20/sum(vT20), ...
            sqrt(sum((vaux20-vaux20'*vT20/sum(vT20)).*(vaux20-vaux20'*vT20/sum(vT20)))/size(stfb.biasST20,1))...
            /sqrt(size(stfb.biasST20,1)), 'color', [0.75 0 0], 'marker', 'o')
        errorbar(j, vaux100'*vT100/sum(vT100), ...
            sqrt(sum((vaux100-vaux100'*vT100/sum(vT100)).*(vaux100-vaux100'*vT100/sum(vT100)))/size(stfb.biasST100,1))...
            /sqrt(size(stfb.biasST100,1)), 'color', [0.5 0 0], 'marker', 'o')
    else
        errorbar(j-1, vaux1'*vT1/sum(vT1), ...
            sqrt(sum((vaux1-vaux1'*vT1/sum(vT1)).*(vaux1-vaux1'*vT1/sum(vT1)))/size(stfb.biasST1,1))...
            /sqrt(size(stfb.biasST1,1)), 'color', [0 0 1], 'marker', 'o')
        errorbar(j-1, vaux20'*vT20/sum(vT20), ...
            sqrt(sum((vaux20-vaux20'*vT20/sum(vT20)).*(vaux20-vaux20'*vT20/sum(vT20)))/size(stfb.biasST20,1))...
            /sqrt(size(stfb.biasST20,1)), 'color', [0 0 0.75], 'marker', 'o')
        errorbar(j-1, vaux100'*vT100/sum(vT100), ...
            sqrt(sum((vaux100-vaux100'*vT100/sum(vT100)).*(vaux100-vaux100'*vT100/sum(vT100)))/size(stfb.biasST100,1))...
            /sqrt(size(stfb.biasST100,1)), 'color', [0 0 0.5], 'marker', 'o')
    end
end
axis([0.5 8.5 0 1])

%% Plot ST max WT
figure,
hold on
lcents = -50:20:1500;
inds = [5 6 7 8 9 10 1 2 3 4];
cmap1 = autumn(15);
cmap2 = winter(15);
for j = 1 : length(pTypes)
    vaux = [];
    a = inds(j);
    for n = 1 : size(stfb.lST,1);
        vaux = vertcat(vaux, ...
            hist(abs(stfb.lST{n,a}), lcents)/sum(hist(stfb.lST{n,a}, lcents)));
    end
    
    if mod(j,2) == 1
        plot(lcents, smooth(mean(vaux)), 'color', cmap1(j,:), 'linewidth', 2)
    else
        plot(lcents, smooth(mean(vaux)), 'color', cmap2(j,:), 'linewidth', 2)
    end
end
axis([0 1500 0 0.05])








%% Plot %FB T4T5
figure,
hold on
% inds = [5 6 7 8 9 10 1 2 3 4];
inds = [3 4 5 6 7 8 1 2];
for j = 1 : length(pTypesT4T5)
    vauxcnt1 = [];
    vTcnt1 = [];
    a = inds(j);
    for n = 1 : size(stfbcntT4T51.pFB,1);
        stfbcntT4T51.pT{n,a}(isnan(stfbcntT4T51.pFB{n,a})) = 1;
        stfbcntT4T51.pFB{n,a}(isnan(stfbcntT4T51.pFB{n,a})) = 0;
        stfbcntT4T51.pT{n,a}(isempty(stfbcntT4T51.pFB{n,a})) = 1;
        stfbcntT4T51.pFB{n,a}(isempty(stfbcntT4T51.pFB{n,a})) = 0;
        vauxcnt1 = vertcat(vauxcnt1, stfbcntT4T51.pFB{n,a}'*stfbcntT4T51.pT{n,a}/sum(stfbcntT4T51.pT{n,a}));
        vTcnt1 = vertcat(vTcnt1, sum(stfbcntT4T51.pT{n,a}));
    end
    vTcnt1(isnan(vauxcnt1)) = [];
    vauxcnt1(isnan(vauxcnt1)) = [];
    vauxcnt2 = [];
    vTcnt2 = [];
    a = inds(j);
    for n = 1 : size(stfbcntT4T52.pFB,1);
        stfbcntT4T52.pT{n,a}(isnan(stfbcntT4T52.pFB{n,a})) = 1;
        stfbcntT4T52.pFB{n,a}(isnan(stfbcntT4T52.pFB{n,a})) = 0;
        stfbcntT4T52.pT{n,a}(isempty(stfbcntT4T52.pFB{n,a})) = 1;
        stfbcntT4T52.pFB{n,a}(isempty(stfbcntT4T52.pFB{n,a})) = 0;
        vauxcnt2 = vertcat(vauxcnt2, stfbcntT4T52.pFB{n,a}'*stfbcntT4T52.pT{n,a}/sum(stfbcntT4T52.pT{n,a}));
        vTcnt2 = vertcat(vTcnt2, sum(stfbcntT4T52.pT{n,a}));
    end
    vTcnt2(isnan(vauxcnt2)) = [];
    vauxcnt2(isnan(vauxcnt2)) = [];
    vauxkir = [];
    vTkir = [];
    a = inds(j);
    for n = 1 : size(stfbkirT4T5.pFB,1);
        stfbkirT4T5.pT{n,a}(isnan(stfbkirT4T5.pFB{n,a})) = 1;
        stfbkirT4T5.pFB{n,a}(isnan(stfbkirT4T5.pFB{n,a})) = 0;
        stfbkirT4T5.pT{n,a}(isempty(stfbkirT4T5.pFB{n,a})) = 1;
        stfbkirT4T5.pFB{n,a}(isempty(stfbkirT4T5.pFB{n,a})) = 0;
        vauxkir = vertcat(vauxkir, stfbkirT4T5.pFB{n,a}'*stfbkirT4T5.pT{n,a}/sum(stfbkirT4T5.pT{n,a}));
        vTkir = vertcat(vTkir, sum(stfbkirT4T5.pT{n,a}));
    end
    vTkir(isnan(vauxkir)) = [];
    vauxkir(isnan(vauxkir)) = [];
    
    
    if mod(j,2) == 1
        errorbar(j, vauxcnt1'*vTcnt1/sum(vTcnt1), ...
            sqrt(sum((vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)).*(vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)))...
            /size(stfbcntT4T51.pFB,1))/sqrt(size(stfbcntT4T51.pFB,1)), 'om')
        errorbar(j, vauxcnt2'*vTcnt2/sum(vTcnt2), ...
            sqrt(sum((vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)).*(vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)))...
            /size(stfbcntT4T52.pFB,1))/sqrt(size(stfbcntT4T52.pFB,1)), 'om')
        
        errorbar(j, vauxkir'*vTkir/sum(vTkir), ...
            sqrt(sum((vauxkir-vauxkir'*vTkir/sum(vTkir)).*(vauxkir-vauxkir'*vTkir/sum(vTkir)))...
            /size(stfbkirT4T5.pFB,1))/sqrt(size(stfbkirT4T5.pFB,1)), 'or')
    else
        errorbar(j-1, vauxcnt1'*vTcnt1/sum(vTcnt1),...
            sqrt(sum((vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)).*(vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)))...
            /size(stfbcntT4T51.pFB,1))/sqrt(size(stfbcntT4T51.pFB,1)), 'oc')
        errorbar(j-1, vauxcnt2'*vTcnt2/sum(vTcnt2),...
            sqrt(sum((vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)).*(vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)))...
            /size(stfbcntT4T52.pFB,1))/sqrt(size(stfbcntT4T52.pFB,1)), 'oc')
        errorbar(j-1, vauxkir'*vTkir/sum(vTkir),...
            sqrt(sum((vauxkir-vauxkir'*vTkir/sum(vTkir)).*(vauxkir-vauxkir'*vTkir/sum(vTkir)))...
            /size(stfbkirT4T5.pFB,1))/sqrt(size(stfbkirT4T5.pFB,1)), 'ob')
    end
end
axis([0.5 8.5 0.1 0.5])
xlabel('Dot Size(º)')
ylabel('% Time FW')

%% Plot %Bias FB T4T5
figure,
hold on
% inds = [5 6 7 8 9 10 1 2 3 4];
inds = [3 4 5 6 7 8 1 2];
for j = 1 : length(pTypesT4T5)
    vauxcnt1 = [];
    vTcnt1 = [];
    a = inds(j);
    for n = 1 : size(stfbcntT4T51.biasFB,1);
        stfbcntT4T51.NbiasFB{n,a}(isnan(stfbcntT4T51.biasFB{n,a})) = 1;
        stfbcntT4T51.biasFB{n,a}(isnan(stfbcntT4T51.biasFB{n,a})) = 0;
        stfbcntT4T51.NbiasFB{n,a}(isempty(stfbcntT4T51.biasFB{n,a})) = 1;
        stfbcntT4T51.biasFB{n,a}(isempty(stfbcntT4T51.biasFB{n,a})) = 0;
        if length(stfbcntT4T51.NbiasFB{n,a}) ~= length(stfbcntT4T51.biasFB{n,a})
            stfbcntT4T51.NbiasFB{n,a} = stfbcntT4T51.NbiasFB{n,a}(1:length(stfbcntT4T51.biasFB{n,a}));
        end
        vauxcnt1 = vertcat(vauxcnt1, max(stfbcntT4T51.biasFB{n,a}'*stfbcntT4T51.NbiasFB{n,a}/sum(stfbcntT4T51.NbiasFB{n,a})));
        vTcnt1 = vertcat(vTcnt1, sum(stfbcntT4T51.NbiasFB{n,a}));
    end
    vTcnt1(isnan(vauxcnt1)) = [];
    vauxcnt1(isnan(vauxcnt1)) = [];
    if length(vauxcnt1) ~= length(vTcnt1)
        vauxcnt1(vauxcnt1==0) = [];
    end
    
    vauxcnt2 = [];
    vTcnt2 = [];
    a = inds(j);
    for n = 1 : size(stfbcntT4T52.biasFB,1);
        stfbcntT4T52.NbiasFB{n,a}(isnan(stfbcntT4T52.biasFB{n,a})) = 1;
        stfbcntT4T52.biasFB{n,a}(isnan(stfbcntT4T52.biasFB{n,a})) = 0;
        stfbcntT4T52.NbiasFB{n,a}(isempty(stfbcntT4T52.biasFB{n,a})) = 1;
        stfbcntT4T52.biasFB{n,a}(isempty(stfbcntT4T52.biasFB{n,a})) = 0;
        vauxcnt2 = vertcat(vauxcnt2, max(stfbcntT4T52.biasFB{n,a}'*stfbcntT4T52.NbiasFB{n,a}/sum(stfbcntT4T52.NbiasFB{n,a})));
        vTcnt2 = vertcat(vTcnt2, sum(stfbcntT4T52.NbiasFB{n,a}));
    end
    vTcnt2(isnan(vauxcnt2)) = [];
    vauxcnt2(isnan(vauxcnt2)) = [];
    if length(vauxcnt2) ~= length(vTcnt2)
        vauxcnt2(vauxcnt2==0) = [];
    end
    
    vauxkir = [];
    vTkir = [];
    a = inds(j);
    for n = 1 : size(stfbkirT4T5.biasFB,1);
        stfbkirT4T5.NbiasFB{n,a}(isnan(stfbkirT4T5.biasFB{n,a})) = 1;
        stfbkirT4T5.biasFB{n,a}(isnan(stfbkirT4T5.biasFB{n,a})) = 0;
        stfbkirT4T5.NbiasFB{n,a}(isempty(stfbkirT4T5.biasFB{n,a})) = 1;
        stfbkirT4T5.biasFB{n,a}(isempty(stfbkirT4T5.biasFB{n,a})) = 0;
        vauxkir = vertcat(vauxkir, max(stfbkirT4T5.biasFB{n,a}'*stfbkirT4T5.NbiasFB{n,a}/sum(stfbkirT4T5.NbiasFB{n,a})));
        vTkir = vertcat(vTkir, sum(stfbkirT4T5.NbiasFB{n,a}));
    end
    vTkir(isnan(vauxkir)) = [];
    vauxkir(isnan(vauxkir)) = [];
    if length(vauxkir) ~= length(vTkir)
        vauxkir(vauxkir==0) = [];
    end
    
    if mod(j,2) == 1
        errorbar(j, vauxcnt1'*vTcnt1/sum(vTcnt1), ...
            sqrt(sum((vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)).*(vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)))...
            /size(stfbcntT4T51.NbiasFB,1))/sqrt(size(stfbcntT4T51.NbiasFB,1)), 'om')
        errorbar(j, vauxcnt2'*vTcnt2/sum(vTcnt2), ...
            sqrt(sum((vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)).*(vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)))...
            /size(stfbcntT4T52.NbiasFB,1))/sqrt(size(stfbcntT4T52.NbiasFB,1)), 'om')
        
        errorbar(j, vauxkir'*vTkir/sum(vTkir), ...
            sqrt(sum((vauxkir-vauxkir'*vTkir/sum(vTkir)).*(vauxkir-vauxkir'*vTkir/sum(vTkir)))...
            /size(stfbkirT4T5.NbiasFB,1))/sqrt(size(stfbkirT4T5.NbiasFB,1)), 'or')
    else
        errorbar(j-1, vauxcnt1'*vTcnt1/sum(vTcnt1),...
            sqrt(sum((vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)).*(vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)))...
            /size(stfbcntT4T51.NbiasFB,1))/sqrt(size(stfbcntT4T51.NbiasFB,1)), 'oc')
        errorbar(j-1, vauxcnt2'*vTcnt2/sum(vTcnt2),...
            sqrt(sum((vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)).*(vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)))...
            /size(stfbcntT4T52.NbiasFB,1))/sqrt(size(stfbcntT4T52.NbiasFB,1)), 'oc')
        errorbar(j-1, vauxkir'*vTkir/sum(vTkir),...
            sqrt(sum((vauxkir-vauxkir'*vTkir/sum(vTkir)).*(vauxkir-vauxkir'*vTkir/sum(vTkir)))...
            /size(stfbkirT4T5.NbiasFB,1))/sqrt(size(stfbkirT4T5.NbiasFB,1)), 'ob')
    end
end
axis([0.5 8.5 0.3 1])
xlabel('Dot Size(º)')
ylabel('% Bias FW')

%% Plot Size Bias FB T4T5
figure,
hold on
% inds = [5 6 7 8 9 10 1 2 3 4];
inds = [3 4 5 6 7 8 1 2];
for j = 1 : length(pTypesT4T5)
    vauxcnt1 = [];
    vTcnt1 = [];
    a = inds(j);
    for n = 1 : size(stfbcntT4T51.valFB,1);
        stfbcntT4T51.NvalFB{n,a}(isnan(stfbcntT4T51.valFB{n,a})) = 1;
        stfbcntT4T51.valFB{n,a}(isnan(stfbcntT4T51.valFB{n,a})) = 0;
        stfbcntT4T51.NvalFB{n,a}(isempty(stfbcntT4T51.valFB{n,a})) = 1;
        stfbcntT4T51.valFB{n,a}(isempty(stfbcntT4T51.valFB{n,a})) = 0;
        if length(stfbcntT4T51.NvalFB{n,a}) ~= length(stfbcntT4T51.valFB{n,a})
            stfbcntT4T51.NvalFB{n,a} = stfbcntT4T51.NvalFB{n,a}(1:length(stfbcntT4T51.valFB{n,a}));
        end
        vauxcnt1 = vertcat(vauxcnt1, max(stfbcntT4T51.valFB{n,a}'*stfbcntT4T51.NvalFB{n,a}/sum(stfbcntT4T51.NvalFB{n,a})));
        vTcnt1 = vertcat(vTcnt1, sum(stfbcntT4T51.NvalFB{n,a}));
    end
    vTcnt1(isnan(vauxcnt1)) = [];
    vauxcnt1(isnan(vauxcnt1)) = [];
    if length(vauxcnt1) ~= length(vTcnt1)
        vauxcnt1(vauxcnt1==0) = [];
    end
    
    vauxcnt2 = [];
    vTcnt2 = [];
    a = inds(j);
    for n = 1 : size(stfbcntT4T52.valFB,1);
        stfbcntT4T52.NvalFB{n,a}(isnan(stfbcntT4T52.valFB{n,a})) = 1;
        stfbcntT4T52.valFB{n,a}(isnan(stfbcntT4T52.valFB{n,a})) = 0;
        stfbcntT4T52.NvalFB{n,a}(isempty(stfbcntT4T52.valFB{n,a})) = 1;
        stfbcntT4T52.valFB{n,a}(isempty(stfbcntT4T52.valFB{n,a})) = 0;
        vauxcnt2 = vertcat(vauxcnt2, max(stfbcntT4T52.valFB{n,a}'*stfbcntT4T52.NvalFB{n,a}/sum(stfbcntT4T52.NvalFB{n,a})));
        vTcnt2 = vertcat(vTcnt2, sum(stfbcntT4T52.NvalFB{n,a}));
    end
    vTcnt2(isnan(vauxcnt2)) = [];
    vauxcnt2(isnan(vauxcnt2)) = [];
    if length(vauxcnt2) ~= length(vTcnt2)
        vauxcnt2(vauxcnt2==0) = [];
    end
    
    vauxkir = [];
    vTkir = [];
    a = inds(j);
    for n = 1 : size(stfbkirT4T5.valFB,1);
        stfbkirT4T5.NvalFB{n,a}(isnan(stfbkirT4T5.valFB{n,a})) = 1;
        stfbkirT4T5.valFB{n,a}(isnan(stfbkirT4T5.valFB{n,a})) = 0;
        stfbkirT4T5.NvalFB{n,a}(isempty(stfbkirT4T5.valFB{n,a})) = 1;
        stfbkirT4T5.valFB{n,a}(isempty(stfbkirT4T5.valFB{n,a})) = 0;
        vauxkir = vertcat(vauxkir, max(stfbkirT4T5.valFB{n,a}'*stfbkirT4T5.NvalFB{n,a}/sum(stfbkirT4T5.NvalFB{n,a})));
        vTkir = vertcat(vTkir, sum(stfbkirT4T5.NvalFB{n,a}));
    end
    vTkir(isnan(vauxkir)) = [];
    vauxkir(isnan(vauxkir)) = [];
    if length(vauxkir) ~= length(vTkir)
        vauxkir(vauxkir==0) = [];
    end
    
    if mod(j,2) == 1
        errorbar(j, vauxcnt1'*vTcnt1/sum(vTcnt1), ...
            sqrt(sum((vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)).*(vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)))...
            /size(stfbcntT4T51.NvalFB,1))/sqrt(size(stfbcntT4T51.NvalFB,1)), 'om')
        errorbar(j, vauxcnt2'*vTcnt2/sum(vTcnt2), ...
            sqrt(sum((vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)).*(vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)))...
            /size(stfbcntT4T52.NvalFB,1))/sqrt(size(stfbcntT4T52.NvalFB,1)), 'om')
        
        errorbar(j, vauxkir'*vTkir/sum(vTkir), ...
            sqrt(sum((vauxkir-vauxkir'*vTkir/sum(vTkir)).*(vauxkir-vauxkir'*vTkir/sum(vTkir)))...
            /size(stfbkirT4T5.NvalFB,1))/sqrt(size(stfbkirT4T5.NvalFB,1)), 'or')
    else
        errorbar(j-1, vauxcnt1'*vTcnt1/sum(vTcnt1),...
            sqrt(sum((vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)).*(vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)))...
            /size(stfbcntT4T51.NvalFB,1))/sqrt(size(stfbcntT4T51.NvalFB,1)), 'oc')
        errorbar(j-1, vauxcnt2'*vTcnt2/sum(vTcnt2),...
            sqrt(sum((vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)).*(vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)))...
            /size(stfbcntT4T52.NvalFB,1))/sqrt(size(stfbcntT4T52.NvalFB,1)), 'oc')
        errorbar(j-1, vauxkir'*vTkir/sum(vTkir),...
            sqrt(sum((vauxkir-vauxkir'*vTkir/sum(vTkir)).*(vauxkir-vauxkir'*vTkir/sum(vTkir)))...
            /size(stfbkirT4T5.NvalFB,1))/sqrt(size(stfbkirT4T5.NvalFB,1)), 'ob')
    end
end
% axis([0.5 8.5 0 120])
axis([0.5 8.5 0 250])
xlabel('Dot Size(º)')
ylabel('Drift (º/s)')

%% Plot Size Bias FB T4T5
figure,
hold on
% inds = [5 6 7 8 9 10 1 2 3 4];
inds = [3 4 5 6 7 8 1 2];
for j = 1 : length(pTypesT4T5)
    vauxcnt1 = [];
    vTcnt1 = [];
    a = inds(j);
    for n = 1 : size(stfbcntT4T51.angDisp,1);
        stfbcntT4T51.NangDisp{n,a}(isnan(stfbcntT4T51.angDisp{n,a})) = 1;
        stfbcntT4T51.angDisp{n,a}(isnan(stfbcntT4T51.angDisp{n,a})) = 0;
        stfbcntT4T51.NangDisp{n,a}(isempty(stfbcntT4T51.angDisp{n,a})) = 1;
        stfbcntT4T51.angDisp{n,a}(isempty(stfbcntT4T51.angDisp{n,a})) = 0;
        vauxcnt1 = vertcat(vauxcnt1, max(stfbcntT4T51.angDisp{n,a}'*stfbcntT4T51.NangDisp{n,a}/sum(stfbcntT4T51.NangDisp{n,a})));
        vTcnt1 = vertcat(vTcnt1, sum(stfbcntT4T51.NangDisp{n,a}));
    end
    vTcnt1(isnan(vauxcnt1)) = [];
    vauxcnt1(isnan(vauxcnt1)) = [];
    if length(vauxcnt1) ~= length(vTcnt1)
        vauxcnt1(vauxcnt1==0) = [];
    end
    
    vauxcnt2 = [];
    vTcnt2 = [];
    a = inds(j);
    for n = 1 : size(stfbcntT4T52.angDisp,1);
        stfbcntT4T52.NangDisp{n,a}(isnan(stfbcntT4T52.angDisp{n,a})) = 1;
        stfbcntT4T52.angDisp{n,a}(isnan(stfbcntT4T52.angDisp{n,a})) = 0;
        stfbcntT4T52.NangDisp{n,a}(isempty(stfbcntT4T52.angDisp{n,a})) = 1;
        stfbcntT4T52.angDisp{n,a}(isempty(stfbcntT4T52.angDisp{n,a})) = 0;
        vauxcnt2 = vertcat(vauxcnt2, max(stfbcntT4T52.angDisp{n,a}'*stfbcntT4T52.NangDisp{n,a}/sum(stfbcntT4T52.NangDisp{n,a})));
        vTcnt2 = vertcat(vTcnt2, sum(stfbcntT4T52.NangDisp{n,a}));
    end
    vTcnt2(isnan(vauxcnt2)) = [];
    vauxcnt2(isnan(vauxcnt2)) = [];
    if length(vauxcnt2) ~= length(vTcnt2)
        vauxcnt2(vauxcnt2==0) = [];
    end
    
    vauxkir = [];
    vTkir = [];
    a = inds(j);
    for n = 1 : size(stfbkirT4T5.angDisp,1);
        stfbkirT4T5.NangDisp{n,a}(isnan(stfbkirT4T5.angDisp{n,a})) = 1;
        stfbkirT4T5.angDisp{n,a}(isnan(stfbkirT4T5.angDisp{n,a})) = 0;
        stfbkirT4T5.NangDisp{n,a}(isempty(stfbkirT4T5.angDisp{n,a})) = 1;
        stfbkirT4T5.angDisp{n,a}(isempty(stfbkirT4T5.angDisp{n,a})) = 0;
        vauxkir = vertcat(vauxkir, max(stfbkirT4T5.angDisp{n,a}'*stfbkirT4T5.NangDisp{n,a}/sum(stfbkirT4T5.NangDisp{n,a})));
        vTkir = vertcat(vTkir, sum(stfbkirT4T5.NangDisp{n,a}));
    end
    vTkir(isnan(vauxkir)) = [];
    vauxkir(isnan(vauxkir)) = [];
    if length(vauxkir) ~= length(vTkir)
        vauxkir(vauxkir==0) = [];
    end
    
    if mod(j,2) == 1
        errorbar(j, vauxcnt1'*vTcnt1/sum(vTcnt1), ...
            sqrt(sum((vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)).*(vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)))...
            /size(stfbcntT4T51.NangDisp,1))/sqrt(size(stfbcntT4T51.NangDisp,1)), 'om')
        errorbar(j, vauxcnt2'*vTcnt2/sum(vTcnt2), ...
            sqrt(sum((vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)).*(vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)))...
            /size(stfbcntT4T52.NangDisp,1))/sqrt(size(stfbcntT4T52.NangDisp,1)), 'om')
        
        errorbar(j, vauxkir'*vTkir/sum(vTkir), ...
            sqrt(sum((vauxkir-vauxkir'*vTkir/sum(vTkir)).*(vauxkir-vauxkir'*vTkir/sum(vTkir)))...
            /size(stfbkirT4T5.NangDisp,1))/sqrt(size(stfbkirT4T5.NangDisp,1)), 'or')
    else
        errorbar(j-1, vauxcnt1'*vTcnt1/sum(vTcnt1),...
            sqrt(sum((vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)).*(vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)))...
            /size(stfbcntT4T51.NangDisp,1))/sqrt(size(stfbcntT4T51.NangDisp,1)), 'oc')
        errorbar(j-1, vauxcnt2'*vTcnt2/sum(vTcnt2),...
            sqrt(sum((vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)).*(vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)))...
            /size(stfbcntT4T52.NangDisp,1))/sqrt(size(stfbcntT4T52.NangDisp,1)), 'oc')
        errorbar(j-1, vauxkir'*vTkir/sum(vTkir),...
            sqrt(sum((vauxkir-vauxkir'*vTkir/sum(vTkir)).*(vauxkir-vauxkir'*vTkir/sum(vTkir)))...
            /size(stfbkirT4T5.NangDisp,1))/sqrt(size(stfbkirT4T5.NangDisp,1)), 'ob')
    end
end
% axis([0.5 8.5 0 120])
axis([0.5 8.5 0 20])
xlabel('Dot Size(º)')
ylabel('AngDisp (º/mm)')


%% Plot Size FB T4T5
figure,
hold on
% inds = [5 6 7 8 9 10 1 2 3 4];
inds = [3 4 5 6 7 8 1 2];
lcents = 1:1:150;
cmap1 = autumn(20);
cmap2 = winter(20);
for j = 1 : length(pTypesT4T5)
    vauxcnt1 = [];
    vauxcnt2 = [];
    vauxkir = [];
    a = inds(j);
    for n = 1 : size(stfbcntT4T51.lFB,1);
        vauxcnt1 = vertcat(vauxcnt1, ...
            hist(stfbcntT4T51.lFB{n,a}, lcents)/sum(hist(stfbcntT4T51.lFB{n,a}, lcents)));
    end
    for n = 1 : size(stfbcntT4T52.lFB,1);
        vauxcnt2 = vertcat(vauxcnt2, ...
            hist(stfbcntT4T52.lFB{n,a}, lcents)/sum(hist(stfbcntT4T52.lFB{n,a}, lcents)));
    end
    for n = 1 : size(stfbkirT4T5.lFB,1);
        vauxkir = vertcat(vauxkir, ...
            hist(stfbkirT4T5.lFB{n,a}, lcents)/sum(hist(stfbkirT4T5.lFB{n,a}, lcents)));
    end
    
    if mod(j,2) == 1
        plot(lcents/60, mean(vauxcnt1), 'c', 'linewidth', 2)
%         vm = mean(vauxcnt1)-std(vauxcnt1);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60, vm, 'color', 'c')
%         vm = mean(vauxcnt1)+std(vauxcnt1);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60, vm, 'color', 'c')
        
        plot(lcents/60, mean(vauxcnt2), 'c', 'linewidth', 2)
%         vm = mean(vauxcnt2)-std(vauxcnt2);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60, vm, 'color', 'c')
%         vm = mean(vauxcnt2)+std(vauxcnt2);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60, vm, 'color', 'c')
        
        plot(lcents/60, mean(vauxkir), 'b', 'linewidth', 2)
%         vm = mean(vauxkir)-std(vauxkir);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60, vm, 'color', 'b')
%         vm = mean(vauxkir)+std(vauxkir);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60, vm, 'color', 'b')
    else
        plot(lcents/60, mean(vauxcnt1), 'm', 'linewidth', 2)
%         vm = mean(vauxcnt1)-std(vauxcnt1);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60,vm, 'color', 'm')
%         vm = mean(vauxcnt1)+std(vauxcnt1);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60, vm, 'color', 'm')
        
        plot(lcents/60, mean(vauxcnt2), 'm', 'linewidth', 2)
%         vm = mean(vauxcnt2)-std(vauxcnt2);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60,vm, 'color', 'm')
%         vm = mean(vauxcnt2)+std(vauxcnt2);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60, vm, 'color', 'm')
        
        plot(lcents/60, mean(vauxkir), 'r', 'linewidth', 2)
%         vm = mean(vauxkir)-std(vauxkir);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60,vm, 'color', 'r')
%         vm = mean(vauxkir)+std(vauxkir);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60, vm, 'color', 'r')
    end
end
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel('Duration FW (s)')


%% Plot %ST T4T5
figure,
hold on
% % inds = [5 6 7 8 9 10 1 2 3 4];
inds = [3 4 5 6 7 8 1 2];
for j = 1 : length(pTypesT4T5)
    vauxcnt1 = [];
    vTcnt1 = [];
    a = inds(j);
    for n = 1 : size(stfbcntT4T51.pST,1);
        stfbcntT4T51.pT{n,a}(isnan(stfbcntT4T51.pST{n,a})) = 1;
        stfbcntT4T51.pST{n,a}(isnan(stfbcntT4T51.pST{n,a})) = 0;
        stfbcntT4T51.pT{n,a}(isempty(stfbcntT4T51.pST{n,a})) = 1;
        stfbcntT4T51.pST{n,a}(isempty(stfbcntT4T51.pST{n,a})) = 0;
        vauxcnt1 = vertcat(vauxcnt1, stfbcntT4T51.pST{n,a}'*stfbcntT4T51.pT{n,a}/sum(stfbcntT4T51.pT{n,a}));
        vTcnt1 = vertcat(vTcnt1, sum(stfbcntT4T51.pT{n,a}));
    end
    vTcnt1(isnan(vauxcnt1)) = [];
    vauxcnt1(isnan(vauxcnt1)) = [];
    vauxcnt2 = [];
    vTcnt2 = [];
    a = inds(j);
    for n = 1 : size(stfbcntT4T52.pST,1);
        stfbcntT4T52.pT{n,a}(isnan(stfbcntT4T52.pST{n,a})) = 1;
        stfbcntT4T52.pST{n,a}(isnan(stfbcntT4T52.pST{n,a})) = 0;
        stfbcntT4T52.pT{n,a}(isempty(stfbcntT4T52.pST{n,a})) = 1;
        stfbcntT4T52.pST{n,a}(isempty(stfbcntT4T52.pST{n,a})) = 0;
        vauxcnt2 = vertcat(vauxcnt2, stfbcntT4T52.pST{n,a}'*stfbcntT4T52.pT{n,a}/sum(stfbcntT4T52.pT{n,a}));
        vTcnt2 = vertcat(vTcnt2, sum(stfbcntT4T52.pT{n,a}));
    end
    vTcnt2(isnan(vauxcnt2)) = [];
    vauxcnt2(isnan(vauxcnt2)) = [];
    vauxkir = [];
    vTkir = [];
    a = inds(j);
    for n = 1 : size(stfbkirT4T5.pST,1);
        stfbkirT4T5.pT{n,a}(isnan(stfbkirT4T5.pST{n,a})) = 1;
        stfbkirT4T5.pST{n,a}(isnan(stfbkirT4T5.pST{n,a})) = 0;
        stfbkirT4T5.pT{n,a}(isempty(stfbkirT4T5.pST{n,a})) = 1;
        stfbkirT4T5.pST{n,a}(isempty(stfbkirT4T5.pST{n,a})) = 0;
        vauxkir = vertcat(vauxkir, stfbkirT4T5.pST{n,a}'*stfbkirT4T5.pT{n,a}/sum(stfbkirT4T5.pT{n,a}));
        vTkir = vertcat(vTkir, sum(stfbkirT4T5.pT{n,a}));
    end
    vTkir(isnan(vauxkir)) = [];
    vauxkir(isnan(vauxkir)) = [];
    
    
    if mod(j,2) == 1
%         errorbar(j, vauxcnt1'*vTcnt1/sum(vTcnt1), ...
%             sqrt(sum((vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)).*(vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)))...
%             /size(stfbcntT4T51.pT,1))/sqrt(size(stfbcntT4T51.pT,1)), 'om')
        errorbar(j, vauxcnt2'*vTcnt2/sum(vTcnt2), ...
            sqrt(sum((vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)).*(vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)))...
            /size(stfbcntT4T52.pT,1))/sqrt(size(stfbcntT4T52.pT,1)), 'om')
        
        errorbar(j, vauxkir'*vTkir/sum(vTkir), ...
            sqrt(sum((vauxkir-vauxkir'*vTkir/sum(vTkir)).*(vauxkir-vauxkir'*vTkir/sum(vTkir)))...
            /size(stfbkirT4T5.pT,1))/sqrt(size(stfbkirT4T5.pT,1)), 'or')
    else
%         errorbar(j-1, vauxcnt1'*vTcnt1/sum(vTcnt1),...
%             sqrt(sum((vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)).*(vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)))...
%             /size(stfbcntT4T51.pT,1))/sqrt(size(stfbcntT4T51.pT,1)), 'oc')
        errorbar(j-1, vauxcnt2'*vTcnt2/sum(vTcnt2),...
            sqrt(sum((vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)).*(vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)))...
            /size(stfbcntT4T52.pT,1))/sqrt(size(stfbcntT4T52.pT,1)), 'oc')
        errorbar(j-1, vauxkir'*vTkir/sum(vTkir),...
            sqrt(sum((vauxkir-vauxkir'*vTkir/sum(vTkir)).*(vauxkir-vauxkir'*vTkir/sum(vTkir)))...
            /size(stfbkirT4T5.pT,1))/sqrt(size(stfbkirT4T5.pT,1)), 'ob')
    end
end
axis([0.5 8.5 0.3 0.7])

%% Plot ST Bias T4T5
figure,
hold on
% inds = [5 6 7 8 9 10 1 2 3 4];
inds = [3 4 5 6 7 8 1 2];
for j = 1 : length(pTypesT4T5)
    vauxcnt11 = [];
    vTcnt11 = [];
    a = inds(j);
    for n = 1 : size(stfbcntT4T51.biasST1,1);
        stfbcntT4T51.NbiasST1{n,a}(isnan(stfbcntT4T51.biasST1{n,a})) = 1;
        stfbcntT4T51.biasST1{n,a}(isnan(stfbcntT4T51.biasST1{n,a})) = 0;
        stfbcntT4T51.NbiasST1{n,a}(isempty(stfbcntT4T51.biasST1{n,a})) = 1;
        stfbcntT4T51.biasST1{n,a}(isempty(stfbcntT4T51.biasST1{n,a})) = 0;
        vauxcnt11 = vertcat(vauxcnt11, stfbcntT4T51.biasST1{n,a}'*stfbcntT4T51.NbiasST1{n,a}/sum(stfbcntT4T51.NbiasST1{n,a}));
        vTcnt11 = vertcat(vTcnt11, sum(stfbcntT4T51.NbiasST1{n,a}));
    end
    vTcnt11(isnan(vauxcnt11)) = [];
    vauxcnt11(isnan(vauxcnt11)) = [];
    vauxcnt120 = [];
    vTcnt120 = [];
    a = inds(j);
    for n = 1 : size(stfbcntT4T51.biasST20,1);
        stfbcntT4T51.NbiasST20{n,a}(isnan(stfbcntT4T51.biasST20{n,a})) = 1;
        stfbcntT4T51.biasST20{n,a}(isnan(stfbcntT4T51.biasST20{n,a})) = 0;
        stfbcntT4T51.NbiasST20{n,a}(isempty(stfbcntT4T51.biasST20{n,a})) = 1;
        stfbcntT4T51.biasST20{n,a}(isempty(stfbcntT4T51.biasST20{n,a})) = 0;
        vauxcnt120 = vertcat(vauxcnt120, stfbcntT4T51.biasST20{n,a}'*stfbcntT4T51.NbiasST20{n,a}/sum(stfbcntT4T51.NbiasST20{n,a}));
        vTcnt120 = vertcat(vTcnt120, sum(stfbcntT4T51.NbiasST20{n,a}));
    end
    vTcnt120(isnan(vauxcnt120)) = [];
    vauxcnt120(isnan(vauxcnt120)) = [];
    vauxcnt1100 = [];
    vTcnt1100 = [];
    a = inds(j);
    for n = 1 : size(stfbcntT4T51.biasST100,1);
        stfbcntT4T51.NbiasST100{n,a}(isnan(stfbcntT4T51.biasST100{n,a})) = 1;
        stfbcntT4T51.biasST100{n,a}(isnan(stfbcntT4T51.biasST100{n,a})) = 0;
        stfbcntT4T51.NbiasST100{n,a}(isempty(stfbcntT4T51.biasST100{n,a})) = 1;
        stfbcntT4T51.biasST100{n,a}(isempty(stfbcntT4T51.biasST100{n,a})) = 0;
        vauxcnt1100 = vertcat(vauxcnt1100, stfbcntT4T51.biasST100{n,a}'*stfbcntT4T51.NbiasST100{n,a}/sum(stfbcntT4T51.NbiasST100{n,a}));
        vTcnt1100 = vertcat(vTcnt1100, sum(stfbcntT4T51.NbiasST100{n,a}));
    end
    vTcnt1100(isnan(vauxcnt1100)) = [];
    vauxcnt1100(isnan(vauxcnt1100)) = [];
    vauxcnt21 = [];
    vTcnt21 = [];
    a = inds(j);
    for n = 1 : size(stfbcntT4T52.biasST1,1);
        stfbcntT4T52.NbiasST1{n,a}(isnan(stfbcntT4T52.biasST1{n,a})) = 1;
        stfbcntT4T52.biasST1{n,a}(isnan(stfbcntT4T52.biasST1{n,a})) = 0;
        stfbcntT4T52.NbiasST1{n,a}(isempty(stfbcntT4T52.biasST1{n,a})) = 1;
        stfbcntT4T52.biasST1{n,a}(isempty(stfbcntT4T52.biasST1{n,a})) = 0;
        vauxcnt21 = vertcat(vauxcnt21, stfbcntT4T52.biasST1{n,a}'*stfbcntT4T52.NbiasST1{n,a}/sum(stfbcntT4T52.NbiasST1{n,a}));
        vTcnt21 = vertcat(vTcnt21, sum(stfbcntT4T52.NbiasST1{n,a}));
    end
    vTcnt21(isnan(vauxcnt21)) = [];
    vauxcnt21(isnan(vauxcnt21)) = [];
    vauxcnt220 = [];
    vTcnt220 = [];
    a = inds(j);
    for n = 1 : size(stfbcntT4T52.biasST20,1);
        stfbcntT4T52.NbiasST20{n,a}(isnan(stfbcntT4T52.biasST20{n,a})) = 1;
        stfbcntT4T52.biasST20{n,a}(isnan(stfbcntT4T52.biasST20{n,a})) = 0;
        stfbcntT4T52.NbiasST20{n,a}(isempty(stfbcntT4T52.biasST20{n,a})) = 1;
        stfbcntT4T52.biasST20{n,a}(isempty(stfbcntT4T52.biasST20{n,a})) = 0;
        vauxcnt220 = vertcat(vauxcnt220, stfbcntT4T52.biasST20{n,a}'*stfbcntT4T52.NbiasST20{n,a}/sum(stfbcntT4T52.NbiasST20{n,a}));
        vTcnt220 = vertcat(vTcnt220, sum(stfbcntT4T52.NbiasST20{n,a}));
    end
    vTcnt220(isnan(vauxcnt220)) = [];
    vauxcnt220(isnan(vauxcnt220)) = [];
    vauxcnt2100 = [];
    vTcnt2100 = [];
    a = inds(j);
    for n = 1 : size(stfbcntT4T52.biasST100,1);
        stfbcntT4T52.NbiasST100{n,a}(isnan(stfbcntT4T52.biasST100{n,a})) = 1;
        stfbcntT4T52.biasST100{n,a}(isnan(stfbcntT4T52.biasST100{n,a})) = 0;
        stfbcntT4T52.NbiasST100{n,a}(isempty(stfbcntT4T52.biasST100{n,a})) = 1;
        stfbcntT4T52.biasST100{n,a}(isempty(stfbcntT4T52.biasST100{n,a})) = 0;
        vauxcnt2100 = vertcat(vauxcnt2100, stfbcntT4T52.biasST100{n,a}'*stfbcntT4T52.NbiasST100{n,a}/sum(stfbcntT4T52.NbiasST100{n,a}));
        vTcnt2100 = vertcat(vTcnt2100, sum(stfbcntT4T52.NbiasST100{n,a}));
    end
    vTcnt2100(isnan(vauxcnt2100)) = [];
    vauxcnt2100(isnan(vauxcnt2100)) = [];
    vauxkir1 = [];
    vTkir1 = [];
    a = inds(j);
    for n = 1 : size(stfbkirT4T5.biasST1,1);
        stfbkirT4T5.NbiasST1{n,a}(isnan(stfbkirT4T5.biasST1{n,a})) = 1;
        stfbkirT4T5.biasST1{n,a}(isnan(stfbkirT4T5.biasST1{n,a})) = 0;
        stfbkirT4T5.NbiasST1{n,a}(isempty(stfbkirT4T5.biasST1{n,a})) = 1;
        stfbkirT4T5.biasST1{n,a}(isempty(stfbkirT4T5.biasST1{n,a})) = 0;
        vauxkir1 = vertcat(vauxkir1, stfbkirT4T5.biasST1{n,a}'*stfbkirT4T5.NbiasST1{n,a}/sum(stfbkirT4T5.NbiasST1{n,a}));
        vTkir1 = vertcat(vTkir1, sum(stfbkirT4T5.NbiasST1{n,a}));
    end
    vTkir1(isnan(vauxkir1)) = [];
    vauxkir1(isnan(vauxkir1)) = [];
    vauxkir20 = [];
    vTkir20 = [];
    a = inds(j);
    for n = 1 : size(stfbkirT4T5.biasST20,1);
        stfbkirT4T5.NbiasST20{n,a}(isnan(stfbkirT4T5.biasST20{n,a})) = 1;
        stfbkirT4T5.biasST20{n,a}(isnan(stfbkirT4T5.biasST20{n,a})) = 0;
        stfbkirT4T5.NbiasST20{n,a}(isempty(stfbkirT4T5.biasST20{n,a})) = 1;
        stfbkirT4T5.biasST20{n,a}(isempty(stfbkirT4T5.biasST20{n,a})) = 0;
        vauxkir20 = vertcat(vauxkir20, stfbkirT4T5.biasST20{n,a}'*stfbkirT4T5.NbiasST20{n,a}/sum(stfbkirT4T5.NbiasST20{n,a}));
        vTkir20 = vertcat(vTkir20, sum(stfbkirT4T5.NbiasST20{n,a}));
    end
    vTkir20(isnan(vauxkir20)) = [];
    vauxkir20(isnan(vauxkir20)) = [];
    vauxkir100 = [];
    vTkir100 = [];
    a = inds(j);
    for n = 1 : size(stfbkirT4T5.biasST100,1);
        stfbkirT4T5.NbiasST100{n,a}(isnan(stfbkirT4T5.biasST100{n,a})) = 1;
        stfbkirT4T5.biasST100{n,a}(isnan(stfbkirT4T5.biasST100{n,a})) = 0;
        stfbkirT4T5.NbiasST100{n,a}(isempty(stfbkirT4T5.biasST100{n,a})) = 1;
        stfbkirT4T5.biasST100{n,a}(isempty(stfbkirT4T5.biasST100{n,a})) = 0;
        vauxkir100 = vertcat(vauxkir100, stfbkirT4T5.biasST100{n,a}'*stfbkirT4T5.NbiasST100{n,a}/sum(stfbkirT4T5.NbiasST100{n,a}));
        vTkir100 = vertcat(vTkir100, sum(stfbkirT4T5.NbiasST100{n,a}));
    end
    vTkir100(isnan(vauxkir100)) = [];
    vauxkir100(isnan(vauxkir100)) = [];
    
    dl = 0.35;
    if mod(j,2) == 1
        errorbar(j-dl, vauxcnt11'*vTcnt11/sum(vTcnt11), ...
            sqrt(sum((vauxcnt11-vauxcnt11'*vTcnt11/sum(vTcnt11)).*(vauxcnt11-vauxcnt11'*vTcnt11/sum(vTcnt11)))...
            /size(stfbcntT4T51.biasST1,1))/sqrt(size(stfbcntT4T51.biasST1,1)), 'color', [1 0 1], 'marker', 'o')
        errorbar(j, vauxcnt120'*vTcnt120/sum(vTcnt120), ...
            sqrt(sum((vauxcnt120-vauxcnt120'*vTcnt120/sum(vTcnt120)).*(vauxcnt120-vauxcnt120'*vTcnt120/sum(vTcnt120)))...
            /size(stfbcntT4T51.biasST20,1))/sqrt(size(stfbcntT4T51.biasST20,1)), 'color', [0.75 0 0.75], 'marker', 'o')
        errorbar(j+dl, vauxcnt1100'*vTcnt1100/sum(vTcnt1100), ...
            sqrt(sum((vauxcnt1100-vauxcnt1100'*vTcnt1100/sum(vTcnt1100)).*(vauxcnt1100-vauxcnt1100'*vTcnt1100/sum(vTcnt1100)))...
            /size(stfbcntT4T51.biasST100,1))/sqrt(size(stfbcntT4T51.biasST100,1)), 'color', [0.5 0 0.5], 'marker', 'o')
                
        errorbar(j-dl, vauxcnt21'*vTcnt21/sum(vTcnt21), ...
            sqrt(sum((vauxcnt21-vauxcnt21'*vTcnt21/sum(vTcnt21)).*(vauxcnt21-vauxcnt21'*vTcnt21/sum(vTcnt21)))...
            /size(stfbcntT4T52.biasST1,1))/sqrt(size(stfbcntT4T52.biasST1,1)), 'color', [1 0 1], 'marker', 'o')
        errorbar(j, vauxcnt220'*vTcnt220/sum(vTcnt220), ...
            sqrt(sum((vauxcnt220-vauxcnt220'*vTcnt220/sum(vTcnt220)).*(vauxcnt220-vauxcnt220'*vTcnt220/sum(vTcnt220)))...
            /size(stfbcntT4T52.biasST20,1))/sqrt(size(stfbcntT4T52.biasST20,1)), 'color', [0.75 0 0.75], 'marker', 'o')
        errorbar(j+dl, vauxcnt2100'*vTcnt2100/sum(vTcnt2100), ...
            sqrt(sum((vauxcnt2100-vauxcnt2100'*vTcnt2100/sum(vTcnt2100)).*(vauxcnt2100-vauxcnt2100'*vTcnt2100/sum(vTcnt2100)))...
            /size(stfbcntT4T52.biasST100,1))/sqrt(size(stfbcntT4T52.biasST100,1)), 'color', [0.5 0 0.5], 'marker', 'o')
                
        errorbar(j-dl, vauxkir1'*vTkir1/sum(vTkir1), ...
            sqrt(sum((vauxkir1-vauxkir1'*vTkir1/sum(vTkir1)).*(vauxkir1-vauxkir1'*vTkir1/sum(vTkir1)))...
            /size(stfbkirT4T5.biasST1,1))/sqrt(size(stfbkirT4T5.biasST1,1)), 'color', [1 0 0], 'marker', 'o')
        errorbar(j, vauxkir20'*vTkir20/sum(vTkir20), ...
            sqrt(sum((vauxkir20-vauxkir20'*vTkir20/sum(vTkir20)).*(vauxkir20-vauxkir20'*vTkir20/sum(vTkir20)))...
            /size(stfbkirT4T5.biasST20,1))/sqrt(size(stfbkirT4T5.biasST20,1)), 'color', [0.75 0 0], 'marker', 'o')
        errorbar(j+dl, vauxkir100'*vTkir100/sum(vTkir100), ...
            sqrt(sum((vauxkir100-vauxkir100'*vTkir100/sum(vTkir100)).*(vauxkir100-vauxkir100'*vTkir100/sum(vTkir100)))...
            /size(stfbkirT4T5.biasST100,1))/sqrt(size(stfbkirT4T5.biasST100,1)), 'color', [0.5 0 0], 'marker', 'o')
    else
        errorbar(j-1-dl, vauxcnt11'*vTcnt11/sum(vTcnt11),...
            sqrt(sum((vauxcnt11-vauxcnt11'*vTcnt11/sum(vTcnt11)).*(vauxcnt11-vauxcnt11'*vTcnt11/sum(vTcnt11)))...
            /size(stfbcntT4T51.biasST1,1))/sqrt(size(stfbcntT4T51.biasST1,1)), 'color', [0 1 1], 'marker', 'o')
        errorbar(j-1, vauxcnt120'*vTcnt120/sum(vTcnt120),...
            sqrt(sum((vauxcnt120-vauxcnt120'*vTcnt120/sum(vTcnt120)).*(vauxcnt120-vauxcnt120'*vTcnt120/sum(vTcnt120)))...
            /size(stfbcntT4T51.biasST20,1))/sqrt(size(stfbcntT4T51.biasST20,1)), 'color', [0 0.75 0.75], 'marker', 'o')      
        errorbar(j-1+dl, vauxcnt1100'*vTcnt1100/sum(vTcnt1100),...
            sqrt(sum((vauxcnt1100-vauxcnt1100'*vTcnt1100/sum(vTcnt1100)).*(vauxcnt1100-vauxcnt1100'*vTcnt1100/sum(vTcnt1100)))...
            /size(stfbcntT4T51.biasST100,1))/sqrt(size(stfbcntT4T51.biasST100,1)), 'color', [0 0.5 0.5], 'marker', 'o')      
        
        errorbar(j-1-dl, vauxcnt21'*vTcnt21/sum(vTcnt21),...
            sqrt(sum((vauxcnt21-vauxcnt21'*vTcnt21/sum(vTcnt21)).*(vauxcnt21-vauxcnt21'*vTcnt21/sum(vTcnt21)))...
            /size(stfbcntT4T52.biasST1,1))/sqrt(size(stfbcntT4T52.biasST1,1)), 'color', [0 1 1], 'marker', 'o')
        errorbar(j-1, vauxcnt220'*vTcnt220/sum(vTcnt220),...
            sqrt(sum((vauxcnt220-vauxcnt220'*vTcnt220/sum(vTcnt220)).*(vauxcnt220-vauxcnt220'*vTcnt220/sum(vTcnt220)))...
            /size(stfbcntT4T52.biasST20,1))/sqrt(size(stfbcntT4T52.biasST20,1)), 'color', [0 0.75 0.75], 'marker', 'o')        
        errorbar(j-1+dl, vauxcnt2100'*vTcnt2100/sum(vTcnt2100),...
            sqrt(sum((vauxcnt2100-vauxcnt2100'*vTcnt2100/sum(vTcnt2100)).*(vauxcnt2100-vauxcnt2100'*vTcnt2100/sum(vTcnt2100)))...
            /size(stfbcntT4T52.biasST100,1))/sqrt(size(stfbcntT4T52.biasST100,1)), 'color', [0 0.5 0.5], 'marker', 'o')      

        errorbar(j-1-dl, vauxkir1'*vTkir1/sum(vTkir1),...
            sqrt(sum((vauxkir1-vauxkir1'*vTkir1/sum(vTkir1)).*(vauxkir1-vauxkir1'*vTkir1/sum(vTkir1)))...
            /size(stfbkirT4T5.biasST1,1))/sqrt(size(stfbkirT4T5.biasST1,1)), 'color', [0 0 1], 'marker', 'o')
        errorbar(j-1, vauxkir20'*vTkir20/sum(vTkir20),...
            sqrt(sum((vauxkir20-vauxkir20'*vTkir20/sum(vTkir20)).*(vauxkir20-vauxkir20'*vTkir20/sum(vTkir20)))...
            /size(stfbkirT4T5.biasST20,1))/sqrt(size(stfbkirT4T5.biasST20,1)), 'color', [0 0 0.75], 'marker', 'o')
        errorbar(j-1+dl, vauxkir100'*vTkir100/sum(vTkir100),...
            sqrt(sum((vauxkir100-vauxkir100'*vTkir100/sum(vTkir100)).*(vauxkir100-vauxkir100'*vTkir100/sum(vTkir100)))...
            /size(stfbkirT4T5.biasST100,1))/sqrt(size(stfbkirT4T5.biasST100,1)), 'color', [0 0 0.5], 'marker', 'o')
    end
end
axis([0.5 8.5 0 1])

%% Plot ST max T4T5
figure,
hold on
lcents = -50:20:1500;
% inds = [5 6 7 8 9 10 1 2 3 4];
inds = [3 4 5 6 7 8 1 2];
cmap1 = autumn(15);
cmap2 = winter(15);
cmap3 = spring(15);
cmap4 = summer(15);
for j = 1 : length(pTypesT4T5)
    vaux1 = [];
    a = inds(j);
    for n = 1 : size(stfbcntT4T51.lST,1);
        vaux1 = vertcat(vaux1, ...
            hist(abs(stfbcntT4T51.lST{n,a}), lcents)/sum(hist(stfbcntT4T51.lST{n,a}, lcents)));
    end
    vaux2 = [];
    a = inds(j);
    for n = 1 : size(stfbcntT4T52.lST,1);
        vaux2 = vertcat(vaux2, ...
            hist(abs(stfbcntT4T52.lST{n,a}), lcents)/sum(hist(stfbcntT4T52.lST{n,a}, lcents)));
    end
    vauxk = [];
    a = inds(j);
    for n = 1 : size(stfbkirT4T5.lST,1);
        vauxk = vertcat(vauxk, ...
            hist(abs(stfbkirT4T5.lST{n,a}), lcents)/sum(hist(stfbkirT4T5.lST{n,a}, lcents)));
    end
    
    if mod(j,2) == 1
        plot(lcents, smooth(mean(vaux1)), 'color', cmap1(j,:), 'linewidth', 2)
        plot(lcents, smooth(mean(vaux2)), 'color', cmap1(j,:), 'linewidth', 2)
        plot(lcents, smooth(mean(vauxk)), 'color', cmap3(j,:), 'linewidth', 2)
    else
        plot(lcents, smooth(mean(vaux1)), 'color', cmap2(j,:), 'linewidth', 2)
        plot(lcents, smooth(mean(vaux2)), 'color', cmap2(j,:), 'linewidth', 2)
        plot(lcents, smooth(mean(vauxk)), 'color', cmap4(j,:), 'linewidth', 2)
    end
end
axis([0 1500 0 0.05])











%% Plot %FB HSVS
figure,
hold on
inds = [3 4 5 6 7 8 1 2];
for j = 1 : length(pTypesHS)
    vauxcnt1 = [];
    vTcnt1 = [];
    a = inds(j);
    for n = 1 : size(stfbcntHS1.pFB,1);
        stfbcntHS1.pT{n,a}(isnan(stfbcntHS1.pFB{n,a})) = 1;
        stfbcntHS1.pFB{n,a}(isnan(stfbcntHS1.pFB{n,a})) = 0;
        stfbcntHS1.pT{n,a}(isempty(stfbcntHS1.pFB{n,a})) = 1;
        stfbcntHS1.pFB{n,a}(isempty(stfbcntHS1.pFB{n,a})) = 0;
        vauxcnt1 = vertcat(vauxcnt1, stfbcntHS1.pFB{n,a}'*stfbcntHS1.pT{n,a}/sum(stfbcntHS1.pT{n,a}));
        vTcnt1 = vertcat(vTcnt1, sum(stfbcntHS1.pT{n,a}));
    end
    vTcnt1(isnan(vauxcnt1)) = [];
    vauxcnt1(isnan(vauxcnt1)) = [];
    vauxcnt2 = [];
    vTcnt2 = [];
    a = inds(j);
    for n = 1 : size(stfbcntHS2.pFB,1);
        stfbcntHS2.pT{n,a}(isnan(stfbcntHS2.pFB{n,a})) = 1;
        stfbcntHS2.pFB{n,a}(isnan(stfbcntHS2.pFB{n,a})) = 0;
        stfbcntHS2.pT{n,a}(isempty(stfbcntHS2.pFB{n,a})) = 1;
        stfbcntHS2.pFB{n,a}(isempty(stfbcntHS2.pFB{n,a})) = 0;
        vauxcnt2 = vertcat(vauxcnt2, stfbcntHS2.pFB{n,a}'*stfbcntHS2.pT{n,a}/sum(stfbcntHS2.pT{n,a}));
        vTcnt2 = vertcat(vTcnt2, sum(stfbcntHS2.pT{n,a}));
    end
    vTcnt2(isnan(vauxcnt2)) = [];
    vauxcnt2(isnan(vauxcnt2)) = [];
    vauxkir = [];
    vTkir = [];
    a = inds(j);
    for n = 1 : size(stfbkirHS.pFB,1);
        stfbkirHS.pT{n,a}(isnan(stfbkirHS.pFB{n,a})) = 1;
        stfbkirHS.pFB{n,a}(isnan(stfbkirHS.pFB{n,a})) = 0;
        stfbkirHS.pT{n,a}(isempty(stfbkirHS.pFB{n,a})) = 1;
        stfbkirHS.pFB{n,a}(isempty(stfbkirHS.pFB{n,a})) = 0;
        vauxkir = vertcat(vauxkir, stfbkirHS.pFB{n,a}'*stfbkirHS.pT{n,a}/sum(stfbkirHS.pT{n,a}));
        vTkir = vertcat(vTkir, sum(stfbkirHS.pT{n,a}));
    end
    vTkir(isnan(vauxkir)) = [];
    vauxkir(isnan(vauxkir)) = [];
    
    
    if mod(j,2) == 1
        errorbar(j, vauxcnt1'*vTcnt1/sum(vTcnt1), ...
            sqrt(sum((vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)).*(vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)))...
            /size(stfbcntHS1.pFB,1))/sqrt(size(stfbcntHS1.pFB,1)), 'om')
        errorbar(j, vauxcnt2'*vTcnt2/sum(vTcnt2), ...
            sqrt(sum((vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)).*(vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)))...
            /size(stfbcntHS2.pFB,1))/sqrt(size(stfbcntHS2.pFB,1)), 'om')
        
        errorbar(j, vauxkir'*vTkir/sum(vTkir), ...
            sqrt(sum((vauxkir-vauxkir'*vTkir/sum(vTkir)).*(vauxkir-vauxkir'*vTkir/sum(vTkir)))...
            /size(stfbkirHS.pFB,1))/sqrt(size(stfbkirHS.pFB,1)), 'or')
    else
        errorbar(j-1, vauxcnt1'*vTcnt1/sum(vTcnt1),...
            sqrt(sum((vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)).*(vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)))...
            /size(stfbcntHS1.pFB,1))/sqrt(size(stfbcntHS1.pFB,1)), 'oc')
        errorbar(j-1, vauxcnt2'*vTcnt2/sum(vTcnt2),...
            sqrt(sum((vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)).*(vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)))...
            /size(stfbcntHS2.pFB,1))/sqrt(size(stfbcntHS2.pFB,1)), 'oc')
        errorbar(j-1, vauxkir'*vTkir/sum(vTkir),...
            sqrt(sum((vauxkir-vauxkir'*vTkir/sum(vTkir)).*(vauxkir-vauxkir'*vTkir/sum(vTkir)))...
            /size(stfbkirHS.pFB,1))/sqrt(size(stfbkirHS.pFB,1)), 'ob')
    end
end
axis([0.5 8.5 0.1 0.5])
xlabel('Dot Size(º)')
ylabel('% Time FW')

%% Plot %Bias FB HSVS
clc
figure,
hold on
inds = [3 4 5 6 7 8 1 2];
thst = 350;
for j = 1 : length(pTypesHS)
    vauxcnt1 = [];
    vTcnt1 = [];
    a = inds(j);
    for n = 1 : size(stfbcntHS1.biasFB,1);
        stfbcntHS1.NbiasFB{n,a}(isnan(stfbcntHS1.biasFB{n,a})) = 1;
        stfbcntHS1.biasFB{n,a}(isnan(stfbcntHS1.biasFB{n,a})) = 0;
        stfbcntHS1.NbiasFB{n,a}(isempty(stfbcntHS1.biasFB{n,a})) = 1;
        stfbcntHS1.biasFB{n,a}(isempty(stfbcntHS1.biasFB{n,a})) = 0;
        vauxcnt1 = vertcat(vauxcnt1, max(stfbcntHS1.biasFB{n,a}'*stfbcntHS1.NbiasFB{n,a}/sum(stfbcntHS1.NbiasFB{n,a})));
        vTcnt1 = vertcat(vTcnt1, sum(stfbcntHS1.NbiasFB{n,a}));
    end
    vTcnt1(isnan(vauxcnt1)) = [];
    vauxcnt1(isnan(vauxcnt1)) = [];
    if length(vauxcnt1) ~= length(vTcnt1)
        vauxcnt1(vauxcnt1==0) = [];
    end
    vauxcnt1(vTcnt1<thst) = [];
    vTcnt1(vTcnt1<thst) = [];
    
    vauxcnt2 = [];
    vTcnt2 = [];
    a = inds(j);
    for n = 1 : size(stfbcntHS2.biasFB,1);
        stfbcntHS2.NbiasFB{n,a}(isnan(stfbcntHS2.biasFB{n,a})) = 1;
        stfbcntHS2.biasFB{n,a}(isnan(stfbcntHS2.biasFB{n,a})) = 0;
        stfbcntHS2.NbiasFB{n,a}(isempty(stfbcntHS2.biasFB{n,a})) = 1;
        stfbcntHS2.biasFB{n,a}(isempty(stfbcntHS2.biasFB{n,a})) = 0;
        vauxcnt2 = vertcat(vauxcnt2, max(stfbcntHS2.biasFB{n,a}'*stfbcntHS2.NbiasFB{n,a}/sum(stfbcntHS2.NbiasFB{n,a})));
        vTcnt2 = vertcat(vTcnt2, sum(stfbcntHS2.NbiasFB{n,a}));
    end
    vTcnt2(isnan(vauxcnt2)) = [];
    vauxcnt2(isnan(vauxcnt2)) = [];
    if length(vauxcnt2) ~= length(vTcnt2)
        vauxcnt2(vauxcnt2==0) = [];
    end
    vauxcnt2(vTcnt2<thst) = [];
    vTcnt2(vTcnt2<thst) = [];
    
    vauxkir = [];
    vTkir = [];
    a = inds(j);
    for n = 1 : size(stfbkirHS.biasFB,1);
        stfbkirHS.NbiasFB{n,a}(isnan(stfbkirHS.biasFB{n,a})) = 1;
        stfbkirHS.biasFB{n,a}(isnan(stfbkirHS.biasFB{n,a})) = 0;
        stfbkirHS.NbiasFB{n,a}(isempty(stfbkirHS.biasFB{n,a})) = 1;
        stfbkirHS.biasFB{n,a}(isempty(stfbkirHS.biasFB{n,a})) = 0;
        vauxkir = vertcat(vauxkir, max(stfbkirHS.biasFB{n,a}'*stfbkirHS.NbiasFB{n,a}/sum(stfbkirHS.NbiasFB{n,a})));
        vTkir = vertcat(vTkir, sum(stfbkirHS.NbiasFB{n,a}));
    end
    vTkir(isnan(vauxkir)) = [];
    vauxkir(isnan(vauxkir)) = [];
    if length(vauxkir) ~= length(vTkir)
        vauxkir(vauxkir==0) = [];
    end
    vauxkir(vTkir<thst) = [];
    vTkir(vTkir<thst) = [];
    
    if mod(j,2) == 1
        errorbar(j, vauxcnt1'*vTcnt1/sum(vTcnt1), ...
            sqrt(sum((vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)).*(vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)))...
            /size(stfbcntHS1.NbiasFB,1))/sqrt(size(stfbcntHS1.NbiasFB,1)), 'om', 'markersize', 7, 'markerfacecolor', 'm')
        errorbar(j, vauxcnt2'*vTcnt2/sum(vTcnt2), ...
            sqrt(sum((vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)).*(vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)))...
            /size(stfbcntHS2.NbiasFB,1))/sqrt(size(stfbcntHS2.NbiasFB,1)), 'om', 'markersize', 7, 'markerfacecolor', 'm')
        errorbar(j, vauxkir'*vTkir/sum(vTkir), ...
            sqrt(sum((vauxkir-vauxkir'*vTkir/sum(vTkir)).*(vauxkir-vauxkir'*vTkir/sum(vTkir)))...
            /size(stfbkirHS.NbiasFB,1))/sqrt(size(stfbkirHS.NbiasFB,1)), 'or', 'markersize', 7, 'markerfacecolor', 'r')
        
        statTestmww(vauxkir, vTkir, vauxcnt1, vTcnt1, thst, pTypesHS{a}, 'HSKIR', 'HSCnt1');
        statTestmww(vauxkir, vTkir, vauxcnt2, vTcnt2, thst, pTypesHS{a}, 'HSKIR', 'HSCnt2');
        statTestmww(vauxcnt2, vTcnt2, vauxcnt1, vTcnt1, thst, pTypesHS{a}, 'HSCnt2', 'HSCnt1');
    else
        errorbar(j-1, vauxcnt1'*vTcnt1/sum(vTcnt1),...
            sqrt(sum((vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)).*(vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)))...
            /size(stfbcntHS1.NbiasFB,1))/sqrt(size(stfbcntHS1.NbiasFB,1)), 'oc', 'markersize', 7, 'markerfacecolor', 'c')
        errorbar(j-1, vauxcnt2'*vTcnt2/sum(vTcnt2),...
            sqrt(sum((vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)).*(vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)))...
            /size(stfbcntHS2.NbiasFB,1))/sqrt(size(stfbcntHS2.NbiasFB,1)), 'oc', 'markersize', 7, 'markerfacecolor', 'c')
        errorbar(j-1, vauxkir'*vTkir/sum(vTkir),...
            sqrt(sum((vauxkir-vauxkir'*vTkir/sum(vTkir)).*(vauxkir-vauxkir'*vTkir/sum(vTkir)))...
            /size(stfbkirHS.NbiasFB,1))/sqrt(size(stfbkirHS.NbiasFB,1)), 'ob', 'markersize', 7, 'markerfacecolor', 'b')
        
        statTestmww(vauxkir, vTkir, vauxcnt1, vTcnt1, thst, pTypesHS{a}, 'HSKIR', 'HSCnt1');
        statTestmww(vauxkir, vTkir, vauxcnt2, vTcnt2, thst, pTypesHS{a}, 'HSKIR', 'HSCnt2');
        statTestmww(vauxcnt2, vTcnt2, vauxcnt1, vTcnt1, thst, pTypesHS{a}, 'HSCnt2', 'HSCnt1');
    end
    disp(' ')
end
axis([0.5 8.5 0.5 1])
xlabel('Dot Size(º)')
ylabel('% Bias FW')


%% Plot Size Bias FB HSVS
clc
figure,
hold on
inds = [3 4 5 6 7 8 1 2];
vc1 = cell(length(pTypesHS),1);
vc2 = cell(length(pTypesHS),1);
vk = cell(length(pTypesHS),1);
thst = 350;
for j = 1 : length(pTypesHS)
    vauxcnt1 = [];
    vTcnt1 = [];
    a = inds(j);
    for n = 1 : size(stfbcntHS1.valFB,1);
        stfbcntHS1.NvalFB{n,a}(isnan(stfbcntHS1.valFB{n,a})) = 1;
        stfbcntHS1.valFB{n,a}(isnan(stfbcntHS1.valFB{n,a})) = 0;
        stfbcntHS1.NvalFB{n,a}(isempty(stfbcntHS1.valFB{n,a})) = 1;
        stfbcntHS1.valFB{n,a}(isempty(stfbcntHS1.valFB{n,a})) = 0;
        vauxcnt1 = vertcat(vauxcnt1, max(stfbcntHS1.valFB{n,a}'*stfbcntHS1.NvalFB{n,a}/sum(stfbcntHS1.NvalFB{n,a})));
        vTcnt1 = vertcat(vTcnt1, sum(stfbcntHS1.NvalFB{n,a}));
    end
    vTcnt1(isnan(vauxcnt1)) = [];
    vauxcnt1(isnan(vauxcnt1)) = [];
    if length(vauxcnt1) ~= length(vTcnt1)
        vauxcnt1(vauxcnt1==0) = [];
    end
    vauxcnt1(vTcnt1<thst) = [];
    vTcnt1(vTcnt1<thst) = [];
    
    vauxcnt2 = [];
    vTcnt2 = [];
    a = inds(j);
    for n = 1 : size(stfbcntHS2.valFB,1);
        stfbcntHS2.NvalFB{n,a}(isnan(stfbcntHS2.valFB{n,a})) = 1;
        stfbcntHS2.valFB{n,a}(isnan(stfbcntHS2.valFB{n,a})) = 0;
        stfbcntHS2.NvalFB{n,a}(isempty(stfbcntHS2.valFB{n,a})) = 1;
        stfbcntHS2.valFB{n,a}(isempty(stfbcntHS2.valFB{n,a})) = 0;
        vauxcnt2 = vertcat(vauxcnt2, max(stfbcntHS2.valFB{n,a}'*stfbcntHS2.NvalFB{n,a}/sum(stfbcntHS2.NvalFB{n,a})));
        vTcnt2 = vertcat(vTcnt2, sum(stfbcntHS2.NvalFB{n,a}));
    end
    vTcnt2(isnan(vauxcnt2)) = [];
    vauxcnt2(isnan(vauxcnt2)) = [];
    if length(vauxcnt2) ~= length(vTcnt2)
        vauxcnt2(vauxcnt2==0) = [];
    end
    vauxcnt2(vTcnt2<thst) = [];
    vTcnt2(vTcnt2<thst) = [];
    
    vauxkir = [];
    vTkir = [];
    a = inds(j);
    for n = 1 : size(stfbkirHS.valFB,1);
        stfbkirHS.NvalFB{n,a}(isnan(stfbkirHS.valFB{n,a})) = 1;
        stfbkirHS.valFB{n,a}(isnan(stfbkirHS.valFB{n,a})) = 0;
        stfbkirHS.NvalFB{n,a}(isempty(stfbkirHS.valFB{n,a})) = 1;
        stfbkirHS.valFB{n,a}(isempty(stfbkirHS.valFB{n,a})) = 0;
        vauxkir = vertcat(vauxkir, max(stfbkirHS.valFB{n,a}'*stfbkirHS.NvalFB{n,a}/sum(stfbkirHS.NvalFB{n,a})));
        vTkir = vertcat(vTkir, sum(stfbkirHS.NvalFB{n,a}));
    end
    vTkir(isnan(vauxkir)) = [];
    vauxkir(isnan(vauxkir)) = [];
    if length(vauxkir) ~= length(vTkir)
        vauxkir(vauxkir==0) = [];
    end
    vauxkir(vTkir<thst) = [];
    vTkir(vTkir<thst) = [];
    
    if mod(j,2) == 1
        errorbar(j, vauxcnt1'*vTcnt1/sum(vTcnt1), ...
            sqrt(sum((vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)).*(vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)))...
            /size(stfbcntHS1.NvalFB,1))/sqrt(size(stfbcntHS1.NvalFB,1)), 'om', 'markersize', 7, 'markerfacecolor', 'm')
        errorbar(j, vauxcnt2'*vTcnt2/sum(vTcnt2), ...
            sqrt(sum((vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)).*(vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)))...
            /size(stfbcntHS2.NvalFB,1))/sqrt(size(stfbcntHS2.NvalFB,1)), 'om', 'markersize', 7, 'markerfacecolor', 'm')
        errorbar(j, vauxkir'*vTkir/sum(vTkir), ...
            sqrt(sum((vauxkir-vauxkir'*vTkir/sum(vTkir)).*(vauxkir-vauxkir'*vTkir/sum(vTkir)))...
            /size(stfbkirHS.NvalFB,1))/sqrt(size(stfbkirHS.NvalFB,1)), 'or', 'markersize', 7, 'markerfacecolor', 'r')
        
        statTestmww(vauxkir, vTkir, vauxcnt1, vTcnt1, thst, pTypesHS{a}, 'HSKIR', 'HSCnt1');
        statTestmww(vauxkir, vTkir, vauxcnt2, vTcnt2, thst, pTypesHS{a}, 'HSKIR', 'HSCnt2');
        statTestmww(vauxcnt2, vTcnt2, vauxcnt1, vTcnt1, thst, pTypesHS{a}, 'HSCnt2', 'HSCnt1');
    else
        errorbar(j-1, vauxcnt1'*vTcnt1/sum(vTcnt1),...
            sqrt(sum((vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)).*(vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)))...
            /size(stfbcntHS1.NvalFB,1))/sqrt(size(stfbcntHS1.NvalFB,1)), 'oc', 'markersize', 7, 'markerfacecolor', 'c')
        errorbar(j-1, vauxcnt2'*vTcnt2/sum(vTcnt2),...
            sqrt(sum((vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)).*(vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)))...
            /size(stfbcntHS2.NvalFB,1))/sqrt(size(stfbcntHS2.NvalFB,1)), 'oc', 'markersize', 7, 'markerfacecolor', 'c')
        errorbar(j-1, vauxkir'*vTkir/sum(vTkir),...
            sqrt(sum((vauxkir-vauxkir'*vTkir/sum(vTkir)).*(vauxkir-vauxkir'*vTkir/sum(vTkir)))...
            /size(stfbkirHS.NvalFB,1))/sqrt(size(stfbkirHS.NvalFB,1)), 'ob', 'markersize', 7, 'markerfacecolor', 'b')
        
        statTestmww(vauxkir, vTkir, vauxcnt1, vTcnt1, thst, pTypesHS{a}, 'HSKIR', 'HSCnt1');
        statTestmww(vauxkir, vTkir, vauxcnt2, vTcnt2, thst, pTypesHS{a}, 'HSKIR', 'HSCnt2');
        statTestmww(vauxcnt2, vTcnt2, vauxcnt1, vTcnt1, thst, pTypesHS{a}, 'HSCnt2', 'HSCnt1');
    end
    disp(' ')
end
axis([0.5 8.5 0 270])
xlabel('Dot Size(º)')
ylabel('Drift (º/s)')


%% Plot Size AngDisp FB HSVS
figure,
hold on
inds = [3 4 5 6 7 8 1 2];
vc1 = cell(length(pTypesHS),1);
vc2 = cell(length(pTypesHS),1);
vk = cell(length(pTypesHS),1);
for j = 1 : length(pTypesHS)
    vauxcnt1 = [];
    vTcnt1 = [];
    a = inds(j);
    for n = 1 : size(stfbcntHS1.angDisp,1);
        stfbcntHS1.NangDisp{n,a}(isnan(stfbcntHS1.angDisp{n,a})) = 1;
        stfbcntHS1.angDisp{n,a}(isnan(stfbcntHS1.angDisp{n,a})) = 0;
        stfbcntHS1.NangDisp{n,a}(isempty(stfbcntHS1.angDisp{n,a})) = 1;
        stfbcntHS1.angDisp{n,a}(isempty(stfbcntHS1.angDisp{n,a})) = 0;
        vauxcnt1 = vertcat(vauxcnt1, max(stfbcntHS1.angDisp{n,a}'*stfbcntHS1.NangDisp{n,a}/sum(stfbcntHS1.NangDisp{n,a})));
        vTcnt1 = vertcat(vTcnt1, sum(stfbcntHS1.NangDisp{n,a}));
    end
    vTcnt1(isnan(vauxcnt1)) = [];
    vauxcnt1(isnan(vauxcnt1)) = [];
    if length(vauxcnt1) ~= length(vTcnt1)
        vauxcnt1(vauxcnt1==0) = [];
    end
    vc1{j} = vauxcnt1(vauxcnt1~=0);
    
    vauxcnt2 = [];
    vTcnt2 = [];
    a = inds(j);
    for n = 1 : size(stfbcntHS2.angDisp,1);
        stfbcntHS2.NangDisp{n,a}(isnan(stfbcntHS2.angDisp{n,a})) = 1;
        stfbcntHS2.angDisp{n,a}(isnan(stfbcntHS2.angDisp{n,a})) = 0;
        stfbcntHS2.NangDisp{n,a}(isempty(stfbcntHS2.angDisp{n,a})) = 1;
        stfbcntHS2.angDisp{n,a}(isempty(stfbcntHS2.angDisp{n,a})) = 0;
        vauxcnt2 = vertcat(vauxcnt2, max(stfbcntHS2.angDisp{n,a}'*stfbcntHS2.NangDisp{n,a}/sum(stfbcntHS2.NangDisp{n,a})));
        vTcnt2 = vertcat(vTcnt2, sum(stfbcntHS2.NangDisp{n,a}));
    end
    vTcnt2(isnan(vauxcnt2)) = [];
    vauxcnt2(isnan(vauxcnt2)) = [];
    if length(vauxcnt2) ~= length(vTcnt2)
        vauxcnt2(vauxcnt2==0) = [];
    end
    vc2{j} = vauxcnt2(vauxcnt2~=0);
    
    vauxkir = [];
    vTkir = [];
    a = inds(j);
    for n = 1 : size(stfbkirHS.angDisp,1);
        stfbkirHS.NangDisp{n,a}(isnan(stfbkirHS.angDisp{n,a})) = 1;
        stfbkirHS.angDisp{n,a}(isnan(stfbkirHS.angDisp{n,a})) = 0;
        stfbkirHS.NangDisp{n,a}(isempty(stfbkirHS.angDisp{n,a})) = 1;
        stfbkirHS.angDisp{n,a}(isempty(stfbkirHS.angDisp{n,a})) = 0;
        vauxkir = vertcat(vauxkir, max(stfbkirHS.angDisp{n,a}'*stfbkirHS.NangDisp{n,a}/sum(stfbkirHS.NangDisp{n,a})));
        vTkir = vertcat(vTkir, sum(stfbkirHS.NangDisp{n,a}));
    end
    vTkir(isnan(vauxkir)) = [];
    vauxkir(isnan(vauxkir)) = [];
    if length(vauxkir) ~= length(vTkir)
        vauxkir(vauxkir==0) = [];
    end
    vk{j} = vauxkir(vauxkir~=0);
    
    if mod(j,2) == 1
        errorbar(j, vauxcnt1'*vTcnt1/sum(vTcnt1), ...
            sqrt(sum((vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)).*(vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)))...
            /size(stfbcntHS1.NangDisp,1))/sqrt(size(stfbcntHS1.NangDisp,1)), 'om')
        errorbar(j, vauxcnt2'*vTcnt2/sum(vTcnt2), ...
            sqrt(sum((vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)).*(vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)))...
            /size(stfbcntHS2.NangDisp,1))/sqrt(size(stfbcntHS2.NangDisp,1)), 'om')
        errorbar(j, vauxkir'*vTkir/sum(vTkir), ...
            sqrt(sum((vauxkir-vauxkir'*vTkir/sum(vTkir)).*(vauxkir-vauxkir'*vTkir/sum(vTkir)))...
            /size(stfbkirHS.NangDisp,1))/sqrt(size(stfbkirHS.NangDisp,1)), 'or')
        
        statTestmww(vauxkir, vTkir, vauxcnt1, vTcnt1, 350, pTypesHS{j}, 'HSKIR', 'HSCnt1');
        statTestmww(vauxkir, vTkir, vauxcnt2, vTcnt2, 350, pTypesHS{j}, 'HSKIR', 'HSCnt2');
        statTestmww(vauxcnt2, vTcnt2, vauxcnt1, vTcnt1, 350, pTypesHS{j}, 'HSCnt2', 'HSCnt1');
    else
        errorbar(j-1, vauxcnt1'*vTcnt1/sum(vTcnt1),...
            sqrt(sum((vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)).*(vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)))...
            /size(stfbcntHS1.NangDisp,1))/sqrt(size(stfbcntHS1.NangDisp,1)), 'oc')
        errorbar(j-1, vauxcnt2'*vTcnt2/sum(vTcnt2),...
            sqrt(sum((vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)).*(vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)))...
            /size(stfbcntHS2.NangDisp,1))/sqrt(size(stfbcntHS2.NangDisp,1)), 'oc')
        errorbar(j-1, vauxkir'*vTkir/sum(vTkir),...
            sqrt(sum((vauxkir-vauxkir'*vTkir/sum(vTkir)).*(vauxkir-vauxkir'*vTkir/sum(vTkir)))...
            /size(stfbkirHS.NangDisp,1))/sqrt(size(stfbkirHS.NangDisp,1)), 'ob')
        
        statTestmww(vauxkir, vTkir, vauxcnt1, vTcnt1, 350, pTypesHS{j}, 'HSKIR', 'HSCnt1');
        statTestmww(vauxkir, vTkir, vauxcnt2, vTcnt2, 350, pTypesHS{j}, 'HSKIR', 'HSCnt2');
        statTestmww(vauxcnt2, vTcnt2, vauxcnt1, vTcnt1, 350, pTypesHS{j}, 'HSCnt2', 'HSCnt1');
    end
end
axis([0.5 8.5 0 20])
xlabel('Dot Size(º)')
ylabel('AngDisp (º/mm)')



%% Plot Size FB HSVS
figure,
hold on
inds = [3 4 5 6 7 8 1 2];
lcents = 1:1:150;
cmap1 = autumn(20);
cmap2 = winter(20);
for j = 1 %: length(pTypesHS)
    vauxcnt1 = [];
    vauxcnt2 = [];
    vauxkir = [];
    a = j;%inds(j);
    for n = 1 : size(stfbcntHS1.lFB,1);
        vauxcnt1 = vertcat(vauxcnt1, ...
            hist(stfbcntHS1.lFB{n,a}, lcents)/sum(hist(stfbcntHS1.lFB{n,a}, lcents)));
    end
    for n = 1 : size(stfbcntHS2.lFB,1);
        vauxcnt2 = vertcat(vauxcnt2, ...
            hist(stfbcntHS2.lFB{n,a}, lcents)/sum(hist(stfbcntHS2.lFB{n,a}, lcents)));
    end
    for n = 1 : size(stfbkirHS.lFB,1);
        vauxkir = vertcat(vauxkir, ...
            hist(stfbkirHS.lFB{n,a}, lcents)/sum(hist(stfbkirHS.lFB{n,a}, lcents)));
    end
    
    if mod(j,2) == 1
        plot(lcents/60, nanmean(vauxcnt1), 'c', 'linewidth', 2)
%         vm = mean(vauxcnt1)-std(vauxcnt1);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60, vm, 'color', 'c')
%         vm = mean(vauxcnt1)+std(vauxcnt1);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60, vm, 'color', 'c')
        
        plot(lcents/60, nanmean(vauxcnt2), 'c', 'linewidth', 2)
%         vm = mean(vauxcnt2)-std(vauxcnt2);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60, vm, 'color', 'c')
%         vm = mean(vauxcnt2)+std(vauxcnt2);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60, vm, 'color', 'c')
        
        plot(lcents/60, nanmean(vauxkir), 'b', 'linewidth', 2)
%         vm = mean(vauxkir)-std(vauxkir);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60, vm, 'color', 'b')
%         vm = mean(vauxkir)+std(vauxkir);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60, vm, 'color', 'b')
    else
        plot(lcents/60, mean(vauxcnt1), 'm', 'linewidth', 2)
%         vm = mean(vauxcnt1)-std(vauxcnt1);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60,vm, 'color', 'm')
%         vm = mean(vauxcnt1)+std(vauxcnt1);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60, vm, 'color', 'm')
        
        plot(lcents/60, mean(vauxcnt2), 'm', 'linewidth', 2)
%         vm = mean(vauxcnt2)-std(vauxcnt2);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60,vm, 'color', 'm')
%         vm = mean(vauxcnt2)+std(vauxcnt2);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60, vm, 'color', 'm')
        
        plot(lcents/60, mean(vauxkir), 'r', 'linewidth', 2)
%         vm = mean(vauxkir)-std(vauxkir);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60,vm, 'color', 'r')
%         vm = mean(vauxkir)+std(vauxkir);
%         vm(vm<=0) = 10^-4;
%         plot(lcents/60, vm, 'color', 'r')
    end
end
% set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel('Duration FW (s)')
ylabel('Average Distribution')
axis([0 2 0 0.2])

vauxcnt1(isnan(vauxcnt1))=0;
vauxcnt2(isnan(vauxcnt2))=0;
p1 = mwwtest(vauxcnt1*lcents', vauxkir*lcents',0);
p2 = mwwtest(vauxcnt2*lcents', vauxkir*lcents',0);
disp(num2str(p1.p(2)))
disp(num2str(p2.p(2)))
%% Plot %ST HSVS
figure,
hold on
inds = [3 4 5 6 7 8 1 2];
for j = 1 : length(pTypesHS)
    vauxcnt1 = [];
    vTcnt1 = [];
    a = inds(j);
    for n = 1 : size(stfbcntHS1.pST,1);
        stfbcntHS1.pT{n,a}(isnan(stfbcntHS1.pST{n,a})) = 1;
        stfbcntHS1.pST{n,a}(isnan(stfbcntHS1.pST{n,a})) = 0;
        stfbcntHS1.pT{n,a}(isempty(stfbcntHS1.pST{n,a})) = 1;
        stfbcntHS1.pST{n,a}(isempty(stfbcntHS1.pST{n,a})) = 0;
        vauxcnt1 = vertcat(vauxcnt1, stfbcntHS1.pST{n,a}'*stfbcntHS1.pT{n,a}/sum(stfbcntHS1.pT{n,a}));
        vTcnt1 = vertcat(vTcnt1, sum(stfbcntHS1.pT{n,a}));
    end
    vTcnt1(isnan(vauxcnt1)) = [];
    vauxcnt1(isnan(vauxcnt1)) = [];
    vauxcnt2 = [];
    vTcnt2 = [];
    a = inds(j);
    for n = 1 : size(stfbcntHS2.pST,1);
        stfbcntHS2.pT{n,a}(isnan(stfbcntHS2.pST{n,a})) = 1;
        stfbcntHS2.pST{n,a}(isnan(stfbcntHS2.pST{n,a})) = 0;
        stfbcntHS2.pT{n,a}(isempty(stfbcntHS2.pST{n,a})) = 1;
        stfbcntHS2.pST{n,a}(isempty(stfbcntHS2.pST{n,a})) = 0;
        vauxcnt2 = vertcat(vauxcnt2, stfbcntHS2.pST{n,a}'*stfbcntHS2.pT{n,a}/sum(stfbcntHS2.pT{n,a}));
        vTcnt2 = vertcat(vTcnt2, sum(stfbcntHS2.pT{n,a}));
    end
    vTcnt2(isnan(vauxcnt2)) = [];
    vauxcnt2(isnan(vauxcnt2)) = [];
    vauxkir = [];
    vTkir = [];
    a = inds(j);
    for n = 1 : size(stfbkirHS.pST,1);
        stfbkirHS.pT{n,a}(isnan(stfbkirHS.pST{n,a})) = 1;
        stfbkirHS.pST{n,a}(isnan(stfbkirHS.pST{n,a})) = 0;
        stfbkirHS.pT{n,a}(isempty(stfbkirHS.pST{n,a})) = 1;
        stfbkirHS.pST{n,a}(isempty(stfbkirHS.pST{n,a})) = 0;
        vauxkir = vertcat(vauxkir, stfbkirHS.pST{n,a}'*stfbkirHS.pT{n,a}/sum(stfbkirHS.pT{n,a}));
        vTkir = vertcat(vTkir, sum(stfbkirHS.pT{n,a}));
    end
    vTkir(isnan(vauxkir)) = [];
    vauxkir(isnan(vauxkir)) = [];
    
    
    if mod(j,2) == 1
        errorbar(j, vauxcnt1'*vTcnt1/sum(vTcnt1), ...
            sqrt(sum((vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)).*(vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)))...
            /size(stfbcntHS1.pT,1))/sqrt(size(stfbcntHS1.pT,1)), 'om')
        errorbar(j, vauxcnt2'*vTcnt2/sum(vTcnt2), ...
            sqrt(sum((vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)).*(vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)))...
            /size(stfbcntHS2.pT,1))/sqrt(size(stfbcntHS2.pT,1)), 'om')
        
        errorbar(j, vauxkir'*vTkir/sum(vTkir), ...
            sqrt(sum((vauxkir-vauxkir'*vTkir/sum(vTkir)).*(vauxkir-vauxkir'*vTkir/sum(vTkir)))...
            /size(stfbkirHS.pT,1))/sqrt(size(stfbkirHS.pT,1)), 'or')
    else
        errorbar(j-1, vauxcnt1'*vTcnt1/sum(vTcnt1),...
            sqrt(sum((vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)).*(vauxcnt1-vauxcnt1'*vTcnt1/sum(vTcnt1)))...
            /size(stfbcntHS1.pT,1))/sqrt(size(stfbcntHS1.pT,1)), 'oc')
        errorbar(j-1, vauxcnt2'*vTcnt2/sum(vTcnt2),...
            sqrt(sum((vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)).*(vauxcnt2-vauxcnt2'*vTcnt2/sum(vTcnt2)))...
            /size(stfbcntHS2.pT,1))/sqrt(size(stfbcntHS2.pT,1)), 'oc')
        errorbar(j-1, vauxkir'*vTkir/sum(vTkir),...
            sqrt(sum((vauxkir-vauxkir'*vTkir/sum(vTkir)).*(vauxkir-vauxkir'*vTkir/sum(vTkir)))...
            /size(stfbkirHS.pT,1))/sqrt(size(stfbkirHS.pT,1)), 'ob')
    end
end
axis([0.5 8.5 0.3 0.7])

%% Plot Bias ST HSVS
thst = 100;
clc
figure,
hold on
inds = [3 4 5 6 7 8 1 2];
for j = 1 : length(pTypesHS)
    vauxcnt11 = [];
    vTcnt11 = [];
    a = inds(j);
    for n = 1 : size(stfbcntHS1.biasST1,1);
        stfbcntHS1.NbiasST1{n,a}(isnan(stfbcntHS1.biasST1{n,a})) = 1;
        stfbcntHS1.biasST1{n,a}(isnan(stfbcntHS1.biasST1{n,a})) = 0;
        stfbcntHS1.NbiasST1{n,a}(isempty(stfbcntHS1.biasST1{n,a})) = 1;
        stfbcntHS1.biasST1{n,a}(isempty(stfbcntHS1.biasST1{n,a})) = 0;
        vauxcnt11 = vertcat(vauxcnt11, stfbcntHS1.biasST1{n,a}'*stfbcntHS1.NbiasST1{n,a}/sum(stfbcntHS1.NbiasST1{n,a}));
        vTcnt11 = vertcat(vTcnt11, sum(stfbcntHS1.NbiasST1{n,a}));
    end
    vTcnt11(isnan(vauxcnt11)) = [];
    vauxcnt11(isnan(vauxcnt11)) = [];
    vauxcnt120 = [];
    vTcnt120 = [];
    a = inds(j);
    for n = 1 : size(stfbcntHS1.biasST20,1);
        stfbcntHS1.NbiasST20{n,a}(isnan(stfbcntHS1.biasST20{n,a})) = 1;
        stfbcntHS1.biasST20{n,a}(isnan(stfbcntHS1.biasST20{n,a})) = 0;
        stfbcntHS1.NbiasST20{n,a}(isempty(stfbcntHS1.biasST20{n,a})) = 1;
        stfbcntHS1.biasST20{n,a}(isempty(stfbcntHS1.biasST20{n,a})) = 0;
        vauxcnt120 = vertcat(vauxcnt120, stfbcntHS1.biasST20{n,a}'*stfbcntHS1.NbiasST20{n,a}/sum(stfbcntHS1.NbiasST20{n,a}));
        vTcnt120 = vertcat(vTcnt120, sum(stfbcntHS1.NbiasST20{n,a}));
    end
    vTcnt120(isnan(vauxcnt120)) = [];
    vauxcnt120(isnan(vauxcnt120)) = [];
    vauxcnt120(vTcnt120<thst) = [];
    vTcnt120(vTcnt120<thst) = [];
    
    vauxcnt1100 = [];
    vTcnt1100 = [];
    a = inds(j);
    for n = 1 : size(stfbcntHS1.biasST100,1);
        stfbcntHS1.NbiasST100{n,a}(isnan(stfbcntHS1.biasST100{n,a})) = 1;
        stfbcntHS1.biasST100{n,a}(isnan(stfbcntHS1.biasST100{n,a})) = 0;
        stfbcntHS1.NbiasST100{n,a}(isempty(stfbcntHS1.biasST100{n,a})) = 1;
        stfbcntHS1.biasST100{n,a}(isempty(stfbcntHS1.biasST100{n,a})) = 0;
        vauxcnt1100 = vertcat(vauxcnt1100, stfbcntHS1.biasST100{n,a}'*stfbcntHS1.NbiasST100{n,a}/sum(stfbcntHS1.NbiasST100{n,a}));
        vTcnt1100 = vertcat(vTcnt1100, sum(stfbcntHS1.NbiasST100{n,a}));
    end
    vTcnt1100(isnan(vauxcnt1100)) = [];
    vauxcnt1100(isnan(vauxcnt1100)) = [];
    vauxcnt21 = [];
    vTcnt21 = [];
    a = inds(j);
    for n = 1 : size(stfbcntHS2.biasST1,1);
        stfbcntHS2.NbiasST1{n,a}(isnan(stfbcntHS2.biasST1{n,a})) = 1;
        stfbcntHS2.biasST1{n,a}(isnan(stfbcntHS2.biasST1{n,a})) = 0;
        stfbcntHS2.NbiasST1{n,a}(isempty(stfbcntHS2.biasST1{n,a})) = 1;
        stfbcntHS2.biasST1{n,a}(isempty(stfbcntHS2.biasST1{n,a})) = 0;
        vauxcnt21 = vertcat(vauxcnt21, stfbcntHS2.biasST1{n,a}'*stfbcntHS2.NbiasST1{n,a}/sum(stfbcntHS2.NbiasST1{n,a}));
        vTcnt21 = vertcat(vTcnt21, sum(stfbcntHS2.NbiasST1{n,a}));
    end
    vTcnt21(isnan(vauxcnt21)) = [];
    vauxcnt21(isnan(vauxcnt21)) = [];
    vauxcnt220 = [];
    vTcnt220 = [];
    a = inds(j);
    for n = 1 : size(stfbcntHS2.biasST20,1);
        stfbcntHS2.NbiasST20{n,a}(isnan(stfbcntHS2.biasST20{n,a})) = 1;
        stfbcntHS2.biasST20{n,a}(isnan(stfbcntHS2.biasST20{n,a})) = 0;
        stfbcntHS2.NbiasST20{n,a}(isempty(stfbcntHS2.biasST20{n,a})) = 1;
        stfbcntHS2.biasST20{n,a}(isempty(stfbcntHS2.biasST20{n,a})) = 0;
        vauxcnt220 = vertcat(vauxcnt220, stfbcntHS2.biasST20{n,a}'*stfbcntHS2.NbiasST20{n,a}/sum(stfbcntHS2.NbiasST20{n,a}));
        vTcnt220 = vertcat(vTcnt220, sum(stfbcntHS2.NbiasST20{n,a}));
    end
    vTcnt220(isnan(vauxcnt220)) = [];
    vauxcnt220(isnan(vauxcnt220)) = [];
    vauxcnt220(vTcnt220<thst) = [];
    vTcnt220(vTcnt220<thst) = [];
    
    vauxcnt2100 = [];
    vTcnt2100 = [];
    a = inds(j);
    for n = 1 : size(stfbcntHS2.biasST100,1);
        stfbcntHS2.NbiasST100{n,a}(isnan(stfbcntHS2.biasST100{n,a})) = 1;
        stfbcntHS2.biasST100{n,a}(isnan(stfbcntHS2.biasST100{n,a})) = 0;
        stfbcntHS2.NbiasST100{n,a}(isempty(stfbcntHS2.biasST100{n,a})) = 1;
        stfbcntHS2.biasST100{n,a}(isempty(stfbcntHS2.biasST100{n,a})) = 0;
        vauxcnt2100 = vertcat(vauxcnt2100, stfbcntHS2.biasST100{n,a}'*stfbcntHS2.NbiasST100{n,a}/sum(stfbcntHS2.NbiasST100{n,a}));
        vTcnt2100 = vertcat(vTcnt2100, sum(stfbcntHS2.NbiasST100{n,a}));
    end
    vTcnt2100(isnan(vauxcnt2100)) = [];
    vauxcnt2100(isnan(vauxcnt2100)) = [];
    vauxkir1 = [];
    vTkir1 = [];
    a = inds(j);
    for n = 1 : size(stfbkirHS.biasST1,1);
        stfbkirHS.NbiasST1{n,a}(isnan(stfbkirHS.biasST1{n,a})) = 1;
        stfbkirHS.biasST1{n,a}(isnan(stfbkirHS.biasST1{n,a})) = 0;
        stfbkirHS.NbiasST1{n,a}(isempty(stfbkirHS.biasST1{n,a})) = 1;
        stfbkirHS.biasST1{n,a}(isempty(stfbkirHS.biasST1{n,a})) = 0;
        vauxkir1 = vertcat(vauxkir1, stfbkirHS.biasST1{n,a}'*stfbkirHS.NbiasST1{n,a}/sum(stfbkirHS.NbiasST1{n,a}));
        vTkir1 = vertcat(vTkir1, sum(stfbkirHS.NbiasST1{n,a}));
    end
    vTkir1(isnan(vauxkir1)) = [];
    vauxkir1(isnan(vauxkir1)) = [];
    vauxkir20 = [];
    vTkir20 = [];
    a = inds(j);
    for n = 1 : size(stfbkirHS.biasST20,1);
        stfbkirHS.NbiasST20{n,a}(isnan(stfbkirHS.biasST20{n,a})) = 1;
        stfbkirHS.biasST20{n,a}(isnan(stfbkirHS.biasST20{n,a})) = 0;
        stfbkirHS.NbiasST20{n,a}(isempty(stfbkirHS.biasST20{n,a})) = 1;
        stfbkirHS.biasST20{n,a}(isempty(stfbkirHS.biasST20{n,a})) = 0;
        vauxkir20 = vertcat(vauxkir20, stfbkirHS.biasST20{n,a}'*stfbkirHS.NbiasST20{n,a}/sum(stfbkirHS.NbiasST20{n,a}));
        vTkir20 = vertcat(vTkir20, sum(stfbkirHS.NbiasST20{n,a}));
    end
    vTkir20(isnan(vauxkir20)) = [];
    vauxkir20(isnan(vauxkir20)) = [];
    vauxkir20(vTkir20<thst) = [];
    vTkir20(vTkir20<thst) = [];
    
    vauxkir100 = [];
    vTkir100 = [];
    a = inds(j);
    for n = 1 : size(stfbkirHS.biasST100,1);
        stfbkirHS.NbiasST100{n,a}(isnan(stfbkirHS.biasST100{n,a})) = 1;
        stfbkirHS.biasST100{n,a}(isnan(stfbkirHS.biasST100{n,a})) = 0;
        stfbkirHS.NbiasST100{n,a}(isempty(stfbkirHS.biasST100{n,a})) = 1;
        stfbkirHS.biasST100{n,a}(isempty(stfbkirHS.biasST100{n,a})) = 0;
        vauxkir100 = vertcat(vauxkir100, stfbkirHS.biasST100{n,a}'*stfbkirHS.NbiasST100{n,a}/sum(stfbkirHS.NbiasST100{n,a}));
        vTkir100 = vertcat(vTkir100, sum(stfbkirHS.NbiasST100{n,a}));
    end
    vTkir100(isnan(vauxkir100)) = [];
    vauxkir100(isnan(vauxkir100)) = [];
    
    dl = 0.35;
    if mod(j,2) == 1
%         errorbar(j-dl, vauxcnt11'*vTcnt11/sum(vTcnt11), ...
%             sqrt(sum((vauxcnt11-vauxcnt11'*vTcnt11/sum(vTcnt11)).*(vauxcnt11-vauxcnt11'*vTcnt11/sum(vTcnt11)))...
%             /size(stfbcntHS1.biasST1,1))/sqrt(size(stfbcntHS1.biasST1,1)), 'color', [1 0 1], 'marker', 'o')
        errorbar(j, vauxcnt120'*vTcnt120/sum(vTcnt120), ...
            sqrt(sum((vauxcnt120-vauxcnt120'*vTcnt120/sum(vTcnt120)).*(vauxcnt120-vauxcnt120'*vTcnt120/sum(vTcnt120)))...
            /size(stfbcntHS1.biasST20,1))/sqrt(size(stfbcntHS1.biasST20,1)), 'color', [0.75 0 0.75], 'marker', 'o', 'markersize', 7, 'markerfacecolor', [0.75 0 0.75])
%         errorbar(j+dl, vauxcnt1100'*vTcnt1100/sum(vTcnt1100), ...
%             sqrt(sum((vauxcnt1100-vauxcnt1100'*vTcnt1100/sum(vTcnt1100)).*(vauxcnt1100-vauxcnt1100'*vTcnt1100/sum(vTcnt1100)))...
%             /size(stfbcntHS1.biasST100,1))/sqrt(size(stfbcntHS1.biasST100,1)), 'color', [0.5 0 0.5], 'marker', 'o')
                
%         errorbar(j-dl, vauxcnt21'*vTcnt21/sum(vTcnt21), ...
%             sqrt(sum((vauxcnt21-vauxcnt21'*vTcnt21/sum(vTcnt21)).*(vauxcnt21-vauxcnt21'*vTcnt21/sum(vTcnt21)))...
%             /size(stfbcntHS2.biasST1,1))/sqrt(size(stfbcntHS2.biasST1,1)), 'color', [1 0 1], 'marker', 'o')
        errorbar(j, vauxcnt220'*vTcnt220/sum(vTcnt220), ...
            sqrt(sum((vauxcnt220-vauxcnt220'*vTcnt220/sum(vTcnt220)).*(vauxcnt220-vauxcnt220'*vTcnt220/sum(vTcnt220)))...
            /size(stfbcntHS2.biasST20,1))/sqrt(size(stfbcntHS2.biasST20,1)), 'color', [0.75 0 0.75], 'marker', 'o', 'markersize', 7, 'markerfacecolor', [0.75 0 0.75])
%         errorbar(j+dl, vauxcnt2100'*vTcnt2100/sum(vTcnt2100), ...
%             sqrt(sum((vauxcnt2100-vauxcnt2100'*vTcnt2100/sum(vTcnt2100)).*(vauxcnt2100-vauxcnt2100'*vTcnt2100/sum(vTcnt2100)))...
%             /size(stfbcntHS2.biasST100,1))/sqrt(size(stfbcntHS2.biasST100,1)), 'color', [0.5 0 0.5], 'marker', 'o')
                
%         errorbar(j-dl, vauxkir1'*vTkir1/sum(vTkir1), ...
%             sqrt(sum((vauxkir1-vauxkir1'*vTkir1/sum(vTkir1)).*(vauxkir1-vauxkir1'*vTkir1/sum(vTkir1)))...
%             /size(stfbkirHS.biasST1,1))/sqrt(size(stfbkirHS.biasST1,1)), 'color', [1 0 0], 'marker', 'o')
        errorbar(j, vauxkir20'*vTkir20/sum(vTkir20), ...
            sqrt(sum((vauxkir20-vauxkir20'*vTkir20/sum(vTkir20)).*(vauxkir20-vauxkir20'*vTkir20/sum(vTkir20)))...
            /size(stfbkirHS.biasST20,1))/sqrt(size(stfbkirHS.biasST20,1)), 'color', [0.75 0 0], 'marker', 'o', 'markersize', 7, 'markerfacecolor', [0.75 0 0])
%         errorbar(j+dl, vauxkir100'*vTkir100/sum(vTkir100), ...
%             sqrt(sum((vauxkir100-vauxkir100'*vTkir100/sum(vTkir100)).*(vauxkir100-vauxkir100'*vTkir100/sum(vTkir100)))...
%             /size(stfbkirHS.biasST100,1))/sqrt(size(stfbkirHS.biasST100,1)), 'color', [0.5 0 0], 'marker', 'o')

        statTestmww(vauxkir20, vTkir20, vauxcnt120, vTcnt120, thst, pTypesHS{a}, 'HSKIR', 'HSCnt1');
        statTestmww(vauxkir20, vTkir20, vauxcnt220, vTcnt220, thst, pTypesHS{a}, 'HSKIR', 'HSCnt2');
        statTestmww(vauxcnt220, vTcnt220, vauxcnt120, vTcnt120, thst, pTypesHS{a}, 'HSCnt2', 'HSCnt1');
    else
%         errorbar(j-1-dl, vauxcnt11'*vTcnt11/sum(vTcnt11),...
%             sqrt(sum((vauxcnt11-vauxcnt11'*vTcnt11/sum(vTcnt11)).*(vauxcnt11-vauxcnt11'*vTcnt11/sum(vTcnt11)))...
%             /size(stfbcntHS1.biasST1,1))/sqrt(size(stfbcntHS1.biasST1,1)), 'color', [0 1 1], 'marker', 'o')
        errorbar(j-1, vauxcnt120'*vTcnt120/sum(vTcnt120),...
            sqrt(sum((vauxcnt120-vauxcnt120'*vTcnt120/sum(vTcnt120)).*(vauxcnt120-vauxcnt120'*vTcnt120/sum(vTcnt120)))...
            /size(stfbcntHS1.biasST20,1))/sqrt(size(stfbcntHS1.biasST20,1)), 'color', [0 0.75 0.75], 'marker', 'o', 'markersize', 7, 'markerfacecolor', [0 0.75 0.75])      
%         errorbar(j-1+dl, vauxcnt1100'*vTcnt1100/sum(vTcnt1100),...
%             sqrt(sum((vauxcnt1100-vauxcnt1100'*vTcnt1100/sum(vTcnt1100)).*(vauxcnt1100-vauxcnt1100'*vTcnt1100/sum(vTcnt1100)))...
%             /size(stfbcntHS1.biasST100,1))/sqrt(size(stfbcntHS1.biasST100,1)), 'color', [0 0.5 0.5], 'marker', 'o')      
        
%         errorbar(j-1-dl, vauxcnt21'*vTcnt21/sum(vTcnt21),...
%             sqrt(sum((vauxcnt21-vauxcnt21'*vTcnt21/sum(vTcnt21)).*(vauxcnt21-vauxcnt21'*vTcnt21/sum(vTcnt21)))...
%             /size(stfbcntHS2.biasST1,1))/sqrt(size(stfbcntHS2.biasST1,1)), 'color', [0 1 1], 'marker', 'o')
        errorbar(j-1, vauxcnt220'*vTcnt220/sum(vTcnt220),...
            sqrt(sum((vauxcnt220-vauxcnt220'*vTcnt220/sum(vTcnt220)).*(vauxcnt220-vauxcnt220'*vTcnt220/sum(vTcnt220)))...
            /size(stfbcntHS2.biasST20,1))/sqrt(size(stfbcntHS2.biasST20,1)), 'color', [0 0.75 0.75], 'marker', 'o', 'markersize', 7, 'markerfacecolor', [0 0.75 0.75])              
%         errorbar(j-1+dl, vauxcnt2100'*vTcnt2100/sum(vTcnt2100),...
%             sqrt(sum((vauxcnt2100-vauxcnt2100'*vTcnt2100/sum(vTcnt2100)).*(vauxcnt2100-vauxcnt2100'*vTcnt2100/sum(vTcnt2100)))...
%             /size(stfbcntHS2.biasST100,1))/sqrt(size(stfbcntHS2.biasST100,1)), 'color', [0 0.5 0.5], 'marker', 'o')      

%         errorbar(j-1-dl, vauxkir1'*vTkir1/sum(vTkir1),...
%             sqrt(sum((vauxkir1-vauxkir1'*vTkir1/sum(vTkir1)).*(vauxkir1-vauxkir1'*vTkir1/sum(vTkir1)))...
%             /size(stfbkirHS.biasST1,1))/sqrt(size(stfbkirHS.biasST1,1)), 'color', [0 0 1], 'marker', 'o')
        errorbar(j-1, vauxkir20'*vTkir20/sum(vTkir20),...
            sqrt(sum((vauxkir20-vauxkir20'*vTkir20/sum(vTkir20)).*(vauxkir20-vauxkir20'*vTkir20/sum(vTkir20)))...
            /size(stfbkirHS.biasST20,1))/sqrt(size(stfbkirHS.biasST20,1)), 'color', [0 0 0.75], 'marker', 'o', 'markersize', 7, 'markerfacecolor', [0 0 0.75])      
%         errorbar(j-1+dl, vauxkir100'*vTkir100/sum(vTkir100),...
%             sqrt(sum((vauxkir100-vauxkir100'*vTkir100/sum(vTkir100)).*(vauxkir100-vauxkir100'*vTkir100/sum(vTkir100)))...
%             /size(stfbkirHS.biasST100,1))/sqrt(size(stfbkirHS.biasST100,1)), 'color', [0 0 0.5], 'marker', 'o')
        statTestmww(vauxkir20, vTkir20, vauxcnt120, vTcnt120, thst, pTypesHS{a}, 'HSKIR', 'HSCnt1');
        statTestmww(vauxkir20, vTkir20, vauxcnt220, vTcnt220, thst, pTypesHS{a}, 'HSKIR', 'HSCnt2');
        statTestmww(vauxcnt220, vTcnt220, vauxcnt120, vTcnt120, thst, pTypesHS{a}, 'HSCnt2', 'HSCnt1');
    end
    disp(' ')
end
axis([0.5 8.5 0 1])




















%% Plot ST max HSVS
figure,
hold on
lcents = -50:20:1500;
inds = [3 4 5 6 7 8 1 2];
% cmap1 = autumn(15);
% cmap2 = winter(15);
% cmap3 = spring(15);
% cmap4 = summer(15);
for j =  [7 8]%1 : length(pTypesHS)
    vaux1 = [];
    vaux2 = [];
    vauxk = [];
    a = inds(j);
    for n = 1 : size(stfbcntHS1.lST,1)
        vaux1 = vertcat(vaux1, ...
            hist(abs(stfbcntHS1.lST{n,a}), lcents)/sum(hist(stfbcntHS1.lST{n,a}, lcents)));
    end

    for n = 1 : size(stfbcntHS2.lST,1)
        vaux2 = vertcat(vaux2, ...
            hist(abs(stfbcntHS2.lST{n,a}), lcents)/sum(hist(stfbcntHS2.lST{n,a}, lcents)));
    end

    for n = 1 : size(stfbkirHS.lST,1)
        vauxk = vertcat(vauxk, ...
            hist(abs(stfbkirHS.lST{n,a}), lcents)/sum(hist(stfbkirHS.lST{n,a}, lcents)));
    end
    
    if mod(j,2) == 1
        plot(lcents, smooth(nanmean(vaux1)), 'color', 'm', 'linewidth', 2)
        plot(lcents, smooth(nanmean(vaux2)), 'color', 'm', 'linewidth', 2)
        plot(lcents, smooth(nanmean(vauxk)), 'color', 'r', 'linewidth', 2)
    else
        plot(lcents, smooth(nanmean(vaux1)), 'color', 'c', 'linewidth', 2)
        plot(lcents, smooth(nanmean(vaux2)), 'color', 'c', 'linewidth', 2)
        plot(lcents, smooth(nanmean(vauxk)), 'color', 'b', 'linewidth', 2)
    end
end
axis([0 1500 0 0.05])








%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
% path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Unilateral Stimulus\WT TB\';
% [stfbUni, pTypesUni] = GetSTandFBUnil(path);

path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Unilateral Stimulus\R39VT05\LexAPOKir\';
[stfbUniCnt1, pTypesUni] = GetSTandFBUnil(path);

path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Unilateral Stimulus\R39VT05\R39VT05Cnt\';
[stfbUniCnt2, pTypesUni] = GetSTandFBUnil(path);

path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Unilateral Stimulus\R39VT05\R39VT05Kir\';
[stfbUniKir, pTypesUni] = GetSTandFBUnil(path);
%%
figure,
hold on
inds = [5 6 7 8 1 2 3 4];
thr = 50;
for j = 1 : length(pTypesUni)-2
    vauxcnt11 = [];
    vauxcnt12 = [];
    vTcnt11 = [];
    vTcnt12 = [];
    a = inds(j);
    for n = 1 : size(stfbUniCnt1.valFB,1);
        for tt = 1 : 2
            stfbUniCnt1.NvalFB{n,a,tt}(isnan(stfbUniCnt1.valFB{n,a,tt})) = 1;
            stfbUniCnt1.valFB{n,a,tt}(isnan(stfbUniCnt1.valFB{n,a,tt})) = 0;
            stfbUniCnt1.NvalFB{n,a,tt}(isempty(stfbUniCnt1.valFB{n,a,tt})) = 1;
            stfbUniCnt1.valFB{n,a,tt}(isempty(stfbUniCnt1.valFB{n,a,tt})) = 0;
        end
        fbv1 = vertcat(abs(stfbUniCnt1.valFB{n,a,1}(stfbUniCnt1.valFB{n,a,1}>thr)),...
            abs(stfbUniCnt1.valFB{n,a,2}(stfbUniCnt1.valFB{n,a,2}<-thr)));
        nfv1 = vertcat(stfbUniCnt1.NvalFB{n,a,1}(stfbUniCnt1.valFB{n,a,1}>thr),...
            stfbUniCnt1.NvalFB{n,a,2}(stfbUniCnt1.valFB{n,a,2}<-thr));
        vauxcnt11 = vertcat(vauxcnt11, max(fbv1'*nfv1/sum(nfv1)));
        vTcnt11 = vertcat(vTcnt11, sum(nfv1));
        fbv2 = vertcat(abs(stfbUniCnt1.valFB{n,a,1}(stfbUniCnt1.valFB{n,a,1}<-thr)),...
            abs(stfbUniCnt1.valFB{n,a,2}(stfbUniCnt1.valFB{n,a,2}>thr)));
        nfv2 = vertcat(stfbUniCnt1.NvalFB{n,a,1}(stfbUniCnt1.valFB{n,a,1}<-thr),...
            stfbUniCnt1.NvalFB{n,a,2}(stfbUniCnt1.valFB{n,a,2}>thr));
        vauxcnt12 = vertcat(vauxcnt12, max(fbv2'*nfv2/sum(nfv2)));
        vTcnt12 = vertcat(vTcnt12, sum(nfv2));
    end
    vTcnt11(isnan(vauxcnt11)) = [];
    vauxcnt11(isnan(vauxcnt11)) = [];
    if length(vauxcnt11) ~= length(vTcnt11)
        vauxcnt11(vauxcnt11==0) = [];
    end
    vTcnt12(isnan(vauxcnt12)) = [];
    vauxcnt12(isnan(vauxcnt12)) = [];
    if length(vauxcnt12) ~= length(vTcnt12)
        vauxcnt12(vauxcnt12==0) = [];
    end
    
    
    vauxcnt21 = [];
    vauxcnt22 = [];
    vTcnt21 = [];
    vTcnt22 = [];
    a = inds(j);
    for n = 1 : size(stfbUniCnt2.valFB,1);
        for tt = 1 : 2
            stfbUniCnt2.NvalFB{n,a,tt}(isnan(stfbUniCnt2.valFB{n,a,tt})) = 1;
            stfbUniCnt2.valFB{n,a,tt}(isnan(stfbUniCnt2.valFB{n,a,tt})) = 0;
            stfbUniCnt2.NvalFB{n,a,tt}(isempty(stfbUniCnt2.valFB{n,a,tt})) = 1;
            stfbUniCnt2.valFB{n,a,tt}(isempty(stfbUniCnt2.valFB{n,a,tt})) = 0;
        end
        fbv1 = vertcat(abs(stfbUniCnt2.valFB{n,a,1}(stfbUniCnt2.valFB{n,a,1}>thr)),...
            abs(stfbUniCnt2.valFB{n,a,2}(stfbUniCnt2.valFB{n,a,2}<-thr)));
        nfv1 = vertcat(stfbUniCnt2.NvalFB{n,a,1}(stfbUniCnt2.valFB{n,a,1}>thr),...
            stfbUniCnt2.NvalFB{n,a,2}(stfbUniCnt2.valFB{n,a,2}<-thr));
        vauxcnt21 = vertcat(vauxcnt21, max(fbv1'*nfv1/sum(nfv1)));
        vTcnt21 = vertcat(vTcnt21, sum(nfv1));
        fbv2 = vertcat(abs(stfbUniCnt2.valFB{n,a,1}(stfbUniCnt2.valFB{n,a,1}<-thr)),...
            abs(stfbUniCnt2.valFB{n,a,2}(stfbUniCnt2.valFB{n,a,2}>thr)));
        nfv2 = vertcat(stfbUniCnt2.NvalFB{n,a,1}(stfbUniCnt2.valFB{n,a,1}<-thr),...
            stfbUniCnt2.NvalFB{n,a,2}(stfbUniCnt2.valFB{n,a,2}>thr));
        vauxcnt22 = vertcat(vauxcnt22, max(fbv2'*nfv2/sum(nfv2)));
        vTcnt22 = vertcat(vTcnt22, sum(nfv2));
    end
    vTcnt21(isnan(vauxcnt21)) = [];
    vauxcnt21(isnan(vauxcnt21)) = [];
    if length(vauxcnt21) ~= length(vTcnt21)
        vauxcnt21(vauxcnt21==0) = [];
    end
    vTcnt22(isnan(vauxcnt22)) = [];
    vauxcnt22(isnan(vauxcnt22)) = [];
    if length(vauxcnt22) ~= length(vTcnt22)
        vauxcnt22(vauxcnt22==0) = [];
    end
    
    vauxkir1 = [];
    vauxkir2 = [];
    vTkir1 = [];
    vTkir2 = [];
    a = inds(j);
    for n = 1 : size(stfbUniKir.valFB,1);
        for tt = 1 : 2
            stfbUniKir.NvalFB{n,a,tt}(isnan(stfbUniKir.valFB{n,a,tt})) = 1;
            stfbUniKir.valFB{n,a,tt}(isnan(stfbUniKir.valFB{n,a,tt})) = 0;
            stfbUniKir.NvalFB{n,a,tt}(isempty(stfbUniKir.valFB{n,a,tt})) = 1;
            stfbUniKir.valFB{n,a,tt}(isempty(stfbUniKir.valFB{n,a,tt})) = 0;
        end
        fbv1 = vertcat(abs(stfbUniKir.valFB{n,a,1}(stfbUniKir.valFB{n,a,1}>thr)),...
            abs(stfbUniKir.valFB{n,a,2}(stfbUniKir.valFB{n,a,2}<-thr)));
        nfv1 = vertcat(stfbUniKir.NvalFB{n,a,1}(stfbUniKir.valFB{n,a,1}>thr),...
            stfbUniKir.NvalFB{n,a,2}(stfbUniKir.valFB{n,a,2}<-thr));
        vauxkir1 = vertcat(vauxkir1, max(fbv1'*nfv1/sum(nfv1)));
        vTkir1 = vertcat(vTkir1, sum(nfv1));
        fbv2 = vertcat(abs(stfbUniKir.valFB{n,a,1}(stfbUniKir.valFB{n,a,1}<-thr)),...
            abs(stfbUniKir.valFB{n,a,2}(stfbUniKir.valFB{n,a,2}>thr)));
        nfv2 = vertcat(stfbUniKir.NvalFB{n,a,1}(stfbUniKir.valFB{n,a,1}<-thr),...
            stfbUniKir.NvalFB{n,a,2}(stfbUniKir.valFB{n,a,2}>thr));
        vauxkir2 = vertcat(vauxkir2, max(fbv2'*nfv2/sum(nfv2)));
        vTkir2 = vertcat(vTkir2, sum(nfv2));
    end
    vTkir1(isnan(vauxkir1)) = [];
    vauxkir1(isnan(vauxkir1)) = [];
    if length(vauxkir1) ~= length(vTkir1)
        vauxkir1(vauxkir1==0) = [];
    end
    vTkir2(isnan(vauxkir2)) = [];
    vauxkir2(isnan(vauxkir2)) = [];
    if length(vauxkir2) ~= length(vTkir2)
        vauxkir2(vauxkir2==0) = [];
    end
  
    if mod(j,2) == 1
        errorbar(j, vauxcnt11'*vTcnt11/sum(vTcnt11), ...
            sqrt(sum((vauxcnt11-vauxcnt11'*vTcnt11/sum(vTcnt11)).*(vauxcnt11-vauxcnt11'*vTcnt11/sum(vTcnt11)))...
            /size(stfbUniCnt1.valFB,1))/sqrt(size(stfbUniCnt1.valFB,1)), 'color', [1 0 0], 'marker', 'o');
        errorbar(j, vauxcnt12'*vTcnt12/sum(vTcnt12), ...
            sqrt(sum((vauxcnt12-vauxcnt12'*vTcnt12/sum(vTcnt12)).*(vauxcnt12-vauxcnt12'*vTcnt12/sum(vTcnt12)))...
            /size(stfbUniCnt1.valFB,1))/sqrt(size(stfbUniCnt1.valFB,1)), 'color', [1 0 1], 'marker', 'o');
        
        errorbar(j, vauxcnt21'*vTcnt21/sum(vTcnt21), ...
            sqrt(sum((vauxcnt21-vauxcnt21'*vTcnt21/sum(vTcnt21)).*(vauxcnt21-vauxcnt21'*vTcnt21/sum(vTcnt21)))...
            /size(stfbUniCnt2.valFB,1))/sqrt(size(stfbUniCnt2.valFB,1)), 'color', [1 0 0], 'marker', 'o');
        errorbar(j, vauxcnt22'*vTcnt22/sum(vTcnt22), ...
            sqrt(sum((vauxcnt22-vauxcnt22'*vTcnt22/sum(vTcnt22)).*(vauxcnt22-vauxcnt22'*vTcnt22/sum(vTcnt22)))...
            /size(stfbUniCnt2.valFB,1))/sqrt(size(stfbUniCnt2.valFB,1)), 'color', [1 0 1], 'marker', 'o');
        
        errorbar(j, vauxkir1'*vTkir1/sum(vTkir1), ...
            sqrt(sum((vauxkir1-vauxkir1'*vTkir1/sum(vTkir1)).*(vauxkir1-vauxkir1'*vTkir1/sum(vTkir1)))...
            /size(stfbUniKir.valFB,1))/sqrt(size(stfbUniKir.valFB,1)), 'color', [0.5 0 0], 'marker', 'o');
        errorbar(j, vauxkir2'*vTkir2/sum(vTkir2), ...
            sqrt(sum((vauxkir2-vauxkir2'*vTkir2/sum(vTkir2)).*(vauxkir2-vauxkir2'*vTkir2/sum(vTkir2)))...
            /size(stfbUniKir.valFB,1))/sqrt(size(stfbUniKir.valFB,1)), 'color', [0.5 0 0.5], 'marker', 'o');
    else
        errorbar(j-1, vauxcnt11'*vTcnt11/sum(vTcnt11), ...
            sqrt(sum((vauxcnt11-vauxcnt11'*vTcnt11/sum(vTcnt11)).*(vauxcnt11-vauxcnt11'*vTcnt11/sum(vTcnt11)))...
            /size(stfbUniCnt1.valFB,1))/sqrt(size(stfbUniCnt1.valFB,1)), 'color', [0 0 1], 'marker', 'o');
        errorbar(j-1, vauxcnt12'*vTcnt12/sum(vTcnt12), ...
            sqrt(sum((vauxcnt12-vauxcnt12'*vTcnt12/sum(vTcnt12)).*(vauxcnt12-vauxcnt12'*vTcnt12/sum(vTcnt12)))...
            /size(stfbUniCnt1.valFB,1))/sqrt(size(stfbUniCnt1.valFB,1)), 'color', [0 1 1], 'marker', 'o');
        
        errorbar(j-1, vauxcnt21'*vTcnt21/sum(vTcnt21), ...
            sqrt(sum((vauxcnt21-vauxcnt21'*vTcnt21/sum(vTcnt21)).*(vauxcnt21-vauxcnt21'*vTcnt21/sum(vTcnt21)))...
            /size(stfbUniCnt2.valFB,1))/sqrt(size(stfbUniCnt2.valFB,1)), 'color', [0 0 1], 'marker', 'o');
        errorbar(j-1, vauxcnt22'*vTcnt22/sum(vTcnt22), ...
            sqrt(sum((vauxcnt22-vauxcnt22'*vTcnt22/sum(vTcnt22)).*(vauxcnt22-vauxcnt22'*vTcnt22/sum(vTcnt22)))...
            /size(stfbUniCnt2.valFB,1))/sqrt(size(stfbUniCnt2.valFB,1)), 'color', [0 1 1], 'marker', 'o');
        
        errorbar(j-1, vauxkir1'*vTkir1/sum(vTkir1), ...
            sqrt(sum((vauxkir1-vauxkir1'*vTkir1/sum(vTkir1)).*(vauxkir1-vauxkir1'*vTkir1/sum(vTkir1)))...
            /size(stfbUniKir.valFB,1))/sqrt(size(stfbUniKir.valFB,1)), 'color', [0 0 0.5], 'marker', 'o');
        errorbar(j-1, vauxkir2'*vTkir2/sum(vTkir2), ...
            sqrt(sum((vauxkir2-vauxkir2'*vTkir2/sum(vTkir2)).*(vauxkir2-vauxkir2'*vTkir2/sum(vTkir2)))...
            /size(stfbUniKir.valFB,1))/sqrt(size(stfbUniKir.valFB,1)), 'color', [0 0.5 0.5], 'marker', 'o');
    end
end
axis([0.5 5.5 50 180])
% axis([0.5 8.5 0 150])
xlabel('Dot Size(º)')
ylabel('Drift (º/s)')


%%
figure,
hold on
inds = [5 6 7 8 1 2 3 4];
thr = 25;
for j = 1 : length(pTypesUni)-2
    vauxcnt11 = [];
    vauxcnt12 = [];
    vTcnt11 = [];
    vTcnt12 = [];
    a = inds(j);
    for n = 1 : size(stfbUni.valFB,1);
        for tt = 1 : 2
            stfbUni.NvalFB{n,a,tt}(isnan(stfbUni.valFB{n,a,tt})) = 1;
            stfbUni.valFB{n,a,tt}(isnan(stfbUni.valFB{n,a,tt})) = 0;
            stfbUni.NvalFB{n,a,tt}(isempty(stfbUni.valFB{n,a,tt})) = 1;
            stfbUni.valFB{n,a,tt}(isempty(stfbUni.valFB{n,a,tt})) = 0;
        end
        fbv1 = vertcat(abs(stfbUni.valFB{n,a,1}(stfbUni.valFB{n,a,1}>thr)),...
            abs(stfbUni.valFB{n,a,2}(stfbUni.valFB{n,a,2}<-thr)));
        nfv1 = vertcat(stfbUni.NvalFB{n,a,1}(stfbUni.valFB{n,a,1}>thr),...
            stfbUni.NvalFB{n,a,2}(stfbUni.valFB{n,a,2}<-thr));
        vauxcnt11 = vertcat(vauxcnt11, max(fbv1'*nfv1/sum(nfv1)));
        vTcnt11 = vertcat(vTcnt11, sum(nfv1));
        fbv2 = vertcat(abs(stfbUni.valFB{n,a,1}(stfbUni.valFB{n,a,1}<-thr)),...
            abs(stfbUni.valFB{n,a,2}(stfbUni.valFB{n,a,2}>thr)));
        nfv2 = vertcat(stfbUni.NvalFB{n,a,1}(stfbUni.valFB{n,a,1}<-thr),...
            stfbUni.NvalFB{n,a,2}(stfbUni.valFB{n,a,2}>thr));
        vauxcnt12 = vertcat(vauxcnt12, max(fbv2'*nfv2/sum(nfv2)));
        vTcnt12 = vertcat(vTcnt12, sum(nfv2));
    end
    vTcnt11(isnan(vauxcnt11)) = [];
    vauxcnt11(isnan(vauxcnt11)) = [];
    if length(vauxcnt11) ~= length(vTcnt11)
        vauxcnt11(vauxcnt11==0) = [];
    end
    vTcnt12(isnan(vauxcnt12)) = [];
    vauxcnt12(isnan(vauxcnt12)) = [];
    if length(vauxcnt12) ~= length(vTcnt12)
        vauxcnt12(vauxcnt12==0) = [];
    end

  
    if mod(j,2) == 1
        errorbar(j, vauxcnt11'*vTcnt11/sum(vTcnt11), ...
            sqrt(sum((vauxcnt11-vauxcnt11'*vTcnt11/sum(vTcnt11)).*(vauxcnt11-vauxcnt11'*vTcnt11/sum(vTcnt11)))...
            /size(stfbUni.valFB,1))/sqrt(size(stfbUni.valFB,1)), 'color', [1 0 0], 'marker', 'o');
        errorbar(j, vauxcnt12'*vTcnt12/sum(vTcnt12), ...
            sqrt(sum((vauxcnt12-vauxcnt12'*vTcnt12/sum(vTcnt12)).*(vauxcnt12-vauxcnt12'*vTcnt12/sum(vTcnt12)))...
            /size(stfbUni.valFB,1))/sqrt(size(stfbUni.valFB,1)), 'color', [1 0 1], 'marker', 'o');
    else
        errorbar(j-1, vauxcnt11'*vTcnt11/sum(vTcnt11), ...
            sqrt(sum((vauxcnt11-vauxcnt11'*vTcnt11/sum(vTcnt11)).*(vauxcnt11-vauxcnt11'*vTcnt11/sum(vTcnt11)))...
            /size(stfbUni.valFB,1))/sqrt(size(stfbUni.valFB,1)), 'color', [0 0 1], 'marker', 'o');
        errorbar(j-1, vauxcnt12'*vTcnt12/sum(vTcnt12), ...
            sqrt(sum((vauxcnt12-vauxcnt12'*vTcnt12/sum(vTcnt12)).*(vauxcnt12-vauxcnt12'*vTcnt12/sum(vTcnt12)))...
            /size(stfbUni.valFB,1))/sqrt(size(stfbUni.valFB,1)), 'color', [0 1 1], 'marker', 'o');
    end
end
axis([0.5 5.5 25 180])
% axis([0.5 8.5 0 150])
xlabel('Dot Size(º)')
ylabel('Drift (º/s)')

%%
figure,
hold on
lcents = -50:20:1500;
inds = [5 6 7 8 1 2 3 4];
thr = 25;
cmap1 = autumn(15);
cmap2 = winter(15);
for j = 1 : length(pTypesUni)-2
    vauxcnt11 = [];
    vauxcnt12 = [];
    a = inds(j);
    for n = 1 : size(stfbUni.lST,1);
        for tt = 1 : 2
            stfbUni.lST{n,a,tt}(isnan(stfbUni.lST{n,a,tt})) = 0;
        end
        sp1 = vertcat(abs(stfbUni.lST{n,a,1}(stfbUni.lST{n,a,1}>0)), ...
            abs(stfbUni.lST{n,a,2}(stfbUni.lST{n,a,2}<0)));
        vauxcnt11 = vertcat(vauxcnt11, hist(abs(sp1), lcents)/sum(hist(abs(sp1), lcents)));
        sp2 = vertcat(abs(stfbUni.lST{n,a,1}(stfbUni.lST{n,a,1}<0)), ...
            abs(stfbUni.lST{n,a,2}(stfbUni.lST{n,a,2}>0)));
        vauxcnt12 = vertcat(vauxcnt12, hist(abs(sp2), lcents)/sum(hist(abs(sp2), lcents)));
    end
    vauxcnt11(isnan(vauxcnt11)) = [];
    vauxcnt12(isnan(vauxcnt12)) = [];

    
    if mod(j,2) == 1
%         plot(lcents, smooth(mean(vauxcnt11)), 'color', cmap1(j,:), 'linewidth', 2)
%         plot(lcents, smooth(mean(vauxcnt12)), 'color', cmap2(j,:), 'linewidth', 2)
    else
        plot(lcents, smooth(mean(vauxcnt11)), 'color', cmap1(j,:), 'linewidth', 2)
        plot(lcents, smooth(mean(vauxcnt12)), 'color', cmap2(j,:), 'linewidth', 2)
    end

end
% axis([0 1500 0 0.05])
% % axis([0.5 5.5 25 180])
% % axis([0.5 8.5 0 150])
% xlabel('Dot Size(º)')
% ylabel('Drift (º/s)')

%%
clear
clc

path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\WT TB\';
[stfb, pTypesV] = GetSTandFB(path);

path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Dark\Dark WTTB\';
[stfd, pTypesD] = GetSTandFB(path);

%%
cmap = hot(10);
figure,
hold on
lcents = 1:1:150;
inds = [5 7 9 1];
vauxv = cell(size(inds));
vauxd = [];
for j = 1 : length(inds)
    for n = 1 : size(stfb.lFB,1);
        vauxv{j} = vertcat(vauxv{j}, ...
            hist(stfb.lFB{n,inds(j)}, lcents)/sum(hist(stfb.lFB{n,inds(j)}, lcents)));
    end
end
for n = 1 : size(stfd.lFB,1);
    vauxd = vertcat(vauxd, ...
        hist(stfd.lFB{n,1}, lcents)/sum(hist(stfd.lFB{n,1}, lcents)));
end
for j = 1 : length(inds)
    plot(lcents/60, (nanmean(vauxv{j})), 'color', cmap(j,:), 'linewidth', 2)
end
plot(lcents/60, (nanmean(vauxd)), 'b', 'linewidth', 2)

% set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel('Duration FW (s)')
ylabel('Cumulative Distribution')
axis([0 2 0 0.15])
%%
p = mwwtest(vauxv{3}*lcents', vauxv{4}*lcents', 0);
disp(num2str(p.p(2)))
