%% Teste to spike detector
clear
clc

MSTRD = [];
MSTRL = [];
SDSTRD = [];
SDSTRL = [];

[params] = GetParams();
vfstdtvals = [0 0.1 0.2 0.5 1 2 3 5 10 20];
for i = 1 : length(vfstdtvals)
    % Test 1: Vf std
    params.vfstdt = vfstdtvals(i);
    path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Dark\Dark WTTB\';
    [STRALL1, DSTALL1, NSTRALL1, PSSALL1, NPSSALL1, pTypes1, ANGDALL1, VTALL1] = GetStrAndVisInf(path, params);
    path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\WT TB\';
    [STRALL2, DSTALL2, NSTRALL2, PSSALL2, NPSSALL2, pTypes2, ANGDALL2, VTALL2] = GetStrAndVisInf(path, params);
    params.pTypes = pTypes1;
    params.lthr = 350;
    params.pTypes = pTypes1;
    [MSTR1, SEMSTR1, MF1, NF1] = GetGMSEM(STRALL1, DSTALL1, NSTRALL1, params);
    params.pTypes = pTypes2;
    [MSTR2, SEMSTR2, MF2, NF2] = GetGMSEM(STRALL2, DSTALL2, NSTRALL2, params);
    
    MSTRD = vertcat(MSTRD, MSTR1(1));
    MSTRL = vertcat(MSTRL, MSTR2(1));
    SDSTRD = vertcat(SDSTRD, SEMSTR1(1));
    SDSTRL = vertcat(SDSTRL, SEMSTR2(1));
end

figure;
errorbar(vfstdtvals, MSTRD, SDSTRD, 'ok')
hold on
errorbar(vfstdtvals, MSTRL, SDSTRL, 'or')
xlabel('min Vf sdt (mm/s)')


%%
clear
clc

MSTRD = [];
MSTRL = [];
SDSTRD = [];
SDSTRL = [];

[params] = GetParams();
cutoffvals = [0 0.01 0.05 0.075 0.1 0.12 0.15 0.175 0.2 0.225 0.25 0.3 0.5];
for i = 1 : length(cutoffvals)
    % Test 1: Vf std
    params.cutoff = cutoffvals(i);
    path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Dark\Dark WTTB\';
    [STRALL1, DSTALL1, NSTRALL1, PSSALL1, NPSSALL1, pTypes1, ANGDALL1, VTALL1] = GetStrAndVisInf(path, params);
    path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\WT TB\';
    [STRALL2, DSTALL2, NSTRALL2, PSSALL2, NPSSALL2, pTypes2, ANGDALL2, VTALL2] = GetStrAndVisInf(path, params);
    params.pTypes = pTypes1;
    params.lthr = 350;
    params.pTypes = pTypes1;
    [MSTR1, SEMSTR1, MF1, NF1] = GetGMSEM(STRALL1, DSTALL1, NSTRALL1, params);
    params.pTypes = pTypes2;
    [MSTR2, SEMSTR2, MF2, NF2] = GetGMSEM(STRALL2, DSTALL2, NSTRALL2, params);
    
    MSTRD = vertcat(MSTRD, MSTR1(1));
    MSTRL = vertcat(MSTRL, MSTR2(1));
    SDSTRD = vertcat(SDSTRD, SEMSTR1(1));
    SDSTRL = vertcat(SDSTRL, SEMSTR2(1));
end
%%
figure;
errorbar(cutoffvals, MSTRD, SDSTRD, 'ok')
hold on
errorbar(cutoffvals, MSTRL, SDSTRL, 'or')
xlabel('cuttuff a.u.')


%%
clear
clc

[params] = GetParams();
% cutoffvals = [0 0.01 0.05 0.075 0.1 0.12 0.15 0.175 0.2 0.225 0.25 0.3 0.5];

minFs = [10^-3 0.01 0.05 0.1 0.5 1 5 10 15 20];
maxFs = [10^-3 0.01 0.05 0.1 0.5 1 5 10 15 20 30];
MSTRD = cell(length(minFs),length(maxFs));
MSTRL = cell(length(minFs),length(maxFs));
SDSTRD = cell(length(minFs),length(maxFs));
SDSTRL = cell(length(minFs),length(maxFs));
for i = 1 : length(minFs)
    for j = 1 : length(maxFs)
        if (minFs(i) < maxFs(j))
            % Test 1: Vf std
            params.maxF = maxFs(j);
            params.minF = minFs(i);
            path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Dark\Dark WTTB\';
            [STRALL1, DSTALL1, NSTRALL1, PSSALL1, NPSSALL1, pTypes1, ANGDALL1, VTALL1] = GetStrAndVisInf(path, params);
            path = 'C:\Users\User\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\WT TB\';
            [STRALL2, DSTALL2, NSTRALL2, PSSALL2, NPSSALL2, pTypes2, ANGDALL2, VTALL2] = GetStrAndVisInf(path, params);
            params.pTypes = pTypes1;
            params.lthr = 350;
            params.pTypes = pTypes1;
            [MSTR1, SEMSTR1, MF1, NF1] = GetGMSEM(STRALL1, DSTALL1, NSTRALL1, params);
            params.pTypes = pTypes2;
            [MSTR2, SEMSTR2, MF2, NF2] = GetGMSEM(STRALL2, DSTALL2, NSTRALL2, params);
            
            MSTRD{i,j} = vertcat(MSTRD{i,j}, MSTR1(1));
            MSTRL{i,j} = vertcat(MSTRL{i,j}, MSTR2(1));
            SDSTRD{i,j} = vertcat(SDSTRD{i,j}, SEMSTR1(1));
            SDSTRL{i,j} = vertcat(SDSTRL{i,j}, SEMSTR2(1));
        end
    end
end
%%
MSTRD2 = MSTRD;
MSTRL2 = MSTRL;
SDSTRD2 = SDSTRD;
SDSTRL2 = SDSTRL;
%%
for i = 1 : size(MSTRD2,1)
    for j = 1 : size(MSTRD2,2)
        if (isempty(MSTRD2{i,j}))
            MSTRD2{i,j} = nan;
            MSTRL2{i,j} = nan;
            SDSTRD2{i,j} = nan;
            SDSTRL2{i,j} = nan;
        end
    end
end




figure,
subplot(2,2,1)
imagesc(minFs, maxFs, cell2mat(MSTRL2)')
set(gca, 'YDir', 'normal');
caxis([20 55])
colormap hot
colorbar
subplot(2,2,2)
imagesc(minFs, maxFs, cell2mat(SDSTRL2)')
colormap hot
caxis([0 3])
set(gca, 'YDir', 'normal');
colorbar
subplot(2,2,3)
imagesc(minFs, maxFs, cell2mat(MSTRD2)')
colormap hot
caxis([20 55])
set(gca, 'YDir', 'normal');
colorbar
subplot(2,2,4)
imagesc(minFs, maxFs, cell2mat(SDSTRD2)')
colormap hot
caxis([0 3])
set(gca, 'YDir', 'normal');
colorbar
