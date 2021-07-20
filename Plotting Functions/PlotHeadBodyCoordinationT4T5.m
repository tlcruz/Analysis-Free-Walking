addpath('Head Analysis')
clc
clear
% load data and calculate head body covariance
params = GetParams();
pathC = 'D:\Dropbox\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\LegTracking\Dark\';
% Compute cross covarianves
[CCVAllC, NCCVAllC, CCVFWC, NCCVFWC, pTypesC, DVRFWC, DHANGC, rdC] = GetCCov(pathC, params);
pathS = 'D:\Dropbox\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\LegTracking\Dark\';
% Compute cross covarianves
[CCVAllS, NCCVAllS, CCVFWS, NCCVFWS, pTypesS, DVRFWS, DHANGS, rdS] = GetCCov(pathS, params);

figure('Name',pathC,'NumberTitle','off')
inds = [1 2 3 4 5 6 7 8];
for i = 1 : length(pTypesC)
    subplot(2,ceil(length(pTypesC)/2), i)
    hold on
    % stack the covariances for all the control flies
    GMCCVC = [];
    GMLC = [];
    for k = 1 : size(CCVAllC,2)
        if isempty(find(isnan(CCVAllC{inds(i),k})==1))
            GMCCVC = horzcat(GMCCVC, CCVAllC{inds(i),k});
            GMLC = vertcat(GMLC, NCCVAllC{inds(i),k});
        end
    end
    GMCCVFC = [];
    GMLFC = [];
    for k = 1 : size(CCVFWC,2)
        if isempty(find(isnan(CCVFWC{inds(i),k})==1))
            GMCCVFC = horzcat(GMCCVFC, CCVFWC{inds(i),k});
            GMLFC = vertcat(GMLFC, NCCVFWC{inds(i),k});
        end
    end
    
    % stack the covariances for all the silenced flies
    GMCCVS = [];
    GMLS = [];
    for k = 1 : size(CCVAllS,2)
        if isempty(find(isnan(CCVAllS{inds(i),k})==1))
            GMCCVS = horzcat(GMCCVS, CCVAllS{inds(i),k});
            GMLS = vertcat(GMLS, NCCVAllS{inds(i),k});
        end
    end
    GMCCVFS = [];
    GMLFS = [];
    for k = 1 : size(CCVFWS,2)
        if isempty(find(isnan(CCVFWS{inds(i),k})==1))
            GMCCVFS = horzcat(GMCCVFS, CCVFWS{inds(i),k});
            GMLFS = vertcat(GMLFS, NCCVFWS{inds(i),k});
        end
    end
    
    % stack the bootstrapped covariances for all the control flies
    mRandAllC = [];
    stdRandAllC = [];
    nRandAllC = [];
    for k = 1 : size(CCVAllC,2)
        if isempty(find(isnan(CCVAllC{inds(i),k})==1))
            if ~isempty(rdC.CCVALLMeanDist{inds(i),k}*rdC.CCVALLNDist{inds(i),k}/sum(rdC.CCVALLNDist{inds(i),k})) & ...
                ~isnan(rdC.CCVALLMeanDist{inds(i),k}*rdC.CCVALLNDist{inds(i),k}/sum(rdC.CCVALLNDist{inds(i),k})) & ...
                max(rdC.CCVALLMeanDist{inds(i),k}*rdC.CCVALLNDist{inds(i),k}/sum(rdC.CCVALLNDist{inds(i),k}))<0.32
                mRandAllC = horzcat(mRandAllC, ...
                    rdC.CCVALLMeanDist{inds(i),k}*rdC.CCVALLNDist{inds(i),k}/sum(rdC.CCVALLNDist{inds(i),k}));
                nRandAllC = vertcat(nRandAllC, sum(rdC.CCVALLNDist{inds(i),k}));
            end
        end
    end
    
    mRandFwC = [];
    stdRandFwC = [];
    nRandFwC = [];
    for k = 1 : size(CCVFWC,2)
        if isempty(find(isnan(CCVFWC{inds(i),k})==1))
            if ~isempty(rdC.CCVFWMeanDist{inds(i),k}*rdC.CCVFWNDist{inds(i),k}/sum(rdC.CCVFWNDist{inds(i),k})) & ...
                ~isnan(rdC.CCVFWMeanDist{inds(i),k}*rdC.CCVFWNDist{inds(i),k}/sum(rdC.CCVFWNDist{inds(i),k}))
                mRandFwC = horzcat(mRandFwC, ...
                    rdC.CCVFWMeanDist{inds(i),k}*rdC.CCVFWNDist{inds(i),k}/sum(rdC.CCVFWNDist{inds(i),k}));
                nRandFwC = vertcat(nRandFwC, sum(rdC.CCVFWNDist{inds(i),k}));
            end
        end
    end
    
    
    % stack the bootstrapped covariances for all the silenced flies
    mRandAllS = [];
    stdRandAllS = [];
    nRandAllS = [];
    for k = 1 : size(CCVAllS,2)
        if isempty(find(isnan(CCVAllS{inds(i),k})==1))
            if ~isempty(rdS.CCVALLMeanDist{inds(i),k}*rdS.CCVALLNDist{inds(i),k}/sum(rdS.CCVALLNDist{inds(i),k})) & ...
                ~isnan(rdS.CCVALLMeanDist{inds(i),k}*rdS.CCVALLNDist{inds(i),k}/sum(rdS.CCVALLNDist{inds(i),k})) & ...
                max(rdS.CCVALLMeanDist{inds(i),k}*rdS.CCVALLNDist{inds(i),k}/sum(rdS.CCVALLNDist{inds(i),k}))<0.32
                mRandAllS = horzcat(mRandAllS, ...
                    rdS.CCVALLMeanDist{inds(i),k}*rdS.CCVALLNDist{inds(i),k}/sum(rdS.CCVALLNDist{inds(i),k}));
                nRandAllS = vertcat(nRandAllS, sum(rdS.CCVALLNDist{inds(i),k}));
            end
        end
    end
    
    mRandFwS = [];
    stdRandFwS = [];
    nRandFwS = [];
    for k = 1 : size(CCVFWS,2)
        if isempty(find(isnan(CCVFWS{inds(i),k})==1))
            if ~isempty(rdS.CCVFWMeanDist{inds(i),k}*rdS.CCVFWNDist{inds(i),k}/sum(rdS.CCVFWNDist{inds(i),k})) & ...
                ~isnan(rdS.CCVFWMeanDist{inds(i),k}*rdS.CCVFWNDist{inds(i),k}/sum(rdS.CCVFWNDist{inds(i),k}))
                mRandFwS = horzcat(mRandFwS, ...
                    rdS.CCVFWMeanDist{inds(i),k}*rdS.CCVFWNDist{inds(i),k}/sum(rdS.CCVFWNDist{inds(i),k}));
                nRandFwS = vertcat(nRandFwS, sum(rdS.CCVFWNDist{inds(i),k}));
            end
        end
    end
    
    
    % compute grand mean and error for controlflies
    gmEDC = GMCCVC*GMLC/sum(GMLC);
    gsdEDC = sqrt(((bsxfun(@minus, GMCCVC, gmEDC).*...
                   bsxfun(@minus, GMCCVC, gmEDC)*GMLC))/sum(GMLC))/sqrt(size(CCVAllC,2));
    
    gmEDRandMC = mRandAllC*nRandAllC/sum(nRandAllC);
    gmEDRandSTDC = sqrt(((bsxfun(@minus, mRandAllC, gmEDRandMC).*...
                         bsxfun(@minus, mRandAllC, gmEDRandMC)*nRandAllC))/sum(nRandAllC));

    gmEDFC = GMCCVFC*GMLFC/sum(GMLFC);
    gsdEDFC = sqrt(((bsxfun(@minus, GMCCVFC, gmEDFC).*...
                    bsxfun(@minus, GMCCVFC, gmEDFC)*GMLFC))/sum(GMLFC))/sqrt(size(CCVFWC,2));

    gmEDRandMFC = mRandFwC*nRandFwC/sum(nRandFwC);
    gmEDRandSTDFC = sqrt(((bsxfun(@minus, mRandFwC, gmEDRandMFC).*...
                          bsxfun(@minus, mRandFwC, gmEDRandMFC)*nRandFwC))/sum(nRandFwC));
           
    % compute grand mean and error for silenced flies
    gmEDS = GMCCVS*GMLS/sum(GMLS);
    gsdEDS = sqrt(((bsxfun(@minus, GMCCVS, gmEDS).*...
                   bsxfun(@minus, GMCCVS, gmEDS)*GMLS))/sum(GMLS))/sqrt(size(CCVAllS,2));
    
    gmEDRandMS = mRandAllS*nRandAllS/sum(nRandAllS);
    gmEDRandSTDS = sqrt(((bsxfun(@minus, mRandAllS, gmEDRandMS).*...
                         bsxfun(@minus, mRandAllS, gmEDRandMS)*nRandAllS))/sum(nRandAllS));

    gmEDFS = GMCCVFS*GMLFS/sum(GMLFS);
    gsdEDFS = sqrt(((bsxfun(@minus, GMCCVFS, gmEDFS).*...
                    bsxfun(@minus, GMCCVFS, gmEDFS)*GMLFS))/sum(GMLFS))/sqrt(size(CCVFWS,2));

    gmEDRandMFS = mRandFwS*nRandFwS/sum(nRandFwS);
    gmEDRandSTDFS = sqrt(((bsxfun(@minus, mRandFwS, gmEDRandMFS).*...
                          bsxfun(@minus, mRandFwS, gmEDRandMFS)*nRandFwS))/sum(nRandFwS));

    % plots control flies           
    plot((-40:40)/60,gmEDC,'r', 'linewidth', 3)
    plot((-40:40)/60,gmEDC+gsdEDC,'r', 'linewidth', 1)
    plot((-40:40)/60,gmEDC-gsdEDC,'r', 'linewidth', 1)
    
    plot((-40:40)/60,gmEDRandMC+gmEDRandSTDC,'color', [0.4 0 0], 'linewidth', 1)
    plot((-40:40)/60,gmEDRandMC-gmEDRandSTDC,'color', [0.4 0 0], 'linewidth', 1)
    
    plot((-20:20)/60,gmEDFC,'m', 'linewidth', 3)
    plot((-20:20)/60,gmEDFC+gsdEDFC,'m', 'linewidth', 1)
    plot((-20:20)/60,gmEDFC-gsdEDFC,'m', 'linewidth', 1)

    plot((-20:20)/60,gmEDRandMFC+gmEDRandSTDFC,'color', [0.4 0 0.4], 'linewidth', 1)
    plot((-20:20)/60,gmEDRandMFC-gmEDRandSTDFC,'color', [0.4 0 0.4], 'linewidth', 1)
    
    % plots silenced flies
    plot((-40:40)/60,gmEDS,'b', 'linewidth', 3)
    plot((-40:40)/60,gmEDS+gsdEDS,'b', 'linewidth', 1)
    plot((-40:40)/60,gmEDS-gsdEDS,'b', 'linewidth', 1)
    
    plot((-40:40)/60,gmEDRandMS+gmEDRandSTDS,'color', [0 0 0.4], 'linewidth', 1)
    plot((-40:40)/60,gmEDRandMS-gmEDRandSTDS,'color', [0 0 0.4], 'linewidth', 1)
    
    plot((-20:20)/60,gmEDFS,'c', 'linewidth', 3)
    plot((-20:20)/60,gmEDFS+gsdEDFS,'c', 'linewidth', 1)
    plot((-20:20)/60,gmEDFS-gsdEDFS,'c', 'linewidth', 1)

    plot((-20:20)/60,gmEDRandMFS+gmEDRandSTDFS,'color', [0 0.4 0.4], 'linewidth', 1)
    plot((-20:20)/60,gmEDRandMFS-gmEDRandSTDFS,'color', [0 0.4 0.4], 'linewidth', 1)
    
    axis([-1 1 -1 1])
    title([pTypesC{inds(i)} '   ' num2str(max(GMCCVC*GMLC/sum(GMLC)))])
end

