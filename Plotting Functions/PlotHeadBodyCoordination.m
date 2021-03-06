addpath('Head Analysis')

clc
clear
params = GetParams();
path = '\';
% get head body covariance durinf forward runs and saccades
[CCVAll, NCCVAll, CCVFW, NCCVFW, pTypes, DVRFW, DHANG, rd] = GetCCov(path, params);

figure('Name',path,'NumberTitle','off')
inds = [5 6 7 8 9 10 1 2 3 4];
for i = 1 : length(pTypes)
    subplot(2,ceil(length(pTypes)/2), i)
    hold on
    % stack the covariances for all the flies
    GMCCV = [];
    GML = [];
    for k = 1 : size(CCVAll,2)
        if isempty(find(isnan(CCVAll{inds(i),k})==1))
            GMCCV = horzcat(GMCCV, CCVAll{inds(i),k});
            GML = vertcat(GML, NCCVAll{inds(i),k});
        end
    end
    GMCCVF = [];
    GMLF = [];
    for k = 1 : size(CCVFW,2)
        if isempty(find(isnan(CCVFW{inds(i),k})==1))
            GMCCVF = horzcat(GMCCVF, CCVFW{inds(i),k});
            GMLF = vertcat(GMLF, NCCVFW{inds(i),k});
        end
    end
    
    % stack the bootstrapped covariances for all the flies
    mRandAll = [];
    stdRandAll = [];
    nRandAll = [];
    for k = 1 : size(CCVAll,2)
        if isempty(find(isnan(CCVAll{inds(i),k})==1))
            if ~isempty(rd.CCVALLMeanDist{inds(i),k}*rd.CCVALLNDist{inds(i),k}/sum(rd.CCVALLNDist{inds(i),k})) & ...
                ~isnan(rd.CCVALLMeanDist{inds(i),k}*rd.CCVALLNDist{inds(i),k}/sum(rd.CCVALLNDist{inds(i),k})) & ...
                max(rd.CCVALLMeanDist{inds(i),k}*rd.CCVALLNDist{inds(i),k}/sum(rd.CCVALLNDist{inds(i),k}))<0.32
                mRandAll = horzcat(mRandAll, ...
                    rd.CCVALLMeanDist{inds(i),k}*rd.CCVALLNDist{inds(i),k}/sum(rd.CCVALLNDist{inds(i),k}));
                nRandAll = vertcat(nRandAll, sum(rd.CCVALLNDist{inds(i),k}));
            end
        end
    end
    
    mRandFw = [];
    stdRandFw = [];
    nRandFw = [];
    for k = 1 : size(CCVFW,2)
        if isempty(find(isnan(CCVFW{inds(i),k})==1))
            if ~isempty(rd.CCVFWMeanDist{inds(i),k}*rd.CCVFWNDist{inds(i),k}/sum(rd.CCVFWNDist{inds(i),k})) & ...
                ~isnan(rd.CCVFWMeanDist{inds(i),k}*rd.CCVFWNDist{inds(i),k}/sum(rd.CCVFWNDist{inds(i),k}))
                mRandFw = horzcat(mRandFw, ...
                    rd.CCVFWMeanDist{inds(i),k}*rd.CCVFWNDist{inds(i),k}/sum(rd.CCVFWNDist{inds(i),k}));
                nRandFw = vertcat(nRandFw, sum(rd.CCVFWNDist{inds(i),k}));
            end
        end
    end
    
    % compute grand mean and error
    gmED = GMCCV*GML/sum(GML);
    gsdED = sqrt(((bsxfun(@minus, GMCCV, gmED).*...
                   bsxfun(@minus, GMCCV, gmED)*GML))/sum(GML))/sqrt(size(CCVAll,2));
    
    gmEDRandM = mRandAll*nRandAll/sum(nRandAll);
    gmEDRandSTD = sqrt(((bsxfun(@minus, mRandAll, gmEDRandM).*...
                         bsxfun(@minus, mRandAll, gmEDRandM)*nRandAll))/sum(nRandAll));

    gmEDF = GMCCVF*GMLF/sum(GMLF);
    gsdEDF = sqrt(((bsxfun(@minus, GMCCVF, gmEDF).*...
                    bsxfun(@minus, GMCCVF, gmEDF)*GMLF))/sum(GMLF))/sqrt(size(CCVFW,2));

    gmEDRandMF = mRandFw*nRandFw/sum(nRandFw);
    gmEDRandSTDF = sqrt(((bsxfun(@minus, mRandFw, gmEDRandMF).*...
                          bsxfun(@minus, mRandFw, gmEDRandMF)*nRandFw))/sum(nRandFw));

                      
    plot((-40:40)/60,gmED,'r', 'linewidth', 3)
    plot((-40:40)/60,gmED+gsdED,'r', 'linewidth', 1)
    plot((-40:40)/60,gmED-gsdED,'r', 'linewidth', 1)
    
    plot((-40:40)/60,gmEDRandM+gmEDRandSTD,'color', [0.4 0 0], 'linewidth', 1)
    plot((-40:40)/60,gmEDRandM-gmEDRandSTD,'color', [0.4 0 0], 'linewidth', 1)
    
    plot((-20:20)/60,gmEDF,'b', 'linewidth', 3)
    plot((-20:20)/60,gmEDF+gsdEDF,'b', 'linewidth', 1)
    plot((-20:20)/60,gmEDF-gsdEDF,'b', 'linewidth', 1)

    plot((-20:20)/60,gmEDRandMF+gmEDRandSTDF,'color', [0 0 0.4], 'linewidth', 1)
    plot((-20:20)/60,gmEDRandMF-gmEDRandSTDF,'color', [0 0 0.4], 'linewidth', 1)
    axis([-1 1 -0.8 0.5])
    title([pTypes{inds(i)} '   ' num2str(max(GMCCV*GML/sum(GML)))])
end
