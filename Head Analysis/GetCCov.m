function [CCVAll, NCCVAll, CCVFW, NCCVFW, pTypes, DVRFW, DHANG, rd] = GetCCov(path, params)
% get the list of flies to analyse
addpath('Saccade Segmenter')

flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowAndHighRes.mat']);
pTypes = dt.Fly.Seq;
pTypes = unique(pTypes);

% initialize variables for each fly and trial type 
CCVFW = cell(length(pTypes), length(flies));
NCCVFW = cell(length(pTypes), length(flies));
CCVAll = cell(length(pTypes), length(flies));
NCCVAll = cell(length(pTypes), length(flies));
DVRFW = cell(length(pTypes), length(flies));
DHANG = cell(length(pTypes), length(flies));
DANGDEV = cell(length(pTypes), length(flies));

rd.CCVFWMeanDist = cell(length(pTypes), length(flies));
rd.CCVFWSTDDist = cell(length(pTypes), length(flies));
rd.CCVFWNDist = cell(length(pTypes), length(flies));
rd.CCVALLMeanDist = cell(length(pTypes), length(flies));
rd.CCVALLSTDDist = cell(length(pTypes), length(flies));
rd.CCVALLNDist = cell(length(pTypes), length(flies));

rd.NSpks = cell(length(pTypes), length(flies));
rd.Nfbts = cell(length(pTypes), length(flies));

for n = 1:length(flies)
    % load data for a specific fly
    if exist([path flies(n).name '\DataLowAndHighRes.mat'], 'file') == 2
        dt = load([path flies(n).name '\DataLowAndHighRes.mat']);
        seq = dt.Fly.Seq;
        dt = dt.Fly.Data;
        % iterate across protocol segments
        for k = 1 : length(dt)
            % iterate across trial types
            for l = 1 : length(pTypes)
                xcovsall = [];
                lbtsall = [];
                xcovsfw = [];
                lbtsfw = [];
                xcovfwRand = [];
                nxcovfwRand = [];
                xcovallRand = [];
                nxcovallRand = [];
                
                switch seq{k}
                    % when trial type matches a protocol segment
                    case pTypes{l}
                        % Load walking parameters for this protocol segment
                        vrb = dt{k}.Vr;
                        vfb = dt{k}.Vf;
                        wd = dt{k}.WallDist;
                        actst = dt{k}.actState;
                        hA = dt{k}.HangDS;
                        actst(wd < params.mDistWall) = 0;
                        vrx = vrb(actst == 1);
                        % segment angular spikes and forward segments
                        [locs, ~, cmhSong, thr] = SortSpikes(vrb, actst, params);
                        [locsF, pksIF,pksEF] = CullSpikesBasedOnVf(vrb, locs, cmhSong, thr, vfb, params);
                        [locsF, pksIF, pksEF] = TemplateCompSpike(vrb,locsF,pksIF, pksEF, cmhSong, thr, params);
                        [FBouts] = GetFbouts(actst, vfb, locsF, pksIF,pksEF, params);
                        
                        
                        if isempty(rd.NSpks{l,n})
                            rd.NSpks{l,n} = length(locsF);
                        else
                            rd.NSpks{l,n} = rd.NSpks{l,n} + length(locsF);
                        end
                        
                        if isempty(rd.Nfbts{l,n})
                            rd.Nfbts{l,n} = 0;
                        end
                        
                        % Get covariances during forward segments
                        for j = 1:length(FBouts)
                            if length(FBouts{j})>params.btthr && FBouts{j}(end) < length(hA)
                                vr = vrb(FBouts{j});
                                hang = hA(FBouts{j});
                                vr(isnan(hang)) = [];
                                hang(isnan(hang)) = [];
                                vf = vfb(FBouts{j});
                                vf(isnan(hang)) = [];
                                angDev = sum(abs(vr))/sum(vf);
                                DANGDEV{l,n} = horzcat(DANGDEV{l,n}, angDev);
                                if angDev < 1.25
                                    if (length(hang)>(params.WindowCCVFW+params.btthr))
                                        [x,~] = xcov(vr,hang,params.WindowCCVFW,'coeff');
                                        xcovsfw = horzcat(xcovsfw,x);
                                        lbtsfw = vertcat(lbtsfw,length(hang));
                                    end
                                end
                                rd.Nfbts{l,n} = rd.Nfbts{l,n}+1;
                                % bootstrap and get noise level covariances
                                for i = 1 : params.nBoots
                                    ind = randi(length(vrx)-length(FBouts{j}));
                                    vrBoot = vrx(ind:(ind+length(FBouts{j})-1));
                                    if (length(hang)>(params.WindowCCVFW+params.btthr)&& length(hang)==length(vrBoot))
                                        [x,~] = xcov(vrBoot,hang,params.WindowCCVFW,'coeff');
                                        xcovfwRand = horzcat(xcovfwRand, x);
                                        nxcovfwRand = vertcat(nxcovfwRand, length(hang));
                                    end
                                end
                                
                            end
                        end
                        
                        if(~isempty(xcovsfw))
                            % get mean covariance per fly per trial type
                            CCVFW{l,n} = horzcat(CCVFW{l,n}, xcovsfw*lbtsfw/sum(lbtsfw));
                            NCCVFW{l,n} = vertcat(NCCVFW{l,n}, sum(lbtsfw));
                            mRandDist = xcovfwRand*nxcovfwRand/sum(nxcovfwRand);
                            sDRandDist = sqrt(((bsxfun(@minus, xcovfwRand, mRandDist).*...
                                bsxfun(@minus, xcovfwRand, mRandDist)*nxcovfwRand))/sum(nxcovfwRand));
                            if ~isempty(mRandDist)
                                rd.CCVFWMeanDist{l,n} = horzcat(rd.CCVFWMeanDist{l,n}, mRandDist);
                                rd.CCVFWSTDDist{l,n} = horzcat(rd.CCVFWSTDDist{l,n}, sDRandDist);
                                rd.CCVFWNDist{l,n} = vertcat(rd.CCVFWNDist{l,n}, sum(nxcovfwRand));
                            end
                        end
                        
                        % Get covariances during walking bouts
                        Bouts = dt{k}.Bouts;
                        for j = 1:length(Bouts)
                            if length(Bouts{j})>params.btthr && Bouts{j}(end) < length(hA)
                                vr = vrb(Bouts{j});
                                hang = hA(Bouts{j});
                                vr(isnan(hang)) = [];
                                hang(isnan(hang)) = [];
                                if (length(hang)>(params.WindowCCVALL+params.btthr))
                                    [x,~] = xcov(vr,hang,params.WindowCCVALL,'coeff');
                                    xcovsall = horzcat(xcovsall,x);
                                    lbtsall = vertcat(lbtsall,length(hang));
                                end
                                % bootstrap and get noise level covariances
                                for i = 1 : params.nBoots
                                    ind = randi(length(vrx)-length(Bouts{j}));
                                    vrBoot = vrx(ind:(ind+length(Bouts{j})-1));
                                    if (length(hang)>(params.WindowCCVALL+params.btthr) && length(hang)==length(vrBoot))
                                        [x,~] = xcov(vrBoot,hang,params.WindowCCVALL,'coeff');
                                        xcovallRand = horzcat(xcovallRand, x);
                                        nxcovallRand = vertcat(nxcovallRand, length(hang));
                                    end
                                end
                                
                            end
                        end
                        
                        % get mean covariance per fly per trial type
                        if(~isempty(xcovsall))
                            CCVAll{l,n} = horzcat(CCVAll{l,n}, xcovsall*lbtsall/sum(lbtsall));
                            NCCVAll{l,n} = vertcat(NCCVAll{l,n}, sum(lbtsall));
                            mRandDist = xcovallRand*nxcovallRand/sum(nxcovallRand);
                            sDRandDist = sqrt(((bsxfun(@minus, xcovallRand, mRandDist).*...
                                bsxfun(@minus, xcovallRand, mRandDist)*nxcovallRand))/sum(nxcovallRand));
                            if ~isempty(mRandDist)
                                rd.CCVALLMeanDist{l,n} = horzcat(rd.CCVALLMeanDist{l,n}, mRandDist);
                                rd.CCVALLSTDDist{l,n} = horzcat(rd.CCVALLSTDDist{l,n}, sDRandDist);
                                rd.CCVALLNDist{l,n} = vertcat(rd.CCVALLNDist{l,n}, sum(nxcovallRand));
                            end
                        end
                        
                        % get average angular speed and head angle
                        dVr = [];
                        dHA = [];
                        for j = 1 : length(FBouts)
                            if length(FBouts{j}) > 20 && FBouts{j}(end) < length(hA) 
                                vr = vrb(FBouts{j});
                                vf = vrb(FBouts{j});
                                angDev = sum(pi*abs(vr)/180)/sum(vf);
                                hang = hA(FBouts{j});
                                if angDev > 2.5
                                    for f = 1 : length(FBouts{j})-20
                                        dVr = vertcat(dVr, mean(vr(f:(f+19))));
                                        dHA = vertcat(dHA, mean(hang(f:(f+19))-hang(f)));
                                    end
                                end
                            end
                        end
                        if(~isempty(dHA))
                            DVRFW{l,n} = vertcat(DVRFW{l,n}, dVr);
                            DHANG{l,n} = vertcat(DHANG{l,n}, dHA);
                        end  
                end
            end
        end
    end
    disp([num2str(n) ' out of ' num2str(length(flies))])
end
end