function [STRALL, DSTALL, NSTRALL, PSSALL, NPSSALL, pTypes, ANGDALL, VTALL] = GetStrAndVisInf(path, params)
if nargin < 1
    path = uigetdir;
    path = [path '\'];
    params  = GetParams();
elseif nargin == 1
    params  = GetParams();
end

flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowRes.mat'], 'Flies');
seq = dt.Flies.Seq;
pTypes = unique(seq);
disp(path)
STRALL = cell(length(pTypes), length(flies));
ANGDALL = cell(length(pTypes), length(flies));
VTALL = cell(length(pTypes), length(flies));
NSTRALL = cell(length(pTypes), length(flies));
PSSALL = cell(length(pTypes), length(flies));
NPSSALL = cell(length(pTypes), length(flies));
DSTALL = cell(length(pTypes), length(flies));
for n = 1 : length(flies)
    dt = load([path flies(n).name '\DataLowRes.mat'], 'Flies');
    seq = dt.Flies.Seq;
    dt = dt.Flies.Data;
    for k = 1 : length(dt)
        for j = 1 : length(pTypes)
            str = [];
            angd = [];
            vt = [];
            dst = [];
            nstr = [];
            npss = [];
            pss = [];
            switch seq{k}
                case pTypes{j}
                    % Load walking params for this bout
                    vrb = dt{k}.Vr;
                    vfb = dt{k}.Vf;
                    vtb = dt{k}.Vt;
                    xtb = dt{k}.X;
                    ytb = dt{k}.Y;
                    flp = dt{k}.flp;
                    wd = dt{k}.WallDist;
                    actst = dt{k}.actState;
                    % Get spike times and forward bouts
                    actst(wd < params.mDistWall) = 0;
                    
                    [locs, ~, cmhSong, thr] = SortSpikes(vrb, actst, params);
                    [locsF, pksIF,pksEF] = CullSpikesBasedOnVf(vrb, locs, cmhSong, thr, vfb, params);
                    [locsF, pksIF, pksEF] = TemplateCompSpike(vrb,locsF,pksIF, pksEF, cmhSong, thr, params);
                    [FBouts] = GetFbouts(actst, vfb, locsF, pksIF,pksEF, params);
                    
                    for pp = 1 : length(FBouts)
                        % Calculate Straightness
                        if length(FBouts{pp}) > max(params.windStr+1, params.minStrB);
                            xaux = xtb(FBouts{pp});
                            yaux = ytb(FBouts{pp});
                            vtaux = vtb(FBouts{pp});
                            vraux = vrb(FBouts{pp});
                            nv = length(xaux) - params.windStr;
                            dm = [];
                            ss = [];
                            for l = 1 : nv
                                pi = [xaux(l),yaux(l),0];
                                pf = [xaux(l+params.windStr),yaux(l+params.windStr),0];
                                dis = sum(vtaux((l):(l+params.windStr)))/60;
                                pt = [xaux(l+floor(params.windStr/2)),...
                                    yaux(l+floor(params.windStr/2)),0];
                                dr = point_to_line(pt,pi,pf);
                                dm = vertcat(dm, dr);
                                ss = vertcat(ss, dis);
                            end
                            str = vertcat(str, sum(ss)/sum(dm));
                            dst = vertcat(dst, sum(ss));
                            nstr = vertcat(nstr, length(xaux));
                            vt = vertcat(vt, mean(vtaux));
                            angd = vertcat(angd, sum(abs(vraux))/sum(vtaux));
                        end
                    end
                    
                    % Calculate Probability Same Side
                    %                     for i = 1 : length(dt{k}.Bouts)
                    %                         if length(dt{k}.Bouts{i}) > params.thr
                    %                             vrb = dt{k}.Vr(dt{k}.Bouts{i});
                    %                             vfb = dt{k}.Vf(dt{k}.Bouts{i});
                    %                             [BoutVr,~,~,BoutP] = GetProbVec(vrb,vfb, params.delta,1,-inf,inf);
                    %                             npss = vertcat(npss,length(BoutP));
                    %                             pss = vertcat(pss,mean(BoutP));
                    % %                             pss = vertcat(pss,mean(abs(BoutVr)));
                    %                         end
                    %                     end
                    
%                     
%                     fc = 1;
%                     fs = 60;
% %                     [b,a] = butter(6,fc/(fs/2));
%                     vrb2 = filtfilt(b,a,vrb);

%                     for i = 1 : length(dt{k}.Bouts)
                    for i = 1 : length(FBouts)
%                         if length(dt{k}.Bouts{i}) > params.thr
                        if length(FBouts{i}) > max(params.windStr+1, params.minStrB);
%                             vrb = dt{k}.Vr(dt{k}.Bouts{i});
%                             vfb = dt{k}.Vf(dt{k}.Bouts{i});
                            vrb = dt{k}.Vr(FBouts{i});
                            vfb = dt{k}.Vf(FBouts{i});
                            [BoutVr,~,~,BoutP] = GetProbVec(vrb,vfb, params.delta,1,-inf,inf);
                            npss = vertcat(npss,length(BoutP));
                            pss = vertcat(pss,mean(BoutP));
%                             pss = vertcat(pss,mean(abs(BoutVr)));
                        end
                    end
                    
                    
                    %                     for i = 1 : length(FBouts)
                    %                         if length(FBouts{pp}) > max(params.windStr+1, params.minStrB);
                    %                             vrbvi = vrb(FBouts{i});
                    %                             vfbvi = vfb(FBouts{i});
                    %                             [BoutVr,~,~,BoutP] = GetProbVec(vrbvi,vfbvi, params.delta,1,-inf,inf);
                    %                             npss = vertcat(npss,length(BoutP));
                    %                             pss = vertcat(pss,mean(BoutP));
                    % %                             pss = vertcat(pss,mean(abs(BoutVr)));
                    %                         end
                    %                     end
                    
                    DSTALL{j,n} = vertcat(DSTALL{j,n}, dst);
                    STRALL{j,n} = vertcat(STRALL{j,n}, str);
                    NSTRALL{j,n} = vertcat(NSTRALL{j,n}, nstr);
                    PSSALL{j,n} = vertcat(PSSALL{j,n}, pss);
                    NPSSALL{j,n} = vertcat(NPSSALL{j,n}, npss);
                    ANGDALL{j,n} = vertcat(ANGDALL{j,n}, angd);
                    VTALL{j,n} = vertcat(VTALL{j,n}, vt);
            end
        end
    end
    disp([num2str(floor(100*n/length(flies))) '% Done' ])
end
end