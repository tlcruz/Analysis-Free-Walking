%% Function to calculate straightness and visual influence
function [STRALL, DSTALL, NSTRALL, PSSALL, VRSSALL, NPSSALL, pTypes, ANGDALL, VTALL] = GetStrAndVisInf(path, params)
% if path is not provided display a path selection window
addpath('Saccade Segmenter')
if nargin < 1
    path = uigetdir;
    path = [path '\'];
    params  = GetParams();
elseif nargin == 1
    params  = GetParams();
end

% get the list of flies to analyse
flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowRes.mat'], 'Flies');
seq = dt.Flies.Seq;
pTypes = unique(seq);
disp(path)

% initialize variables for each fly and trial type 
STRALL = cell(length(pTypes), length(flies));
ANGDALL = cell(length(pTypes), length(flies));
VTALL = cell(length(pTypes), length(flies));
NSTRALL = cell(length(pTypes), length(flies));
PSSALL = cell(length(pTypes), length(flies));
VRSSALL = cell(length(pTypes), length(flies));
NPSSALL = cell(length(pTypes), length(flies));
DSTALL = cell(length(pTypes), length(flies));

for n = 1 : length(flies)
    % load data for a specific fly
    dt = load([path flies(n).name '\DataLowRes.mat'], 'Flies');
    seq = dt.Flies.Seq;
    dt = dt.Flies.Data;
    % iterate across protocol segments
    for k = 1 : length(dt)
        % iterate across trial types
        for j = 1 : length(pTypes)
            str = [];
            angd = [];
            vt = [];
            dst = [];
            nstr = [];
            npss = [];
            pss = [];
            vrpss = [];
            switch seq{k}
                % when trial type matches a protocol segment
                case pTypes{j}
                    % Load walking parameters for this protocol segment
                    vrb = dt{k}.Vr;
                    vfb = dt{k}.Vf;
                    vtb = dt{k}.Vt;
                    xtb = dt{k}.X;
                    ytb = dt{k}.Y;
                    wd = dt{k}.WallDist;
                    actst = dt{k}.actState;
                    % Get spike times and forward bouts
                    actst(wd < params.mDistWall) = 0;
                    
                    % segment angular spikes and forward segments
                    [locs, ~, cmhSong, thr] = SortSpikes(vrb, actst, params);
                    [locsF, pksIF,pksEF] = CullSpikesBasedOnVf(vrb, locs, cmhSong, thr, vfb, params);
                    [locsF, pksIF, pksEF] = TemplateCompSpike(vrb,locsF,pksIF, pksEF, cmhSong, thr, params);
                    [FBouts] = GetFbouts(actst, vfb, locsF, pksIF,pksEF, params);
                    
                    % for each forward segment 
                    for pp = 1 : length(FBouts)
                        % Calculate Straightness for a specific forward
                        % segment
                        if length(FBouts{pp}) > max(params.windStr+1, params.minStrB)
                            xaux = xtb(FBouts{pp});
                            yaux = ytb(FBouts{pp});
                            vtaux = vtb(FBouts{pp});
                            vraux = vrb(FBouts{pp});
                            nv = length(xaux) - params.windStr;
                            dm = [];
                            ss = [];
                            % calculate difference of a point from straight line
                            % within a defined window centered on the point 
                            for l = 1 : nv
                                % initial and last point of the window
                                pi = [xaux(l),yaux(l),0];
                                pf = [xaux(l+params.windStr),yaux(l+params.windStr),0];
                                % distance covered within that window
                                dis = sum(vtaux((l):(l+params.windStr)))/60;
                                % get central point
                                pt = [xaux(l+floor(params.windStr/2)),...
                                    yaux(l+floor(params.windStr/2)),0];
                                % get distance from the central point to
                                % the line defined by the initial and last
                                % point
                                dr = point_to_line(pt,pi,pf);
                                dm = vertcat(dm, dr);
                                ss = vertcat(ss, dis);
                            end
                            % save straightness by dividing the total
                            % distance by the distance to straight line
                            str = vertcat(str, sum(ss)/sum(dm));
                            % save total distance
                            dst = vertcat(dst, sum(ss));
                            % save forward segment duration
                            nstr = vertcat(nstr, length(xaux));
                            % save average speed
                            vt = vertcat(vt, mean(vtaux));
                            % save angular deflection
                            angd = vertcat(angd, sum(abs(vraux))/sum(vtaux));
                        end
                    end
                    
                    % for each walking bout 
                    for i = 1 : length(dt{k}.Bouts)
                        if length(dt{k}.Bouts{i}) > params.thr
                            vrb = dt{k}.Vr(dt{k}.Bouts{i});
                            vfb = dt{k}.Vf(dt{k}.Bouts{i});
                            % calculate probability same side with a time
                            % difference delta
                            [BoutVr,~,~,BoutP] = GetProbVec(vrb,vfb, params.delta,1,-inf,inf);
                            npss = vertcat(npss,length(BoutP));
                            pss = vertcat(pss,mean(BoutP));
                            vrpss = vertcat(vrpss,mean(abs(BoutVr)));
                        end
                    end
                    
                    DSTALL{j,n} = vertcat(DSTALL{j,n}, dst);
                    STRALL{j,n} = vertcat(STRALL{j,n}, str);
                    NSTRALL{j,n} = vertcat(NSTRALL{j,n}, nstr);
                    PSSALL{j,n} = vertcat(PSSALL{j,n}, pss);
                    VRSSALL{j,n} = vertcat(VRSSALL{j,n}, vrpss);
                    NPSSALL{j,n} = vertcat(NPSSALL{j,n}, npss);
                    ANGDALL{j,n} = vertcat(ANGDALL{j,n}, angd);
                    VTALL{j,n} = vertcat(VTALL{j,n}, vt);
            end
        end
    end
    disp([num2str(floor(100*n/length(flies))) '% Done' ])
end
end