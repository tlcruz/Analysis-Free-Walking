function [ot] = UnilPSSValsVrForwBouts(path, stL)
flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowRes.mat'], 'Flies');
seq = dt.Flies.Seq;
pTypes = unique(seq);

PSSF = cell(length(pTypes), 2, 2);
NPF = cell(length(pTypes), 2, 2);
SPKPROM = cell(length(pTypes), 2, 2);
SPKSZ = cell(length(pTypes), 2, 2);
NSPK = cell(length(pTypes), 2, 2);

STR = cell(length(pTypes), length(flies));
N = cell(length(pTypes), length(flies));


delta = 15;
stR=1+mod(stL,2);
params  = GetParams();
for n = 1 : length(flies)
    dt = load([path flies(n).name '\DataLowRes.mat'], 'Flies');
    seq = dt.Flies.Seq;
    dt = dt.Flies.Data;
    for k = 1 : length(dt)
        for j = 1 : length(pTypes)
            switch seq{k}
                case pTypes{j}
                    vrb = dt{k}.Vr;
                    vfb = dt{k}.Vf;
                    vtb = dt{k}.Vt;
                    xx = dt{k}.X;
                    yy = dt{k}.Y;
                    flp = dt{k}.flp;
                    wd = dt{k}.WallDist;
                    actst = dt{k}.actState;
                    % Get spike times and forward bouts
                    actst(wd < params.mDistWall) = 0;
                    
                    [locs, ~, cmhSong, thr] = SortSpikes(vrb, actst, params);
                    [locsF, pksIF,pksEF] = CullSpikesBasedOnVf(vrb, locs, cmhSong, thr, vfb, params);
                    [locsF, pksIF, pksEF] = TemplateCompSpike(vrb,locsF,pksIF, pksEF, cmhSong, thr, params);
                    [FBouts] = GetFbouts(actst, vfb, locsF, pksIF,pksEF, params);
                    
                    str = [];
                    angd = [];
                    dist = [];
                    nn = [];
                    for pp = 1 : length(FBouts)
                        % Calculate Straightness
                        if length(FBouts{pp}) > max(params.windStr+1, params.minStrB);
                            xaux = xx(FBouts{pp});
                            yaux = yy(FBouts{pp});
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
                            dist = vertcat(dist, sum(ss));
                            nn = vertcat(nn, length(xaux));
                            %         ANGD = vertcat(ANGD, sum(abs(vraux))/sum(vtaux));
                            angd = vertcat(angd, mean(abs(vraux)));
                        end
                    end
                    
                    if ~isempty(str)
                        STR{j,n} = str;
                        N{j,n} = nn;
                    end
                    
                    spksrr = [];
                    spkslr = [];
                    spksrl = [];
                    spksll = [];
                    for i = 1 : length(locsF)
                        if flp(locsF(i)) == 0
                            if vrb(locsF(i))>0
                                spksrr = vertcat(spksrr, vrb(locsF(i))...
                                    -((vrb(locsF(i)-pksIF(i))+vrb(locsF(i)+pksEF(i)))/2));
                            else
                                spkslr = vertcat(spkslr, vrb(locsF(i))...
                                    -((vrb(locsF(i)-pksIF(i))+vrb(locsF(i)+pksEF(i)))/2));
                            end
                        else
                            if vrb(locsF(i))>0
                                spksrl = vertcat(spksrl, vrb(locsF(i))...
                                    -((vrb(locsF(i)-pksIF(i))+vrb(locsF(i)+pksEF(i)))/2));
                            else
                                spksll = vertcat(spksll, vrb(locsF(i))...
                                    -((vrb(locsF(i)-pksIF(i))+vrb(locsF(i)+pksEF(i)))/2));
                            end
                        end
                    end
                    
                    pssrr = [];
                    psslr = [];
                    nprr = [];
                    nplr = [];
                    pssrl = [];
                    pssll = [];
                    nprl = [];
                    npll = [];
                    for i = 1 : length(FBouts)
                        if length(FBouts{i}) > 20
                            vrbx = vrb(FBouts{i});
                            if mean(flp(FBouts{i})) == 0
                                vr1 = vrbx(1:end-delta+1);
                                vr2 = vrbx(delta:end);
                                pssrr = vertcat(pssrr, mean(abs(vr2(vr1>0))));
                                nprr = vertcat(nprr, length(abs(vr2(vr1>0))));
                                psslr = vertcat(psslr, mean(abs(vr2(vr1<=0))));
                                nplr = vertcat(nplr, length(abs(vr2(vr1<=0))));
                            elseif mean(flp(FBouts{i})) == 180
                                vr1 = vrbx(1:end-delta+1);
                                vr2 = vrbx(delta:end);
                                pssrl = vertcat(pssrl, mean(abs(vr2(vr1>0))));
                                nprl = vertcat(nprl, length(abs(vr2(vr1>0))));
                                pssll = vertcat(pssll, mean(abs(vr2(vr1<=0))));
                                npll = vertcat(npll, length(abs(vr2(vr1<=0))));
                            else
                            end
                        end
                    end
                    if ~isempty(pssrr)
                        PSSF{j,stR,1} = vertcat(PSSF{j,stR,1}, pssrr'*nprr/sum(nprr));
                        PSSF{j,stR,2} = vertcat(PSSF{j,stR,2}, psslr'*nplr/sum(nplr));
                        
                        SPKPROM{j,stR,1} = vertcat(SPKPROM{j,stR,1}, ...
                            nanmean(abs(spksrr))/(pssrr'*nprr/sum(nprr)));
                        SPKPROM{j,stR,2} = vertcat(SPKPROM{j,stR,2}, ...
                            nanmean(abs(spkslr))/(psslr'*nplr/sum(nplr)));
                        
                        SPKSZ{j,stR,1} = vertcat(SPKSZ{j,stR,1}, ...
                            nanmean(abs(spksrr)));
                        SPKSZ{j,stR,2} = vertcat(SPKSZ{j,stR,2}, ...
                            nanmean(abs(spkslr)));
                        
                        NSPK{j,stR,1} = vertcat(NSPK{j,stR,1}, length(spksrr));
                        NSPK{j,stR,2} = vertcat(NSPK{j,stR,2}, length(spkslr));
                        
                        NPF{j,stR,1} = vertcat(NPF{j,stR,1}, sum(nprr));
                        NPF{j,stR,2} = vertcat(NPF{j,stR,2}, sum(nplr));
                    end
                    if ~isempty(pssrl)
                        PSSF{j,stL,1} = vertcat(PSSF{j,stL,1}, pssrl'*nprl/sum(nprl));
                        PSSF{j,stL,2} = vertcat(PSSF{j,stL,2}, pssll'*npll/sum(npll));
                        
                        SPKPROM{j,stL,1} = vertcat(SPKPROM{j,stL,1}, ...
                            mean(abs(spksrl))/(pssrl'*nprl/sum(nprl)));
                        SPKPROM{j,stL,2} = vertcat(SPKPROM{j,stL,2}, ...
                            mean(abs(spksll))/(pssll'*npll/sum(npll)));
                        
                        SPKSZ{j,stL,1} = vertcat(SPKSZ{j,stL,1}, ...
                            mean(abs(spksrl)));
                        SPKSZ{j,stL,2} = vertcat(SPKSZ{j,stL,2}, ...
                            mean(abs(spksll)));
                        
                        NSPK{j,stL,1} = vertcat(NSPK{j,stL,1}, length(spksrl));
                        NSPK{j,stL,2} = vertcat(NSPK{j,stL,2}, length(spksll));
                        
                        NPF{j,stL,1} = vertcat(NPF{j,stL,1}, sum(nprl));
                        NPF{j,stL,2} = vertcat(NPF{j,stL,2}, sum(npll));
                    end
            end
        end
    end
end

ot.PSSF = PSSF;
ot.NPF = NPF;
ot.SPKPROM = SPKPROM;
ot.SPKSZ = SPKSZ;
ot.NSPK = NSPK;
ot.pTypes = pTypes;
ot.STR = STR;
ot.N = N;

end