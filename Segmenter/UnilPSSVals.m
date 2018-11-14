function [ot] = UnilPSSVals(path, stL)
flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowRes.mat'], 'Flies');
seq = dt.Flies.Seq;
pTypes = unique(seq);
params  = GetParams();
PSS = cell(length(pTypes), 2);
PSSF = cell(length(pTypes), 2, 2);
PS = cell(length(pTypes), 2);
NP = cell(length(pTypes), 2);
NPF = cell(length(pTypes), 2, 2);

SSPK = cell(length(pTypes), 2, 2);
NSPK = cell(length(pTypes), 2, 2);
delta = 15;
stR=1+mod(stL,2);

for n = 1 : length(flies)
    dt = load([path flies(n).name '\DataLowRes.mat'], 'Flies');
    seq = dt.Flies.Seq;
    dt = dt.Flies.Data;
    for k = 1 : length(dt)
        for j = 1 : length(pTypes)
            switch seq{k}
                case pTypes{j}
                    vr = dt{k}.Vr;
                    flp = dt{k}.flp;
                    psr = [];
                    pssr = [];
                    npsr = [];
                    pssrr = [];
                    psslr = [];
                    nprr = [];
                    nplr = [];
                    psl = [];
                    pssl = [];
                    npsl = [];
                    pssrl = [];
                    pssll = [];
                    nprl = [];
                    npll = [];
                    
                    nspkslr = [];
                    sspkslr = [];
                    nspksrr = [];
                    sspksrr = [];
                    nspksll = [];
                    sspksll = [];
                    nspksrl = [];
                    sspksrl = [];
                    for i = 1 : length(dt{k}.Bouts)
                        if length(dt{k}.Bouts{i}) > 100
                            vrb = vr(dt{k}.Bouts{i});
                            inds = vrb;
                            inds(inds>0) = 1;
                            inds(inds<=0) = -1;
                            if mean(flp(dt{k}.Bouts{i})) == 0
                                vrb = dt{k}.Vr(dt{k}.Bouts{i});
                                vfb = dt{k}.Vf(dt{k}.Bouts{i});
                                actst = dt{k}.actState(dt{k}.Bouts{i});
                                % Get spike times and forward bouts
                                [locs, ~, cmhSong, thr] = SortSpikes(vrb, actst, params);
                                [locsF, pksIF,pksEF] = CullSpikesBasedOnVf(vrb, locs, cmhSong, thr, vfb, params);
                                [locsF, pksIF, pksEF] = TemplateCompSpike(vrb,locsF,pksIF, pksEF, cmhSong, thr, params);
                                [FBouts] = GetFbouts(actst, vfb, locsF, pksIF,pksEF, params);
                                
                                pksIF(locsF<delta) = [];
                                pksEF(locsF<delta) = [];
                                locsF(locsF<delta) = [];
                                
                                vrps = vrb(locsF-delta+1);
                                locsfr = locsF(vrps>0);
                                locsfl = locsF(vrps<=0);
                                nspkslr = vertcat(nspkslr, length(locsfl));
                                %                                 sspkslr = vertcat(sspkslr, mean(abs(vrb(locsF(vrps<=0))...
                                %                                     -(vrb(pksIF(vrps<=0))+vrb(pksEF(vrps<=0)))/2)));
                                alr = mean(abs(vrb(locsF(vrps<=0)))); %vertcat(sspkslr, ));
                                nspksrr = vertcat(nspksrr, length(locsfr));
                                %                                 sspksrr = vertcat(sspksrr, mean(abs(vrb(locsF(vrps>0))...
                                %                                     -(vrb(pksIF(vrps>0))+vrb(pksEF(vrps>0)))/2)));
                                %                                 sspksrr = vertcat(sspksrr, mean(abs(vrb(locsF(vrps>0)))));
                                arr = mean(abs(vrb(locsF(vrps>0))));
                                
                                vrrr = [];
                                vrlr = [];
                                for ij = 1 : length(FBouts)
                                    vrbx = vrb(FBouts{ij});
                                    vr1 = vrbx(1:end-delta+1);
                                    vr2 = vrbx(delta:end);
                                    vrrr = vertcat(vrrr, mean(abs(vr2(vr1>0))));
                                    vrlr = vertcat(vrlr, mean(abs(vr2(vr1<=0))));
                                end
                                
%                                 sspksrr = vertcat(sspksrr, arr/nanmean(vrrr));
%                                 sspkslr = vertcat(sspkslr, alr/nanmean(vrlr));
                                sspksrr = vertcat(sspksrr, nanmean(vrrr));
                                sspkslr = vertcat(sspkslr, nanmean(vrlr));
                                
                                psr = vertcat(psr, length(inds(inds>0))/length(inds));
                                npsr = vertcat(npsr, length(inds));
                                aux1 = inds(1:end-delta+1);
                                aux2 = inds(delta:end);
                                aux3 = abs(aux1+aux2)/2;
                                pssr = vertcat(pssr, mean(aux3));
                                pssrr = vertcat(pssrr, ...
                                    length(find((aux1+aux2)/2 == 1))/(length(find(aux1 == 1))+1));
                                nprr = vertcat(nprr, length(find(aux1 == 1)));
                                psslr = vertcat(psslr, ...
                                    length(find((aux1+aux2)/2 == -1))/(length(find(aux1 == -1))+1));
                                nplr = vertcat(nplr, length(find(aux1 == -1)));
                            elseif mean(flp(dt{k}.Bouts{i})) == 180
                                vrb = dt{k}.Vr(dt{k}.Bouts{i});
                                vfb = dt{k}.Vf(dt{k}.Bouts{i});
                                actst = dt{k}.actState(dt{k}.Bouts{i});
                                % Get spike times and forward bouts
                                [locs, ~, cmhSong, thr] = SortSpikes(vrb, actst, params);
                                [locsF, pksIF,pksEF] = CullSpikesBasedOnVf(vrb, locs, cmhSong, thr, vfb, params);
                                [locsF, pksIF, pksEF] = TemplateCompSpike(vrb,locsF,pksIF, pksEF, cmhSong, thr, params);
                                [FBouts] = GetFbouts(actst, vfb, locsF, pksIF,pksEF, params);
                                
                                pksIF(locsF<delta) = [];
                                pksEF(locsF<delta) = [];
                                locsF(locsF<delta) = [];
                                
                                vrps = vrb(locsF-delta+1);
                                locsfr = locsF(vrps>0);
                                locsfl = locsF(vrps<=0);
                                nspksll = vertcat(nspksll, length(locsfl));
                                %                                 sspksll = vertcat(sspksll, mean(abs(vrb(locsF(vrps<=0))...
                                %                                     -(vrb(pksIF(vrps<=0))+vrb(pksEF(vrps<=0)))/2)));
                                %                                 sspksll = vertcat(sspksll, mean(abs(vrb(locsF(vrps<=0)))));
                                all = mean(abs(vrb(locsF(vrps<=0))));
                                nspksrl = vertcat(nspksrl, length(locsfr));
                                %                                 sspksrl = vertcat(sspksrl, mean(abs(vrb(locsF(vrps>0))...
                                %                                     -(vrb(pksIF(vrps>0))+vrb(pksEF(vrps>0)))/2)));
                                %                                 sspksrl = vertcat(sspksrl, mean(abs(vrb(locsF(vrps>0)))));
                                
                                arl = mean(abs(vrb(locsF(vrps>0))));
                                
                                vrrl = [];
                                vrll = [];
                                for ij = 1 : length(FBouts)
                                    vrbx = vrb(FBouts{ij});
                                    vr1 = vrbx(1:end-delta+1);
                                    vr2 = vrbx(delta:end);
                                    vrrl = vertcat(vrrl, mean(abs(vr2(vr1>0))));
                                    vrll = vertcat(vrll, mean(abs(vr2(vr1<=0))));
                                end
                                
%                                 sspksrl = vertcat(sspksrl, arl/nanmean(vrrl));
%                                 sspksll = vertcat(sspksll, all/nanmean(vrll));
                                sspksrl = vertcat(sspksrl, nanmean(vrrl));
                                sspksll = vertcat(sspksll, nanmean(vrll));
                                
                                psl = vertcat(psl, length(inds(inds<0))/length(inds));
                                npsl = vertcat(npsl, length(inds));
                                aux1 = inds(1:end-delta+1);
                                aux2 = inds(delta:end);
                                aux3 = abs(aux1+aux2)/2;
                                pssl = vertcat(pssl, mean(aux3));
                                pssrl = vertcat(pssrl, ...
                                    length(find((aux1+aux2)/2 ==1))/(length(find(aux1 == 1))+1));
                                nprl = vertcat(nprl, length(find(aux1 == 1)));
                                pssll = vertcat(pssll, ...
                                    length(find((aux1+aux2)/2 ==-1))/(length(find(aux1 == -1))+1));
                                npll = vertcat(npll, length(find(aux1 == -1)));
                            else
                            end
                        end
                    end
                    if ~isempty(psr)
                        SSPK{j,stR,1} = vertcat(SSPK{j,stR,1}, sspksrr'*nspksrr/sum(nspksrr));
                        SSPK{j,stR,2} = vertcat(SSPK{j,stR,2}, sspkslr'*nspkslr/sum(nspkslr));
                        NSPK{j,stR,1} = vertcat(NSPK{j,stR,1}, sum(nspksrr));
                        NSPK{j,stR,2} = vertcat(NSPK{j,stR,2}, sum(nspkslr));
                        PS{j,stR} = vertcat(PS{j,stR}, psr'*npsr/sum(npsr));
                        PSS{j,stR} = vertcat(PSS{j,stR}, pssr'*npsr/sum(npsr));
                        PSSF{j,stR,1} = vertcat(PSSF{j,stR,1}, pssrr'*nprr/sum(nprr));
                        PSSF{j,stR,2} = vertcat(PSSF{j,stR,2}, psslr'*nplr/sum(nplr));
                        NP{j,stR} = vertcat(NP{j,stR}, sum(npsr));
                        NPF{j,stR,1} = vertcat(NPF{j,stR,1}, sum(nprr));
                        NPF{j,stR,2} = vertcat(NPF{j,stR,2}, sum(nplr));
                    end
                    if ~isempty(psl)
                        SSPK{j,stL,1} = vertcat(SSPK{j,stL,1}, sspksrl'*nspksrl/sum(nspksrl));
                        SSPK{j,stL,2} = vertcat(SSPK{j,stL,2}, sspksll'*nspksll/sum(nspksll));
                        NSPK{j,stL,1} = vertcat(NSPK{j,stL,1}, sum(nspksrl));
                        NSPK{j,stL,2} = vertcat(NSPK{j,stL,2}, sum(nspksll));
                        PS{j,stL} = vertcat(PS{j,stL}, psl'*npsl/sum(npsl));
                        PSS{j,stL} = vertcat(PSS{j,stL}, pssl'*npsl/sum(npsl));
                        PSSF{j,stL,1} = vertcat(PSSF{j,stL,1}, pssrl'*nprl/sum(nprl));
                        PSSF{j,stL,2} = vertcat(PSSF{j,stL,2}, pssll'*npll/sum(npll));
                        NP{j,stL} = vertcat(NP{j,stL}, sum(npsl));
                        NPF{j,stL,1} = vertcat(NPF{j,stL,1}, sum(nprl));
                        NPF{j,stL,2} = vertcat(NPF{j,stL,2}, sum(npll));
                    end
            end
        end
    end
end

ot.PSS = PSS;
ot.NP = NP;
ot.PSSF = PSSF;
ot.NPF = NPF;
ot.pTypes = pTypes;
ot.SSPK = SSPK;
ot.NSPK = NSPK;

end