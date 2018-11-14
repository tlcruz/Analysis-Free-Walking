function [ot] = UnilPSSValsVr(path, stL)
flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowRes.mat'], 'Flies');
seq = dt.Flies.Seq;
pTypes = unique(seq);

PSSF = cell(length(pTypes), 2, 2);
NPF = cell(length(pTypes), 2, 2);
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

                    pssrl = [];
                    pssll = [];
                    nprl = [];
                    npll = [];
                    for i = 1 : length(dt{k}.Bouts)
                        if length(dt{k}.Bouts{i}) > 100
                            vrb = vr(dt{k}.Bouts{i});
                            inds = vrb;
                            inds(inds>0) = 1;
                            inds(inds<=0) = -1;
                            if mean(flp(dt{k}.Bouts{i})) == 0
                                vr1 = vrb(1:end-delta+1);
                                vr2 = vrb(delta:end);
                                pssrr = vertcat(pssrr, mean(abs(vr2(vr1>0))));
                                nprr = vertcat(nprr, length(abs(vr2(vr1>0))));
                                psslr = vertcat(psslr, mean(abs(vr2(vr1<=0))));
                                nplr = vertcat(nplr, length(abs(vr2(vr1<=0))));
                            elseif mean(flp(dt{k}.Bouts{i})) == 180
                                vr1 = vrb(1:end-delta+1);
                                vr2 = vrb(delta:end);
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
                        NPF{j,stR,1} = vertcat(NPF{j,stR,1}, sum(nprr));
                        NPF{j,stR,2} = vertcat(NPF{j,stR,2}, sum(nplr));
                    end
                    if ~isempty(pssrl)
                        PSSF{j,stL,1} = vertcat(PSSF{j,stL,1}, pssrl'*nprl/sum(nprl));
                        PSSF{j,stL,2} = vertcat(PSSF{j,stL,2}, pssll'*npll/sum(npll));
                        NPF{j,stL,1} = vertcat(NPF{j,stL,1}, sum(nprl));
                        NPF{j,stL,2} = vertcat(NPF{j,stL,2}, sum(npll));
                    end
            end
        end
    end
end

ot.PSSF = PSSF;
ot.NPF = NPF;
ot.pTypes = pTypes;

end