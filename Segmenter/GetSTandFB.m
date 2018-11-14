function [stfb, pTypes] = GetSTandFB(path)
flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowRes.mat'], 'Flies');
seq = dt.Flies.Seq;
pTypes = unique(seq);

pST = cell(length(flies),length(pTypes));
pFB = cell(length(flies),length(pTypes));
angDisp = cell(length(flies),length(pTypes));
NangDisp = cell(length(flies),length(pTypes));
pT = cell(length(flies),length(pTypes));
lFB = cell(length(flies),length(pTypes));
lST = cell(length(flies),length(pTypes));
biasFB = cell(length(flies),length(pTypes));
NbiasFB = cell(length(flies),length(pTypes));
valFB = cell(length(flies),length(pTypes));
NvalFB = cell(length(flies),length(pTypes));
biasST1 = cell(length(flies),length(pTypes));
NbiasST1 = cell(length(flies),length(pTypes));
biasST20 = cell(length(flies),length(pTypes));
NbiasST20 = cell(length(flies),length(pTypes));
biasST100 = cell(length(flies),length(pTypes));
NbiasST100 = cell(length(flies),length(pTypes));

disp(path)
for n = 1 : length(flies)
    
    dt = load([path flies(n).name '\DataLowRes.mat'], 'Flies');
    seq = dt.Flies.Seq;
    dt = dt.Flies.Data;
    for k = 1 : length(dt)
        for j = 1 : length(pTypes)
            
            pfb = [];
            angD = [];
            nangD = [];
            pst = [];
            pt = [];
            lfb = [];
            lst = [];
            valfb = [];
            nvalfb = [];
            bfb = [];
            nbfb = [];
            bst1 = [];
            bst20 = [];
            bst100 = [];
            
            switch seq{k}
                case pTypes{j}
                    Vf = dt{k}.Vf;
                    Vs = dt{k}.Vs;
                    Vr = dt{k}.Vr;
                    Vt = dt{k}.Vt;
                    acSt = dt{k}.actState;
                    
                    [params] = GetParams();
                    if(strcmp(pTypes{j}(end-1:end),'NG'))
                        params.cutoff = 0.15;
                    else
                        params.cutoff = 0.24;
                    end
                    
                    
                    [locs, ~, cmhSong, thr, dx] = SortSpikes(Vr, acSt, params);
                    [locsF, pksIF, pksEF, ~] = CullSpikesBasedOnVf(Vr, locs, cmhSong, thr, Vf, params);
                    [locsF, pksIF, pksEF, scores] = TemplateCompSpike(Vr,locsF,pksIF, pksEF, cmhSong, thr, params);
                    [FWBoutStr, FWBout] = GetFbouts(acSt, Vf, locsF, pksIF,pksEF, params);
                    
                    FVect = zeros(length(Vr),1);
                    for i = 1 : length(FWBout)
                        FVect(FWBout{i}) = 1;
                        lfb = vertcat(lfb, length(FWBout{i}));
                    end
                    vect = zeros(length(Vr),1);
                    for i = 1 : length(locsF)
                        vect((locsF(i)-pksIF(i)):(locsF(i)+pksEF(i))) = 1;
                        lst = vertcat(lst, Vr(locsF(i)));
                    end
                    npss = [];
                    pss = [];
                    for i = 1 : length(FWBout)
                        if length(FWBout{i}) > 10
                            npss = vertcat(npss,length(FWBout{i}));
                            pss = vertcat(pss,abs(mean(Vr(FWBout{i})))/mean(abs(Vr(FWBout{i}))));
                            valfb = vertcat(valfb, abs(mean(Vr(FWBout{i}))));
                            nvalfb = vertcat(nvalfb, length(FWBout{i}));
                            
                            angD = vertcat(angD, sum(abs(Vr(FWBout{i})))/sum(Vt(FWBout{i})));
                            nangD = vertcat(nangD, length(FWBout{i}));
                        end
                    end
                    tw = 2;
                    w = 6;
                    for i = 1 : length(locsF)
                        if (locsF(i)-pksIF(i)) > w+3
                            auxVr = Vr((locsF(i)-pksIF(i) - w - tw):(locsF(i)-pksIF(i) - w + tw));
                            auxVf = Vf((locsF(i)-pksIF(i) - w - tw):(locsF(i)-pksIF(i) - w + tw));
                            if abs(auxVr) > 1
                                if sign(mean(auxVr)) == sign(Vr(locsF(i)))
                                    bst1 = vertcat(bst1, 1);
                                else
                                    bst1 = vertcat(bst1, 0);
                                end
                            end
                            if abs(auxVr) > 20
                                if sign(mean(auxVr)) == sign(Vr(locsF(i)))
                                    bst20 = vertcat(bst20, 1);
                                else
                                    bst20 = vertcat(bst20, 0);
                                end
                            end
                            if abs(auxVr) > 100
                                if sign(mean(auxVr)) == sign(Vr(locsF(i)))
                                    bst100 = vertcat(bst100, 1);
                                else
                                    bst100 = vertcat(bst100, 0);
                                end
                            end
                        end
                    end
                    
                    pfb = vertcat(pfb, sum(abs(FVect))/sum(abs(acSt)));
                    pst = vertcat(pst, sum(abs(vect))/sum(abs(acSt)));
                    pt = vertcat(pt, sum(abs(acSt)));
                    bfb = vertcat(bfb, pss'*npss/sum(npss));
                    nbfb = vertcat(nbfb, sum(npss));
                    
            end
            
            pFB{n,j} = vertcat(pFB{n,j}, pfb);
            pST{n,j} = vertcat(pST{n,j}, pst);
            angDisp{n,j} = vertcat(angDisp{n,j}, angD);
            NangDisp{n,j} = vertcat(NangDisp{n,j}, nangD);
            pT{n,j} = vertcat(pT{n,j}, pt);
            lFB{n,j} = vertcat(lFB{n,j}, lfb);
            lST{n,j} = vertcat(lST{n,j}, lst);
            valFB{n,j} = vertcat(valFB{n,j}, valfb);
            NvalFB{n,j} = vertcat(NvalFB{n,j}, nvalfb);
            biasFB{n,j} = vertcat(biasFB{n,j}, bfb);
            NbiasFB{n,j} = vertcat(NbiasFB{n,j}, nbfb);
            biasST1{n,j} = vertcat(biasST1{n,j}, mean(bst1));
            NbiasST1{n,j} = vertcat(NbiasST1{n,j}, length(bst1));
            biasST20{n,j} = vertcat(biasST20{n,j}, mean(bst20));
            NbiasST20{n,j} = vertcat(NbiasST20{n,j}, length(bst20));
            biasST100{n,j} = vertcat(biasST100{n,j}, mean(bst100));
            NbiasST100{n,j} = vertcat(NbiasST100{n,j}, length(bst100));
            
        end
    end
    disp([num2str(floor(100*n/length(flies))) '% Done' ])
end

stfb.pST = pST;
stfb.angDisp = angDisp;
stfb.NangDisp = NangDisp;
stfb.pFB = pFB;
stfb.pT = pT;
stfb.lFB = lFB;
stfb.lST = lST;
stfb.valFB = valFB;
stfb.NvalFB = NvalFB;
stfb.biasFB = biasFB;
stfb.NbiasFB = NbiasFB;
stfb.biasST1 = biasST1;
stfb.NbiasST1 = NbiasST1;
stfb.biasST20 = biasST20;
stfb.NbiasST20 = NbiasST20;
stfb.biasST100 = biasST100;
stfb.NbiasST100 = NbiasST100;

end