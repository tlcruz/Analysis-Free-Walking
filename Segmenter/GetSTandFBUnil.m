function [stfb, pTypes] = GetSTandFBUnil(path)
flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowRes.mat'], 'Flies');
seq = dt.Flies.Seq;
pTypes = unique(seq);

pST = cell(length(flies),length(pTypes),2);
pFB = cell(length(flies),length(pTypes),2);
pT = cell(length(flies),length(pTypes),2);
lFB = cell(length(flies),length(pTypes),2);
lST = cell(length(flies),length(pTypes),2);
biasFB = cell(length(flies),length(pTypes),2);
NbiasFB = cell(length(flies),length(pTypes),2);
valFB = cell(length(flies),length(pTypes),2);
NvalFB = cell(length(flies),length(pTypes),2);
biasST1 = cell(length(flies),length(pTypes),2);
NbiasST1 = cell(length(flies),length(pTypes),2);
biasST20 = cell(length(flies),length(pTypes),2);
NbiasST20 = cell(length(flies),length(pTypes),2);
biasST100 = cell(length(flies),length(pTypes),2);
NbiasST100 = cell(length(flies),length(pTypes),2);

flpv = [0 180];

for tt = 1 : 2
    for n = 1 : length(flies)
        disp(num2str(n))
        dt = load([path flies(n).name '\DataLowRes.mat'], 'Flies');
        seq = dt.Flies.Seq;
        dt = dt.Flies.Data;
        for k = 1 : length(dt)
            for j = 1 : length(pTypes)
                pfb = [];
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
                        Vr = dt{k}.Vr;
                        flp = dt{k}.flp;
                        acSt = dt{k}.actState;
                        % Vr = Vr(flp == flpv(tt));
                        % Vf = Vf(flp == flpv(tt));
                        % acSt = acSt(flp == flpv(tt));
                        
                        
                        [params] = GetParams();
                        if(strcmp(pTypes{j}(end-1:end),'NG'))
                            params.cutoff = 0.15;
                            params.RG = 0;
                        else
                            params.cutoff = 0.24;
                            params.RG = 1;
                        end
                        
                        [locs, ~, cmhSong, thr, dx] = SortSpikes(Vr, acSt, params);
                        [locsF, pksIF, pksEF, ~] = CullSpikesBasedOnVf(Vr, locs, cmhSong, thr, Vf, params);
                        [locsF, pksIF, pksEF, scores] = TemplateCompSpike(Vr,locsF,pksIF, pksEF, cmhSong, thr, params);
                        [FWBoutStr, FWBout] = GetFbouts(acSt, Vf, locsF, pksIF,pksEF, params);
                        
                        FVect = zeros(length(Vr),1);
                        for i = 1 : length(FWBout)
                            FVect(FWBout{i}) = 1;
                            if(flp(FWBout{i}(1)) == flpv(tt))
                                lfb = vertcat(lfb, length(FWBout{i}));
                            end
                        end
                        FVect(flp ~= flpv(tt)) = [];
                        acSt(flp ~= flpv(tt)) = [];
                        
                        vect = zeros(length(Vr),1);
                        for i = 1 : length(locsF)
                            vect((locsF(i)-pksIF(i)):(locsF(i)+pksEF(i))) = 1;
                            if(flp(locsF(i)) == flpv(tt))
                                lst = vertcat(lst, Vr(locsF(i)));
                            end
                        end
                        vect(flp~= flpv(tt)) = [];
                        
                        npss = [];
                        pss = [];
                        for i = 1 : length(FWBout)
                            if(flp(FWBout{i}(1)) == flpv(tt))
                                if length(FWBout{i}) > 10
                                    npss = vertcat(npss,length(FWBout{i}));
                                    pss = vertcat(pss,abs(median(Vr(FWBout{i})))/median(abs(Vr(FWBout{i}))));
                                    valfb = vertcat(valfb, (median(Vr(FWBout{i}))));
                                    nvalfb = vertcat(nvalfb, length(FWBout{i}));
                                end
                            end
                        end
                        
                        tw = 2;
                        w = 6;
                        for i = 1 : length(locsF)
                            if(flp(locsF(i)) == flpv(tt))
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
                        end
                        
                        pfb = vertcat(pfb, sum(abs(FVect))/sum(abs(acSt)));
                        pst = vertcat(pst, sum(abs(vect))/sum(abs(acSt)));
                        pt = vertcat(pt, sum(abs(acSt)));
                        bfb = vertcat(bfb, pss'*npss/sum(npss));
                        nbfb = vertcat(nbfb, sum(npss));
                        
                end
                
                pFB{n,j,tt} = vertcat(pFB{n,j,tt}, pfb);
                pST{n,j,tt} = vertcat(pST{n,j,tt}, pst);
                pT{n,j,tt} = vertcat(pT{n,j,tt}, pt);
                lFB{n,j,tt} = vertcat(lFB{n,j,tt}, lfb);
                lST{n,j,tt} = vertcat(lST{n,j,tt}, lst);
                valFB{n,j,tt} = vertcat(valFB{n,j,tt}, valfb);
                NvalFB{n,j,tt} = vertcat(NvalFB{n,j,tt}, nvalfb);
                biasFB{n,j,tt} = vertcat(biasFB{n,j,tt}, bfb);
                NbiasFB{n,j,tt} = vertcat(NbiasFB{n,j,tt}, nbfb);
                biasST1{n,j,tt} = vertcat(biasST1{n,j,tt}, mean(bst1));
                NbiasST1{n,j,tt} = vertcat(NbiasST1{n,j,tt}, length(bst1));
                biasST20{n,j,tt} = vertcat(biasST20{n,j,tt}, mean(bst20));
                NbiasST20{n,j,tt} = vertcat(NbiasST20{n,j,tt}, length(bst20));
                biasST100{n,j,tt} = vertcat(biasST100{n,j,tt}, mean(bst100));
                NbiasST100{n,j,tt} = vertcat(NbiasST100{n,j,tt}, length(bst100));
                
            end
        end
    end
end
stfb.pST = pST;
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