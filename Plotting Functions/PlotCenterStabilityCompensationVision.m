%% Plot center of stability displacement
clear
clc
% load leg parameters during forward runs
path = '\';
params = GetParams();
[LegMovForwSeg, pTypes] = GetForwSegLegData(path, params);

thr = 0.02;
DeltaFLD1 = cell(size(LegMovForwSeg, 1),size(LegMovForwSeg, 2));
DeltaFLD2 = cell(size(LegMovForwSeg, 1),size(LegMovForwSeg, 2));
DeltaSC1 = cell(size(LegMovForwSeg, 1),size(LegMovForwSeg, 2));
DeltaSC2 = cell(size(LegMovForwSeg, 1),size(LegMovForwSeg, 2));
for j = 1 : 2 : size(LegMovForwSeg, 1)
    for n = 1 : size(LegMovForwSeg, 2)
        % clean the leg position traces
        [lD] = GetCleanLegData(LegMovForwSeg{j,n}, thr);
        for i = 1 : length(lD.VR)
            % get the parameters of the triangle in stance at each step
            [TrSt, TrSw, indsT1, indsT2] = GetStepTriangle(lD.FLLY{i}, lD.FRLY{i}, lD.MLLY{i}, lD.MRLY{i}, lD.HLLY{i}, ...
                lD.HRLY{i}, lD.FLLX{i}, lD.FRLX{i}, lD.MLLX{i}, lD.MRLX{i}, lD.HLLX{i}, lD.HRLX{i}, lD.ThX{i}, lD.ThY{i});
            for t = 1: size(TrSt.Tr1X,2)
                if ~isempty(TrSt.Tr1X) && ~isempty(TrSw.Tr1X) && ~isempty(TrSt.Tr2X) && ~isempty(TrSw.Tr2X) &&...
                        TrSt.Tr1X(1,t)~=0 && TrSw.Tr1X(1,t)~=0
                    FLX = TrSt.Tr1X(1,t);
                    FLDX = TrSt.Tr1X(1,t)-TrSw.Tr1X(1,t);
                    SCX = TrSt.cTr1X(t);
                    SCDX = TrSt.cTr1X(t)-TrSw.cTr1X(t);
                    indsT2b = indsT2;
                    indNSt = find((indsT2b-indsT1(t)) <= 5 & (indsT2b-indsT1(t))>0);
                    
                    if ~isempty(indNSt)
                        if TrSt.Tr2X(1,indNSt)~=0 && TrSw.Tr2X(1,indNSt)~=0
                            FL2X =  TrSt.Tr2X(1,indNSt);
                            FLD2X =  TrSt.Tr2X(1,indNSt)-TrSw.Tr2X(1,indNSt);
                            SC2X = TrSt.cTr2X(indNSt);
                            SCD2X = TrSt.cTr2X(indNSt)-TrSw.cTr2X(indNSt);
                            DeltaFLD1{j,n} = vertcat(DeltaFLD1{j,n}, FLDX);
                            DeltaFLD2{j,n} = vertcat(DeltaFLD2{j,n}, FLD2X);
                            DeltaSC1{j,n} = vertcat(DeltaSC1{j,n}, SCX);
                            DeltaSC2{j,n} = vertcat(DeltaSC2{j,n}, SC2X);
                        end
                    end
                end
            end
            
            for t = 1: size(TrSt.Tr2X,2)
                if ~isempty(TrSt.Tr2X) && ~isempty(TrSw.Tr2X) && ~isempty(TrSt.Tr1X) && ~isempty(TrSw.Tr1X) &&...
                        TrSt.Tr2X(1,t)~=0 && TrSw.Tr2X(1,t)~=0
                    FLX = TrSt.Tr2X(1,t);
                    FLDX = TrSt.Tr2X(1,t)-TrSw.Tr2X(1,t);
                    SCX = TrSt.cTr2X(t);
                    SCDX = TrSt.cTr2X(t)-TrSw.cTr2X(t);
                    indsT1b = indsT1;
                    indNSt = find((indsT1b-indsT2(t)) <= 5 & (indsT1b-indsT2(t))>0);
                    
                    if ~isempty(indNSt)
                        if TrSt.Tr1X(1,indNSt)~=0 && TrSw.Tr1X(1,indNSt)~=0
                            FL2X =  TrSt.Tr1X(1,indNSt);
                            FLD2X =  TrSt.Tr1X(1,indNSt)-TrSw.Tr1X(1,indNSt);
                            SC2X = TrSt.cTr1X(indNSt);
                            SCD2X = TrSt.cTr1X(indNSt)-TrSw.cTr1X(indNSt);
                            DeltaFLD1{j,n} = vertcat(DeltaFLD1{j,n}, FLDX);
                            DeltaFLD2{j,n} = vertcat(DeltaFLD2{j,n}, FLD2X);
                            DeltaSC1{j,n} = vertcat(DeltaSC1{j,n}, SCX);
                            DeltaSC2{j,n} = vertcat(DeltaSC2{j,n}, SC2X);
                        end
                    end
                end
            end
        end
    end
end



figure,
subplot(1,2,1)
hold on
plot([-30 30], [-30 30], '--k', 'color', [.8 .8 .8])
sc1 = [];
sc2 = [];
title(pTypes{1})
DeltaSC = cell(size(DeltaFLD1));
for j = 1: size(DeltaFLD1,1)
    for n = 1 : size(DeltaFLD1,2)
        if ~isempty(DeltaFLD1{j,n})
            FL1 = DeltaFLD1{j,n};
            FL2 = DeltaFLD2{j,n};
            SC1 = DeltaSC1{j,n};
            SC2 = DeltaSC2{j,n};
            % given an initial displacement
            indR = find(FL1>5 & SC1>5);
            indB = find(FL1<-5 & SC1<-5);
            sc1 = vertcat(sc1, abs(SC1(indR)));
            sc2 = vertcat(sc2, abs(SC2(indR)));
            sc1 = vertcat(sc1, abs(SC1(indB)));
            sc2 = vertcat(sc2, abs(SC2(indB)));
            % get the percentage compensation
            DeltaSC{j,n} = vertcat(DeltaSC{j,n}, (SC1(indR)-SC2(indR))./abs(SC1(indR)));
            DeltaSC{j,n} = vertcat(DeltaSC{j,n}, -(SC1(indB)-SC2(indB))./abs(SC1(indB)));

            subplot(1,2,1)
            hold on
            scatter(SC1(indR), SC2(indR), 100,'r', 'filled')
            scatter(SC1(indB), SC2(indB), 100,'b', 'filled')

            subplot(1,2,2)
            for i = 1 : length(indR)
                hold on
                plot([1 2], abs([SC1(indR(i)) SC2(indR(i))]), 'color', [0.5 0.5 0.5])
                scatter(1, abs(SC1(indR(i))), 70, [0.5 0.5 0.5], 'filled')
                scatter(2, abs(SC2(indR(i))), 70, [0.8 0.8 0.8], 'filled')
                axis([0 3 0 30])
            end
            for i = 1 : length(indB)
                hold on
                plot([1 2], abs(-[SC1(indB(i)) SC2(indB(i))]))
                scatter(1, abs(-SC1(indB(i))), 70, [0.5 0.5 0.5], 'filled')
                scatter(2, abs(-SC2(indB(i))), 70, [0.8 0.8 0.8], 'filled')
                axis([0 3 0 30])
            end
        end
    end
end
bar([1 2], [mean(sc1) mean(sc2)])
errorbar([1 2], [mean(sc1) mean(sc2)], [std(sc1) std(sc2)]/sqrt(25), 'ok')
axis([0 3 0 30])
title(num2str(signrank(sc1,sc2)))
set(gcf,'renderer','Painters')

% plot average percentage compensation for every protocol
GM = [];
SEM = [];
for j = 1 : 2: size(DeltaSC,1)
    SCaux = [];
    nSCaux = [];
    for n = 1 : size(DeltaSC,2)
        if ~isempty(DeltaSC{j,n})
            SCaux = vertcat(SCaux, nanmean(DeltaSC{j,n}));
            nSCaux = vertcat(nSCaux, length(DeltaSC{j,n}));
        end
    end
    SCaux = SCaux(nSCaux>2);
    nSCaux = nSCaux(nSCaux>2);
    GM = vertcat(GM, SCaux'*nSCaux/sum(nSCaux));
    SEM = vertcat(SEM, sqrt(((SCaux-SCaux'*nSCaux/sum(nSCaux)).*(SCaux-SCaux'*nSCaux/sum(nSCaux)))'*nSCaux/sum(nSCaux))/sqrt(length(SCaux)));
    
    disp([pTypes{j} ' ' num2str(gm) '+/-' num2str(sem)])
    
end
inds = [2 3 4 1];
figure,
hold on
errorbar(1:ceil(size(DeltaSC,1)/2), GM(inds), SEM(inds), 'ok')
plot(1:4, nanmean(DeltaSCBS)*ones(4,1), 'g')
axis([0.5 4.5 0 0.6])

