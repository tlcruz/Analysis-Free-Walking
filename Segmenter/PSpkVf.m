path = 'C:\Users\tomas\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\WT TB\';
flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowRes.mat'], 'Flies');
seq = dt.Flies.Seq;
pTypes = unique(seq);

vFcents = -12:5:43;

PSSKF = cell(length(pTypes), length(flies), length(vFcents)-1);
NPKF = cell(length(pTypes), length(flies), length(vFcents)-1);
delta = 15;
params = GetParams();
for n = 1 : length(flies)
    dt = load([path flies(n).name '\DataLowRes.mat'], 'Flies');
    seq = dt.Flies.Seq;
    dt = dt.Flies.Data;
    for k = 1 : length(dt)
        for j = 1 : length(pTypes)
            switch seq{k}
                case pTypes{j}
                    pssk = cell(length(vFcents)-1, 1);
                    npssk = cell(length(vFcents)-1, 1);
                    vr = dt{k}.Vr;
                    vf = dt{k}.Vf;
                    acSt = dt{k}.actState;
                    for i = 1 : length(dt{k}.Bouts)
                        if length(dt{k}.Bouts{i}) > 100
                              vrb = vr(dt{k}.Bouts{i});
                              vfb = vf(dt{k}.Bouts{i});
                              acStb = acSt(dt{k}.Bouts{i});
%                             vrb = vr((dt{k}.Bouts{i}(1)-10):(dt{k}.Bouts{i}(1)+10));
%                             vfb = vf((dt{k}.Bouts{i}(1)-10):(dt{k}.Bouts{i}(1)+10));
%                             acStb = acSt((dt{k}.Bouts{i}(1)-10):(dt{k}.Bouts{i}(1)+10));
%                             
                            [locs, ~, cmhSong, thr] = SortSpikes(vrb, acStb, params);
                            [locsF, pksIF,pksEF, ~] = CullSpikesBasedOnVf(vrb, locs, cmhSong, thr, vfb, params);
                            [locsF, pksIF, pksEF] = TemplateCompSpike(vrb,locsF,pksIF, pksEF, cmhSong, thr, params);
                            
                            spkk = zeros(length(vrb), 1);
                            for kk = 1 : length(locsF)
                                spkk((locsF(kk)-pksIF(kk)):(locsF(kk)+pksEF(kk))) = 1;
                            end
                            for kk = 1 : length(vrb)
                                for kj = 1 : length(vFcents)-1
                                    if vfb(kk) > vFcents(kj) && vfb(kk) < vFcents(kj+1)
                                        pssk{kj} = vertcat(pssk{kj}, spkk(kk));
                                        npssk{kj} = vertcat(npssk{kj}, 1);
                                    end
                                end
                            end
                        end
                    end
                    for kj = 1 : length(vFcents)-1
                        if ~isempty(PSSKF{j, n, kj})
                            if~isempty(pssk{kj})
                                PSSKF{j, n, kj} = (NPKF{j, n, kj}*PSSKF{j, n, kj} ...
                                    + nanmean(pssk{kj})*nansum(npssk{kj}))/(nansum(npssk{kj})+NPKF{j, n, kj});
                                NPKF{j, n, kj} = NPKF{j, n, kj} + nansum(npssk{kj});
                            end
                        else
                            if~isempty(pssk{kj})
                                PSSKF{j, n, kj} = nanmean(pssk{kj});
                                NPKF{j, n, kj} = nansum(npssk{kj});
                            else
                                PSSKF{j, n, kj} = 0;
                                NPKF{j, n, kj} = 0;
                            end
                        end
                    end
            end
            
            clearvars pss npss
        end
    end
end


%%

figure,
cmap = autumn(5);
for j = 1 : length(pTypes)/2
    subplot(1,5,j)
    plot([-1 50], [0 0], '--g', 'linewidth', 2)
    hold on
    VI = [];
    for kj = 1 : length(vFcents)-1
        vrg = [];
        nvrg = [];
        for n = 1 : length(flies)
            vrg = vertcat(vrg, PSSKF{2*j, n, kj});
            nvrg = vertcat(nvrg, NPKF{2*j, n, kj});
        end
        vmrg = vrg'*nvrg/(sum(nvrg)+1);
        semvmrg = sqrt(((vrg-vmrg).*(vrg-vmrg))'*nvrg/sum(nvrg))/sqrt(length(flies));
        
        vng = [];
        nvng = [];
        for n = 1 : length(flies)
            vng = vertcat(vng, PSSKF{2*j-1, n, kj});
            nvng = vertcat(nvng, NPKF{2*j-1, n, kj});
        end
        vmng = vng'*nvng/(sum(nvng)+1);
        semvmng = sqrt(((vng-vmng).*(vng-vmng))'*nvng/sum(nvng))/sqrt(length(flies));
        
        vii = (vmrg-vmng)/(vmrg+vmng);
        semvii = sqrt(((2*vmng/((vmrg+vmng)^2))^2)*semvmrg^2 + ((2*vmrg/((vmrg+vmng)^2))^2)*semvmng^2);
        
        errorbar(vFcents(kj), vmng, semvmng, 'o', 'color', cmap(j,:))
        VI = vertcat(VI, vii);
    end
%     plot(vFcents(2:9), VI(2:9), 'color', cmap(j,:), 'linewidth', 2)
%     axis([-5 31 -0.1 0.4])
    title(pTypes{2*j-1})
end

%%

figure,
cmap = autumn(5);

%     subplot(1,5,j)
plot([-1 50], [0 0], '--g', 'linewidth', 2)
hold on
VI = [];
for kj = 1 : length(vFcents)-1
    vng = [];
    nvng = [];
    for j = 1 : length(pTypes)/2
        for n = 1 : length(flies)
            vng = vertcat(vng, PSSKF{2*j-1, n, kj});
            nvng = vertcat(nvng, NPKF{2*j-1, n, kj});
        end
    end
    vmng = vng'*nvng/(sum(nvng)+1);
    semvmng = sqrt(((vng-vmng).*(vng-vmng))'*nvng/sum(nvng))/sqrt(length(flies));
    
    errorbar(vFcents(kj), vmng, semvmng, 'o', 'color', 'k')
end
axis([-15 40 0 1])
xlabel('Forward Speed (mm/s)')
ylabel('Probability Spike Turn')








