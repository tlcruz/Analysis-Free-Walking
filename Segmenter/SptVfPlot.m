clear
clc
[params] = GetParams();
path = 'C:\Users\tomas\Dropbox (Sensorimotor)\ChiappeLabNew\data\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\WT TB\';
flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowRes.mat'], 'Flies');
seq = dt.Flies.Seq;
pTypes = unique(seq);
disp(path)

vFALL = cell(length(pTypes),1);
STRALL = cell(length(pTypes),1);
NSTRALL = cell(length(pTypes),1);
for n = 1 : length(flies)
    dt = load([path flies(n).name '\DataLowRes.mat'], 'Flies');
    seq = dt.Flies.Seq;
    dt = dt.Flies.Data;
    for k = 1 : length(dt)
        for j = 1 : length(pTypes)
            str = [];
            vt = [];
            nstr = [];
            switch seq{k}
                case pTypes{j}
                    vrb = dt{k}.Vr;
                    vfb = dt{k}.Vf;
                    vtb = dt{k}.Vt;
                    xtb = dt{k}.X;
                    ytb = dt{k}.Y;
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
%                             dst = vertcat(dst, sum(ss));
                            nstr = vertcat(nstr, length(xaux));
                            vt = vertcat(vt, mean(vtaux));
%                             angd = vertcat(angd, sum(abs(vraux))/sum(vtaux));
                        end
                    end
                    
                    vFALL{j} = vertcat(vFALL{j}, vt);
                    STRALL{j} = vertcat(STRALL{j}, str);
                    NSTRALL{j} = vertcat(NSTRALL{j}, nstr);
            end
        end
    end
    disp([num2str(floor(100*n/length(flies))) '% Done' ])
end

%%

vfCents = -5:0.1:50;
cmap = autumn(5);
ds = [5 7 9 1 3];
figure,
hold on
vf = [];
str = [];
for i = 1 : 5
    vf = vertcat(vf, vFALL{ds(i)});
    str = vertcat(str, STRALL{ds(i)});
%     scatter(vFALL{ds(i)}, STRALL{ds(i)}, 100, cmap(i,:), 'filled')
    plot(vfCents, cumsum(smooth(hist(vFALL{ds(i)}, vfCents)/sum(hist(vFALL{ds(i)}, vfCents)))), ...
        'color', cmap(i,:))
    plot([mean(vFALL{2*i-1}) mean(vFALL{2*i-1})], [0 1],'color', cmap(i,:))
    plot([-5 50], [0.2 0.2], 'b')
    plot([-5 50], [0.8 0.8], 'b')
    
    plot([-5 50], [0.1 0.1], 'g')
    plot([-5 50], [0.9 0.9], 'g')
    axis([-5 50 0 1])
end
%%
clc

figure,
minb = 15;
maxb = 27;
v = 120;
hold on
cmap = autumn(4);
for dsize = 1:4;

    vf = vFALL{ds(dsize)};
    str = STRALL{ds(dsize)};
    nstr = NSTRALL{ds(dsize)};
    indsL = find(vf<minb & str<v);
    lv = str(indsL)'*nstr(indsL)/sum(nstr(indsL));
    slv = sqrt(((str(indsL)-lv).*(str(indsL)-lv))'*nstr(indsL)/sum(nstr(indsL)))/sqrt(33);
    indsM = find(vf>minb & vf<maxb & str<v);
    mv = str(indsM)'*nstr(indsM)/sum(nstr(indsM));
    indsH = find(vf>maxb & str<v);
    hv = str(indsH)'*nstr(indsH)/sum(nstr(indsH));
    shv = sqrt(((str(indsH)-hv).*(str(indsH)-hv))'*nstr(indsH)/sum(nstr(indsH)))/sqrt(33);
    disp([num2str(dsize) ':  ' num2str(lv) '+/-' num2str(slv) '    ' num2str(hv) '+/- ' num2str(shv)])
    
    errorbar([1 2], [lv hv], [slv shv], '-o', 'color', cmap(dsize,:))
end
axis([0.5 2.5 30 65])

%%

clear
clc
[params] = GetParams();
path = 'C:\Users\tomas\Dropbox (Sensorimotor)\ChiappeLabNew\data\TOMÁS\Free Walking VR\Data\Bilateral Stimulus Size\WT TB\';
flies = dir(path);
flies = flies(3:end);
dt = load([path flies(1).name '\DataLowRes.mat'], 'Flies');
seq = dt.Flies.Seq;
pTypes = unique(seq);
disp(path)

vFSPALL = cell(length(pTypes),length(flies));
vRSPALL = cell(length(pTypes),length(flies));
NSPALL = cell(length(pTypes),length(flies));

vFALL = cell(length(pTypes),length(flies));
vRALL = cell(length(pTypes),length(flies));
NALL = cell(length(pTypes),length(flies));
for n = 1 : length(flies)
    dt = load([path flies(n).name '\DataLowRes.mat'], 'Flies');
    seq = dt.Flies.Seq;
    dt = dt.Flies.Data;
    for k = 1 : length(dt)
        for j = 1 : length(pTypes)
            str = [];
            vt = [];
            nstr = [];
            switch seq{k}
                case pTypes{j}
                    vrb = dt{k}.Vr;
                    vfb = dt{k}.Vf;
                    vtb = dt{k}.Vt;
                    xtb = dt{k}.X;
                    ytb = dt{k}.Y;
                    wd = dt{k}.WallDist;
                    actst = dt{k}.actState;
                    % Get spike times and forward bouts
                    actst(wd < params.mDistWall) = 0;
                    
                    [locs, ~, cmhSong, thr] = SortSpikes(vrb, actst, params);
                    [locsF, pksIF,pksEF] = CullSpikesBasedOnVf(vrb, locs, cmhSong, thr, vfb, params);
                    [locsF, pksIF, pksEF] = TemplateCompSpike(vrb,locsF,pksIF, pksEF, cmhSong, thr, params);
                    [FBouts] = GetFbouts(actst, vfb, locsF, pksIF,pksEF, params);
                    
                    for pp = 1 : length(FBouts)
                        vfsp = vfb(FBouts{pp});
                        vrsp = vrb(FBouts{pp});
                        if length(FBouts{pp}) > 11
                            for i = 1: floor((length(FBouts{pp})-10)/10)
                                indi = (5+(10*(i-1)));
                                indf = (10*i+5);
                                vFALL{j,n} = vertcat(vFALL{j,n}, min(vfsp(indi:indf)));
                                vRALL{j,n} = vertcat(vRALL{j,n}, vrsp(indf-5));
                            end
                        end
                    end
                    
%                     for pp = 1 : length(FBouts)
%                         % Calculate Straightness
%                         if length(FBouts{pp}) > max(params.windStr+1, params.minStrB);
%                             xaux = xtb(FBouts{pp});
%                             yaux = ytb(FBouts{pp});
%                             vtaux = vtb(FBouts{pp});
%                             vraux = vrb(FBouts{pp});
%                             nv = length(xaux) - params.windStr;
%                             dm = [];
%                             ss = [];
%                             for l = 1 : nv
%                                 pi = [xaux(l),yaux(l),0];
%                                 pf = [xaux(l+params.windStr),yaux(l+params.windStr),0];
%                                 dis = sum(vtaux((l):(l+params.windStr)))/60;
%                                 pt = [xaux(l+floor(params.windStr/2)),...
%                                     yaux(l+floor(params.windStr/2)),0];
%                                 dr = point_to_line(pt,pi,pf);
%                                 dm = vertcat(dm, dr);
%                                 ss = vertcat(ss, dis);
%                             end
%                             str = vertcat(str, sum(ss)/sum(dm));
% %                             dst = vertcat(dst, sum(ss));
%                             nstr = vertcat(nstr, length(xaux));
%                             vt = vertcat(vt, mean(vtaux));
% %                             angd = vertcat(angd, sum(abs(vraux))/sum(vtaux));
%                         end
%                     end
                    for i = 1 : length(locsF)
                        vfsp = vfb((locsF(i)-pksIF(i)):(locsF(i)+pksIF(i)));
                        vrsp = vrb((locsF(i)-pksIF(i)):(locsF(i)+pksIF(i)));
                        
                        vFSPALL{j,n} = vertcat(vFSPALL{j,n}, min(vfsp));
                        vRSPALL{j,n} = vertcat(vRSPALL{j,n}, vrb(locsF(i)));
                    end

%                     NSPALL{j} = vertcat(NSPALL{j}, nstr);
            end
        end
    end
    disp([num2str(floor(100*n/length(flies))) '% Done' ])
end

%%

figure,
for i = 1 : 5
   subplot(1,5,i)
   scatter(vRALL{2*i-1},vFALL{2*i-1}, 20, '.k')
   hold on
   scatter(vRSPALL{2*i-1},vFSPALL{2*i-1}, 20, '.r')
   title(pTypes{2*i-1})
   axis([-1200 1200 -20 40])
end
%%
mvfsp = cell(5,1);
nsp = cell(5,1);
mvf = cell(5,1);
semvfsp = cell(5,1);
semvf = cell(5,1);
nf = cell(5,1);
figure,
for n = 1 : length(flies)
    for i = 1 : 5
        subplot(1,5,i)
        hold on
        vfsp = vFSPALL{2*i-1,n};
        vf = vFALL{2*i-1,n};
        mvfsp{i} = vertcat(mvfsp{i}, mean(vfsp));
        semvfsp{i} = vertcat(semvfsp{i}, std(vfsp));
        nsp{i} = vertcat(nsp{i}, length(vfsp));
        
        mvf{i} = vertcat(mvf{i}, mean(vf));
        semvf{i} = vertcat(semvf{i}, std(vf));
        nf{i} = vertcat(nf{i}, length(vf));
        
        
        errorbar([1 2], [mean(vfsp) mean(vf)], [std(vfsp) std(vf)], '-o', 'color', [0.85 0.85 0.85])
        title([pTypes{2*i-1}])
        axis([0.5 2.5 -5 30])
    end
end

for i = 1 : 5
    subplot(1,5,i)
    hold on
    gmsp = mvfsp{i}' * nsp{i} / sum(nsp{i});
    semsp = sqrt(((mvfsp{i}-gmsp).*(mvfsp{i}-gmsp))'* nsp{i} / sum(nsp{i}))/sqrt(33);
    
    gmf = mvf{i}' * nf{i} / sum(nf{i});
    semf = sqrt(((mvf{i}-gmf).*(mvf{i}-gmf))'* nf{i} / sum(nf{i}))/sqrt(33);
    errorbar([1 2], [gmsp gmf], [semsp semf], '-o', 'color', [0 0 0], 'linewidth', 2)
    title([pTypes{2*i-1}])
    axis([0.5 2.5 -5 30])
end























            
            
            
            

