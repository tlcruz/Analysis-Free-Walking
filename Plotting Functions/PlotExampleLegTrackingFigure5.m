clear
clc 
% get list of videos for a fly
flyID = 'Fly 3539LegTrack';
flyVidID = 'Fly3539CrpVids';
pathVids = ['\' flyVidID];
vids = dir(pathVids);
vids = vids(3:end);

% load low and high resolution data
path = ['\' flyID];
pathi = [path '\DataLowHighRes.mat'];
dt = load(pathi);
dt = dt.Flies;
patho = [path '\DataLowAndHighResHH.mat'];
dto = load(patho);
dto = dto.Fly;

%% Crop
% select the window to plot
i = 3;
delta = 6;
pType = dt.Seq{i};
ap = 79;
ti = ap;
tf = ap+2;
frameI = dt.Data{i}.FramesC2(floor(ti*60));
frameF = dt.Data{i}.FramesC2(floor(tf*60));
vidID = ceil(frameI/17999);
frameI = mod(frameI, 17999);
frameF = mod(frameF, 17999);

% load the video and select the first frame
vRead = VideoReader([pathVids '\' vids(vidID).name]);
vRead.CurrentTime = (frameI-delta)/vRead.FrameRate;
frameI =  readFrame(vRead);
% load the frames
Mat = zeros(vRead.Height, vRead.Width, floor((tf-ti)*120));
for in = 1 : floor((tf-ti)*120)
    Mat(:,:,in) = (readFrame(vRead));
end

% high and low resolution time window
thr = 1;
tH = (1:length(dt.Data{i}.TopHeadX))/120;
tL = (1:length(dt.Data{i}.Vr))/60;

% remove frames with more unreliable tracking
vr = dt.Data{i}.Vr;
vf = -1000 + 10*dt.Data{i}.Vf;
ha = 800+10*dto.Data{i}.Hang;
lfl = dt.Data{i}.LeftFrontLegY;
lfl(dt.Data{i}.LeftFrontLegErr > thr) = nan;
rfl = dt.Data{i}.RightFrontLegY;
rfl(dt.Data{i}.RightFrontLegErr > thr) = nan;
lml = dt.Data{i}.LeftMiddleLegY;
lml(dt.Data{i}.LeftMiddleLegErr > thr) = nan;
rml = dt.Data{i}.RightMiddleLegY;
rml(dt.Data{i}.RightMiddleLegErr > thr) = nan;
lhl = dt.Data{i}.LeftHindLegY;
lhl(dt.Data{i}.LeftHindLegErr > thr) = nan;
rhl = dt.Data{i}.RightHindLegY;
rhl(dt.Data{i}.RightHindLegErr > thr) = nan;
lflx = dt.Data{i}.LeftFrontLegX;
lflx(dt.Data{i}.LeftFrontLegErr > thr) = nan;
lflx2 = lflx - nanmean(lflx) + 100;
rflx = dt.Data{i}.RightFrontLegX;
rflx(dt.Data{i}.RightFrontLegErr > thr) = nan;
rflx2 = rflx - nanmean(rflx) + 100;
lmlx = dt.Data{i}.LeftMiddleLegX;
lmlx(dt.Data{i}.LeftMiddleLegErr > thr) = nan;
lmlx2 = lmlx - nanmean(lmlx) + 200;
rmlx = dt.Data{i}.RightMiddleLegX;
rmlx(dt.Data{i}.RightMiddleLegErr > thr) = nan;
rmlx2 = rmlx - nanmean(rmlx) + 200;
lhlx = dt.Data{i}.LeftHindLegX;
lhlx(dt.Data{i}.LeftHindLegErr > thr) = nan;
lhlx2 = lhlx - nanmean(lhlx) + 300;
rhlx = dt.Data{i}.RightHindLegX;
rhlx(dt.Data{i}.RightHindLegErr > thr) = nan;
rhlx2 = rhlx - nanmean(rhlx) + 300;

% leg color code
aB = 0.2;
cmap = [.2 .5 .7 aB; .8 .3 .2 aB; 1 .5 0 aB; .5 0 .5 aB; 0 .5 0 aB; 0 .8 .8 aB];

% Plot First Frame
tcurr = ti + 1/120;
indIL = find(tL>=tcurr,1);
indIH = find(tH>=tcurr,1);
figure,
set(gcf, 'Position', [0, 0, 1000, 600])
% low resolution data
subplot(4,5,[1 2 3 6 7 8])
hold on
plot(tL, zeros(length(vr),1), '--g', 'linewidth', 2)
plot(tL, vr, 'color', [0 0 1 aB], 'linewidth', 2)
plot(tL(1:indIL), vr(1:indIL), 'color', [0 0 1], 'linewidth', 2)
scatter(tL(indIL), vr(indIL), 80, [0 0 1], 'filled', 'linewidth', 2)
plot(tL, vf, 'color', [0 0 0.5 aB], 'linewidth', 2)
plot(tL(1:indIL), vf(1:indIL), 'color', [0 0 0.5], 'linewidth', 2)
scatter(tL(indIL), vf(indIL), 80, [0 0 0.5], 'filled', 'linewidth', 2)
plot(tH, 2*ha-700, 'color', [1 0.5 0 aB], 'linewidth', 2)
plot(tH(1:indIH), 2*ha(1:indIH)-700, 'color', [1 0.5 0], 'linewidth', 2)
scatter(tH(indIH), 2*ha(indIH)-700, 80, [1 0.5 0], 'filled', 'linewidth', 2)
axis([ti tf -1200 1000])
ylabel('Head and Body Speeds')
set(gca,'Box','off')
% high resolution data Y direction
subplot(4,5,[11 12 13])
hold on
plot(tH, lfl, 'color', cmap(1,:), 'linewidth', 1.5)
plot(tH(1:indIH), lfl(1:indIH), 'color', cmap(1,1:3), 'linewidth', 1.5)
scatter(tH(indIH), lfl(indIH), 80, cmap(1,1:3), 'filled', 'linewidth', 1.5)
plot(tH, rfl, 'color', cmap(2,:), 'linewidth', 1.5)
plot(tH(1:indIH), rfl(1:indIH), 'color', cmap(2,1:3), 'linewidth', 1.5)
scatter(tH(indIH), rfl(indIH), 80, cmap(2,1:3), 'filled', 'linewidth', 1.5)
plot(tH, lml, 'color', cmap(3,:), 'linewidth', 1.5)
plot(tH(1:indIH), lml(1:indIH), 'color', cmap(3,1:3), 'linewidth', 1.5)
scatter(tH(indIH), lml(indIH), 80, cmap(3,1:3), 'filled', 'linewidth', 1.5)
plot(tH, rml, 'color', cmap(4,:), 'linewidth', 1.5)
plot(tH(1:indIH), rml(1:indIH), 'color', cmap(4,1:3), 'linewidth', 1.5)
scatter(tH(indIH), rml(indIH), 80, cmap(4,1:3), 'filled', 'linewidth', 1.5)
plot(tH, lhl, 'color', cmap(5,:), 'linewidth', 1.5)
plot(tH(1:indIH), lhl(1:indIH), 'color', cmap(5,1:3), 'linewidth', 1.5)
scatter(tH(indIH), lhl(indIH), 80, cmap(5,1:3), 'filled', 'linewidth', 1.5)
plot(tH, rhl, 'color', cmap(6,:), 'linewidth', 1.5)
plot(tH(1:indIH), rhl(1:indIH), 'color', cmap(6,1:3), 'linewidth', 1.5)
scatter(tH(indIH), rhl(indIH), 80, cmap(6,1:3), 'filled', 'linewidth', 1.5)
axis([ti tf 0 360])
ylabel('Leg Y Position')
set(gca,'Box','off')
% high resolution data X direction
subplot(4,5,[16 17 18])
hold on
plot(tH, lflx2, 'color', cmap(1,:), 'linewidth', 1.5)
plot(tH(1:indIH), lflx2(1:indIH), 'color', cmap(1,1:3), 'linewidth', 1.5)
scatter(tH(indIH), lflx2(indIH), 80, cmap(1,1:3), 'filled', 'linewidth', 1.5)
plot(tH, rflx2, 'color', cmap(2,:), 'linewidth', 1.5)
plot(tH(1:indIH), rflx2(1:indIH), 'color', cmap(2,1:3), 'linewidth', 1.5)
scatter(tH(indIH), rflx2(indIH), 80, cmap(2,1:3), 'filled', 'linewidth', 1.5)
plot(tH, lmlx2, 'color', cmap(3,:), 'linewidth', 1.5)
plot(tH(1:indIH), lmlx2(1:indIH), 'color', cmap(3,1:3), 'linewidth', 1.5)
scatter(tH(indIH), lmlx2(indIH), 80, cmap(3,1:3), 'filled', 'linewidth', 1.5)
plot(tH, rmlx2, 'color', cmap(4,:), 'linewidth', 1.5)
plot(tH(1:indIH), rmlx2(1:indIH), 'color', cmap(4,1:3), 'linewidth', 1.5)
scatter(tH(indIH), rmlx2(indIH), 80, cmap(4,1:3), 'filled', 'linewidth', 1.5)
plot(tH, lhlx2, 'color', cmap(5,:), 'linewidth', 1.5)
plot(tH(1:indIH), lhlx2(1:indIH), 'color', cmap(5,1:3), 'linewidth', 1.5)
scatter(tH(indIH), lhlx2(indIH), 80, cmap(5,1:3), 'filled', 'linewidth', 1.5)
plot(tH, rhlx2, 'color', cmap(6,:), 'linewidth', 1.5)
plot(tH(1:indIH), rhlx2(1:indIH), 'color', cmap(6,1:3), 'linewidth', 1.5)
scatter(tH(indIH), rhlx2(indIH), 80, cmap(6,1:3), 'filled', 'linewidth', 1.5)
axis([ti tf 0 360])
ylabel('Leg X Position')
set(gca,'Box','off')
subplot(4,5,[4 5 9 10 14 15 19 20])
imagesc(Mat(:,:,1))
caxis([0 255])
hold on
colormap gray
scatter(lflx(indIH), lfl(indIH), 100, cmap(1,1:3), 'filled', 'linewidth', 1.5)
scatter(rflx(indIH), rfl(indIH), 100, cmap(2,1:3), 'filled', 'linewidth', 1.5)
scatter(lmlx(indIH), lml(indIH), 100, cmap(3,1:3), 'filled', 'linewidth', 1.5)
scatter(rmlx(indIH), rml(indIH), 100, cmap(4,1:3), 'filled', 'linewidth', 1.5)
scatter(lhlx(indIH), lhl(indIH), 100, cmap(5,1:3), 'filled', 'linewidth', 1.5)
scatter(rhlx(indIH), rhl(indIH), 100, cmap(6,1:3), 'filled', 'linewidth', 1.5)
box off
axis off
axis square


%% Make Video
v2 = VideoWriter('C:\Users\tomas\Desktop\newfile.avi','Motion JPEG AVI');
v2.FrameRate = 12;
open(v2);
% iterate across all frames to make a video
for fn = 1: floor((tf-ti)*120)
    tcurr = ti + fn/120;
    indIL = find(tL>=tcurr,1);
%     indIH = find(tH>=tcurr,1);
    indIH = find(tH>=ti,1) + fn;
    figure,
    set(gcf, 'Position', [0, 0, 1000, 600])
    subplot(4,5,[1 2 3 6 7 8])
    hold on
    plot(tL, zeros(length(vr),1), '--g', 'linewidth', 2)
    plot(tL, vr, 'color', [0 0 1 aB], 'linewidth', 2)
    plot(tL(1:indIL), vr(1:indIL), 'color', [0 0 1], 'linewidth', 2)
    scatter(tL(indIL), vr(indIL), 80, [0 0 1], 'filled', 'linewidth', 2)
    plot(tL, vf, 'color', [0 0 0.5 aB], 'linewidth', 2)
    plot(tL(1:indIL), vf(1:indIL), 'color', [0 0 0.5], 'linewidth', 2)
    scatter(tL(indIL), vf(indIL), 80, [0 0 0.5], 'filled', 'linewidth', 2)
    plot(tH, 2*ha-700, 'color', [1 0.5 0 aB], 'linewidth', 2)
    plot(tH(1:indIH), 2*ha(1:indIH)-700, 'color', [1 0.5 0], 'linewidth', 2)
    scatter(tH(indIH), 2*ha(indIH)-700, 80, [1 0.5 0], 'filled', 'linewidth', 2)
    axis([ti tf -1200 1000])
    ylabel('Head and Body Speeds')
    set(gca,'Box','off')
    subplot(4,5,[11 12 13])
    hold on
    plot(tH, lfl, 'color', cmap(1,:), 'linewidth', 1.5)
    plot(tH(1:indIH), lfl(1:indIH), 'color', cmap(1,1:3), 'linewidth', 1.5)
    scatter(tH(indIH), lfl(indIH), 80, cmap(1,1:3), 'filled', 'linewidth', 1.5)
    plot(tH, rfl, 'color', cmap(2,:), 'linewidth', 1.5)
    plot(tH(1:indIH), rfl(1:indIH), 'color', cmap(2,1:3), 'linewidth', 1.5)
    scatter(tH(indIH), rfl(indIH), 80, cmap(2,1:3), 'filled', 'linewidth', 1.5)
    plot(tH, lml, 'color', cmap(3,:), 'linewidth', 1.5)
    plot(tH(1:indIH), lml(1:indIH), 'color', cmap(3,1:3), 'linewidth', 1.5)
    scatter(tH(indIH), lml(indIH), 80, cmap(3,1:3), 'filled', 'linewidth', 1.5)
    plot(tH, rml, 'color', cmap(4,:), 'linewidth', 1.5)
    plot(tH(1:indIH), rml(1:indIH), 'color', cmap(4,1:3), 'linewidth', 1.5)
    scatter(tH(indIH), rml(indIH), 80, cmap(4,1:3), 'filled', 'linewidth', 1.5)
    plot(tH, lhl, 'color', cmap(5,:), 'linewidth', 1.5)
    plot(tH(1:indIH), lhl(1:indIH), 'color', cmap(5,1:3), 'linewidth', 1.5)
    scatter(tH(indIH), lhl(indIH), 80, cmap(5,1:3), 'filled', 'linewidth', 1.5)
    plot(tH, rhl, 'color', cmap(6,:), 'linewidth', 1.5)
    plot(tH(1:indIH), rhl(1:indIH), 'color', cmap(6,1:3), 'linewidth', 1.5)
    scatter(tH(indIH), rhl(indIH), 80, cmap(6,1:3), 'filled', 'linewidth', 1.5)
    axis([ti tf 0 360])
    ylabel('Leg Y Position')
    set(gca,'Box','off')
    subplot(4,5,[16 17 18])
    hold on
    plot(tH, lflx2, 'color', cmap(1,:), 'linewidth', 1.5)
    plot(tH(1:indIH), lflx2(1:indIH), 'color', cmap(1,1:3), 'linewidth', 1.5)
    scatter(tH(indIH), lflx2(indIH), 80, cmap(1,1:3), 'filled', 'linewidth', 1.5)
    plot(tH, rflx2, 'color', cmap(2,:), 'linewidth', 1.5)
    plot(tH(1:indIH), rflx2(1:indIH), 'color', cmap(2,1:3), 'linewidth', 1.5)
    scatter(tH(indIH), rflx2(indIH), 80, cmap(2,1:3), 'filled', 'linewidth', 1.5)
    plot(tH, lmlx2, 'color', cmap(3,:), 'linewidth', 1.5)
    plot(tH(1:indIH), lmlx2(1:indIH), 'color', cmap(3,1:3), 'linewidth', 1.5)
    scatter(tH(indIH), lmlx2(indIH), 80, cmap(3,1:3), 'filled', 'linewidth', 1.5)
    plot(tH, rmlx2, 'color', cmap(4,:), 'linewidth', 1.5)
    plot(tH(1:indIH), rmlx2(1:indIH), 'color', cmap(4,1:3), 'linewidth', 1.5)
    scatter(tH(indIH), rmlx2(indIH), 80, cmap(4,1:3), 'filled', 'linewidth', 1.5)
    plot(tH, lhlx2, 'color', cmap(5,:), 'linewidth', 1.5)
    plot(tH(1:indIH), lhlx2(1:indIH), 'color', cmap(5,1:3), 'linewidth', 1.5)
    scatter(tH(indIH), lhlx2(indIH), 80, cmap(5,1:3), 'filled', 'linewidth', 1.5)
    plot(tH, rhlx2, 'color', cmap(6,:), 'linewidth', 1.5)
    plot(tH(1:indIH), rhlx2(1:indIH), 'color', cmap(6,1:3), 'linewidth', 1.5)
    scatter(tH(indIH), rhlx2(indIH), 80, cmap(6,1:3), 'filled', 'linewidth', 1.5)
    axis([ti tf 0 360])
    ylabel('Leg X Position')
    set(gca,'Box','off')
    subplot(4,5,[4 5 9 10 14 15 19 20])
    imagesc(Mat(:,:,fn))
    hold on
    colormap gray
    axis square
    scatter(lflx(indIH), lfl(indIH), 100, cmap(1,1:3), 'filled', 'linewidth', 1.5)
    scatter(rflx(indIH), rfl(indIH), 100, cmap(2,1:3), 'filled', 'linewidth', 1.5)
    scatter(lmlx(indIH), lml(indIH), 100, cmap(3,1:3), 'filled', 'linewidth', 1.5)
    scatter(rmlx(indIH), rml(indIH), 100, cmap(4,1:3), 'filled', 'linewidth', 1.5)
    scatter(lhlx(indIH), lhl(indIH), 100, cmap(5,1:3), 'filled', 'linewidth', 1.5)
    scatter(rhlx(indIH), rhl(indIH), 100, cmap(6,1:3), 'filled', 'linewidth', 1.5)
    box off
    axis off
    set(gcf, 'Visible', 'off')
    F = getframe(gcf);
    writeVideo(v2,F);
    disp(['Frame ' num2str(fn) ' out of ' num2str(floor((tf-ti)*120))])
    clf
end
close(v2)
