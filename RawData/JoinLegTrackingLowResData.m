clear
clc
% Load sturcture raw leg tracking and low resolution tracking data
path = '\Fly 060';
pathLRfile = [path '\DataLowAndHighRes.mat'];
pathHRfile = [path '\StructTrackData.mat'];
fLR = load(pathLRfile);
fHR = load(pathHRfile);
fLR.Flies = fLR.Fly;
% iterate for all the trials
for i = 1 : length(fLR.Flies.Data)
    % identify the high resolution frames that correspond to the trial
    iniVal = fLR.Flies.Data{i}.FramesC2(1);
    endVal = min(length(fHR.TrackData.frameNo), fLR.Flies.Data{i}.FramesC2(end));
    % store leg tracking data for that trial 
    fLR.Flies.Data{i}.TopHeadX = fHR.TrackData.topHeadX(iniVal:endVal);
    fLR.Flies.Data{i}.TopHeadY = fHR.TrackData.topHeadY(iniVal:endVal);
    fLR.Flies.Data{i}.TopHeadErr = 1-fHR.TrackData.topHeadLikelihood(iniVal:endVal);
    fLR.Flies.Data{i}.NeckX = fHR.TrackData.NeckX(iniVal:endVal);
    fLR.Flies.Data{i}.NeckY = fHR.TrackData.NeckY(iniVal:endVal);
    fLR.Flies.Data{i}.NeckErr = 1-fHR.TrackData.NeckLikelihood(iniVal:endVal);
    fLR.Flies.Data{i}.LeftEyeX = fHR.TrackData.LeftEyeX(iniVal:endVal);
    fLR.Flies.Data{i}.LeftEyeY = fHR.TrackData.LeftEyeY(iniVal:endVal);
    fLR.Flies.Data{i}.LeftEyeErr = 1-fHR.TrackData.LeftEyeLikelihood(iniVal:endVal);
    fLR.Flies.Data{i}.RightEyeX = fHR.TrackData.RightEyeX(iniVal:endVal);
    fLR.Flies.Data{i}.RightEyeY = fHR.TrackData.RightEyeY(iniVal:endVal);
    fLR.Flies.Data{i}.RightEyeErr = 1-fHR.TrackData.RightEyeLikelihood(iniVal:endVal);
    fLR.Flies.Data{i}.ThoraxX = fHR.TrackData.ThoraxX(iniVal:endVal);
    fLR.Flies.Data{i}.ThoraxY = fHR.TrackData.ThoraxY(iniVal:endVal);
    fLR.Flies.Data{i}.ThoraxErr = 1-fHR.TrackData.ThoraxLikelihood(iniVal:endVal);
    fLR.Flies.Data{i}.AbdomenX = fHR.TrackData.AbdomenX(iniVal:endVal);
    fLR.Flies.Data{i}.AbdomenY = fHR.TrackData.AbdomenY(iniVal:endVal);
    fLR.Flies.Data{i}.AbdomenErr = 1-fHR.TrackData.AbdomenLikelihood(iniVal:endVal);
    fLR.Flies.Data{i}.LeftFrontLegX = fHR.TrackData.LeftFrontLegX(iniVal:endVal);
    fLR.Flies.Data{i}.LeftFrontLegY = fHR.TrackData.LeftFrontLegY(iniVal:endVal);
    fLR.Flies.Data{i}.LeftFrontLegErr = 1-fHR.TrackData.LeftFrontLegLikelihood(iniVal:endVal);
    fLR.Flies.Data{i}.RightFrontLegX = fHR.TrackData.RightFrontLegX(iniVal:endVal);
    fLR.Flies.Data{i}.RightFrontLegY = fHR.TrackData.RightFrontLegY(iniVal:endVal);
    fLR.Flies.Data{i}.RightFrontLegErr = 1-fHR.TrackData.RightFrontLegLikelihood(iniVal:endVal);
    fLR.Flies.Data{i}.LeftMiddleLegX = fHR.TrackData.LeftMiddleLegX(iniVal:endVal);
    fLR.Flies.Data{i}.LeftMiddleLegY = fHR.TrackData.LeftMiddleLegY(iniVal:endVal);
    fLR.Flies.Data{i}.LeftMiddleLegErr = 1-fHR.TrackData.LeftMiddleLegLikelihood(iniVal:endVal);
    fLR.Flies.Data{i}.RightMiddleLegX = fHR.TrackData.RightMiddleLegX(iniVal:endVal);
    fLR.Flies.Data{i}.RightMiddleLegY = fHR.TrackData.RightMiddleLegY(iniVal:endVal);
    fLR.Flies.Data{i}.RightMiddleLegErr = 1-fHR.TrackData.RightMiddleLegLikelihood(iniVal:endVal);
    fLR.Flies.Data{i}.LeftHindLegX = fHR.TrackData.LeftHindLegX(iniVal:endVal);
    fLR.Flies.Data{i}.LeftHindLegY = fHR.TrackData.LeftHindLegY(iniVal:endVal);
    fLR.Flies.Data{i}.LeftHindLegErr = 1-fHR.TrackData.LeftHindLegLikelihood(iniVal:endVal);
    fLR.Flies.Data{i}.RightHindLegX = fHR.TrackData.RightHindLegX(iniVal:endVal);
    fLR.Flies.Data{i}.RightHindLegY = fHR.TrackData.RightHindLegY(iniVal:endVal);
    fLR.Flies.Data{i}.RightHindLegErr = 1-fHR.TrackData.RightHindLegLikelihood(iniVal:endVal);
end
Flies = fLR.Flies;
% save in a new matlab file
save([path '\DataLowHighRes.mat'], 'Flies')
disp(path)