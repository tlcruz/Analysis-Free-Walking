clc
path = '\Fly 060';
fileList = dir(fullfile(path, '*.csv'));

% Initialize variable list to store tracking data
TrackData.frameNo = [];
TrackData.topHeadX = [];
TrackData.topHeadY = [];
TrackData.topHeadLikelihood = [];
TrackData.LeftEyeX = [];
TrackData.LeftEyeY = [];
TrackData.LeftEyeLikelihood = [];
TrackData.RightEyeX = [];
TrackData.RightEyeY = [];
TrackData.RightEyeLikelihood = [];
TrackData.NeckX = [];
TrackData.NeckY = [];
TrackData.NeckLikelihood = [];
TrackData.ThoraxX = [];
TrackData.ThoraxY = [];
TrackData.ThoraxLikelihood = [];
TrackData.AbdomenX = [];
TrackData.AbdomenY = []; 
TrackData.AbdomenLikelihood = [];
TrackData.LeftFrontLegX = [];
TrackData.LeftFrontLegY = []; 
TrackData.LeftFrontLegLikelihood = [];
TrackData.RightFrontLegX = [];
TrackData.RightFrontLegY = [];
TrackData.RightFrontLegLikelihood = [];
TrackData.LeftMiddleLegX = [];
TrackData.LeftMiddleLegY = [];
TrackData.LeftMiddleLegLikelihood = [];
TrackData.RightMiddleLegX = [];
TrackData.RightMiddleLegY = [];
TrackData.RightMiddleLegLikelihood = [];
TrackData.LeftHindLegX = [];
TrackData.LeftHindLegY = [];
TrackData.LeftHindLegLikelihood = [];
TrackData.RightHindLegX = [];
TrackData.RightHindLegY = [];
TrackData.RightHindLegLikelihood = [];

for n = 1 : length(fileList)
    % load tracking data
    pathFile = [path '\' fileList(n).name];
    M = csvread(pathFile, 3, 0);
    M = vertcat(M, M(end,:));
    M(end,1) = M(end-1,1)+1;
    % Frame Number
    if isempty(TrackData.frameNo)
        TrackData.frameNo = 1+M(:,1); 
    else
        TrackData.frameNo = vertcat(TrackData.frameNo, TrackData.frameNo(end) + 1 + M(:,1)); 
    end
    % Top Head
    TrackData.topHeadX = vertcat(TrackData.topHeadX, M(:,2));
    TrackData.topHeadY = vertcat(TrackData.topHeadY, M(:,3));
    TrackData.topHeadLikelihood = vertcat(TrackData.topHeadLikelihood, M(:,4));
    % Left Eye
    TrackData.LeftEyeX = vertcat(TrackData.LeftEyeX, M(:,5));
    TrackData.LeftEyeY = vertcat(TrackData.LeftEyeY, M(:,6));
    TrackData.LeftEyeLikelihood = vertcat(TrackData.LeftEyeLikelihood, M(:,7));
    % Right Eye
    TrackData.RightEyeX = vertcat(TrackData.RightEyeX, M(:,8));
    TrackData.RightEyeY = vertcat(TrackData.RightEyeY, M(:,9));
    TrackData.RightEyeLikelihood = vertcat(TrackData.RightEyeLikelihood, M(:,10));
    % Neck
    TrackData.NeckX = vertcat(TrackData.NeckX, M(:,11));
    TrackData.NeckY = vertcat(TrackData.NeckY, M(:,12));
    TrackData.NeckLikelihood = vertcat(TrackData.NeckLikelihood, M(:,13));
    % Thorax
    TrackData.ThoraxX = vertcat(TrackData.ThoraxX, M(:,14));
    TrackData.ThoraxY = vertcat(TrackData.ThoraxY, M(:,15));
    TrackData.ThoraxLikelihood = vertcat(TrackData.ThoraxLikelihood, M(:,16));
    % Abdomen
    TrackData.AbdomenX = vertcat(TrackData.AbdomenX, M(:,17));
    TrackData.AbdomenY = vertcat(TrackData.AbdomenY, M(:,18));
    TrackData.AbdomenLikelihood = vertcat(TrackData.AbdomenLikelihood, M(:,19));
    % Left Front Leg
    TrackData.LeftFrontLegX = vertcat(TrackData.LeftFrontLegX, M(:,20));
    TrackData.LeftFrontLegY = vertcat(TrackData.LeftFrontLegY, M(:,21));
    TrackData.LeftFrontLegLikelihood = vertcat(TrackData.LeftFrontLegLikelihood, M(:,22));
    % Right Front Leg
    TrackData.RightFrontLegX = vertcat(TrackData.RightFrontLegX, M(:,23));
    TrackData.RightFrontLegY = vertcat(TrackData.RightFrontLegY, M(:,24));
    TrackData.RightFrontLegLikelihood = vertcat(TrackData.RightFrontLegLikelihood, M(:,25));
    % Left Middle Leg
    TrackData.LeftMiddleLegX = vertcat(TrackData.LeftMiddleLegX, M(:,26));
    TrackData.LeftMiddleLegY = vertcat(TrackData.LeftMiddleLegY, M(:,27));
    TrackData.LeftMiddleLegLikelihood = vertcat(TrackData.LeftMiddleLegLikelihood, M(:,28));
    % Right Middle Leg
    TrackData.RightMiddleLegX = vertcat(TrackData.RightMiddleLegX, M(:,29));
    TrackData.RightMiddleLegY = vertcat(TrackData.RightMiddleLegY, M(:,30));
    TrackData.RightMiddleLegLikelihood = vertcat(TrackData.RightMiddleLegLikelihood, M(:,31));
    % Left Hind Leg
    TrackData.LeftHindLegX = vertcat(TrackData.LeftHindLegX, M(:,32));
    TrackData.LeftHindLegY = vertcat(TrackData.LeftHindLegY, M(:,33));
    TrackData.LeftHindLegLikelihood = vertcat(TrackData.LeftHindLegLikelihood, M(:,34));
    % Right Hind Leg
    TrackData.RightHindLegX = vertcat(TrackData.RightHindLegX, M(:,35));
    TrackData.RightHindLegY = vertcat(TrackData.RightHindLegY, M(:,36));
    TrackData.RightHindLegLikelihood = vertcat(TrackData.RightHindLegLikelihood, M(:,37));
end
% save matlab structure with tracking data
save([path '\StructTrackData.mat'], 'TrackData')
disp('Done')