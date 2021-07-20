function [Flies] = GetRawDataUniS(protocol, path, pathSave, seq, iOr, name)
% Load Data
Flies.Protocol = protocol;
Flies.Seq = seq;
Flies.Type = name;
T = readtable(path, 'Delimiter',' ','ReadVariableNames',false);
% Divide Data Into Protocol
startI = [];
endI = [];
Framei1 = T.Var1(1);
Framei2 = T.Var3(1);
FrameC1 = T.Var1-Framei1+1;
FrameC2 = T.Var3-Framei2+1;

for k = 1 : length(protocol)-1
    if protocol(k) < max(FrameC1)
        auxStart = find(FrameC1 == protocol(k));
        if protocol(k+1) < max(FrameC1)
            auxEnd = find(FrameC1 == protocol(k+1));
        else
            auxEnd = max(FrameC1);
        end
        startI = vertcat(startI, auxStart);
        endI = vertcat(endI, auxEnd);
    else
        startI = vertcat(startI, []);
        endI = vertcat(endI, []);
    end
end
% Get Tracking Traces
Angle = T.Var7(1:endI(end));
X = -80*(T.Var6(1:endI(end))+0.9)/9;
Y = 80*T.Var5(1:endI(end))/9;
clearvars T
VfRaw = 60*((X(2:end)-X(1:end-1)).*cos(pi*Angle(1:end-1)/180) + (Y(2:end)-Y(1:end-1)).*sin(pi*Angle(1:end-1)/180));
VsRaw = 60*(-(X(2:end)-X(1:end-1)).*sin(pi*Angle(1:end-1)/180) + (Y(2:end)-Y(1:end-1)).*cos(pi*Angle(1:end-1)/180));
VrRaw = 60*diff(Angle);
Angle = Angle(1:end);
VtRaw = sqrt(VfRaw.*VfRaw + VsRaw.*VsRaw);
WallDist = 45 - sqrt(X.*X + Y.*Y); 
Flies.VfRaw = VfRaw;
Flies.VsRaw = VsRaw;
Flies.VrRaw = VrRaw;
Flies.AngRaw = Angle;
nJb = 0;
nJ = 1;
flp = iOr * ones(length(VfRaw),1);
while (nJ ~= nJb)
    nJb = nJ;
    [VfRaw, VsRaw, VrRaw, VtRaw, Angle, nJ, flp] = ClearJumps(VfRaw, VsRaw, VrRaw, VtRaw, Angle, flp, iOr);
end
flp = vertcat(flp, flp(end));
% Smooth Data
nPoints = 6;
Vr = smooth(vertcat(VrRaw,0), nPoints/length(VrRaw), 'lowess');
Vf = smooth(vertcat(VfRaw,0), nPoints/length(VfRaw), 'lowess');
Vt = smooth(vertcat(VtRaw,0), nPoints/length(VtRaw), 'lowess');
Vs = smooth(vertcat(VsRaw,0), nPoints/length(VsRaw), 'lowess');

% Divide the data into protocol parts and bouts
Flies.Data = cell(length(startI), 1);
Flies.Paths = cell(length(startI), 1);
percWalk = 0;
for k = 1 : length(startI)
    Flies.Data{k}.Vf = Vf(startI(k):endI(k));
    Flies.Data{k}.Vr = Vr(startI(k):endI(k));
    Flies.Data{k}.Vs = Vs(startI(k):endI(k));
    Flies.Data{k}.Vt = Vt(startI(k):endI(k));
    Flies.Data{k}.X = X(startI(k):endI(k));
    Flies.Data{k}.Y = Y(startI(k):endI(k));
    Flies.Data{k}.Angle = Angle(startI(k):endI(k));
    Flies.Data{k}.flp = flp(startI(k):endI(k));
    Flies.Data{k}.FramesC1 = FrameC1(startI(k):endI(k));
    Flies.Data{k}.FramesC2 = FrameC2(startI(k):endI(k));
    Flies.Data{k}.WallDist = WallDist(startI(k):endI(k));
    [Flies.Data{k}.actState, Flies.Data{k}.Bouts] = GetBouts(Flies.Data{k}.Vr, Flies.Data{k}.Vf, Flies.Data{k}.Vs);
    Flies.Data{k}.StartBout = zeros(size(Flies.Data{k}.Bouts,1),1);
    Flies.Data{k}.EndBout = zeros(size(Flies.Data{k}.Bouts,1),1);
    Flies.Data{k}.StartBoutExt = zeros(size(Flies.Data{k}.Bouts,1),1);
    Flies.Data{k}.EndBoutExt = zeros(size(Flies.Data{k}.Bouts,1),1);
    Flies.Data{k}.StartBoutC2 = zeros(size(Flies.Data{k}.Bouts,1),1);
    Flies.Data{k}.EndBoutC2 = zeros(size(Flies.Data{k}.Bouts,1),1);
    for n = 1 : size(Flies.Data{k}.Bouts,1)
        percWalk = percWalk + length(Flies.Data{k}.Bouts{n});
        if length(Flies.Data{k}.Bouts{n}) > 0*60 && ...
                mean(Flies.Data{k}.Vf(Flies.Data{k}.Bouts{n})) < 0
            Flies.Data{k}.flp(Flies.Data{k}.Bouts{n}) = 180;
            Flies.Data{k}.Vf(Flies.Data{k}.Bouts{n}) = -Flies.Data{k}.Vf(Flies.Data{k}.Bouts{n});
            Flies.Data{k}.Vs(Flies.Data{k}.Bouts{n}) = -Flies.Data{k}.Vs(Flies.Data{k}.Bouts{n});
            Flies.Data{k}.Angle(Flies.Data{k}.Bouts{n}) = Flies.Data{k}.Angle(Flies.Data{k}.Bouts{n})-180;
        end
        Flies.Data{k}.StartBout(n) = Flies.Data{k}.Bouts{n}(1);
        Flies.Data{k}.EndBout(n) = Flies.Data{k}.Bouts{n}(end);
        if Flies.Data{k}.EndBout(n) < length(Flies.Data{k}.FramesC1) - 30
            Flies.Data{k}.EndBoutExt(n) = Flies.Data{k}.Bouts{n}(end)+30;
            Flies.Data{k}.EndBoutC2(n) = Flies.Data{k}.FramesC2(Flies.Data{k}.EndBoutExt(n));
        else
            Flies.Data{k}.EndBoutExt(n) = Flies.Data{k}.Bouts{n}(end);
            Flies.Data{k}.EndBoutC2(n) = Flies.Data{k}.FramesC2(Flies.Data{k}.EndBoutExt(n));
        end
        if Flies.Data{k}.StartBout(n) > 30
            Flies.Data{k}.StartBoutExt(n) = Flies.Data{k}.Bouts{n}(1)-30;
            Flies.Data{k}.StartBoutC2(n) = Flies.Data{k}.FramesC2(Flies.Data{k}.StartBoutExt(n));
        else
            Flies.Data{k}.StartBoutExt(n) = Flies.Data{k}.Bouts{n}(1);
            Flies.Data{k}.StartBoutC2(n) = Flies.Data{k}.FramesC2(Flies.Data{k}.StartBoutExt(n));
        end
    end

end
percWalk = percWalk / protocol(end);
Flies.percWalk = percWalk;
if nargin > 2
    disp(num2str(percWalk))
   save([pathSave '\DataLowRes.mat'], 'Flies'); 
end

end