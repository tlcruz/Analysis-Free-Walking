function [params] = GetParams()
    % params initial detection
    params.Fs = 60;
    params.minF = 10;
    params.maxF = 15;
    params.nstds = 2; % 
    params.MinPeakDist = 0;
    params.thrDist = 5;
    params.MinPeakPromF = 0;
    params.MinPeakPromV = 20;
    params.maxP = 9;
    params.minP = 7;
    
    params.pkvt = 200;
    params.promt = 150;
    params.vfstdt = 3;
    params.vfmt = 6;
    
    params.vft = 6;
    
    params.spkTempPath = 'SpikeTemplateL.mat';
    params.cutoff = 0.15;
    params.RG = 0;
    params.thr = 2*60;
    params.thrspk = 100;
    params.delta = 15;
    params.sep = 20;
    params.windStr = 20;
    params.maxvf = 0.4;
    params.minvf = 5;
    params.minStrB = 20; %
    params.mDistWall = 0;
    
    params.WindowCCVFW = 20;
    params.WindowCCVALL = 40;
    params.btthr = 21;
end