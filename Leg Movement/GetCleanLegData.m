%% function to clean leg data by thresolding in the tracking quality parameter
function [lD] = GetCleanLegData(LegMovForwSeg, thr)
lD.ThX = cell(1,1);
lD.ThY = cell(1,1);
lD.FLLX = cell(1,1);
lD.FLLY = cell(1,1);
lD.FRLX = cell(1,1);
lD.FRLY = cell(1,1);
lD.MLLX = cell(1,1);
lD.MLLY = cell(1,1);
lD.MRLX = cell(1,1);
lD.MRLY = cell(1,1);
lD.HLLX = cell(1,1);
lD.HLLY = cell(1,1);
lD.HRLX = cell(1,1);
lD.HRLY = cell(1,1);
lD.VR = cell(1,1);
lD.VF = cell(1,1);
lD.HAng = cell(1,1);
lD.Time = cell(1,1);
jk = 1;
for i = 1 : length(LegMovForwSeg)
    if ~isempty(LegMovForwSeg{i})
        if (length(LegMovForwSeg{i}.VrLR) > 3)
            % resize low res variables to match the high res data
            vr = LegMovForwSeg{i}.VrLR;
            vf = LegMovForwSeg{i}.VfLR;
            inds = LegMovForwSeg{i}.Inds;
            
            lD.VR{jk} = resample(vr, length(LegMovForwSeg{i}.FLX), length(vr));
            lD.VF{jk} = resample(vf, length(LegMovForwSeg{i}.FLX), length(vf));
            
            lD.Time{jk} = linspace(inds(1)/60,inds(end)/60, length(LegMovForwSeg{i}.FLX));
            % clean leg position traces
            fflx = LegMovForwSeg{i}.FLX;
            ffly = LegMovForwSeg{i}.FLY;
            fflerr = LegMovForwSeg{i}.FLErr;
            fflx(fflerr > thr) = nan;
            ffly(fflerr > thr) = nan;
            lD.FLLX{jk} = fflx;
            lD.FLLY{jk} = ffly;
            
            ffrx = LegMovForwSeg{i}.FRX;
            ffry = LegMovForwSeg{i}.FRY;
            ffrerr = LegMovForwSeg{i}.FRErr;
            ffrx(ffrerr > thr) = nan;
            ffry(ffrerr > thr) = nan;
            lD.FRLX{jk} = ffrx;
            lD.FRLY{jk} = ffry;
            
            mmlx = LegMovForwSeg{i}.MLX;
            mmly = LegMovForwSeg{i}.MLY;
            mmlerr = LegMovForwSeg{i}.MLErr;
            mmlx(mmlerr > thr) = nan;
            mmly(mmlerr > thr) = nan;
            lD.MLLX{jk} = mmlx;
            lD.MLLY{jk} = mmly;
            
            mmrx = LegMovForwSeg{i}.MRX;
            mmry = LegMovForwSeg{i}.MRY;
            mmrerr = LegMovForwSeg{i}.MRErr;
            mmrx(mmrerr > thr) = nan;
            mmry(mmrerr > thr) = nan;
            lD.MRLX{jk} = mmrx;
            lD.MRLY{jk} = mmry;
            
            hhlx = LegMovForwSeg{i}.HLX;
            hhly = LegMovForwSeg{i}.HLY;
            hhlerr = LegMovForwSeg{i}.HLErr;
            hhlx(hhlerr > thr) = nan;
            hhly(hhlerr > thr) = nan;
            lD.HLLX{jk} = hhlx;
            lD.HLLY{jk} = hhly;
            
            hhrx = LegMovForwSeg{i}.HRX;
            hhry = LegMovForwSeg{i}.HRY;
            hhrerr = LegMovForwSeg{i}.HRErr;
            hhrx(hhrerr > thr) = nan;
            hhry(hhrerr > thr) = nan;
            lD.HRLX{jk} = hhrx;
            lD.HRLY{jk} = hhry;
            
            lD.ThX{jk} = LegMovForwSeg{i}.ThX;
            lD.ThY{jk} = LegMovForwSeg{i}.ThY;
            jk = jk + 1;
        end
    end
end
end