%% Function to create a 2D histogram for variables X and Y with bins CentX and CentY
function [mapLength] = Occupancy(X, Y, CentX, CentY)
for x = 1 : length(CentX)
    if x < length(CentX)
        xbin = CentX(x+1)-CentX(x);
    else
        xbin = -CentX(x-1)+CentX(x);
    end
    tmpa=find(X >= CentX(1)+(x-1)*xbin);
    tmpb=find(X < CentX(1)+x*xbin);
    tmpx=tmpa(ismembc(tmpa,tmpb));
    for y = 1 : length(CentY)
        if y < length(CentY)
            ybin = CentY(y+1)-CentY(y);
        else
            ybin = -CentY(y-1)+CentY(y);
        end
        tmpa=find(Y >= CentY(1)+(y-1)*ybin);
        tmpb=find(Y < CentY(1)+y*ybin);
        tmpy=tmpa(ismembc(tmpa,tmpb));
        tmp=tmpx(ismembc(tmpx,tmpy));
        if ~isempty(tmp)
            MPMDVl{x,y} = length(tmp);
        else
            MPMDVl{x,y} = 0;
        end
    end
end
mapLength = cell2mat(MPMDVl);
end

