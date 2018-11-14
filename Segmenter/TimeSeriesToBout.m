function [ActBouts] = TimeSeriesToBout(actState, dtw)
%% Get Bouts
ind1 = find(actState == 1);
ind0 = find(actState == 0);
b = 1;
for j = 1 : (length(ind1) - 1)
    if (ind1(j+1) - ind1(j) == 1)
    else
        b = b + 1;
    end
end
c = 1;
a = 1;
ActBouts = cell(b, 1);
for j = 1 : (length(ind1) - 1)
    if (ind1(j+1) - ind1(j) == 1)
        inds(a) = ind1(j);
        a = a + 1;
    else
        inds(a) = ind1(j);
        if length(inds) > dtw
            ActBouts{c} = inds;
            c = c + 1;
        end
        a = 1;
        clearvars inds;
    end
end
if length(ind1)~= 0 && length(ind0)~= 0
    if exist('inds')
        if length(inds) > dtw
            ActBouts{c} = inds;
        end
    end
elseif length(ind1)~= 0 && length(ind0)== 0
    if exist('inds')
        if length(inds) > dtw
            ActBouts{c} = inds;
        end
    end
end
clearvars inds ind1 ind2;
while(isempty(ActBouts{end}))
    ActBoutsT = ActBouts;
    clearvars ActBouts1
    ActBouts = cell(length(ActBoutsT)-1,1);
    for i = 1 : length(ActBoutsT)-1
        ActBouts(i) = ActBoutsT(i);
    end
    clearvars ActBoutsT
    if length(ActBouts) < 1
        break;
    end
end
end

