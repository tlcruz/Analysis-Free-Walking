function [Bouts, actState] = Vec2Bout(actState, t1, t0, exp)
if nargin < 4
    exp = 0;
end
a = 1;
thisBout = [];
for i = 1 : length(actState) - 1
    if actState(i) == 1 && actState(i + 1) == 1
        thisBout = vertcat(thisBout, i);
    elseif actState(i) == 1 && actState(i + 1) == 0
        thisBout = vertcat(thisBout, i);
        if length(thisBout) < t1
            actState(thisBout) = 0;
        end
        thisBout = [];
    end
end
thisBout = [];
for i = 1 : length(actState) - 1
    if actState(i) == 0 && actState(i + 1) == 0
        thisBout = vertcat(thisBout, i);
    elseif actState(i) == 0 && actState(i + 1) == 1
        thisBout = vertcat(thisBout, i);
        if length(thisBout) < t0
            if exp ~= 0
                if thisBout(1) - exp > 0
                    thisBout = vertcat((thisBout(1)-(1:exp)'), thisBout);
                end
                if thisBout(end)+ exp < length(actState)
                    thisBout = vertcat(thisBout, thisBout(end)+(1:exp)');
                end
            end
            actState(thisBout) = 1;
        end
        thisBout = [];
    end
end
% Get Bouts
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
Bouts = cell(b, 1);
for j = 1 : (length(ind1) - 1)
    if (ind1(j+1) - ind1(j) == 1)
        inds(a) = ind1(j);
        a = a + 1;
    else
        inds(a) = ind1(j);
        if length(inds) > t1
            Bouts{c} = inds;
            c = c + 1;
        end
        a = 1;
        clearvars inds;
    end
end
if length(ind1)~= 0 && length(ind0)~= 0
    if exist('inds')
        if length(inds) > t1
            Bouts{c} = inds;
        end
    end
elseif length(ind1)~= 0 && length(ind0)== 0
    if exist('inds')
        if length(inds) > t1
            Bouts{c} = inds;
        end
    end
end
clearvars inds ind1 ind2;
while(isempty(Bouts{end}))
    ActBoutsT = Bouts;
    clearvars ActBouts1
    Bouts = cell(length(ActBoutsT)-1,1);
    for i = 1 : length(ActBoutsT)-1
        Bouts(i) = ActBoutsT(i);
    end
    clearvars ActBoutsT
    if length(Bouts) < 1
        break;
    end
end
end