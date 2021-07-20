function [M, SEM, MF, NF] = GetGMSEM(STR, DST, N, params)
% calculate weighted grand mean and standard error
% DST works as weights and STR is the variable to be averaged
M = zeros(length(params.pTypes),1);
SEM = zeros(length(params.pTypes),1);
MF = zeros(length(params.pTypes),size(STR,2));
NF = zeros(length(params.pTypes),size(STR,2));
a = 0;
b = 0;
v = 120;
for j = 1 : length(params.pTypes)
    st = [];
    nst = [];
    for n = 1 : size(STR,2)
        if~isempty(STR{j,n})
            b = b + length(STR{j,n});
            if ~isempty(find(STR{j,n}>v))
                a=a+length(DST{j,n}(STR{j,n}>v));
            end
            DST{j,n}(STR{j,n}>v) = 0;
            st = vertcat(st, STR{j,n}'*DST{j,n}/sum(DST{j,n}));
            nst = vertcat(nst, sum(N{j,n}));
            MF(j,n) = STR{j,n}'*DST{j,n}/sum(DST{j,n});
            NF(j,n) = sum(N{j,n});
        end
    end
    st = st(nst > params.lthr);
    nst = nst(nst > params.lthr);
    M(j) = st'*nst/sum(nst);
    SEM(j) = sqrt(((st-M(j)).*(st-M(j)))'*nst/sum(nst))/sqrt(length(st));
end
end