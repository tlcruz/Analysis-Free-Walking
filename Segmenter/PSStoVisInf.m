function [VI, NVI, MVI, SEMVI] = PSStoVisInf(PSSALL, NPSSALL, params)
VI = [];
NVI = [];
for n = 1 : size(PSSALL,2) 
    pssng = [];
    npssng = [];
    pssrg = [];
    npssrg = [];
    
    
    if length(params.pTypes) > 1
        for j = 1 : length(params.pTypes)
            if mod(j,2)==0
                if~isempty(find((NPSSALL{j,n}>params.vthr)==1))
                    npssa = NPSSALL{j,n}(~isnan(PSSALL{j,n}));
                    pssa = PSSALL{j,n}(~isnan(PSSALL{j,n}));
                    pssa = pssa(npssa>params.vthr);
                    npssa = npssa(npssa>params.vthr);
                    pssrg = vertcat(pssrg,pssa'*npssa/sum(npssa));
                    npssrg = vertcat(npssrg, sum(npssa));
                else
                    pssrg = vertcat(pssrg,nan);
                    npssrg = vertcat(npssrg,nan);
                end
            else
                if~isempty(find((NPSSALL{j,n}>params.vthr)==1))
                    npssa = NPSSALL{j,n}(~isnan(PSSALL{j,n}));
                    pssa = PSSALL{j,n}(~isnan(PSSALL{j,n}));
                    pssa = pssa(npssa>params.vthr);
                    npssa = npssa(npssa>params.vthr);
                    pssng = vertcat(pssng,pssa'*npssa/sum(npssa));
                    npssng = vertcat(npssng, sum(npssa));
                else
                    pssng = vertcat(pssng,nan);
                    npssng = vertcat(npssng,nan);
                end
            end
        end
    else
        if~isempty(find((NPSSALL{1,n}>params.vthr)==1))
            npssa = NPSSALL{1,n}(~isnan(PSSALL{1,n}));
            pssa = PSSALL{1,n}(~isnan(PSSALL{1,n}));
            pssa = pssa(npssa>params.vthr);
            npssa = npssa(npssa>params.vthr);
            if length(npssa) > 1
                pssrg = vertcat(pssrg,pssa(1:floor(end/2))'*npssa(1:floor(end/2))/sum(npssa(1:floor(end/2))));
                npssrg = vertcat(npssrg, sum(npssa(1:floor(end/2))));
                pssng = vertcat(pssng,pssa(floor(end/2):end)'*npssa(floor(end/2):end)/sum(npssa(floor(end/2):end)));
                npssng = vertcat(npssng, sum(npssa(floor(end/2):end)));
            end

        end
    end
    
    if(isempty(find(isnan(pssrg)==1,1)) && isempty(find(isnan(pssng)==1,1)))
        VI = horzcat(VI,(pssrg-pssng)./(pssrg+pssng));
        NVI = horzcat(NVI,npssrg+npssng);
%         VI = horzcat(VI,(pssrg));
%         NVI = horzcat(NVI,npssrg);
    else
        aux = (pssrg-pssng)./(pssrg+pssng);
%         aux = (pssrg);
        aux(isnan(aux)) = 0;
        VI = horzcat(VI,aux);
        aux = npssrg+npssng;
%         aux = npssrg;
        aux(isnan(aux)) = 0;
        NVI = horzcat(NVI,aux);
    end
end

MVI = [];%zeros(length(params.pTypes)/2,1);
SEMVI = [];%zeros(length(params.pTypes)/2,1);
for i = 1 : size(VI,1)
    vi = VI(i,:)';
    nvi = NVI(i,:)';
    nvi(vi>250)=0;
    MVI = vertcat(MVI, vi'*nvi/sum(nvi));
    SEMVI = vertcat(SEMVI, sqrt(((vi-MVI(i)).*(vi-MVI(i)))'*nvi/sum(nvi))/sqrt(length(vi)));
end

end