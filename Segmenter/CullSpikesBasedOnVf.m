function [locsF, pksIF,pksEF, pV] = CullSpikesBasedOnVf(Vr, locs, cmhSong, thr, Vf, params)
spksVect = zeros(length(Vr),1);
spksIni = zeros(length(locs),1);
spksEnd = zeros(length(locs),1);
for i = 1 : length(locs)
   if (locs(i) <= params.minP)
       vraux = Vr(locs(i):(locs(i)+params.maxP));
       aaux = diff(vraux);
       vraux = vraux(1:end-1).*vraux(2:end);
       env = find(vraux(5:end)<0);
       aaux = aaux(1:end-1).*aaux(2:end);
       enp = find(aaux(5:end)<0);
       if(isempty(env) && isempty(enp))
           en = params.maxP;
       elseif(isempty(env) && ~isempty(enp))
           en = enp(1)+2;
       elseif(~isempty(env) && isempty(enp))
           en = env(1)+1;
       else
           en = min(enp(1)+2, env(1)+1);
       end
       spksEnd(i) = en;
       spksIni(i) = locs(i);
       spksVect(1:(locs(i)+spksEnd(i))) = 1;
   elseif (locs(i)+params.maxP >= length(Vr))
       vraux = Vr(locs(i)-params.minP:locs(i));
       aaux = diff(vraux);
       vraux = vraux(1:end-1).*vraux(2:end);
       inv = find(vraux(1:5)<0);
       aaux = aaux(1:end-1).*aaux(2:end);
       inp = find(aaux(1:5)<0);
       if(isempty(inv) && isempty(inp))
           in = 0;
       elseif(isempty(inv) && ~isempty(inp))
           in = inp(1);
       elseif(~isempty(inv) && isempty(inp))
           in = inv(1);
       else
           in = min(inp(1), inv(1));
       end
       spksEnd(i) = length(Vr)-locs(i);
       spksIni(i) = params.minP-in;
       spksVect((locs(i)-spksIni(i)):end) = 1;
   else
       vraux = Vr(locs(i):(locs(i)+params.maxP));
       aaux = diff(vraux);
       vraux = vraux(1:end-1).*vraux(2:end);
       env = find(vraux(4:end)<0);
       aaux = aaux(1:end-1).*aaux(2:end);
       enp = find(aaux(5:end)<0);
       if(isempty(env) && isempty(enp))
           en = params.maxP;
       elseif(isempty(env) && ~isempty(enp))
           en = enp(1)+2;
       elseif(~isempty(env) && isempty(enp))
           en = env(1)+2;
       else
           en = min(enp(1)+2, env(1)+2);
       end
       vraux = Vr(locs(i)-params.minP:locs(i));
       aaux = diff(vraux);
       vraux = vraux(1:end-1).*vraux(2:end);
       inv = find(vraux(1:5)<0);
       aaux = aaux(1:end-1).*aaux(2:end);
       inp = find(aaux(1:5)<0);
       if(isempty(inv) && isempty(inp))
           in = 0;
       elseif(isempty(inv) && ~isempty(inp))
           in = inp(1);
       elseif(~isempty(inv) && isempty(inp))
           in = inv(1);
       else
           in = min(inp(1), inv(1));
       end
       spksEnd(i) = en;
       spksIni(i) = params.minP-in;
       spksVect((locs(i)-spksIni(i)):(locs(i)+spksEnd(i))) = 1;
   end
end

locsF = [];
pksIF= [];
pksEF = [];

pV.peakVal = [];
pV.peakProm = [];
pV.vfstd = [];
pV.vfmeam = [];
pV.vfmin = [];
for i = 1 : length(locs)
    if(locs(i)>spksIni(i) && locs(i)+spksEnd(i)<length(Vr))
        vrx = Vr((locs(i)-spksIni(i)):(locs(i)+spksEnd(i)));
        vfx = Vf((locs(i)-spksIni(i)):(locs(i)+spksEnd(i)));
        cwtx = cmhSong((locs(i)-spksIni(i)):(locs(i)+spksEnd(i)));
        
        peakVal = max(abs(vrx));
        peakProm = max(abs(vrx))/((vrx(1)+vrx(end))/2);
        vfstd = std(vfx);
        vfmeam = mean(vfx);
        
            pV.peakVal = vertcat(pV.peakVal, peakVal);            
            pV.peakProm = vertcat(pV.peakProm, peakProm);           
            pV.vfstd = vertcat(pV.vfstd, vfstd);            
            pV.vfmeam = vertcat(pV.vfmeam, vfmeam);
            pV.vfmin = vertcat(pV.vfmin, min(vfx));
        
            inds = (locs(i)-2):(locs(i)+2);
            cmh = cmhSong(inds);
            cmh = max(cmh)/thr;
        if (cmh < 1.8 && peakVal < params.pkvt && peakProm < params.promt && vfstd < params.vfstdt)
        else
            locsF = vertcat(locsF, locs(i));
            pksIF = vertcat(pksIF, spksIni(i));
            pksEF = vertcat(pksEF, spksEnd(i));


        end
    end
end































end
