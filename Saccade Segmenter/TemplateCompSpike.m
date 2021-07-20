function [locsF, pksIF, pksEF, score] = TemplateCompSpike(Vr,locsF,pksIF, pksEF, cmhSong, thr, params)
% load spike template
templ = load(params.spkTempPath);
templ = templ.templ;
score = ones(size(locsF));
for i = 1 : length(locsF)
   if (locsF(i) > 20) && (locsF(i) + 20 < length(Vr))
      % pad the spike so the size matches the template
      aux = sign(Vr(locsF(i)))*(Vr((locsF(i)-20):(locsF(i)+20)));
      inds = (locsF(i)-2):(locsF(i)+2);
      if(pksIF(i) < 20 && params.RG)
          aux(1:(20-max(pksIF(i),9))) = 0;
      end
      if(pksEF(i) < 20 && params.RG)
          aux((20+max(pksEF(i),11)):end) = 0;
      end
      % compute the distance from the template
      score(i) = (aux'*templ)/sqrt(aux'*aux);
      cmh = cmhSong(inds);
      cmh = max(cmh)/thr;
      if cmh > 1.8
          score(i) = 0.3;
      end
   end
end

% keep the spikes that pass the simmilarity threshold
locsF = locsF(score>params.cutoff);
pksIF = pksIF(score>params.cutoff);
pksEF = pksEF(score>params.cutoff);
score = score(score>params.cutoff);
