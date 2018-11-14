function [BoutVr, BoutVf, BoutP, BoutP2] = GetProbVec(Vr, Vf, delta, mode, min, max)
    inds = Vr;
    inds(inds>0) = 1;
    inds(inds<=0) = -1;

    aux1 = inds(1:end-delta+1);
    Vr1 = Vr(1:end-delta+1);
    aux2 = inds(delta:end);
    aux3 = abs(aux1+aux2)/2;

    BoutVr = Vr(1:end-delta+1);
    BoutVf = Vf(1:end-delta+1);
    BoutP2 = aux3;
    if mode == 1
        aux3 = aux3(Vr1 > min & Vr1 < max);
        BoutP = aux3;
    elseif mode == 2
        aux2 = aux2(Vr1 > min & Vr1 < max);
        BoutP = aux2;
    elseif mode == 3
        aux3 = aux3(abs(Vr1) > min & abs(Vr1) < max);
        BoutP = aux3;
    else
        BoutP = aux2;
    end
end