function [traceFreq] = GetMainFreqs(sig, params)
fc = params.minF:0.1:params.maxF;             % wavelet scales, in Hertz
DoGwvlt = [2 3];               % Derivative of Gaussian wavelets examined
ngw = numel(DoGwvlt);
wvlt = cell(1,ngw);
for i = 1:ngw
    wvlt{i} = ['gaus' num2str(DoGwvlt(i))];
end
sc = zeros(ngw,numel(fc));
for i = 1:numel(wvlt)
    K = scal2frq(1,wvlt{i},1/params.Fs);
    sc(i,:) = K./fc;
end
traceFreq = single(zeros(1,numel(sig))); 
cmh_dog = int8(zeros(1,length(sig)));
cmh_sc = int8(zeros(1,length(sig)));
for i= 1:numel(wvlt)
    for j = 1:size(sc,2)
        temp = single(abs(cwt(sig,sc(i,j),wvlt{i})));
        cmh_sc(temp>traceFreq) = j;
        cmh_dog(temp>traceFreq) = i;
        traceFreq(temp>traceFreq) = temp(temp>traceFreq);
        clear temp;
    end
end
end