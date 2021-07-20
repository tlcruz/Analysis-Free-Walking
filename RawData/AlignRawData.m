% Align Raw Data
path = uigetdir;
name = ''; % genotype of the flies
IsDark = 0; % mark 1 if the protocol was run in darkness

flies = dir(path);
protocol = 1:5400:(12*5400+1);
seq = cell(length(protocol)-1,1);
if IsDark == 1
    for i = 1 : length(seq)
        seq{i} = 'Dark';
    end
else
    % define the sequence of trials in the protocol
    seq{1} = '10NG';
    seq{2} = '10RG';
    seq{3} = '1NG';
    seq{4} = '1RG';
    seq{5} = '5NG';
    seq{6} = '5RG';
    seq{7} = '10NG';
    seq{8} = '10RG';
    seq{9} = '1NG';
    seq{10} = '1RG';
    seq{11} = '5NG';
    seq{12} = '5RG';
end
%
for i = 3 : length(flies)
    pathFly = [path '\' flies(i).name];
    flies(i).name
    d = dir(pathFly);
    if exist([pathFly '\DataLowRes.mat'], 'file') == 2
        delete([pathFly '\DataLowRes.mat'])
    end
    pathData = [pathFly '\Cameras.txt'];
    if exist(pathData, 'file') == 2
        [Flies] = GetRawDataUniS(protocol, pathData, [pathFly '\'], seq,0, name);
    else
        disp('Data File not found.')
    end
end