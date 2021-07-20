function GuiPlotRawData()
%% Get Directory and Individual Fly Folders
path = uigetdir;
path = [path '\'];
flies = dir([path, '\*']);
flies = flies(3:end);
flies(~[flies.isdir]) = [];
NFlies = length(flies);
currentFly = 1;
%% Load DataFile Current Fly
dt = load([path flies(currentFly).name '\DataLowRes.mat']);
dt = dt.Flies;
nParts = length(dt.Data);
% nParts = 10;
currentPart = 1;
%% Define GUI elements
f = figure('Visible', 'off', 'Position', [360, 300, 500, 320]);
f.Name = 'NavigateVideo';
movegui(f, 'center')
f.Visible = 'on';

% FlyNo // PartNo
NFlyplus = uicontrol('Style', 'pushbutton', 'String', '+',...
    'Position', [70, 298, 10, 7], 'Callback', {@NFlyplus_Callback});
NFlyplus.Units = 'normalized';
NFlyminus = uicontrol('Style', 'pushbutton', 'String', '-',...
    'Position', [70, 291, 10, 7], 'Callback', {@NFlyminus_Callback});
NFlyminus.Units = 'normalized';
NFlyT = uicontrol(f,'Style','edit','String','FlyNo',...
    'Position',[50, 291, 20, 14], 'Callback', {@NFlyT_CallBack});
NFlyT.Units = 'normalized';
NPartplus = uicontrol('Style', 'pushbutton', 'String', '+',...
    'Position', [70, 282, 10, 7],'Callback', {@NPartplus_Callback});
NPartplus.Units = 'normalized';
NPartminus = uicontrol('Style', 'pushbutton', 'String', '-', ...
    'Position', [70, 275, 10, 7],'Callback', {@NPartminus_Callback});
NPartminus.Units = 'normalized';
NPartT = uicontrol(f,'Style','edit','String','PartNo',...
    'Position',[50, 275, 20, 14],'Callback', {@NPartT_CallBack});
NPartT.Units = 'normalized';

% ZoomPlot
textMin = uicontrol(f,'Style','text', 'String','Start:',...
    'Position',[85, 285, 15, 15]);
textMin.Units = 'normalized';
sliderMin = uicontrol(f,'Style','slider', 'Min',0,...
    'Max',length(dt.Data{currentPart}.Vr)/60,'Value',0,...
    'Position',[100, 290, 60, 15], 'Callback',{@sumin});
sliderMin.Units = 'normalized';
txtboxMin = uicontrol(f,'Style','edit','String','TMin',...
    'Position',[160, 290, 15, 15], 'Callback', {@txtbxmin});
txtboxMin.Units = 'normalized';

textMax = uicontrol(f,'Style','text', 'String','End:',...
    'Position',[85, 270, 15, 15]);
textMax.Units = 'normalized';
sliderMax = uicontrol(f,'Style','slider', 'Min',0,...
    'Max',length(dt.Data{currentPart}.Vr)/60,'Value', ...
    length(dt.Data{currentPart}.Vr)/60,...
    'Position',[100, 275, 60, 15], 'Callback',{@sumax});
sliderMax.Units = 'normalized';
txtboxMax = uicontrol(f,'Style','edit','String','TMax',...
    'Position',[160, 275, 15, 15], 'Callback', {@txtbxmax});
txtboxMax.Units = 'normalized';


zPlotMin = 0;
zPlotMax = length(dt.Data{currentPart}.Vr)/60;


ssmin = -100;
ssmax = 100;
ssmode = 1;

for n = 1 : nParts
    seq = dt.Seq;
    if n == floor(currentPart)
        nseq{n} = uicontrol('Style','text','position',[50+20*n, 290-30, 25, 10],'String', seq{n}, 'BackgroundColor', 'r');
    else
        nseq{n} = uicontrol('Style','text','position',[50+20*n, 290-30, 25, 10],'String', seq{n}, 'BackgroundColor', 'w');
    end
    nseq{n}.Units = 'normalized';
end

ha = axes('Units', 'pixels', 'Position', [250, 50, 220, 200]);
ha.Units = 'normalized';
hb = axes('Units', 'pixels', 'Position', [50, 50, 200, 200]);
hb.Units = 'normalized';
axes(hb);
axis off;
axes(ha);
axis off;
%% CallBacks
% FlyNo // PartNo
    function NFlyplus_Callback(source, eventdata)
        if currentFly <= NFlies
            currentFly = currentFly + 1;
            updateAxes(currentFly,currentPart);
        end
    end
    function NFlyminus_Callback(source, eventdata)
        if currentFly > 1
            currentFly = currentFly - 1;
            updateAxes(currentFly,currentPart);
        end
    end
    function NFlyT_CallBack(hObject, eventdata, handles)
        input = str2double(get(hObject,'String'));
        if isnan(input)
            errordlg('You must enter a numeric value','Invalid Input','modal')
            uicontrol(hObject)
            return
        else
            input = floor(input);
            if input <= NFlies && input >= 1
                currentFly = input;
                updateAxes(currentFly,currentPart);
            else
                currentFly = 1;
                updateAxes(currentFly,currentPart);
            end
        end
    end
    function NPartplus_Callback(source, eventdata)
        if currentPart <= nParts
            currentPart = currentPart + 1;
            zPlotMin = 0;
            zPlotMax = length(dt.Data{currentPart}.Vr)/60;
            updateAxes(currentFly,currentPart);
        end
    end
    function NPartminus_Callback(source, eventdata)
        if currentPart > 1
            currentPart = currentPart - 1;
            zPlotMin = 0;
            zPlotMax = length(dt.Data{currentPart}.Vr)/60;
            updateAxes(currentFly,currentPart);
        end
    end
    function NPartT_CallBack(hObject, eventdata, handles)
        input = str2double(get(hObject,'String'));
        if isnan(input)
            errordlg('You must enter a numeric value','Invalid Input','modal')
            uicontrol(hObject)
            return
        else
            input = floor(input);
            if input <= nParts && input >= 1
                currentPart = input;
                zPlotMin = 0;
                zPlotMax = length(dt.Data{currentPart}.Vr)/60;
                updateAxes(currentFly,currentFly);
            else
                currentPart = 1;
                zPlotMin = 0;
                zPlotMax = length(dt.Data{currentPart}.Vr)/60;
                updateAxes(currentFly,currentPart);
            end
        end
    end

% ZoomPlot
    function sumin(h,event)
        val = 0.01*floor(100*get(h,'Value'));
        if((val) < zPlotMax)
            zPlotMin=(val);
        else
            set(h,'Value',0);
            zPlotMin = 0;
        end
        txtboxMin.String = num2str(zPlotMin);
        updateAxes(currentFly,currentPart);
    end
    function txtbxmin(hObject, eventdata, handles)
        input = str2double(get(hObject,'String'));
        if isnan(input)
            errordlg('You must enter a numeric value','Invalid Input','modal')
            uicontrol(hObject)
            return
        else
            if(input < zPlotMax)
                zPlotMin=input;
            else
                zPlotMin = 0;
                txtboxMin.String = num2str(zPlotMin);
            end
            updateAxes(currentFly,currentPart);
        end
    end
    function sumax(h,event)
        val = 0.01*floor(100*get(h,'Value'));
        if((val) > zPlotMin)
            zPlotMax=(val);
        else
            set(h,'Value',length(dt.Data{currentPart}.Vr)/60);
            zPlotMax = length(dt.Data{currentPart}.Vr)/60;
        end
        txtboxMax.String = num2str(zPlotMax);
        updateAxes(currentFly,currentPart);
    end
    function txtbxmax(hObject, eventdata, handles)
        input = str2double(get(hObject,'String'));
        if isnan(input)
            errordlg('You must enter a numeric value','Invalid Input','modal')
            uicontrol(hObject)
            return
        else
            if(input > zPlotMin)
                zPlotMax = input;
            else
                zPlotMax = length(dt.Data{currentPart}.Vr)/60;
                txtboxMax.String = num2str(zPlotMax);
            end
            updateAxes(currentFly,currentPart);
        end
    end

%% UpdateAxes
    function updateAxes(cFly, cpart)
        dt = load([path flies(cFly).name '\DataLowRes.mat']);
        dt = dt.Flies;
        nParts = length(dt.Data);
        if cpart > nParts
            cpart = 1;
            currentPart = 1;
        end
        NFlyT.String = num2str(cFly);
        NPartT.String = num2str(cpart);
        for k = 1 : nParts
            seq = dt.Seq;
            if k == floor(currentPart)
                nseq{k}.String = seq{k};
                set(nseq{k},'BackgroundColor',[1 0 0]);
            else
                nseq{k}.String = seq{k};
                set(nseq{k},'BackgroundColor',[1 1 1]);
            end
            nseq{k}.Units = 'normalized';
        end
        
        t = ((1+ceil(zPlotMin*60)):(floor(zPlotMax*60)-1));
        
        Vrout = dt.Data{currentPart}.Vr(t);
        Vfout = dt.Data{currentPart}.Vf(t);
        Vsout = dt.Data{currentPart}.Vs(t);
        Vtout = dt.Data{currentPart}.Vt(t);
        Xout = dt.Data{currentPart}.X(t);
        Yout = dt.Data{currentPart}.Y(t);
        flp = dt.Data{currentPart}.flp(t);
        actstout = dt.Data{currentPart}.actState(t);

        axes(hb);
        cla reset;
        hold on
        an = -0.1:0.1:6.3;
        plot(40*sin(an),40*cos(an), 'k') % Plot outline of the arena
        plot(Xout, Yout,'k', 'linewidth', 1.5) % Plot fly trajectory
        axis square
        sz = 40;
        axis([-sz sz -sz sz])
        axes(ha);
        cla reset;

        plot(t/60, zeros(length(t),1), 'color',[0.8 0.8 0.8], 'linewidth', 2)
        hold on
        plot(t/60,Vrout, 'k', 'linewidth', 1); % Plot fly rotation
        hold on
        plot(t/60, zeros(length(t),1)-1300, 'color',[0.8 0.8 0.8], 'linewidth', 1)
        hold on
        plot(t/60, 10*Vfout-1300, 'k', 'linewidth', 1) % Plot fly forward
        hold on
        plot(t/69, 1*flp+1200, 'k', 'linewidth',0.5)
        plot(t/60, (10*22)*ones(length(t),1)+1000, 'color', [0.8 0.8 0.8])
        hold on
        plot(t/60, 10*dt.Data{currentPart}.WallDist(t)+1000, 'k', 'linewidth', 2) % Plot fly distance to wall
        xlabel('Time (s)')
        axis([t(1)/60 t(end)/60 -1500 1500])
    end

end

