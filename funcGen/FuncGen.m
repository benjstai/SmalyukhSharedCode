%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphical User Interface
% 
% Things to do:

function FuncGen
% [GUI] Brief description:
%   Run to generate the Graphical User Interface (GUI).
%   
% setappdata(hgui,'fieldname',value) creates or replaces a structure field
% data = getappdata(hgui) retrieves the GUI information structure
% field = getappdata(hgui,'fieldname') retrieves a specific field
%   
daqreset

% Variables needed in creation of GUI
gui_vars.dev_popupmenu = {'device 1','device 2','device 3','device 4',};
gui_vars.cha_popupmenu = {'channel 0','channel 1','channel 2','channel 3'};
gui_vars.car_popupmenu = {'Square wave','Sine wave','test'};
gui_vars.mod_popupmenu = {'Square wave','Sine wave'};
gui_vars.defalt_parameters = [0,0,60,2,50,0];


% enforces singleton (only one GUI at a time)
hgui = findall(0, 'Name', mfilename);
if isempty(hgui) % if no GUI existing, run make_fig_fcn
    % Create and hide the GUI figure as it is being constructed.
    hgui = make_fig_fcn(gui_vars);
end
connect_button = findobj('Tag','connect_button');
set(connect_button,'Interruptible','off')
start_button = findobj('Tag','start_button');
set(start_button,'BusyAction','cancel')

% Variables sent to app data container
setappdata(hgui,'gui_vars',gui_vars)
setappdata(hgui,'update',0)
setappdata(hgui,'amplitude_edit',gui_vars.defalt_parameters(1));
setappdata(hgui,'phase_edit',gui_vars.defalt_parameters(2));
setappdata(hgui,'carfreq_edit',gui_vars.defalt_parameters(3));
setappdata(hgui,'modfreq_edit',gui_vars.defalt_parameters(4));
setappdata(hgui,'duty_edit',gui_vars.defalt_parameters(5));
setappdata(hgui,'fill_edit',gui_vars.defalt_parameters(6));
setappdata(hgui,'carwave',1)
setappdata(hgui,'modwave',1)

% raise GUI and make visible to user
figure(hgui)
hgui.Visible = 'on';


%% %%%%%%%%%%%%%%%%%%%% Utility functions for GUI %%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function signal=dynamic_data(session,~)
hgui = findall(0, 'Name', mfilename);

    amp = getappdata(hgui,'amplitude_edit');
    pha = getappdata(hgui,'phase_edit')*pi/180;
    f_c = getappdata(hgui,'carfreq_edit');
    f_m = getappdata(hgui,'modfreq_edit');
    dc = getappdata(hgui,'duty_edit')/100;
    fill = getappdata(hgui,'fill_edit')/100;

    switch getappdata(hgui,'carwave')
        case 1
            % Sin wave with f_c frequency (one second)
            sps = floor(session.Rate); % number of samples per second
            data = sin(linspace(0,2*pi*f_c,sps+1) + pha)';
            data(end) = [];

            % Square wave with f_c frequency (one second)
            data(data>=0) = 1;
            data(data<0) = -1;
        case 2
            % Sin wave with f_c frequency (one second)
            sps = floor(session.Rate); % number of samples per second
            data = sin(linspace(0,2*pi*f_c,sps+1) + pha)';
            data(end) = [];
        case 3
            % test dc... Must use DC coupling!!!!!!!!!
            sps = floor(session.Rate);
            data = ones(sps,1);
            data(end) = [];
    end
    
    switch getappdata(hgui,'modwave')
        case 1 % modulate with a square wave
            % Generate one period of the modulated signal
            data = repmat(data,ceil(2/f_m),1); % make data long enough according to f_m
            signal = data(1:round(sps/f_m)); % one period of the modulated signal
            signal(1:round(end*(1-dc))) = signal(1:round(end*(1-dc)))*fill;
        case 2 % modulate with a sine wave
            % Generate one period of the modulated signal
            data = repmat(data,ceil(2/f_m),1); % make data long enough according to f_m
            signal = data(1:round(sps/f_m)); % one period of the modulated signal
            vec = linspace(0,2*pi-2*pi/length(signal),length(signal))'-pi/2;
            signal = signal.*((sin(vec)+1)/2+(cos(vec+pi/2)+1)/2*fill);
    end
    
    % Generate more than one second of output for the listener
    times = 1;
    while session.Rate>length(signal)*times
        times = times+1;
    end
    signal = abs(amp)*repmat(signal,times,1);
    session.queueOutputData(signal)
    
    h_ax = findobj(hgui,'Tag','axis1');
    cla(h_ax) % -------------------------------------------------------------------------------
    plot(linspace(0,times/f_m,length(signal)),signal) % -------------------------------------------------------------------------------
    grid on



%% %%%%%%%%%%%%%%%%%%%% Callback functions for GUI %%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These functions are called when the user interacts with the GUI

function amplitude_edit(~,~,hgui)

disp('amplitude_edit')

setappdata(hgui,'update',1)
h_start = findobj('Tag','start_button');
set(h_start,'String','Update')

function phase_edit(~,~,hgui)

disp('phase_edit')

setappdata(hgui,'update',1)
h_start = findobj('Tag','start_button');
set(h_start,'String','Update')

function carfreq_edit(~,~,hgui)

disp('carfreq_edit')

setappdata(hgui,'update',1)
h_start = findobj('Tag','start_button');
set(h_start,'String','Update')

function modfreq_edit(~,~,hgui)

disp('modfreq_edit')

setappdata(hgui,'update',1)
h_start = findobj('Tag','start_button');
set(h_start,'String','Update')

function duty_edit(~,~,hgui)

disp('duty_edit')

setappdata(hgui,'update',1)
h_start = findobj('Tag','start_button');
set(h_start,'String','Update')

function fill_edit(~,~,hgui)

disp('fill_edit')

setappdata(hgui,'update',1)
h_start = findobj('Tag','start_button');
set(h_start,'String','Update')



function dev_popupmenu(~,~)
disp('dev_popupmenu')

function cha_popupmenu(~,~)
disp('cha_popupmenu')

function car_popupmenu(hpop,~,hgui)
setappdata(hgui,'carwave',get(hpop,'Value'))
disp('car_popupmenu')

function mod_popupmenu(hpop,~,hgui)
setappdata(hgui,'modwave',get(hpop,'Value'))
disp('mod_popupmenu')

function connect_button(h_button,~,hgui)
set(h_button,'BackgroundColor',[.7,.7,.7])
h_dev = findobj('Tag','dev_popupmenu');
h_cha = findobj('Tag','cha_popupmenu');
dev_num = get(h_dev,'Value');
cha_num = get(h_cha,'Value');
    try
        devices = daq.getDevices;
        if isempty(devices)
            errordlg('No data acquisition devices available.')
            set(h_button,'BackgroundColor',[.7,.7,.7])
            return
        end
    catch
        errordlg('daq.getDevices (ERROR)')
        set(h_button,'BackgroundColor',[.7,.7,.7])
        return
    end
    try
        dev_id = devices(dev_num).ID;
        setappdata(hgui,'dev_id',dev_id)
    catch
        errordlg('Device ID not valid')
        set(h_button,'BackgroundColor',[.7,.7,.7])
        return
    end
    try
        ao = devices(dev_num).Subsystems(2);
        if ~strcmp(devices(dev_num).Subsystems(2).SubsystemType,'AnalogOutput')
            errordlg('NOT an Analog Output')
            set(h_button,'BackgroundColor',[.7,.7,.7])
            return
        end
    catch
        errordlg('Analog Output not valid')
        set(h_button,'BackgroundColor',[.7,.7,.7])
        return
    end
    try
        channel_name = ao.ChannelNames{cha_num};
        setappdata(hgui,'channel_name',channel_name)
    catch
        errordlg('Channel not valid')
        set(h_button,'BackgroundColor',[.7,.7,.7])
        return
    end
    
    % Create data acquisition session for specific vendor hardware
    session = daq.createSession('ni');
    
    % Add continuous Analog Output channel with listener
    addAnalogOutputChannel(session,dev_id,channel_name,'Voltage');    
    addlistener(session,'DataRequired',@dynamic_data);
    session.IsContinuous = true;
    
    max_rate = session.RateLimit(2);
    setappdata(hgui,'max_rate',max_rate)
    session.Rate = floor(max_rate);
    
setappdata(hgui,'session',session)
set(h_button,'BackgroundColor',[.7,.93,.7])
disp('connect')


function start_button(~,~,hgui)
    disp('start')
    session = getappdata(hgui,'session');
    if isempty(session)
        errordlg('Device not connected!')
        return
    end
    
    % if session is running warn...
    if session.IsRunning && ~getappdata(hgui,'update')
        disp('Session is running.')
    else
        if getappdata(hgui,'update')
            % get updated parameters and put in app container
            setappdata(hgui,'amplitude_edit',str2double(get(findobj('Tag','amplitude'),'String')))
            setappdata(hgui,'phase_edit',str2double(get(findobj('Tag','phase'),'String')))
            setappdata(hgui,'carfreq_edit',str2double(get(findobj('Tag','carfreq'),'String')))
            setappdata(hgui,'modfreq_edit',str2double(get(findobj('Tag','modfreq'),'String')))
            setappdata(hgui,'duty_edit',str2double(get(findobj('Tag','duty'),'String')))
            setappdata(hgui,'fill_edit',str2double(get(findobj('Tag','fill'),'String')))
            % reset start button
            setappdata(hgui,'update',0)
            h_start = findobj('Tag','start_button');
            set(h_start,'String','Start')
        else
            queueOutputData(session,dynamic_data(session, 0))
            session.startBackground()
        end
    end



function stop_button(~,~,hgui)
    disp('stop')
    session = getappdata(hgui,'session');
    if isempty(session)
        return
    end
    session.stop()
    h_ax = findobj(hgui,'Tag','axis1');
    cla(h_ax) % -------------------------------------------------------------------------------
    plot([0 1],[0 0])% -------------------------------------------------------------------------------
    grid on
    try
        % Generate a Single Scan to "Zero" the voltage
        outputSingleScan(session,0);
    catch
        return
    end



function quit_button(~,~,hgui)
user_closereq(hgui)


%% %%%%%%%%%%%%%%%%%%%%%%%%% Generate the GUI %%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These functions generate the GUI when the function is first run

function hgui = make_fig_fcn(gui_vars)
% make_fig_fcn Brief description:
%   Generates the GUI axis and the GUI organization. Initializes the data
%   structure in the aplication data. Outputs the GUI handle

% figure and figure properties
hgui = figure('Visible','off'); % generate an invisible figure
hgui.Name = mfilename; % set the GUI name to file name
hgui.IntegerHandle = 'off'; % turn off integer handle to prevent overwrite
hgui.Color = [0.87 0.9 0.93]; % background color
hgui.DockControls = 'off'; % prevents GUI from "docking" into matlab
hgui.MenuBar = 'none';
hgui.ToolBar = 'none';
hgui.Units = 'normalized';
hgui.OuterPosition = [0.5 0.1 0.5 0.85];
hgui.Resize = 'on';

elmt = 0;
% elements with callback functions

    % Device listbox
    elmt = elmt+1;
    ui(elmt).tag = 'dev_popupmenu';
    ui(elmt).style = 'popupmenu';
    ui(elmt).str = gui_vars.dev_popupmenu;
    ui(elmt).fontsize = 12;
    ui(elmt).bgc = [1 1 1];
    ui(elmt).rcwh = [.05 .85 .1 .1];
    ui(elmt).callback = {str2func('dev_popupmenu')};
    ui(elmt).value = 1;
    
    % Channel listbox
    elmt = elmt+1;
    ui(elmt).tag = 'cha_popupmenu';
    ui(elmt).style = 'popupmenu';
    ui(elmt).str = gui_vars.cha_popupmenu;
    ui(elmt).fontsize = 12;
    ui(elmt).bgc = [1 1 1];
    ui(elmt).rcwh = [.05 .80 .1 .1];
    ui(elmt).callback = {str2func('cha_popupmenu')};
    ui(elmt).value = 1;
    
    % Carrier waveform listbox
    elmt = elmt+1;
    ui(elmt).tag = 'car_popupmenu';
    ui(elmt).style = 'popupmenu';
    ui(elmt).str = gui_vars.car_popupmenu;
    ui(elmt).fontsize = 12;
    ui(elmt).bgc = [1 1 1];
    ui(elmt).rcwh = [.37 .85 .2 .1];
    ui(elmt).callback = {str2func('car_popupmenu'),hgui};
    ui(elmt).value = 1;
    
    % Modulation waveform listbox
    elmt = elmt+1;
    ui(elmt).tag = 'mod_popupmenu';
    ui(elmt).style = 'popupmenu';
    ui(elmt).str = gui_vars.mod_popupmenu;
    ui(elmt).fontsize = 12;
    ui(elmt).bgc = [1 1 1];
    ui(elmt).rcwh = [.37 .80 .2 .1];
    ui(elmt).callback = {str2func('mod_popupmenu'),hgui};
    ui(elmt).value = 1;

    % Connect button
    elmt = elmt+1;
    ui(elmt).tag = 'connect_button';
    ui(elmt).style = 'pushbutton';
    ui(elmt).str = 'Connect';
    ui(elmt).fontsize = 16;
    ui(elmt).bgc = [0.7 0.7 0.7];
    ui(elmt).rcwh = [.17 .87 .13 .08];
    ui(elmt).callback = {str2func('connect_button'),hgui};
    ui(elmt).value = [];
    
    % START button
    elmt = elmt+1;
    ui(elmt).tag = 'start_button';
    ui(elmt).style = 'pushbutton';
    ui(elmt).str = 'START';
    ui(elmt).fontsize = 30;
    ui(elmt).bgc = [0.3 0.93 0.3];
    ui(elmt).rcwh = [.05 .07 .2 .15];
    ui(elmt).callback = {str2func('start_button'),hgui};
    ui(elmt).value = [];

    % STOP button
    elmt = elmt+1;
    ui(elmt).tag = 'stop_button';
    ui(elmt).style = 'pushbutton';
    ui(elmt).str = 'STOP';
    ui(elmt).fontsize = 30;
    ui(elmt).bgc = [0.3 0.3 0.93];
    ui(elmt).rcwh = [.4 .07 .2 .15];
    ui(elmt).callback = {str2func('stop_button'),hgui};
    ui(elmt).value = [];

    % QUIT button
    elmt = elmt+1;
    ui(elmt).tag = 'quit_button';
    ui(elmt).style = 'pushbutton';
    ui(elmt).str = 'QUIT';
    ui(elmt).fontsize = 30;
    ui(elmt).bgc = [0.93 0.3 0.3];
    ui(elmt).rcwh = [.75 .07 .2 .15];
    ui(elmt).callback = {str2func('quit_button'),hgui};
    ui(elmt).value = [];
    
    % Amplitude
    elmt = elmt+1;
    ui(elmt).tag = 'amplitude';
    ui(elmt).style = 'edit';
    ui(elmt).str = num2str(gui_vars.defalt_parameters(1));
    ui(elmt).fontsize = 12;
    ui(elmt).bgc = [1 1 1];
    ui(elmt).rcwh = [.04 .76 .06 .03];
    ui(elmt).callback = {str2func('amplitude_edit'),hgui};
    ui(elmt).value = [];

    % Phase
    elmt = elmt+1;
    ui(elmt).tag = 'phase';
    ui(elmt).style = 'edit';
    ui(elmt).str = num2str(gui_vars.defalt_parameters(2));
    ui(elmt).fontsize = 12;
    ui(elmt).bgc = [1 1 1];
    ui(elmt).rcwh = [.04 .72 .06 .03];
    ui(elmt).callback = {str2func('phase_edit'),hgui};
    ui(elmt).value = [];
    
    % Carrier Frequency
    elmt = elmt+1;
    ui(elmt).tag = 'carfreq';
    ui(elmt).style = 'edit';
    ui(elmt).str = num2str(gui_vars.defalt_parameters(3));
    ui(elmt).fontsize = 12;
    ui(elmt).bgc = [1 1 1];
    ui(elmt).rcwh = [.04 .68 .06 .03];
    ui(elmt).callback = {str2func('carfreq_edit'),hgui};
    ui(elmt).value = [];
    
    % Modulation Frequency
    elmt = elmt+1;
    ui(elmt).tag = 'modfreq';
    ui(elmt).style = 'edit';
    ui(elmt).str = num2str(gui_vars.defalt_parameters(4));
    ui(elmt).fontsize = 12;
    ui(elmt).bgc = [1 1 1];
    ui(elmt).rcwh = [.04 .64 .06 .03];
    ui(elmt).callback = {str2func('modfreq_edit'),hgui};
    ui(elmt).value = [];
    
    % Modulation Duty Cycle
    elmt = elmt+1;
    ui(elmt).tag = 'duty';
    ui(elmt).style = 'edit';
    ui(elmt).str = num2str(gui_vars.defalt_parameters(5));
    ui(elmt).fontsize = 12;
    ui(elmt).bgc = [1 1 1];
    ui(elmt).rcwh = [.04 .60 .06 .03];
    ui(elmt).callback = {str2func('duty_edit'),hgui};
    ui(elmt).value = [];
    
    % Modulation Fill Percent
    elmt = elmt+1;
    ui(elmt).tag = 'fill';
    ui(elmt).style = 'edit';
    ui(elmt).str = num2str(gui_vars.defalt_parameters(6));
    ui(elmt).fontsize = 12;
    ui(elmt).bgc = [1 1 1];
    ui(elmt).rcwh = [.04 .56 .06 .03];
    ui(elmt).callback = {str2func('fill_edit'),hgui};
    ui(elmt).value = [];
    
% elements with no user interaction
    
    % Amplitude text box
    elmt = elmt+1;
    ui(elmt).tag = 'amplitude_text';
    ui(elmt).style = 'text';
    ui(elmt).str = 'Amplitude [Vrms]              ';
    ui(elmt).fontsize = 12;
    ui(elmt).bgc = hgui.Color;
    ui(elmt).rcwh = [.1 .76 .2 .03];
    ui(elmt).value = [];

    % Phase text box
    elmt = elmt+1;
    ui(elmt).tag = 'phase_text';
    ui(elmt).style = 'text';
    ui(elmt).str = 'Phase [deg]                       ';
    ui(elmt).fontsize = 12;
    ui(elmt).bgc = hgui.Color;
    ui(elmt).rcwh = [.1 .72 .2 .03];
    ui(elmt).value = [];

    % Carrier Fequency text box
    elmt = elmt+1;
    ui(elmt).tag = 'carfreq_text';
    ui(elmt).style = 'text';
    ui(elmt).str = 'Carrier Fequency [Hz]      ';
    ui(elmt).fontsize = 12;
    ui(elmt).bgc = hgui.Color;
    ui(elmt).rcwh = [.1 .68 .2 .03];
    ui(elmt).value = [];

    % Modulation Fequency text box
    elmt = elmt+1;
    ui(elmt).tag = 'modfreq_text';
    ui(elmt).style = 'text';
    ui(elmt).str = 'Modulation Fequency [Hz]';
    ui(elmt).fontsize = 12;
    ui(elmt).bgc = hgui.Color;
    ui(elmt).rcwh = [.1 .64 .2 .03];
    ui(elmt).value = [];

    % Duty Cycle text box
    elmt = elmt+1;
    ui(elmt).tag = 'dc_text';
    ui(elmt).style = 'text';
    ui(elmt).str = 'Duty Cycle [%]                   ';
    ui(elmt).fontsize = 12;
    ui(elmt).bgc = hgui.Color;
    ui(elmt).rcwh = [.1 .60 .2 .03];
    ui(elmt).value = [];

    % Fill text box
    elmt = elmt+1;
    ui(elmt).tag = 'fill_text';
    ui(elmt).style = 'text';
    ui(elmt).str = 'Fill [%]                                ';
    ui(elmt).fontsize = 12;
    ui(elmt).bgc = hgui.Color;
    ui(elmt).rcwh = [.1 .56 .2 .03];
    ui(elmt).value = [];

    % Carrier waveform text box
    elmt = elmt+1;
    ui(elmt).tag = 'carwave_text';
    ui(elmt).style = 'text';
    ui(elmt).str = 'Carrier      ';
    ui(elmt).fontsize = 12;
    ui(elmt).bgc = hgui.Color;
    ui(elmt).rcwh = [.58 .92 .1 .03]; %[.37 .80 .2 .1]
    ui(elmt).value = [];

    % Modulation waveform text box
    elmt = elmt+1;
    ui(elmt).tag = 'modwave_text';
    ui(elmt).style = 'text';
    ui(elmt).str = 'Modulation';
    ui(elmt).fontsize = 12;
    ui(elmt).bgc = hgui.Color;
    ui(elmt).rcwh = [.58 .87 .1 .03];
    ui(elmt).value = [];

make_uicontrols(hgui,ui)

function make_uicontrols(hgui,ui)
% [make_uicontrols] Brief description:
%   Generates the GUI UIcontrols from make_fig_fcn organization

% Create button groups
uibuttongroup('BackgroundColor',hgui.Color,'Position',[.02 .85 .3 .15],...
    'Title','Connect to DAQ','FontSize',16);

uibuttongroup('BackgroundColor',hgui.Color,'Position',[.34 .85 .35 .15],...
    'Title','Waveforms','FontSize',16);

uibuttongroup('BackgroundColor',hgui.Color,'Position',[.02 .54 .3 .3],...
    'Title','Parameters','FontSize',16);

% Create axis
axis1 = axes('Parent',hgui,...
    'OuterPosition',[0 .25 1 .3],...
    'tag','axis1',...
    'box','on');
setappdata(hgui,'axis1',axis1) % axis1 information

for n=1:length(ui)
    pos = [ui(n).rcwh(1) ui(n).rcwh(2) ui(n).rcwh(3) ui(n).rcwh(4)];
    uicontrol(hgui,...
        'Units','normalized',...
        'FontSize',ui(n).fontsize,...
        'tag',ui(n).tag,...
        'Style',ui(n).style,...
        'String',ui(n).str,...
        'BackgroundColor',ui(n).bgc,...
        'OuterPosition',pos,...
        'Callback',ui(n).callback,...
        'Value',ui(n).value);
end
set(hgui,'CloseRequestFcn',@user_closereq)

function user_closereq(hgui,~)
% [user_closereq] Brief description:
%   Close request function, runs when GUI is closed.
%   questdlg options: 'Yes', 'No'. 
stop_button(0,0,hgui)

selection = questdlg('QUIT?',...
  'Close Request',...
  'Yes','No','Yes');

switch selection 
  case 'Yes' %%%%%%%%%%%%%%%%% Quit
      delete(hgui) % close GUI
  case 'No' %%%%%%%%%%%%%%%%%% cancel 
  return 
end






