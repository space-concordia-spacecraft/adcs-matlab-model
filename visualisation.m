function visualisation(block)
%MSFUNTMPL_BASIC A Template for a Level-2 MATLAB S-Function
%   The MATLAB S-function is written as a MATLAB function with the
%   same name as the S-function. Replace 'msfuntmpl_basic' with the 
%   name of your S-function.
%
%   It should be noted that the MATLAB S-function is very similar
%   to Level-2 C-Mex S-functions. You should be able to get more
%   information for each of the block methods by referring to the
%   documentation for C-Mex S-functions.
%
%   Copyright 2003-2010 The MathWorks, Inc.

%%
%% The setup method is used to set up the basic attributes of the
%% S-function such as ports, parameters, etc. Do not add any other
%% calls to the main body of the function.
%%
setup(block);

%endfunction

%% Function: setup ===================================================
%% Abstract:
%%   Set up the basic characteristics of the S-function block such as:
%%   - Input ports
%%   - Output ports
%%   - Dialog parameters
%%   - Options
%%
%%   Required         : Yes
%%   C-Mex counterpart: mdlInitializeSizes
%%
function setup(block)

% Register number of ports
block.NumInputPorts  = 16;
block.NumOutputPorts = 0;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
block.InputPort(1).Dimensions        = 1;
block.InputPort(1).DatatypeID  = 0;  % double
block.InputPort(1).Complexity  = 'Real';
block.InputPort(1).DirectFeedthrough = true;

block.InputPort(2).Dimensions        = 3;
block.InputPort(2).DatatypeID  = 0;  % double
block.InputPort(2).Complexity  = 'Real';
block.InputPort(2).DirectFeedthrough = true;

block.InputPort(3).Dimensions        = 3;

block.InputPort(4).Dimensions        = 3;
block.InputPort(4).DatatypeID  = 0;  % double
block.InputPort(4).Complexity  = 'Real';
block.InputPort(4).DirectFeedthrough = true;

block.InputPort(5).Dimensions        = 3;
block.InputPort(6).Dimensions        = 3;
block.InputPort(7).Dimensions        = 3;

block.InputPort(8).Dimensions        = 1;
block.InputPort(8).DatatypeID  = 0;  % double
block.InputPort(8).Complexity  = 'Real';
block.InputPort(8).DirectFeedthrough = true;

block.InputPort(9).Dimensions        = 3;

block.InputPort(10).Dimensions        = 1;
block.InputPort(10).DatatypeID  = 0;  % double
block.InputPort(10).Complexity  = 'Real';
block.InputPort(10).DirectFeedthrough = true;

block.InputPort(11).Dimensions        = 3;
block.InputPort(12).Dimensions        = 3;
block.InputPort(13).Dimensions        = 3;
block.InputPort(14).Dimensions        = 3;

block.InputPort(15).Dimensions        = 2;
block.InputPort(15).DatatypeID  = 0;  % double
block.InputPort(15).Complexity  = 'Real';
block.InputPort(15).DirectFeedthrough = true;
block.InputPort(16).Dimensions        = 1;
% Override output port properties
% block.OutputPort(1).Dimensions       = 1;
% block.OutputPort(1).DatatypeID  = 0; % double
% block.OutputPort(1).Complexity  = 'Real';

% Register parameters
block.NumDialogPrms     = 1;

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [-1, 0];%[1/block.DialogPrm(1).Data 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'CustomSimState',  < Has GetSimState and SetSimState methods
%    'DisallowSimState' < Error out when saving or restoring the model sim state
block.SimStateCompliance = 'DefaultSimState';

%% -----------------------------------------------------------------
%% The MATLAB S-function uses an internal registry for all
%% block methods. You should register all relevant methods
%% (optional and required) as illustrated below. You may choose
%% any suitable name for the methods and implement these methods
%% as local functions within the same file. See comments
%% provided for each function for more information.
%% -----------------------------------------------------------------

block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
block.RegBlockMethod('InitializeConditions', @InitializeConditions);
block.RegBlockMethod('Start', @Start);
block.RegBlockMethod('Outputs', @Outputs);     % Required
block.RegBlockMethod('Update', @Update);
block.RegBlockMethod('Derivatives', @Derivatives);
block.RegBlockMethod('Terminate', @Terminate); % Required

%end setup

%%
%% PostPropagationSetup:
%%   Functionality    : Setup work areas and state variables. Can
%%                      also register run-time methods here
%%   Required         : No
%%   C-Mex counterpart: mdlSetWorkWidths
%%
function DoPostPropSetup(block)
block.NumDworks = 1;
  
  block.Dwork(1).Name            = 'flag';
  block.Dwork(1).Dimensions      = 1;
  block.Dwork(1).DatatypeID      = 0; 
  block.Dwork(1).Complexity      = 'Real';
  block.Dwork(1).UsedAsDiscState = false;


%%
%% InitializeConditions:
%%   Functionality    : Called at the start of simulation and if it is 
%%                      present in an enabled subsystem configured to reset 
%%                      states, it will be called when the enabled subsystem
%%                      restarts execution to reset the states.
%%   Required         : No
%%   C-MEX counterpart: mdlInitializeConditions
%%
function InitializeConditions(block)

%end InitializeConditions


%%
%% Start:
%%   Functionality    : Called once at start of model execution. If you
%%                      have states that should be initialized once, this 
%%                      is the place to do it.
%%   Required         : No
%%   C-MEX counterpart: mdlStart
%%
function Start(block)
%% VISUALISATION (initialization)

 %% 3d Figure
%add the satelitte
figure('Name','Visualisation of ADCS','NumberTitle','off')
subplot(3,3,[1,6])
 [F, V, C] = rndread('mariofirstcad.slp');
 rotate3d; %clf; %opens figure area
 sat(1) = patch('faces', F, 'vertices' ,V-[200 150 200],'facec', 'flat','FaceVertexCData', C,'EdgeColor','none');
 rotate(sat(1), [0 1 -1], 180);
%sat(1)=patch(0,0,0,0);%replace sat(1) to remove the cad     

%set the axis
% axis manual %prevent automatic resizing of the plot
xlim([-175 175])
ylim([-175 175])
zlim([-175 175])
set(gcf,'Position',get(gcf,'Position').*[0.5 0.2 1.7 1.7],'Resize','Off'); %prevent automatic resizing of the window
set(gca,'Position',get(gca,'Position').*[0.65 1.03 1 1],'visible','off'); %moves the graph in the window and removes the background
light('Position', [1 0 0],'Visible','on','tag','light_sun'); % add a default light
daspect([1 1 1])                    % Setting the aspect ratio
view([1,1,1]);
camproj('perspective')
camva('manual')
camva(10);

% draw the three coloured frame vectors
sat(2) = line([0 150],[0,0],[0,0],'color',[1,0,0]);
sat(3) = line([0 0],[0,150],[0,0],'color',[0,1,0]);
sat(4) = line([0 0],[0,0],[0,230],'color',[0,0,1]);
set([sat(2),sat(3),sat(4)],'linewidth',3)

% add text for x,y and z for each frame vector
sat(5) = text('position',[170 0 0],'string','x','fontw','b');
sat(6) = text('position',[ 0 170 0],'string','y','fontw','b');
sat(7) = text('position',[0 0 250],'string','z','fontw','b');

%combine the satellite and its vectors to 1 object easy to rotate (hgtransform)
set(sat,'Clipping','off');
ax = gca;
combinedobject = hgtransform('parent',ax);
set(sat, 'parent', combinedobject)
set(combinedobject,'tag','satfigure')

% draw the three axis vectors for reference (note the origin offset)
line([-250,250],[0,0],[0,0],'color',[0,0,0],'linewidth',1,'Clipping','off');
line([0,0],[-250,250],[0,0],'color',[0,0,0],'linewidth',1,'Clipping','off');
line([0,0],[0,0],[-250,250],'color',[0,0,0],'linewidth',1,'Clipping','off');
	
% add text for x0,y0 and z0 for each base frame line
text('position',[275 0 0],'string','Velocity','fontw','b','HorizontalAlignment', 'center');
text('position',[0 275 0],'string','~Normal','fontw','b','HorizontalAlignment', 'center');
text('position',[0 0 275],'string','~Position','fontw','b','HorizontalAlignment', 'center');

%% Draw Vectors
Bline=line([0 150],[0,0],[0,0], 'color',[1 0 1],'linewidth',3);
set(Bline,'Clipping','off');
set(Bline,'tag','satfigureBline')
Bline_text = text('position',[170 0 0],'string','B','fontw','b'); %add letter
set(Bline_text,'tag','satfigureBlinetext')


%Draw Sun Vector
Sline=line([0 150],[0,0],[0,0], 'color',[1 1 0],'linewidth',3);
set(Sline,'Clipping','off');
set(Sline,'tag','satfigureSline')
Sline_text = text('position',[170 0 0],'string','S','fontw','b'); %add letter
set(Sline_text,'tag','satfigureSlinetext')


%Draw ground station line
Gline=line([0 150],[0,0],[0,0], 'color',[1 0.5 0],'linewidth',3);
set(Gline,'Clipping','off');
set(Gline,'tag','satfigureGline')
Gline_text = text('position',[170 0 0],'string','G','fontw','b'); %add letter
set(Gline_text,'tag','satfigureGlinetext')

%Draw Argentina line
Aline=line([0 150],[0,0],[0,0], 'color',[0 0.5 0.5],'linewidth',3);
set(Aline,'Clipping','off');
set(Aline,'tag','satfigureAline')
Aline_text = text('position',[170 0 0],'string','A','fontw','b'); %add letter
set(Aline_text,'tag','satfigureAlinetext')

%Draw Namibia line
Nline=line([0 150],[0,0],[0,0], 'color',[0 1 1],'linewidth',3);
set(Nline,'Clipping','off');
set(Nline,'tag','satfigureNline')
Nline_text = text('position',[170 0 0],'string','N','fontw','b'); %add letter
set(Nline_text,'tag','satfigureNlinetext')

%% add UICONTROL
uicontrol('Style','popupmenu','Position',[850 420+250 100 35],'String',{'Auto','Sun Pointing','Ground','Imaging-A','Imaging-N','Nadir Pointing','Detumbling','OFF','Calibration-M'},'Callback',@(src,evt)modeoverdrive, 'tag', 'modeOverdrive');
uicontrol('Position',[890, 55, 60, 20],'String','Minimap','Value',1,'Callback',@(src,evt,lat,lon) Minimap);
uicontrol('Position',[890, 30, 60, 20],'String','Pause','Value',1,'Callback',@(src,evt) Buttonpause, 'tag','Button');
uicontrol('Position',[890, 5, 60, 20],'String','Stop','Value',1,'Callback',@(src,evt) Buttonstop);
uicontrol('Style','edit','Position',[35, 25, 60, 20],'String','1','Callback',@(src,evt) Sampletime, 'tag','Sampletimebutton');
uicontrol('Position',[7, 25, 20, 20],'String','-','Callback',@(src,evt) Sampletimedown);
uicontrol('Position',[103, 25, 20, 20],'String','+','Callback',@(src,evt) Sampletimeup);
uicontrol('Style','edit','Position',[35, 70, 60, 20],'String','inf','Callback',@(src,evt) Pacer, 'tag','Pacerbutton');
uicontrol('Position',[7, 70, 20, 20],'String','-','Callback',@(src,evt) Pacerdown);
uicontrol('Position',[103, 70, 20, 20],'String','+','Callback',@(src,evt) Pacerup);
uicontrol('Style','togglebutton','Position',[7, 47.5, 20, 20],'String',compose("\xD83D\xDD13"),'Callback',@(src,evt) Lock, 'tag','Lockbutton');

%% add labels
uicontrol('Style','text','String','ADCS Mode:','Position',[10 425+250 200 35],'FontSize',25,'ForegroundColor',[0 0 1],'HorizontalAlignment', 'left')

uicontrol('Style','text','String','Euler Angles to Target:','Position',[10 370+250 290 35],'FontSize',20,'ForegroundColor',[1 0 0],'HorizontalAlignment', 'left')
uicontrol('Style','text','String',[compose("\x03D5") compose("\x03B8") compose("\x03C8")],'Position',[15 270+250 100 100],'FontSize',20,'ForegroundColor',[1 0 0],'HorizontalAlignment', 'left')
uicontrol('Style','text','String',['=' newline '=' newline '='],'Position',[40 270+250 50 100],'FontSize',20,'ForegroundColor',[1 0 0],'HorizontalAlignment', 'left')

uicontrol('Style','text','String','Reaction Wheels Usage:','Position',[645 370+250 390 35],'FontSize',20,'ForegroundColor',[0 0.5 1],'HorizontalAlignment', 'left')
uicontrol('Style','text','String',['w   =' newline 'w   =' newline 'w   ='],'Position',[650 270+250 100 100],'FontSize',20,'ForegroundColor',[0 0.5 1],'HorizontalAlignment', 'left')
uicontrol('Style','text','String',['x' newline 'y' newline 'z'],'Position',[675 260+250 20 100],'FontSize',20,'ForegroundColor',[0 0.5 1],'HorizontalAlignment', 'left')

uicontrol('Style','text','String','Angular Rates:','Position',[10 210+250 200 35],'FontSize',20,'ForegroundColor',[0 0.8 0],'HorizontalAlignment', 'left')
uicontrol('Style','text','String',[char(969) '   =' newline char(969) '   =' newline char(969) '   ='],'Position',[10 110+250 100 100],'FontSize',20,'ForegroundColor',[0 0.8 0],'HorizontalAlignment', 'left')
uicontrol('Style','text','String',['x' newline 'y' newline 'z'],'Position',[35 100+250 20 100],'FontSize',20,'ForegroundColor',[0 0.8 0],'HorizontalAlignment', 'left')

uicontrol('Style','text','String','Magnetorquers Usage:','Position',[655 210+250 300 35],'FontSize',20,'ForegroundColor',[1 0.7 0],'HorizontalAlignment', 'left')
uicontrol('Style','text','String',['M   =' newline 'M   =' newline 'M   ='],'Position',[650 110+250 100 100],'FontSize',20,'ForegroundColor',[1 0.7 0],'HorizontalAlignment', 'left')
uicontrol('Style','text','String',['x' newline 'y' newline 'z'],'Position',[675 100+250 20 100],'FontSize',20,'ForegroundColor',[1 0.7 0],'HorizontalAlignment', 'left')

uicontrol('Style','text','String','Momentum:','Position',[10 50+250 200 35],'FontSize',20,'ForegroundColor',[1 0 1],'HorizontalAlignment', 'left')

uicontrol('Style','text','String','Time:','Position',[430 50+250 200 35],'FontSize',20,'ForegroundColor',[1 0.5 0],'HorizontalAlignment', 'left')

uicontrol('Style','text','String','Overdrive:','Position',[715 425+250 130 35],'FontSize',20,'ForegroundColor',[0 0.5 0.5],'HorizontalAlignment', 'left')

uicontrol('Style','text','String','Timestep:','Position',[35 45 60 20],'FontSize',10,'ForegroundColor',[0 0 0],'HorizontalAlignment', 'left')
uicontrol('Style','text','String','RT Pacer:','Position',[35 90 60 20],'FontSize',10,'ForegroundColor',[0 0 0],'HorizontalAlignment', 'left')

%add dashboard variables
MODE = uicontrol('Style','text','String','','Position',[200 420+250 310 40],'FontSize',25,'ForegroundColor',[0 0 0],'HorizontalAlignment', 'left');
EUL  = uicontrol('Style','text','String',{''},'Position',[55 270+250 105 100],'FontSize',20,'ForegroundColor',[0 0 0],'HorizontalAlignment', 'right');
ARAT = uicontrol('Style','text','String',{''},'Position',[80 110+250 100 100],'FontSize',20,'ForegroundColor',[0 0 0],'HorizontalAlignment', 'right'); 
MOME = uicontrol('Style','text','String','','Position',[160 50+250 250 35],'FontSize',20,'ForegroundColor',[0 0 0],'HorizontalAlignment', 'left');
TIME = uicontrol('Style','text','String','','Position',[510 50+250 450 35],'FontSize',20,'ForegroundColor',[0 0 0],'HorizontalAlignment', 'left');
ECLIPSE = uicontrol('Style','text','String','','Position',[555 425+250 130 35],'FontSize',20,'ForegroundColor',[0 0 0.8],'HorizontalAlignment', 'left','String','ECLIPSE','Visible',0);
RPM  = uicontrol('Style','text','String',{''},'Position',[720 270+250 220 100],'FontSize',20,'ForegroundColor',[0 0 0],'HorizontalAlignment', 'right');
M = uicontrol('Style','text','String',{''},'Position',[710 100+250 240 110],'FontSize',20,'ForegroundColor',[0 0 0],'HorizontalAlignment', 'right'); 
performance = uicontrol('Style','text','String',{''},'Position',[5 5 90 15],'FontSize',10,'ForegroundColor',[0 0 0],'HorizontalAlignment', 'left'); 

%make tags to retrieve them in the loop
set(MODE,'tag','satfigureMODE')
set(TIME,'tag','satfigureTIME')
set(EUL,'tag','satfigureEUL')
set(ARAT,'tag','satfigureARAT')
set(MOME,'tag','satfigureMOME')
set(ECLIPSE,'tag','satfigureECLIPSE')
set(RPM,'tag','satfigureRPM')
set(M,'tag','satfigureMo')
set(performance,'tag','performance')

 block.Dwork(1).data = 0;
 
 %% 2D Figures
orbitplan = block.DialogPrm(1).data;
jumps2 = diff(orbitplan(:,3));
jumps1 = abs(jumps2)>200;
jumps = find(jumps1); %jumps from longitude 180 to -180

%Big Map
 
 subplot(3,3,[7,9])

if (size(jumps,1)>1)
     geoplot(orbitplan(1:jumps(1),2),orbitplan(1:jumps(1),3),'DisplayName','Orbit 0');
     hold on 
     for i=1:1:size(jumps)-1
     geoplot(orbitplan(jumps(i)+1:jumps(i+1),2),orbitplan(jumps(i)+1:jumps(i+1),3),'DisplayName',['Orbit ' num2str(i)]);
     end
     geoplot(orbitplan(jumps(i+1)+1:end,2),orbitplan(jumps(i+1)+1:end,3),'DisplayName',['Orbit ' num2str(i+1)]);
elseif (size(jumps)==1)
    geoplot(orbitplan(1:jumps(1),2),orbitplan(1:jumps(1),3),'DisplayName','Orbit 0');
    hold on
    geoplot(orbitplan(jumps(1)+1:end,2),orbitplan(jumps(1)+1:end,3),'DisplayName','Orbit 1');
else
    geoplot(orbitplan(:,2),orbitplan(:,3),'DisplayName','Orbit 0');
end %plot data
%put legend
olegend = legend('Fontsize',12);
olegendposition=get(olegend,'Position');%to get the legend height
set(olegend,'Position',[0.92 0.4-0.5*olegendposition(4) 0 0],'Box','Off');

%put scodin text
 SCODIN = text(45,-73,[char(176) 'SCODIN']);
 set(SCODIN,'tag','MapSCODIN2')
 set(gca,'Position',[0.15 0 0.7 0.4],'clipping','off'); %moves the graph in the window
 geolimits([-60 60],[-180 180])
 
 %small map
 axes('position',[0.05 .16 .25 .25]);
 box on % put box around new pair of axes
if (size(jumps,1)>1)
     geoplot(orbitplan(1:jumps(1),2),orbitplan(1:jumps(1),3));
     hold on 
     for i=1:1:size(jumps)-1
     geoplot(orbitplan(jumps(i)+1:jumps(i+1),2),orbitplan(jumps(i)+1:jumps(i+1),3));
     end
     geoplot(orbitplan(jumps(i+1)+1:end,2),orbitplan(jumps(i+1)+1:end,3));
elseif (size(jumps)==1)
    geoplot(orbitplan(1:jumps(1),2),orbitplan(1:jumps(1),3));
    hold on
    geoplot(orbitplan(jumps(1)+1:end,2),orbitplan(jumps(1)+1:end,3));
else
    geoplot(orbitplan(:,2),orbitplan(:,3));
end %plot data

SCODIN = text(45,-73,compose("\xD83D\xDEF0"),'FontSize',30,'HorizontalAlignment', 'center','VerticalAlignment', 'middle');
set(SCODIN,'tag','MapSCODIN1','Visible','off') 
set(gca,'tag','minimap','Zoomlevel',4,'GridAlphaMode','manual','Visible','Off')
geolimits('manual')
allCurves = findobj(gca,'type','line');
set(allCurves,'Visible','off')

%%
%block.Dwork(1).data = TIME;
%% Function to get the CAD
function [fout, vout, cout] = rndread(filename)
% Reads CAD STL ASCII files, which most CAD programs can export.
% Used to create Matlab patches of CAD 3D data.
% Returns a vertex list and face list, for Matlab patch command.
%
fid=fopen(filename, 'r'); %Open the file, assumes STL ASCII format.
if fid == -1 
    error('File could not be opened, check name or path.')
end

vnum=0;       %Vertex number counter.
report_num=0; %Report the status as we go.
VColor = 0;
%
while feof(fid) == 0                    % test for end of file, if not then do stuff
    tline = fgetl(fid);                 % reads a line of data from file.
    fword = sscanf(tline, '%s ');       % make the line a character string
% Check for color
    if strncmpi(fword, 'c',1) == 1    % Checking if a "C"olor line, as "C" is 1st char.
       VColor = sscanf(tline, '%*s %f %f %f'); % & if a C, get the RGB color data of the face.
    end                                % Keep this color, until the next color is used.
    if strncmpi(fword, 'v',1) == 1   % Checking if a "V"ertex line, as "V" is 1st char.
       vnum = vnum + 1;                % If a V we count the # of V's
       report_num = report_num + 1;    % Report a counter, so long files show status
       if report_num > 249
           %fprintf('Reading vertix num: %d.\n',vnum);
           report_num = 0;
       end
       v(:,vnum) = sscanf(tline, '%*s %f %f %f'); % & if a V, get the XYZ data of it.
       c(:,vnum) = VColor;              % A color for each vertex, which will color the faces.
    end                                 % we "*s" skip the name "color" and get the data.                                          
end
%   Build face list; The vertices are in order, so just number them.
%
fnum = vnum/3;      %Number of faces, vnum is number of vertices.  STL is triangles.
flist = 1:vnum;     %Face list of vertices, all in order.
F = reshape(flist, 3,fnum); %Make a "3 by fnum" matrix of face list data.
%
%   Return the faces and vertexs.
%
fout = F';  %Orients the array for direct use in patch.
vout = v';  % "
cout = c';
%
fclose(fid);

%end Start

%%
%% Outputs:
%%   Functionality    : Called to generate block outputs in
%%                      simulation step
%%   Required         : Yes
%%   C-MEX counterpart: mdlOutputs
%%
function Outputs(block)

%get values from simulink
    vis_t = block.InputPort(1).Data; %time
    vis_b = block.InputPort(2).Data; %magnetic field
    vis_o = block.InputPort(3).Data; %sat attitude
    vis_w = block.InputPort(4).Data; %angular velocity
    vis_l = block.InputPort(5).Data; %angular momentum
    vis_s = block.InputPort(6).Data; %sun
    vis_G = block.InputPort(7).Data; %Ground
    vis_J = block.InputPort(8).Data; %julian date
    vis_e = block.InputPort(9).Data; %euler angles to target
    vis_M = block.InputPort(10).Data;%Mode
    vis_A = block.InputPort(11).Data;%Argentina
    vis_N = block.InputPort(12).Data;%Namibia
    vis_rpm = block.InputPort(13).Data;%Wheel rpm
    vis_Mo = block.InputPort(14).Data;%Magnetorquers moment
    latlon = block.InputPort(15).Data;%SCODIN Position on Map
%     
%Initially find the visual features
    persistent combinedobject Bline Blinetext TIME EUL ARAT MOME Sline Slinetext light_sun...
        Gline Glinetext MODE modes Aline Alinetext  Nline Nlinetext ECLIPSE RPM Mo SCODIN1...
        SCODIN2 minimap performance elapsed fps timestep

if  block.Dwork(1).data ~=1
    
    combinedobject=findobj('tag','satfigure');
    Bline=findobj('tag','satfigureBline');
    Blinetext=findobj('tag','satfigureBlinetext');
    Sline=findobj('tag','satfigureSline');
    Slinetext=findobj('tag','satfigureSlinetext');
    Gline=findobj('tag','satfigureGline');
    Glinetext=findobj('tag','satfigureGlinetext');
    Aline=findobj('tag','satfigureAline');
    Alinetext=findobj('tag','satfigureAlinetext');
    Nline=findobj('tag','satfigureNline');
    Nlinetext=findobj('tag','satfigureNlinetext');
    TIME=findobj('tag','satfigureTIME');
    EUL=findobj('tag','satfigureEUL');
    ARAT=findobj('tag','satfigureARAT');
    MOME=findobj('tag','satfigureMOME');
    MODE=findobj('tag','satfigureMODE');
    light_sun=findobj('tag','light_sun');
    ECLIPSE=findobj('tag','satfigureECLIPSE');
    RPM=findobj('tag','satfigureRPM');
    Mo=findobj('tag','satfigureMo');
    %map2d
    SCODIN1=findobj('tag','MapSCODIN1');
    SCODIN2=findobj('tag','MapSCODIN2');
    minimap = findobj('tag','minimap');
    performance = findobj('tag','performance');
    elapsed = 0;
    fps = 0;
    timestep = 0;
    
    modes = ["Sun Pointing","Ground","Imaging-Argentina","Imaging-Namibia","Nadir Pointing","Detumbling","OFF","Calibration-Compass"];
    tic
    block.Dwork(1).data = 1; %flag to search only once
end

try 
%Rotates satellite
Rot = makehgtform('zrotate',vis_o(1),'yrotate',vis_o(2),'xrotate',vis_o(3));
set(combinedobject,'Matrix',Rot)



%Updates Dashboard variables
TIME.String = [datestr(vis_J-1.7210585e6) ' (+' datestr(seconds(vis_t),'HH:MM:SS') ')'];
EUL.String = {[num2str(vis_e(3), '%.2f') char(176)] [num2str(vis_e(2), '%.2f') char(176)] [num2str(vis_e(1), '%.2f') char(176)]};
ARAT.String = {[num2str(vis_w(1), '%.2f') char(176) '/s'] [num2str(vis_w(2), '%.2f') char(176) '/s'] [num2str(vis_w(3), '%.2f') char(176) '/s']};
MOME.String = [num2str(norm(vis_l)*1000, '%1.2e') ' mNms'];
MODE.String = modes(vis_M);
RPM.String = {[num2str(vis_rpm(1), '%.0f') ' rpm (' num2str(abs(vis_rpm(1))/80, '%02.f') '%)'] [num2str(vis_rpm(2), '%.0f') ' rpm (' num2str(abs(vis_rpm(2))/80, '%02.f') '%)'] [num2str(vis_rpm(3), '%.0f') ' rpm (' num2str(abs(vis_rpm(3))/80, '%02.f') '%)']};
Mo.String = {[num2str(vis_Mo(1), '%.3f') ' Am' char(178) ' (' num2str(abs(vis_Mo(1))/0.24*100, '%03.f') '%)'] [num2str(vis_Mo(2), '%.3f') ' Am' char(178) ' (' num2str(abs(vis_Mo(2))/0.24*100, '%03.f') '%)'] [num2str(vis_Mo(3), '%.3f') ' Am' char(178) ' (' num2str(abs(vis_Mo(3))/0.24*100, '%03.f') '%)']};

% %change color for saturation warning
if norm(vis_l) > 1.77e-3
    MOME.ForegroundColor = [1 0 0];
    else
    MOME.ForegroundColor = [0 0 0];
end



%Rotates B-line
Bdirection = 200*[vis_b(1) vis_b(2) vis_b(3)]/norm([vis_b(1),vis_b(2),vis_b(3)]);
set(Bline,'ZData',[0 Bdirection(3)],'YData',[0 Bdirection(2)],'XData',[0 Bdirection(1)]);%changes the purple Bline
set(Blinetext,'position',1.1*Bdirection);

%Rotates S-line
if acos(vis_s(3)/norm(vis_s))<1.9%rad
%     set(Sline,'ZData',[0 0],'YData',[0 0],'XData',[0 0]);%changes the yellow Sline
%     light_sun.Visible = 'off';
set(ECLIPSE,'Visible',0)
else
set(ECLIPSE,'Visible',1)
end

Sdirection = 200*[vis_s(1) vis_s(2) vis_s(3)]/norm([vis_s(1),vis_s(2),vis_s(3)]);
set(Sline,'ZData',[0 Sdirection(3)],'YData',[0 Sdirection(2)],'XData',[0 Sdirection(1)]);%changes the yellow Sline
light_sun.Position = vis_s;%moves light
set(Slinetext,'position',1.1*Sdirection);

%Rotates G-line
Gdirection = 200*[vis_G(1) vis_G(2) vis_G(3)]/norm([vis_G(1),vis_G(2),vis_G(3)]);
set(Gline,'ZData',[0 Gdirection(3)],'YData',[0 Gdirection(2)],'XData',[0 Gdirection(1)]);%changes the Orange Tline
set(Glinetext,'position',1.1*Gdirection);

%Rotates A-line
Adirection = 200*[vis_A(1) vis_A(2) vis_A(3)]/norm([vis_A(1),vis_A(2),vis_A(3)]);
set(Aline,'ZData',[0 Adirection(3)],'YData',[0 Adirection(2)],'XData',[0 Adirection(1)]);%changes the light blue Aline
set(Alinetext,'position',1.1*Adirection);

%Rotates N-line
Ndirection = 200*[vis_N(1) vis_N(2) vis_N(3)]/norm([vis_N(1),vis_N(2),vis_N(3)]);
set(Nline,'ZData',[0 Ndirection(3)],'YData',[0 Ndirection(2)],'XData',[0 Ndirection(1)]);%changes the light blue Nline
set(Nlinetext,'position',1.1*Ndirection);

% if j > 5 %the map doesnt update as often
%Moves SCODIN on map
set(SCODIN2,'position',[latlon(1),latlon(2)]);

if block.InputPort(16).Data == 1
set(minimap,'mapcenter',[latlon(1),latlon(2)]);
set(SCODIN1,'position',[latlon(1),latlon(2)]);
end


%Update Performance data 
fps = ((1/toc-fps)*0.00002)*elapsed+fps;%filters the output
timestep = (((vis_t-elapsed)/toc-timestep)*0.00002)*elapsed+timestep; %filters the output
performance.String = [num2str(fps,'%2.0f') 'fps ' num2str(timestep, '%2.1f') 'x' ];
elapsed = vis_t;
tic
%update drawing
drawnow limitrate;

catch
warning('An error occurred while running the simulation and the simulation was terminated')
end

%end Outputs

%%
%% Update:
%%   Functionality    : Called to update discrete states
%%                      during simulation step
%%   Required         : No
%%   C-MEX counterpart: mdlUpdate
%%
function Update(block)

%end Update

%%
%% Derivatives:
%%   Functionality    : Called to update derivatives of
%%                      continuous states during simulation step
%%   Required         : No
%%   C-MEX counterpart: mdlDerivatives
%%
function Derivatives(block)

%end Derivatives

%%
%% Terminate:
%%   Functionality    : Called at the end of simulation for cleanup
%%   Required         : Yes
%%   C-MEX counterpart: mdlTerminate
%%
function Terminate(block)
clear ini
%end Terminate

function Buttonpause(src,evt)

button = findobj('tag','Button');

if get_param(bdroot,'SimulationStatus') == "paused"
set_param(bdroot, 'SimulationCommand', 'continue');    
button.String = 'Pause';

else 
set_param(bdroot, 'SimulationCommand', 'pause');
button.String = 'Continue';

end

function Buttonstop(src,evt)
set_param(bdroot, 'SimulationCommand', 'Stop');    

function Minimap(src,evt)
SCODIN1=findobj('tag','MapSCODIN1');
SCODIN2=findobj('tag','MapSCODIN2');
minimap = findobj('tag','minimap');
allCurves = findobj(minimap,'type','line');

if get_param('ADCS_Outline/VISUALISATION/Minimap','Value')=='1'
set_param('ADCS_Outline/VISUALISATION/Minimap','Value','0')
set(minimap,'Visible',0);
set(SCODIN1,'Visible',0);
set(allCurves,'visible',0);
else 
latlon = get(SCODIN2,'position');
set_param('ADCS_Outline/VISUALISATION/Minimap','Value','1')
set(minimap,'mapcenter',[latlon(1),latlon(2)]);
set(SCODIN1,'position',[latlon(1),latlon(2)]);
set(minimap,'Visible',1);
set(SCODIN1,'Visible',1);
set(allCurves,'visible',1);
end

function modeoverdrive(src,evt)
menu = findobj('tag','modeOverdrive');
set_param('ADCS_Outline/Guidance&Determination/Mode Selection/Constant3','Value',num2str(menu.Value-1))

function Sampletime(src,evt)
Pacerbutton = findobj('tag','Pacerbutton');
Lockbutton = findobj('tag','Lockbutton');
Timebutton = findobj('tag','Sampletimebutton');

Pacertime = str2double(Pacerbutton.String);
sampletime = str2double(get_param('ADCS_Outline/Pulse Generator','Period'));
ratio = Pacertime/sampletime;

if str2double(Timebutton.String)>=0
set_param('ADCS_Outline/Pulse Generator','Period',Timebutton.String)
else
Timebutton.String=0;
end

if Lockbutton.Value == 1
    sampletime = str2double(Timebutton.String);
    Pacerbutton.String=num2str(sampletime*ratio);
    set_param('ADCS_Outline/Real-Time Pacer','MaskValueString',num2str(sampletime*ratio))
end 

function Sampletimedown(src,evt)
Pacerbutton = findobj('tag','Pacerbutton');
Lockbutton = findobj('tag','Lockbutton');
Timebutton = findobj('tag','Sampletimebutton');

Pacertime = str2double(Pacerbutton.String);
sampletime = str2double(Timebutton.String);
ratio = Pacertime/sampletime;

if sampletime>1 && sampletime<2
    Timebutton.String='1';
elseif sampletime>1
    Timebutton.String=num2str(sampletime-1);
elseif sampletime>0.1
    Timebutton.String=num2str(sampletime-0.1); 
end
set_param('ADCS_Outline/Pulse Generator','Period',Timebutton.String)

if Lockbutton.Value == 1
    sampletime = str2double(Timebutton.String);
    Pacerbutton.String=num2str(sampletime*ratio);
    set_param('ADCS_Outline/Real-Time Pacer','MaskValueString',num2str(sampletime*ratio))
end 

function Sampletimeup(src,evt)
Pacerbutton = findobj('tag','Pacerbutton');
Lockbutton = findobj('tag','Lockbutton');
Timebutton = findobj('tag','Sampletimebutton');

Pacertime = str2double(Pacerbutton.String);
sampletime = str2double(Timebutton.String);
ratio = Pacertime/sampletime;

if sampletime>=1
Timebutton.String=num2str(sampletime+1);
elseif sampletime<1
Timebutton.String=num2str(sampletime+0.1); 
end
set_param('ADCS_Outline/Pulse Generator','Period',Timebutton.String)

if Lockbutton.Value == 1
    sampletime = str2double(Timebutton.String);
    Pacerbutton.String=num2str(sampletime*ratio);
    set_param('ADCS_Outline/Real-Time Pacer','MaskValueString',num2str(sampletime*ratio))
end 

function Pacer(src,evt)
Pacerbutton = findobj('tag','Pacerbutton');
Lockbutton = findobj('tag','Lockbutton');
Timebutton = findobj('tag','Sampletimebutton');

Pacertime = str2double(get_param('ADCS_Outline/Real-Time Pacer','MaskValueString'));
sampletime = str2double(Timebutton.String);
ratio = sampletime/Pacertime;
flag = 0;

if str2double(Pacerbutton.String)>=0
set_param('ADCS_Outline/Real-Time Pacer','MaskValueString',Pacerbutton.String)
else
set_param('ADCS_Outline/Real-Time Pacer','MaskValueString','inf')
Pacerbutton.String='inf';
flag = 1;
end

if Lockbutton.Value == 1 && flag ~= 1
    Pacertime = str2double(Pacerbutton.String);
    Timebutton.String=num2str(Pacertime*ratio);
    set_param('ADCS_Outline/Pulse Generator','Period',num2str(Pacertime*ratio))
end 

function Pacerdown(src,evt)
Pacerbutton = findobj('tag','Pacerbutton');
Lockbutton = findobj('tag','Lockbutton');
Timebutton = findobj('tag','Sampletimebutton');

Pacertime = str2double(Pacerbutton.String);
sampletime = str2double(Timebutton.String);
ratio = sampletime/Pacertime;
flag = 0;

if isinf(Pacertime)
Pacerbutton.String='512';
flag = 1;
elseif Pacertime > 0.25
Pacerbutton.String=num2str(Pacertime/2);  
else
Pacerbutton.String='0.25';
end
set_param('ADCS_Outline/Real-Time Pacer','MaskValueString',Pacerbutton.String)

if Lockbutton.Value == 1 && flag ~= 1
    Pacertime = str2double(Pacerbutton.String);
    Timebutton.String=num2str(Pacertime*ratio);
    set_param('ADCS_Outline/Pulse Generator','Period',num2str(Pacertime*ratio))
end 

function Pacerup(src,evt)
Pacerbutton = findobj('tag','Pacerbutton');
Lockbutton = findobj('tag','Lockbutton');
Timebutton = findobj('tag','Sampletimebutton');
flag = 0;

Pacertime = str2double(Pacerbutton.String);
sampletime = str2double(Timebutton.String);
ratio = sampletime/Pacertime;

if Pacertime<512
Pacerbutton.String=num2str(Pacertime*2);    
else
Pacerbutton.String='inf';
flag = 1;
end

set_param('ADCS_Outline/Real-Time Pacer','MaskValueString',Pacerbutton.String)

if Lockbutton.Value == 1 && flag ~= 1
    Pacertime = str2double(Pacerbutton.String);
    Timebutton.String=num2str(Pacertime*ratio);
    set_param('ADCS_Outline/Pulse Generator','Period',num2str(Pacertime*ratio))
end 

function Lock(src,evt)
Lockbutton = findobj('tag','Lockbutton');
if Lockbutton.Value == 1
    Lockbutton.String = compose("\xD83D\xDD12");
else
    Lockbutton.String = compose("\xD83D\xDD13");
end



