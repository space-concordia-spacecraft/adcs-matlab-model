function MagnetometerCalibration(block)
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
block.NumInputPorts  = 6;
block.NumOutputPorts = 0;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
block.InputPort(1).Dimensions        = 3;
block.InputPort(1).DatatypeID  = 0;  % double
block.InputPort(1).Complexity  = 'Real';
block.InputPort(1).DirectFeedthrough = true;

block.InputPort(2).Dimensions  = [3,3];
block.InputPort(2).DatatypeID  = 0;  % double
block.InputPort(2).Complexity  = 'Real';
block.InputPort(2).DirectFeedthrough = true;

block.InputPort(3).Dimensions        = 3;
block.InputPort(3).DatatypeID  = 0;  % double
block.InputPort(3).Complexity  = 'Real';
block.InputPort(3).DirectFeedthrough = true;

block.InputPort(4).Dimensions        = 3;
block.InputPort(4).DatatypeID  = 0;  % double
block.InputPort(4).Complexity  = 'Real';
block.InputPort(4).DirectFeedthrough = true;

block.InputPort(5).Dimensions        = 1;
block.InputPort(5).DatatypeID  = 0;  % double
block.InputPort(5).Complexity  = 'Real';
block.InputPort(5).DirectFeedthrough = true;

block.InputPort(6).Dimensions        = 3;
block.InputPort(6).DatatypeID  = 0;  % double
block.InputPort(6).Complexity  = 'Real';
block.InputPort(6).DirectFeedthrough = true;
% Override output port properties
% block.OutputPort(1).Dimensions       = 1;
% block.OutputPort(1).DatatypeID  = 0; % double
% block.OutputPort(1).Complexity  = 'Real';

% Register parameters
block.NumDialogPrms     = 0;

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
 block.Dwork(1).data = 0;
%end Start

%%
%% Outputs:
%%   Functionality    : Called to generate block outputs in
%%                      simulation step
%%   Required         : Yes
%%   C-MEX counterpart: mdlOutputs
%%
function Outputs(block)


mpoints = block.InputPort(1).Data;
invA_mag = block.InputPort(2).Data;
b_mag = block.InputPort(3).Data;
radii_ref=block.InputPort(4).Data;
Mode=block.InputPort(5).Data;
B_E=block.InputPort(6).Data;
persistent j magpoints numofpoints trigger magearth

if  block.Dwork(1).data == 0
    numofpoints=600;
    j=numofpoints+1;
    magpoints = zeros(numofpoints,3);
    magearth = zeros(numofpoints,3);
    block.Dwork(1).data = 1;
end

if Mode-trigger > 0
    if Mode==8 %trigger calibration
    j=1;
    end
end
trigger=Mode;

if j<numofpoints
    magpoints(j,:)=mpoints;
    magearth(j,:)=B_E;
elseif j==numofpoints
    magpoints(j,:)=mpoints;
    magearth(j,:)=B_E;
%[ center, radii, evecs, v, chi2 ] = ellipsoid_fit( magpoints, '' );
figure('Name','Magnetometer Calibration','NumberTitle','off')
plot3(magpoints(:,1),magpoints(:,2),magpoints(:,3),'o','Color','k','MarkerSize',5,'MarkerFaceColor','r');
hold on
grid on
xlabel('x') 
ylabel('y')
zlabel('z') 
camva('manual')
camva(10)
daspect([1 1 1]) 
axis equal
set(gca,'Position',get(gca,'Position').*[2 1.5 1 1]);

[centerE, radiic, evecs, v, chi2 ] = ellipsoid_fit(magpoints,'');
% [theta,phi] = ndgrid(linspace(0,pi,20),linspace(0,2*pi,40));
% x2 = radiic(1)*sin(theta).*cos(phi);
% y2 = radiic(2)*sin(theta).*sin(phi);
% z2 = radiic(3)*cos(theta);
% P = bsxfun(@plus,evecs*[x2(:),y2(:),z2(:)]',centerE);
% X = reshape(P(1,:),size(x2,1),size(x2,2));
% Y = reshape(P(2,:),size(y2,1),size(y2,2));
% Z = reshape(P(3,:),size(z2,1),size(z2,2));
% surf(X,Y,Z,'FaceColor',[0.8500 0.3250 0.0980],'FaceAlpha',0.2);
% 
% %get the calibrated values
scale=[radiic(1),0,0;
       0,radiic(2),0;
       0,0,radiic(3)]/mean(radiic);
Amag = evecs/scale*evecs';
bmag=centerE;
%pointc=zeros(numofpoints,3);
point_scaled=magpoints;
pointc=(Amag*(point_scaled-bmag.').').';%calibrated B

for i = 1:1:5

    for count=1:1:numofpoints
        pointc(count,:)=pointc(count,:)*norm(magearth(1,:))/norm(magearth(count,:));
        point_scaled(count,:)=((Amag\pointc(count,:).')+bmag).';%calibrated B
    end

    [centerE, radiic, evecs, v, chi2 ] = ellipsoid_fit(point_scaled,'');
    %get the calibrated values
    scale=[radiic(1),0,0;
           0,radiic(2),0;
           0,0,radiic(3)]/mean(radiic);
    Amag = evecs/scale*evecs';
    bmag=centerE;

    pointc=(Amag*(magpoints-bmag.').').';%calibrated B
end

[theta,phi] = ndgrid(linspace(0,pi,20),linspace(0,2*pi,40));
x2 = radiic(1)*sin(theta).*cos(phi);
y2 = radiic(2)*sin(theta).*sin(phi);
z2 = radiic(3)*cos(theta);
P = bsxfun(@plus,evecs*[x2(:),y2(:),z2(:)]',centerE);
X = reshape(P(1,:),size(x2,1),size(x2,2));
Y = reshape(P(2,:),size(y2,1),size(y2,2));
Z = reshape(P(3,:),size(z2,1),size(z2,2));
surf(X,Y,Z,'FaceColor',[0.8500 0.3250 0.0980],'FaceAlpha',0.2);

plot3(pointc(:,1),pointc(:,2),pointc(:,3),'o','Color','k','MarkerSize',5,'MarkerFaceColor','b');
%plot3(cpoints(:,1),cpoints(:,2),cpoints(:,3),'o','Color','k','MarkerSize',10,'MarkerFaceColor','b');

[centerE, radii, evecs, v, chi2 ] = ellipsoid_fit(pointc,'');
x2 = radii(1)*sin(theta).*cos(phi);
y2 = radii(2)*sin(theta).*sin(phi);
z2 = radii(3)*cos(theta);
P = bsxfun(@plus,evecs*[x2(:),y2(:),z2(:)]',centerE);
X = reshape(P(1,:),size(x2,1),size(x2,2));
Y = reshape(P(2,:),size(y2,1),size(y2,2));
Z = reshape(P(3,:),size(z2,1),size(z2,2));
surf(X,Y,Z,'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.2);

legend([{'Raw Measurements'},{'Raw Ellipsoid fit'},{'Calibrated Values'},{'Calibrated Ellipsoid fit'}],'FontSize',10,'Position',[0.74 0 0.2 0.2])

save('MagnetometerBias.mat','Amag','bmag')
A_mag = inv(invA_mag);
uicontrol('Style','text','String',[' Calculated Bias = ' newline...
                                    ' x = ' num2str(bmag(1),'%0.3u') 'T' newline...
                                    ' y = ' num2str(bmag(2),'%0.3u') 'T' newline...
                                    ' z = ' num2str(bmag(3),'%0.3u') 'T' newline...
                                    newline...
                                    ' Real Bias = ' newline...
                                    ' x = ' num2str(b_mag(1),'%0.3u') 'T' newline...
                                    ' y = ' num2str(b_mag(2),'%0.3u') 'T' newline...
                                    ' z = ' num2str(b_mag(3),'%0.3u') 'T' newline...
                                    newline...
                                    ' Error = ' num2str(abs(norm(b_mag-bmag)),'%0.2u') 'T' newline...
                                    newline...
                                    ' Calculated A = ' newline...
                                    ' ' num2str(Amag(1,1),'%0.3f') ', ' num2str(Amag(1,2),'%0.3f') ', ' num2str(Amag(1,3),'%0.3f') newline...
                                    ' ' num2str(Amag(2,1),'%0.3f') ', ' num2str(Amag(2,2),'%0.3f') ', ' num2str(Amag(2,3),'%0.3f') newline...
                                    ' ' num2str(Amag(3,1),'%0.3f') ', ' num2str(Amag(3,2),'%0.3f') ', ' num2str(Amag(3,3),'%0.3f') newline...
                                    newline...   
                                    ' Real A = ' newline...
                                    ' ' num2str(A_mag(1,1),'%0.3f') ', ' num2str(A_mag(1,2),'%0.3f') ', ' num2str(A_mag(1,3),'%0.3f') newline...
                                    ' ' num2str(A_mag(2,1),'%0.3f') ', ' num2str(A_mag(2,2),'%0.3f') ', ' num2str(A_mag(2,3),'%0.3f') newline...
                                    ' ' num2str(A_mag(3,1),'%0.3f') ', ' num2str(A_mag(3,2),'%0.3f') ', ' num2str(A_mag(3,3),'%0.3f') newline...
                                    newline... 
                                    ' Error = ' num2str(abs(( sum(abs(sort(radiic,'descend')/mean(radiic)-sort(radii_ref,'descend')/mean(radii_ref))./sort(radii_ref,'descend')))*100),'%0.2f') '%'...
                                    ],'Position',[0 0 160 420],'FontSize',11,'ForegroundColor',[0 0 0],'HorizontalAlignment', 'left')

set_param('ADCS_Outline/Guidance&Determination/Mode Selection/Constant3','Value','0') %reset modeOverdrive
menu = findobj('tag','modeOverdrive');
menu.Value=1;
assignin('base','Amag',Amag)
assignin('base','bmag',bmag)
end
j=j+1;

function [ center, radii, evecs, v, chi2 ] = ellipsoid_fit( X, equals )
%
% Fit an ellispoid/sphere/paraboloid/hyperboloid to a set of xyz data points:
%
%   [center, radii, evecs, pars ] = ellipsoid_fit( X )
%   [center, radii, evecs, pars ] = ellipsoid_fit( [x y z] );
%   [center, radii, evecs, pars ] = ellipsoid_fit( X, 1 );
%   [center, radii, evecs, pars ] = ellipsoid_fit( X, 2, 'xz' );
%   [center, radii, evecs, pars ] = ellipsoid_fit( X, 3 );
%
% Parameters:
% * X, [x y z]   - Cartesian data, n x 3 matrix or three n x 1 vectors
% * flag         - '' or empty fits an arbitrary ellipsoid (default),
%                - 'xy' fits a spheroid with x- and y- radii equal
%                - 'xz' fits a spheroid with x- and z- radii equal
%                - 'xyz' fits a sphere
%                - '0' fits an ellipsoid with its axes aligned along [x y z] axes
%                - '0xy' the same with x- and y- radii equal
%                - '0xz' the same with x- and z- radii equal
%
% Output:
% * center    -  ellispoid or other conic center coordinates [xc; yc; zc]
% * radii     -  ellipsoid or other conic radii [a; b; c]
% * evecs     -  the radii directions as columns of the 3x3 matrix
% * v         -  the 10 parameters describing the ellipsoid / conic algebraically: 
%                Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz + J = 0
% * chi2      -  residual sum of squared errors (chi^2), this chi2 is in the 
%                coordinate frame in which the ellipsoid is a unit sphere.
%
% Author:
% Yury Petrov, Oculus VR
% Date:
% September, 2015
%

narginchk( 1, 3 ) ;  % check input arguments
if nargin == 1
    equals = ''; % no constraints by default
end
    
if size( X, 2 ) ~= 3
    error( 'Input data must have three columns!' );
else
    x = X( :, 1 );
    y = X( :, 2 );
    z = X( :, 3 );
end

% need nine or more data points
if length( x ) < 9 && strcmp( equals, '' ) 
   error( 'Must have at least 9 points to fit a unique ellipsoid' );
end
if length( x ) < 8 && ( strcmp( equals, 'xy' ) || strcmp( equals, 'xz' ) )
   error( 'Must have at least 8 points to fit a unique ellipsoid with two equal radii' );
end
if length( x ) < 6 && strcmp( equals, '0' )
   error( 'Must have at least 6 points to fit a unique oriented ellipsoid' );
end
if length( x ) < 5 && ( strcmp( equals, '0xy' ) || strcmp( equals, '0xz' ) )
   error( 'Must have at least 5 points to fit a unique oriented ellipsoid with two equal radii' );
end
if length( x ) < 4 && strcmp( equals, 'xyz' );
   error( 'Must have at least 4 points to fit a unique sphere' );
end

% fit ellipsoid in the form Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx +
% 2Hy + 2Iz + J = 0 and A + B + C = 3 constraint removing one extra
% parameter
if strcmp( equals, '' )
    D = [ x .* x + y .* y - 2 * z .* z, ...
        x .* x + z .* z - 2 * y .* y, ...
        2 * x .* y, ...
        2 * x .* z, ...
        2 * y .* z, ...
        2 * x, ...
        2 * y, ...
        2 * z, ...
        1 + 0 * x ];  % ndatapoints x 9 ellipsoid parameters
elseif strcmp( equals, 'xy' )
    D = [ x .* x + y .* y - 2 * z .* z, ...
        2 * x .* y, ...
        2 * x .* z, ...
        2 * y .* z, ...
        2 * x, ...
        2 * y, ...
        2 * z, ...
        1 + 0 * x ];  % ndatapoints x 8 ellipsoid parameters
elseif strcmp( equals, 'xz' )
    D = [ x .* x + z .* z - 2 * y .* y, ...
        2 * x .* y, ...
        2 * x .* z, ...
        2 * y .* z, ...
        2 * x, ...
        2 * y, ...
        2 * z, ...
        1 + 0 * x ];  % ndatapoints x 8 ellipsoid parameters
    % fit ellipsoid in the form Ax^2 + By^2 + Cz^2 + 2Gx + 2Hy + 2Iz = 1
elseif strcmp( equals, '0' )
    D = [ x .* x + y .* y - 2 * z .* z, ...
          x .* x + z .* z - 2 * y .* y, ...
          2 * x, ...
          2 * y, ... 
          2 * z, ... 
          1 + 0 * x ];  % ndatapoints x 6 ellipsoid parameters
    % fit ellipsoid in the form Ax^2 + By^2 + Cz^2 + 2Gx + 2Hy + 2Iz = 1,
    % where A = B or B = C or A = C
elseif strcmp( equals, '0xy' )
    D = [ x .* x + y .* y - 2 * z .* z, ...
          2 * x, ...
          2 * y, ... 
          2 * z, ... 
          1 + 0 * x ];  % ndatapoints x 5 ellipsoid parameters
elseif strcmp( equals, '0xz' )
    D = [ x .* x + z .* z - 2 * y .* y, ...
          2 * x, ...
          2 * y, ... 
          2 * z, ... 
          1 + 0 * x ];  % ndatapoints x 5 ellipsoid parameters
     % fit sphere in the form A(x^2 + y^2 + z^2) + 2Gx + 2Hy + 2Iz = 1
elseif strcmp( equals, 'xyz' )
    D = [ 2 * x, ...
          2 * y, ... 
          2 * z, ... 
          1 + 0 * x ];  % ndatapoints x 4 ellipsoid parameters
else
    error( [ 'Unknown parameter value ' equals '!' ] );
end

% solve the normal system of equations
d2 = x .* x + y .* y + z .* z; % the RHS of the llsq problem (y's)
u = ( D' * D ) \ ( D' * d2 );  % solution to the normal equations

% find the residual sum of errors
% chi2 = sum( ( 1 - ( D * u ) ./ d2 ).^2 ); % this chi2 is in the coordinate frame in which the ellipsoid is a unit sphere.

% find the ellipsoid parameters
% convert back to the conventional algebraic form
v=zeros(1,10);
if strcmp( equals, '' )
    v(1) = u(1) +     u(2) - 1;
    v(2) = u(1) - 2 * u(2) - 1;
    v(3) = u(2) - 2 * u(1) - 1;
    v( 4 : 10 ) = u( 3 : 9 );
elseif strcmp( equals, 'xy' )
    v(1) = u(1) - 1;
    v(2) = u(1) - 1;
    v(3) = -2 * u(1) - 1;
    v( 4 : 10 ) = u( 2 : 8 );
elseif strcmp( equals, 'xz' )
    v(1) = u(1) - 1;
    v(2) = -2 * u(1) - 1;
    v(3) = u(1) - 1;
    v( 4 : 10 ) = u( 2 : 8 );
elseif strcmp( equals, '0' )
    v(1) = u(1) +     u(2) - 1;
    v(2) = u(1) - 2 * u(2) - 1;
    v(3) = u(2) - 2 * u(1) - 1;
    v = [ v(1) v(2) v(3) 0 0 0 u( 3 : 6 )' ];
elseif strcmp( equals, '0xy' )
    v(1) = u(1) - 1;
    v(2) = u(1) - 1;
    v(3) = -2 * u(1) - 1;
    v = [ v(1) v(2) v(3) 0 0 0 u( 2 : 5 )' ];
elseif strcmp( equals, '0xz' )
    v(1) = u(1) - 1;
    v(2) = -2 * u(1) - 1;
    v(3) = u(1) - 1;
    v = [ v(1) v(2) v(3) 0 0 0 u( 2 : 5 )' ];
elseif strcmp( equals, 'xyz' )
    v = [ -1 -1 -1 0 0 0 u( 1 : 4 )' ];
end
v = v';

% form the algebraic form of the ellipsoid
A = [ v(1) v(4) v(5) v(7); ...
      v(4) v(2) v(6) v(8); ...
      v(5) v(6) v(3) v(9); ...
      v(7) v(8) v(9) v(10) ];
% find the center of the ellipsoid
center = -A( 1:3, 1:3 ) \ v( 7:9 );
% form the corresponding translation matrix
T = eye( 4 );
T( 4, 1:3 ) = center';
% translate to the center
R = T * A * T';
% solve the eigenproblem
[ evecs, evals ] = eig( R( 1:3, 1:3 ) / -R( 4, 4 ) );
radii = sqrt( 1 ./ diag( abs( evals ) ) );
sgns = sign( diag( evals ) );
radii = radii .* sgns;

% calculate difference of the fitted points from the actual data normalized by the conic radii
d = [ x - center(1), y - center(2), z - center(3) ]; % shift data to origin
d = d * evecs; % rotate to cardinal axes of the conic;
d = [ d(:,1) / radii(1), d(:,2) / radii(2), d(:,3) / radii(3) ]; % normalize to the conic radii
chi2 = sum( abs( 1 - sum( d.^2 .* repmat( sgns', size( d, 1 ), 1 ), 2 ) ) );

if abs( v(end) ) > 1e-6
    v = -v / v(end); % normalize to the more conventional form with constant term = -1
else
    v = -sign( v(end) ) * v;
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