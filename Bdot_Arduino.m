function M = Bdot_Arduino(B)

persistent icount a %integrator

if (isempty(icount))
   try
      a = arduino;
      disp(a)
      %integrator = [0;0;0];
      B=[1;0;0];
   catch
      warning('Unable to connect, user input required')
      disp('For Windows:')
      disp('  Open device manager, select "Ports (COM & LPT)"')
      disp('  Look for COM port of Arduino such as COM4')
      disp('For MacOS:')
      disp('  Open terminal and type: ls /dev/*.')
      disp('  Search for /dev/tty.usbmodem* or /dev/tty.usbserial*. The port number is *.')
      disp('For Linux')
      disp('  Open terminal and type: ls /dev/tty*')
      disp('  Search for /dev/ttyUSB* or /dev/ttyACM*. The port number is *.')
      disp('')
      com_port = input('Specify port (e.g. COM4 for Windows or /dev/ttyUSB0 for Linux): ','s');
      a = arduino(com_port,'nano');
      disp(a)
   end
   icount = 0;
end  
%% code starts here

% B = B/norm(B);
% 
% Bdot =0.2*(B-integrator);
% 
% scale = max([abs(Bdot(1)/0.24);abs(Bdot(2)/0.24);abs(Bdot(3)/0.13)]);
% 
% M = -Bdot/scale;
% 
% integrator = integrator + Bdot*.01;


%% for 2 independent arduinos

%functions to write and read from second arduino
bx = @(level) writePWMDutyCycle(a,'D11',max(0,min(1,(level+80000)*0.00003125/5)));
by = @(level) writePWMDutyCycle(a,'D10',max(0,min(1,(level+80000)*0.00003125/5)));
bz = @(level) writePWMDutyCycle(a,'D9' ,max(0,min(1,(level+80000)*0.00003125/5)));

mx = @() readVoltage(a, 'A3');
my = @() readVoltage(a, 'A2');
mz = @() readVoltage(a, 'A1');

%send simulink data to 2nd arduino for B
bx(B(1))
by(B(2))
bz(B(3))

%send arduino data to simulink for M

M(1) = (mx()-2.5)*0.1;
M(2) = (my()-2.5)*0.1;
M(3) = (mz()-2.5)*0.1;

end







