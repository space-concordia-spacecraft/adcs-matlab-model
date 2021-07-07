clear
close all
gx = geoaxes;
hold on

latA=15*pi/180;
lonA=-150*pi/180;
latC=14*pi/180;
lonC=179*pi/180;
latB=16.5*pi/180;
lonB=17.4*pi/180;

greatcircleN(latA,lonA,latB,lonB,gx)
greatcircleN(latB,lonB,latC,lonC,gx)


function greatcircleN(latA,lonA,latB,lonB,gx)
N = lonB-lonA;
a = pi/2-latB;
b = pi/2-latA;
n = acos(cos(a)*cos(b)+sin(a)*sin(b)*cos(N));
A = asin(sin(N)/sin(n)*sin(a));
B = asin(sin(N)/sin(n)*sin(b));
lat_old=latA;
lon_old=lonA;
if abs(N)>pi/2
    B1 = acos(sin(A)*cos(b));
    n1 = asin(sin(b)/sin(B1));
    step=n1/50;
    for ni=step:step:n1
    ai =acos(cos(b)*cos(ni)+sin(b)*sin(ni)*cos(A)); 
    lat = pi/2-ai;
    Ni = real(asin(sin(A)/sin(ai)*sin(ni)));
    lon = lonA+Ni;
    geoplot(gx,[lat_old*180/pi lat*180/pi],[lon_old*180/pi lon*180/pi],'color','blue');
    lat_old = lat;
    lon_old = lon;
    end
    A2 = pi-B1;
    N2 = N-pi/2;
    n2 = asin(sin(a)/sin(A2)*sin(N2));
    b2 = asin(sin(a)/sin(A2)*sin(B));
    step = n2/50;
    for ni=step:step:n2
    ai =acos(cos(b2)*cos(ni)+sin(b2)*sin(ni)*cos(A2)); 
    lat = pi/2-ai;
    Ni = asin(sin(A2)/sin(ai)*sin(ni));
    lon = lonA+pi/2+Ni;
    geoplot(gx,[lat_old*180/pi lat*180/pi],[lon_old*180/pi lon*180/pi],'color','blue');
    lat_old = lat;
    lon_old = lon;
    end
else
    step=n/100;
    for ni=step:step:n
    ai =acos(cos(b)*cos(ni)+sin(b)*sin(ni)*cos(A)); 
    lat = pi/2-ai;
    Ni = asin(sin(A)/sin(ai)*sin(ni));
    lon = lonA+Ni;
    geoplot(gx,[lat_old*180/pi lat*180/pi],[lon_old*180/pi lon*180/pi],'color','blue');
    lat_old = lat;
    lon_old = lon;
    end
end   
    
geoplot(gx,[latA*180/pi latB*180/pi],[lonA*180/pi lonB*180/pi]);
end



%geolimits(gx,[-90 90],[-180 180]);
% Plot your Patch data here
% ax2 = axes;
% x = [0.25 0.6 0.6 0.25]; % Modify x coordinates of the polygon 
% y = [0.25 0.25 0.4 0.4]; % Modify y coordinates of the polygon
% 
% t = linspace(0, 2*pi);
% r = 0.5;
% x = r*cos(t)+0.5;
% y = r*sin(t)+0.5;
% 
% patch(ax2, x, y,'black','FaceAlpha',.4,'EdgeColor','none'); % Modify patch color and transparency 
% axis([0 1 0 1]);
% % Set ax2 visibility to 'off'
% ax2.Visible = 'off'; 
% ax2.XTick = []; 
% ax2.YTick = []; 
