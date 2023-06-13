clear all

Ic1x = 0.0099039; Ic1y = 0.0076632; Ic1z = 0.001701;
Ic2x = 0.0010639; Ic2y = 0.00062812; Ic2z = 0.00052489;
Igax = 0.00073868; Igtr = 0.0005303;
Im = 0.00164;
C1damping = 0.02*180/pi;
C2damping = 0.01*180/pi;
spinspeed = 3000;

syms time

% FORCING DEFINITION
TankRoll = (time^3/(0.001+time^3))*(1e-16*pi/180)*sin(2*pi*10*time);
TankPitch = (time^3/(0.001+time^3))*(10*pi/180)*sin(2*pi*10*time);
TankYaw = (time^3/(0.001+time^3))*(10*pi/180)*cos(2*pi*15*time);
% END FORCING DEFINITION
TankRolldot = diff(TankRoll, time);
TankPitchdot = diff(TankPitch, time);
TankYawdot = diff(TankYaw, time);
TankRollddot = diff(TankRolldot, time);
TankPitchddot = diff(TankPitchdot, time);
TankYawddot = diff(TankYawdot, time);

matlabFunction([TankRoll; TankPitch; TankYaw; TankRolldot; TankPitchdot; ...
    TankYawdot; TankRollddot; TankPitchddot; TankYawddot], 'File', 'forcing',...
    'Vars', time);

systemparameters = [Ic1x, Ic1y, Ic1z, Ic2x, Ic2y, Ic2z, Igax, ...
    Igtr, Im, C1damping, C2damping, spinspeed];

simulationtime = 0.5;

initialconditions = [0;0;0;0] %[thy; thz; thydot; thzdot]

tspan = linspace(0, simulationtime, 10000);

tolerance = 1e-10;
options = odeset('RelTol', tolerance, 'AbsTol', tolerance);

odehandle = @(t, z) myrhs(t, z, systemparameters);

[t, zvals] = ode45(odehandle, tspan, initialconditions, options);

forcingvector = forcing(t');

xExcitation = forcingvector(1,:);
yExcitation = forcingvector(2,:);
zExcitation = forcingvector(3,:);

error = [zvals(:,1) + yExcitation', zvals(:,2) + zExcitation'];

clf

tankrollplot = plot(t, (180/pi)*xExcitation, ':k');
hold on;
tankpitchplot = plot(t, (180/pi)*yExcitation, '--k');
tankyawplot = plot(t, (180/pi)*zExcitation, '-.k');
lospitchplot = plot(t, (180/pi)*error(:,1), '--k');
losyawplot = plot(t, (180/pi)*error(:,2), '-.k');
lospitchplot.LineWidth = 1.5;
losyawplot.LineWidth = 1.5;
lgnd = legend('Tank Roll', 'Tank Pitch', 'Tank Yaw', ...
    'L.o.S. Pitch', 'L.o.S. Yaw', ...
    'Location', 'southoutside');
lgnd.NumColumns = 3;
xlabel('Time (sec)');
ylabel('Angle (deg)');

set(gcf,'units','centimeters','position',[10,10,15,10]);
shg

function zdot = myrhs(t, z, params)

openedparams = num2cell(params);
[Ic1x, Ic1y, Ic1z, Ic2x, Ic2y, Ic2z, Igax, Igtr, Im, ...
    C1damping, C2damping, spinspeed] = deal(openedparams{:});

currentforcing = forcing(t);

thxT = currentforcing(1);
thyT = currentforcing(2);
thzT = currentforcing(3);
thxTdot = currentforcing(4);
thyTdot = currentforcing(5);
thzTdot = currentforcing(6);
thxTddot = currentforcing(7);
thyTddot = currentforcing(8);
thzTddot = currentforcing(9);

thy = z(1);
thz = z(2);

thydot = z(3);
thzdot = z(4);

MatrixA = gyromatrixA(thxT, thyT, thzT, thxTdot, thyTdot, ...
    thzTdot, thxTddot, thyTddot, thzTddot, thy, thz, thydot, thzdot, ...
    Ic1x, Ic1y, Ic1z, Ic2x, Ic2y, Ic2z, Igax, Igtr, Im, ...
    C1damping, C2damping, spinspeed);

MatrixB = gyromatrixB(thxT, thyT, thzT, thxTdot, thyTdot, ...
    thzTdot, thxTddot, thyTddot, thzTddot, thy, thz, thydot, thzdot, ...
    Ic1x, Ic1y, Ic1z, Ic2x, Ic2y, Ic2z, Igax, Igtr, Im, ...
    C1damping, C2damping, spinspeed);

thddots = MatrixA\MatrixB;

thyddot = thddots(1);
thzddot = thddots(2);

zdot = [thydot; thzdot; thyddot; thzddot];

end