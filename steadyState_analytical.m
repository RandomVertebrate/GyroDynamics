clear all

syms t spinspeed
syms Ic1x Ic1y Ic1z Ic2x Ic2y Ic2z Igax Igtr Im
syms thy thz thydot thzdot thyddot thzddot
syms thxT thyT thzT thxTdot thyTdot thzTdot
syms thxTddot thyTddot thzTddot
syms A1 A2 A3 A4 B1 B2 B3 B4 exf sinwave coswave

spinspeedvalue = 3000;

EoM1 = Ic2y*thyddot + Ic2y*thyTddot + Igtr*thyddot + Igtr*thyTddot + (Im*thyddot)/4 + (Im*thyTddot)/2 + Igax*spinspeed*thzdot + Igax*spinspeed*thzTdot == 0;
EoM2 = Ic1z*thzddot + Ic2z*thzddot + Ic1z*thzTddot + Ic2z*thzTddot + Igtr*thzddot + Igtr*thzTddot + (Im*thzddot)/4 + (Im*thzTddot)/2 - Igax*spinspeed*thydot - Igax*spinspeed*thyTdot == 0;

thetayT = A1*sinwave
thetayTdot = A1*exf*coswave
thetayTddot = -A1*exf^2*sinwave

thetazT = A1*sinwave
thetazTdot = A1*exf*coswave
thetazTddot = -A1*exf^2*sinwave

errory = A2*sinwave+B2*coswave
errorydot = A2*exf*coswave-B2*exf*sinwave
erroryddot = -A2*exf^2*sinwave-B2*exf^2*coswave

errorz = A3*sinwave+B3*coswave
errorzdot = A3*exf*coswave-B3*exf*sinwave
errorzddot = -A3*exf^2*sinwave-B3*exf^2*coswave

EoMsubbed = subs([EoM1;EoM2], [thy, thz, thydot, thzdot, thyddot, thzddot, ...
    thyT, thzT, thyTdot, thzTdot, thyTddot, thzTddot], [errory-thetayT, ...
    errorz-thetazT, errorydot-thetayTdot, errorzdot-thetazTdot, ...
    erroryddot-thetayTddot, errorzddot-thetazTddot, thetayT, thetazT, ...
    thetayTdot, thetazTdot, thetayTddot, thetazTddot])

AmplitudeEquations = [diff(EoMsubbed, coswave); diff(EoMsubbed, sinwave)]

soln = solve(AmplitudeEquations, [A2, B2, A3, B3]);

diary("SteadyStateResponse.txt")

netAmplitudey = sqrt(simplify(soln.A2^2+soln.B2^2))

netAmplitudez = sqrt(simplify(soln.A3^2+soln.B3^2))

netPhasey = atan(soln.B2/soln.A2)

netPhasez = atan(soln.B3/soln.A3)

diary off

anony = @(Ic1x, Ic1y, Ic1z, Ic2x, Ic2y, Ic2z, Igax, Igtr, Im, A1, exf) simplify((1-netAmplitudey)*100)

anonz = @(Ic1x, Ic1y, Ic1z, Ic2x, Ic2y, Ic2z, Igax, Igtr, Im, A1, exf) simplify((1-netAmplitudez)*100)

totalAmplitude = sqrt(netAmplitudey^2+netAmplitudez^2)