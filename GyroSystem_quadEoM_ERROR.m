clear all

syms spinspeed
syms Ic1x Ic1y Ic1z Ic2x Ic2y Ic2z Igax Igtr Im
syms thy thz thydot thzdot thyddot thzddot
syms thxT thyT thzT thxTdot thyTdot thzTdot
syms thxTddot thyTddot thzTddot
syms errory errorz errorydot errorzdot erroryddot errorzddot

EoM = [Ic2y*thyddot + Ic2y*thyTddot + Igtr*thyddot + Igtr*thyTddot + (Im*thyddot)/4 + (Im*thyTddot)/2 + Igax*spinspeed*thzdot + Igax*spinspeed*thzTdot - Ic2y*thxTddot*thz - Ic2y*thxTddot*thzT + Ic2x*thxTdot*thzdot - Ic2y*thxTdot*thzdot - Ic2z*thxTdot*thzdot + Ic2x*thxTdot*thzTdot - Ic2y*thxTdot*thzTdot - Ic2z*thxTdot*thzTdot + Igax*thxTdot*thzdot + Igax*thxTdot*thzTdot - Igtr*thxTddot*thz - Igtr*thxTddot*thzT - 2*Igtr*thxTdot*thzdot - 2*Igtr*thxTdot*thzTdot - (Im*thxTddot*thz)/4 - (Im*thxTddot*thzT)/2 - (Im*thxTdot*thzdot)/4 - (Im*thxTdot*thzTdot)/2 + Igax*spinspeed*thy*thxTdot + Igax*spinspeed*thyT*thxTdot == 0; ...
    Ic1z*thzddot + Ic2z*thzddot + Ic1z*thzTddot + Ic2z*thzTddot + Igtr*thzddot + Igtr*thzTddot + (Im*thzddot)/4 + (Im*thzTddot)/2 - Igax*spinspeed*thydot - Igax*spinspeed*thyTdot - Ic2x*thy*thxTddot + Ic2z*thy*thxTddot + Ic1z*thyT*thxTddot + Ic2z*thyT*thxTddot - Ic2x*thxTdot*thydot + Ic2y*thxTdot*thydot + Ic2z*thxTdot*thydot - Ic1x*thxTdot*thyTdot + Ic1y*thxTdot*thyTdot - Ic2x*thxTdot*thyTdot + Ic1z*thxTdot*thyTdot + Ic2y*thxTdot*thyTdot + Ic2z*thxTdot*thyTdot - Igax*thy*thxTddot - Igax*thxTdot*thydot - Igax*thxTdot*thyTdot + Igtr*thy*thxTddot + Igtr*thyT*thxTddot + 2*Igtr*thxTdot*thydot + 2*Igtr*thxTdot*thyTdot + (Im*thyT*thxTddot)/2 + (Im*thxTdot*thydot)/4 + (Im*thxTdot*thyTdot)/2 + Igax*spinspeed*thxTdot*thz + Igax*spinspeed*thxTdot*thzT == 0]

EoMsubbed = subs(EoM, [thy, thz, thydot, thzdot, thyddot, thzddot, ...
    thyT, thzT, thyTdot, thzTdot, thyTddot, thzTddot], [errory-thyT, ...
    errorz-thzT, errorydot-thyTdot, errorzdot-thzTdot, ...
    erroryddot-thyTddot, errorzddot-thzTddot, thyT, thzT, ...
    thyTdot, thzTdot, thyTddot, thzTddot])

soln = solve(EoMsubbed, [erroryddot, errorzddot]);

vars = [errory;errorz;errorydot;errorzdot]

erryddot1 = soln.erroryddot
errzddot1 = soln.errorzddot

derivatives = [errorydot;errorzdot;erryddot1;errzddot1]

pretty(derivatives)

H = jacobian(derivatives, vars)

pretty(H)