clear all

syms spinspeed
syms Ic1x Ic1y Ic1z Ic2x Ic2y Ic2z Igax Igtr Im
syms thy thz thydot thzdot thyddot thzddot
syms thxT thyT thzT thxTdot thyTdot thzTdot
syms thxTddot thyTddot thzTddot
syms errory errorz errorydot errorzdot erroryddot errorzddot

EoM = [Ic2y*thyddot + Ic2y*thyTddot + Igtr*thyddot + Igtr*thyTddot + (Im*thyddot)/4 + (Im*thyTddot)/2 + Igax*spinspeed*thzdot + Igax*spinspeed*thzTdot == 0; ...
    Ic1z*thzddot + Ic2z*thzddot + Ic1z*thzTddot + Ic2z*thzTddot + Igtr*thzddot + Igtr*thzTddot + (Im*thzddot)/4 + (Im*thzTddot)/2 - Igax*spinspeed*thydot - Igax*spinspeed*thyTdot == 0]

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