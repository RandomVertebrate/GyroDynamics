clear all

ModelOrder = 0

syms t thetaxtank(t) thetaytank(t) thetaztank(t)
syms thetaygyro(t) thetazgyro(t) thetaxgyro(t) spinspeed
syms Ic1x Ic1y Ic1z Ic2x Ic2y Ic2z Igax Igtr Im
syms epsilon thy thz thydot thzdot thyddot thzddot
syms thxT thyT thzT thxTdot thyTdot thzTdot
syms thxTddot thyTddot thzTddot
syms C1damping C2damping

R1tank = [1 0 0; 0 cos(thetaxtank) -sin(thetaxtank); 0 sin(thetaxtank) cos(thetaxtank)];
R2tank = [cos(thetaytank) 0 sin(thetaytank); 0 1 0; -sin(thetaytank) 0 cos(thetaytank)];
R3tank = [cos(thetaztank) -sin(thetaztank) 0; sin(thetaztank) cos(thetaztank) 0; 0 0 1];
Rtank = R1tank*R2tank*R3tank

Rthetaz = [cos(thetazgyro) -sin(thetazgyro) 0; sin(thetazgyro) cos(thetazgyro) 0; 0 0 1];
Rthetay = [cos(thetaygyro) 0 sin(thetaygyro); 0 1 0; -sin(thetaygyro) 0 cos(thetaygyro)];
Rthetax = [1 0 0; 0 cos(thetaxgyro) -sin(thetaxgyro); 0 sin(thetaxgyro) cos(thetaxgyro)];
Rhalfthetaz = [cos((1/2)*thetazgyro) -sin((1/2)*thetazgyro) 0; sin((1/2)*thetazgyro) cos((1/2)*thetazgyro) 0; 0 0 1];
Rhalfthetay = [cos((1/2)*thetaygyro) 0 sin((1/2)*thetaygyro); 0 1 0; -sin((1/2)*thetaygyro) 0 cos((1/2)*thetaygyro)];

Rchassis1 = simplify(Rtank*Rthetaz)
Rchassis2 = simplify(Rchassis1*Rthetay)
Rgyro = simplify(Rchassis2*Rthetax)
Rmirror = simplify(Rtank*Rhalfthetaz*Rhalfthetay)

R = Rtank;
S = diff(R(t),t)*transpose(R(t));
OmegaTank = simplify(expand([S(3,2); S(1,3); S(2,1)]))

R = Rchassis1;
S = diff(R(t),t)*transpose(R(t));
OmegaChassis1 = simplify(expand([S(3,2); S(1,3); S(2,1)]))

R = Rchassis2;
S = diff(R(t),t)*transpose(R(t));
OmegaChassis2 = simplify(expand([S(3,2); S(1,3); S(2,1)]))

R = Rgyro;
S = diff(R(t),t)*transpose(R(t));
OmegaGyro = simplify(expand([S(3,2); S(1,3); S(2,1)]))

R = Rmirror;
S = diff(R(t),t)*transpose(R(t));
OmegaMirror = simplify(expand([S(3,2); S(1,3); S(2,1)]))

I0chassis1 = [Ic1x 0 0; 0 Ic1y 0; 0 0 Ic1z];
Ichassis1 = simplify(expand(Rchassis1*I0chassis1*transpose(Rchassis1)))

I0chassis2 = [Ic2x 0 0; 0 Ic2y 0; 0 0 Ic2z];
Ichassis2 = simplify(expand(Rchassis2*I0chassis2*transpose(Rchassis2)))

I0gyro = [Igax 0 0; 0 Igtr 0; 0 0 Igtr];
Igyro = simplify(expand(Rgyro*I0gyro*transpose(Rgyro)))

Imirror = [Im 0 0; 0 Im 0; 0 0 Im];

KEchassis1 = (1/2)*transpose(OmegaChassis1)*Ichassis1*OmegaChassis1;
KEchassis1 = simplify(expand(KEchassis1))
KEchassis2 = (1/2)*transpose(OmegaChassis2)*Ichassis2*OmegaChassis2;
KEchassis2 = simplify(expand(KEchassis2))
KEgyro = (1/2)*transpose(OmegaGyro)*Igyro*OmegaGyro;
KEgyro = simplify(expand(KEgyro))
KEmirror = (1/2)*transpose(OmegaMirror)*Imirror*OmegaMirror;
KEmirror = simplify(expand(KEmirror))

KEtotal = simplify(KEchassis1 + KEchassis2 + KEgyro + KEmirror);
KEtotal = subs(KEtotal, thetaxgyro, spinspeed*t)

KEsmallangle = subs(KEtotal, [thetaxtank, thetaytank, thetaztank, thetaygyro, thetazgyro], epsilon.*[thetaxtank, thetaytank, thetaztank, thetaygyro, thetazgyro])

KEapprox = subs(taylor(KEsmallangle, epsilon, 'Order', ModelOrder+2), epsilon, 1)

if ModelOrder>0
    disp(['BUILDING MODEL OF ORDER  ', ModelOrder]);
    L = KEapprox;
else
    disp('BUILDING FULL NONLINEAR MODEL')
    L = KEtotal;
end

L = subs(L, [thetaygyro, thetazgyro, diff(thetaygyro,t), diff(thetazgyro,t)], [thy, thz, thydot, thzdot])

delLdelq = [diff(L, thy); diff(L, thz)]

delLdelqdot = subs([diff(L, thydot); diff(L, thzdot)], [thy, thz, thydot, thzdot], [thetaygyro, thetazgyro, diff(thetaygyro,t), diff(thetazgyro,t)])

dbydtofdelLdelqdot = subs(subs(subs(diff(delLdelqdot, t), [diff(thetaygyro, t, 2), diff(thetazgyro, t, 2)], [thyddot, thzddot]), [diff(thetaygyro, t), diff(thetazgyro, t)], [thydot, thzdot]), [thetaygyro, thetazgyro], [thy, thz])

gendampingforces = [thydot*C2damping; thzdot*C1damping];

EoM = dbydtofdelLdelqdot - delLdelq == - gendampingforces;

EoM = subs(EoM, [diff(thetaxtank, t, t), diff(thetaytank, t, t), diff(thetaztank, t, t)], [thxTddot thyTddot thzTddot]);
EoM = subs(EoM, [diff(thetaxtank, t), diff(thetaytank, t), diff(thetaztank, t)], [thxTdot thyTdot thzTdot]);
EoM = subs(EoM, [thetaxtank, thetaytank, thetaztank], [thxT thyT thzT])

[MatrixA, MatrixB] = equationsToMatrix(EoM, [thyddot, thzddot]);

matlabFunction(MatrixA, 'File', 'gyromatrixA', 'Vars', ...
    [thxT, thyT, thzT, thxTdot, thyTdot, thzTdot, ...
    thxTddot, thyTddot, thzTddot, thy, thz, thydot, thzdot, ...
    Ic1x, Ic1y, Ic1z, Ic2x, Ic2y, Ic2z, Igax, Igtr, Im, ...
    C1damping, C2damping, spinspeed])

matlabFunction(MatrixB, 'File', 'gyromatrixB', 'Vars', ...
    [thxT, thyT, thzT, thxTdot, thyTdot, thzTdot, ...
    thxTddot, thyTddot, thzTddot, thy, thz, thydot, thzdot, ...
    Ic1x, Ic1y, Ic1z, Ic2x, Ic2y, Ic2z, Igax, Igtr, Im, ...
    C1damping, C2damping, spinspeed])