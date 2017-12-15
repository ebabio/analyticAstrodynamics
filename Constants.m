%% AAE 532 Table of Constants

%% Global constants

% as returned by google search
G = 6.67408e-11;
AU2meters = 149597870700;
J2000 = juliandate(datetime('2000-01-01 12:00:00'));

%% Orbit characteristics

% as listed in http://www.met.rdg.ac.uk/~ross/Astronomy/Planets.html
% except for gravitational parameters as in https://en.wikipedia.org/wiki/Standard_gravitational_parameter

Sun.gravP   = 1.32712440018e20;
Sun.a       = 0;
Sun.radius = 695700e3;

Mercury.gravP   = 2.2032e13;
Mercury.a   = 0.38709893*AU2meters;
Mercury.e = .2056;
Mercury.lonAscNode = deg2rad(48.33167);
Mercury.inclination = deg2rad(7.00487);
Mercury.lonPeriapsis = deg2rad(77.45645);
Mercury.meanLon = deg2rad(252.25084);

Venus.gravP     = 3.24859e14;
Venus.a     = 0.72333199*AU2meters;
Venus.e = .0068;
Venus.lonAscNode = deg2rad(76.68069);
Venus.inclination = deg2rad(3.39471);
Venus.lonPeriapsis = deg2rad(131.53298);
Venus.meanLon = deg2rad(181.97973);

Earth.gravP     = 3.986004418e14;
Earth.a     = 1.00000011*AU2meters;
Earth.e = 0.0167;
Earth.radius = earthRadius;
Earth.lonAscNode = deg2rad(-11.26064);
Earth.inclination = deg2rad(0.00005);
Earth.lonPeriapsis = deg2rad(102.94719);
Earth.meanLon = deg2rad(100.46435);

Mars.gravP  = 4.282837e13;
Mars.a      = 1.52366231*AU2meters;
Mars.e = 0.0934;
Mars.radius = 3390e3;
Mars.lonAscNode = deg2rad(49.57854);
Mars.inclination = deg2rad(1.85061);
Mars.lonPeriapsis = deg2rad(336.04084);
Mars.meanLon = deg2rad(355.45332);

Jupiter.gravP	= 1.26686534e17;
Jupiter.a	= 5.20336301*AU2meters;
Jupiter.e = .0484;
Jupiter.lonAscNode = deg2rad(100.55615);
Jupiter.inclination = deg2rad(1.30530);
Jupiter.lonPeriapsis = deg2rad(14.75385);
Jupiter.meanLon = deg2rad(34.40438);

Saturn.gravP	= 3.7931187e16;
Saturn.a	= 9.53707032*AU2meters;
Saturn.e = .0542;

Uranus.gravP	= 5.793939e15;
Uranus.a	= 19.19126393*AU2meters;
Uranus.e = .0472;

Neptune.gravP	= 6.836529e15;
Neptune.a	= 30.06896348*AU2meters;
Neptune.e = .0086;

Pluto.gravP     = 8.71e11;
Pluto.a     = 39.48168677*AU2meters;
Pluto.e = 0.2488;

Moon.gravP      = 4.9048695e12;
Moon.radius = 1737e3;
