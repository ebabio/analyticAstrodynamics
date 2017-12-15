%% Test script

%% Setup workspace
clear
clear orbit state
clc
Constants

% Add path
path = [pwd,'/library']; 
addpath(path);

%% Test 1: orbit visualization
% using predefined solar system body orbits to visualize the 3d orbits
mercuryOrbit = orbit(Mercury, Sun);
venusOrbit = orbit(Venus, Sun);
earthOrbit = orbit(Earth, Sun);
marsOrbit = orbit(Mars, Sun);
jupiterOrbit = orbit(Jupiter, Sun);

mercury0 = mercuryOrbit.computeState(0);
venus0 = venusOrbit.computeState(0);
earth0 = earthOrbit.computeState(0);
mars0 = marsOrbit.computeState(0);
jupiter0 = jupiterOrbit.computeState(0);

earthOrbit = orbit(Earth.lonAscNode, Earth.inclination, Earth.lonPeriapsis,...
    Earth.a, Earth.e, ...
    Earth.meanLon - Earth.lonAscNode - Earth.lonPeriapsis,...
    Sun);

figure(1)
clf reset
quiverScale = 1e7;
mercuryOrbit.plot()
mercury0.plot(quiverScale)
venusOrbit.plot()
venus0.plot(quiverScale)
earthOrbit.plot()
earth0.plot(quiverScale)
marsOrbit.plot()
mars0.plot(quiverScale)
jupiterOrbit.plot()
jupiter0.plot(quiverScale)
title('Inner planets and Jupiter orbits')

%% Test 2: porkchop plots
% compute porkchop plot for a 2-impulse transfer from Earth to Mars using
% lambert's theorem

figure(2)
transfer = transferSolver(earthOrbit, marsOrbit);
transfer = transfer.computePorkChop(200:1:600,150:1:360, 1);

figure(3)
clf reset
transfer.plot()
title('Porkchop plot for Earth to Mars')
saveas(3, 'results\Earth2MarsPorkchop', 'png')

%% Test 3: patched conics fly-by
% fly-by from Earth to Mars using a predefined transfer and the resulting
% orbit

transfer = transferSolver(earthOrbit, marsOrbit);
transfer = transfer.solveFixedTransfer(455,215);

flyby = flyBySolver(Mars, marsOrbit);
flyby = flyby.solveFlyBy(transfer.arrivingState, [0; 0; 6e6]);

figure(4)
transfer.plot()
title('Chosen transfer to Mars vicinity')

figure(5)
transfer.departureOrbit.plot()
transfer.transferOrbit.plot()
transfer.initialState.plot(quiverScale)
transfer.leavingState.plot(quiverScale)
transfer.arrivingState.plot(quiverScale)
flyby.outboundOrbit.plot()
flyby.outboundState.plot(quiverScale)
title('Resulting orbit from a low Mars fly-by')
saveas(5, 'results\Earth2MarsFlybyTrajectory', 'png')

figure(6)
flyby.localOrbit.plot()
title('Local orbit for Mars fly-by')
saveas(6, 'results\MarsFlyBy', 'png')