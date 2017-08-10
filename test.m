%% Test script

%% Setup workspace
clear
clear orbit state
clc
close all

Constants

%% Test1

venusOrbit = orbit(Venus, Sun);
earthOrbit = orbit(Earth, Sun);
marsOrbit = orbit(Mars, Sun);
jupiterOrbit = orbit(Jupiter, Sun);

venus0 = venusOrbit.computeState(0);
earth0 = earthOrbit.computeState(0);
mars0 = marsOrbit.computeState(0);
jupiter0 = jupiterOrbit.computeState(0);

venusOrbit = orbit(venus0, 0, Sun);
earthOrbit = orbit(Earth.lonAscNode, Earth.inclination, Earth.lonPeriapsis,...
    Earth.a, Earth.e, ...
    Earth.meanLon - Earth.lonAscNode - Earth.lonPeriapsis,...
    Sun);

figure
quiverScale = 1e7;
venusOrbit.plot()
venus0.plot(quiverScale)
earthOrbit.plot()
earth0.plot(quiverScale)
marsOrbit.plot()
mars0.plot(quiverScale)
jupiterOrbit.plot()
jupiter0.plot(quiverScale);




%% Test2

venusOrbit = orbit(venus0, 0, Sun);
earthOrbit = orbit(Earth.lonAscNode, Earth.inclination, Earth.lonPeriapsis,...
    Earth.a, Earth.e, ...
    Earth.meanLon - Earth.lonAscNode - Earth.lonPeriapsis,...
    earthOrbit);

figure
quiverScale = 3e6;
venusOrbit.plot()
venus0.plot(quiverScale)
earthOrbit.plot()
earth0.plot(quiverScale)
marsOrbit.plot()
mars0.plot(quiverScale)
jupiterOrbit.plot()
jupiter0.plot(quiverScale);

%% Test3

transfer = transferSolver(earthOrbit, marsOrbit);
transfer = transfer.computePorkChop(200:1:600,150:1:360, 1);

figure
transfer.plot()

%% Test4

transfer = transferSolver(earthOrbit, marsOrbit);
transfer = transfer.solveFixedTransfer(455,215);

flyby = flyBySolver(Mars, marsOrbit);
flyby = flyby.solveFlyBy(transfer.arrivingState, [0; 0; 6e6]);

figure
transfer.plot()

figure
transfer.departureOrbit.plot()
transfer.transferOrbit.plot()
transfer.initialState.plot(quiverScale)
transfer.leavingState.plot(quiverScale)
transfer.arrivingState.plot(quiverScale)
flyby.outboundOrbit.plot()
flyby.outboundState.plot(quiverScale)

figure
flyby.localOrbit.plot()