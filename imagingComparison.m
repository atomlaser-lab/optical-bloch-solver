%% Example script
%
%   This demonstrates how to use the optical-bloch-solver classes to solve
%   optical pumping problems
%
dt = 10e-9;
T = 100e-6;
op = opticalSystem('Rb87','D2');
th = pi*0;ph = 0;
I = 1*16;
B = 10;

%% Pump from |1,-1> to F = 2 manifold
op.laser1.setIntensity(I)...
    .setPolarization([1,0,0],'linear')...
    .setStates([1,-1],[2,-2],0);

op.setMagneticField(B,[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]);
op.initGroundPop = zeros(op.numGroundStates,1);
op.initGroundPop(3) = 1;
op.intConstField(dt,T);

P = op.getPopulations('ground');
Pf = P(:,end);

%% Imaging atoms already in |2,0> state
op.laser1.setStates([2,-2],[3,-3],0)...
    .setPolarization([1,0,0],'spherical');
op.initGroundPop = zeros(op.numGroundStates,1);
op.initGroundPop(6) = 1;
op.refresh.intConstField(dt,T);

figure(1);clf;
subplot(1,2,1);
op.plotPopulations('ground');
subplot(1,2,2);
op.plotScatteringRates;
grid on;
plot_format('Time [s]','Photon Scattering Rate [s^{-1}]','',10);

Photons = trapz(op.t,sum(op.getScatteringRates,1));

%% Image atoms after optical pumping from F = 1 manifold
op.initGroundPop = Pf;
op.refresh.intConstField(dt,T);

figure(2);clf;
subplot(1,2,1);
op.plotPopulations('ground');
subplot(1,2,2);
op.plotScatteringRates;
grid on;
plot_format('Time [s]','Photon Scattering Rate [s^{-1}]','',10);

Photons(2) = trapz(op.t,sum(op.getScatteringRates,1));

%%
fprintf(1,'Total Photons = %.2e\n',Photons);
fprintf(1,'Ratio = %.3f\n',Photons(1)/Photons(2));