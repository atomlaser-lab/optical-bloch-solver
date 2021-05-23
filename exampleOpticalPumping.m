%% Example script
%
%   This demonstrates how to use the optical-bloch-solver classes to solve
%   optical pumping problems
%
op = opticalSystem('Rb87','D2');
op.laser1.setIntensity(1e-2*16)...
    .setPolarization([1,0,0],'spherical')...
    .setStates([2,-1],[3,-1],0);

th = 0;ph = 0;
op.setMagneticField(3,[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]);
op.initGroundPop = zeros(op.numGroundStates,1);
op.initGroundPop(4:5) = [0.53,0.327];
op.intConstField(1e-9,100e-6);

%%
figure(1);clf;
op.plotPopulations('ground');

figure(2);clf;
op.plotScatteringRates;
grid on;
plot_format('Time [s]','Photon Scattering Rate [s^{-1}]','',10);

fprintf(1,'Total Photons = %.2e\n',trapz(op.t,sum(op.getScatteringRates,1)));