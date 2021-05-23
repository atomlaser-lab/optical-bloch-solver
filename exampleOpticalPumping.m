%% Example script
%
%   This demonstrates how to use the optical-bloch-solver classes to solve
%   optical pumping problems
%
op = opticalSystem('Rb87','D2');
op.laser1.setIntensity(1e0*16)...
    .setPolarization([1,0,0],'spherical')...
    .setStates([2,-1],[3,-1],0);

th = 0;ph = 0;
op.setMagneticField(3,[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]);
op.initGroundPop = zeros(op.numGroundStates,1);
op.initGroundPop(4:5) = [0.53,0.327];
op.intConstField(10e-9,100e-6);

%%
figure(1);clf;
op.plotPopulations('ground');

figure(2);clf;
op.plotScatteringRates;
grid on;
plot_format('Time [s]','Photon Scattering Rate [s^{-1}]','',10);

fprintf(1,'Total Photons = %.2e\n',trapz(op.t,sum(op.getScatteringRates,1)));

figure(3);clf;
basis = 'spherical';
if strcmpi(basis,'spherical')
    str = {'-1','0','+1'};
else
    str = {'x','y','z'};
end
P = op.getPolarization(basis);
subplot(2,1,1);
plot(op.t,real(P),'-','linewidth',2);
legend(str);
plot_format('Time [s]','Absorbance','Absorbance',10);
subplot(2,1,2);
plot(op.t,imag(P),'-','linewidth',2);
legend(str);
plot_format('Time [s]','Reactance','Reactance',10);