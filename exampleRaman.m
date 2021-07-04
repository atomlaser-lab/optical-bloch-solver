%% Example script
%
%   This demonstrates how to use the optical-bloch-solver classes to solve
%   optical pumping problems
%
op = opticalSystem('Rb87','D2');
op.laser1.setGaussBeam(160e-3,10e-3)...
    .setPolarization([1,0,0],'spherical')...
    .setStates([2,0],[2,0],-1e3);
op.laser2.setGaussBeam(160e-3,10e-3)...
    .setPolarization([0,0,1],'spherical')...
    .setStates([1,-1],[2,0],-1e3+0.075);

th = 90*pi/180;
ph = 0;
op.setMagneticField(1,[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]);
op.initGroundPop = zeros(op.numGroundStates,1);
op.initGroundPop(3) = 1;

tmp = 2*pi*1e2*ones(8,8)*0;
tmp = tmp - diag(diag(tmp));
op.decay(1:8,1:8) = tmp;

op.intConstField(0.1e-6,50e-6);

t0 = op.t;
%%
op.t = t0*1e6;
figure(1);clf;
subplot(2,1,1);
op.plotPopulations('ground');
xlabel('Time [us]');
xlim([0,Inf]);
subplot(2,1,2);
P = op.getPopulations('ground');
plot(op.t,sum(P(1:3,:),1),'k-','linewidth',2);
hold on;
plot(op.t,sum(P(4:8,:),1),'r-','linewidth',2);
plot_format('Time [us]','Populations','',10);
legend('F = 1','F = 2');
xlim([0,Inf]);

% figure(2);clf;
% op.plotScatteringRates;
% grid on;
% plot_format('Time [s]','Photon Scattering Rate [s^{-1}]','',10);
% 
% fprintf(1,'Total Photons = %.2e\n',trapz(op.t,sum(op.getScatteringRates,1)));