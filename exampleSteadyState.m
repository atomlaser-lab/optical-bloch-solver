%% Example script
%
%   This demonstrates how to use the optical-bloch-solver classes to solve
%   for steady-state spectra
%
op = opticalSystem('Rb87','D2');
op.laser1.setIntensity(1e-2*16)...
    .setPolarization([1,0,0],'linear')...
    .setStates([2,0],[3,0],0);
op.laser2.setIntensity(1e-2*16)...
    .setPolarization([1,0,0],'linear')...
    .setStates([1,0],[2,0],0);

th = 0;ph = 0;
op.setMagneticField(1e-3,[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]);

D = linspace(-1e3,1e3,100);
Pg = zeros(op.ground.numStates,numel(D));
Pe = zeros(op.excited.numStates,numel(D));
for nn = 1:numel(D)
    op.laser1.detuning = D(nn);
    op.refresh.solveSteadyState;
    Pg(:,nn) = op.getPopulations('ground');
    Pe(:,nn) = op.getPopulations('excited');
end

%%
figure(1);clf;
ax = gca;
grid on;
set(gca,'NextPlot','ReplaceChildren','LineStyleOrder',{'-','--','.-'})
plot(D,Pg);
legend(op.getPopLegend('ground'));