%% Example script
%
%   This demonstrates how to use the optical-bloch-solver classes to solve
%   optical pumping problems
%
op = opticalSystem('Rb87','D2');
op.laser1.setGaussBeam(100e-3,10e-3)...
    .setPolarization([1,0,0],'spherical')...
    .setStates([2,0],[2,0],-1e3);
op.laser2.setGaussBeam(100e-3,10e-3)...
    .setPolarization([0,0,1],'spherical')...
    .setStates([1,0],[2,0],-1e3+0.075);

th = 90*pi/180;
ph = 0;
op.setMagneticField(1,[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]);
op.initGroundPop = zeros(op.numGroundStates,1);
op.initGroundPop(2) = 1;

tmp = 2*pi*1e2*ones(8,8)*0;
tmp = tmp - diag(diag(tmp));
op.decay(1:8,1:8) = tmp;


%% Solve for each detuning
f2 = (-5000:20:5000)*1e-3;
P = zeros(op.numGroundStates,numel(f2));
for nn = 1:numel(f2)
    op.laser2.detuning = -1e3 + f2(nn);
    op.calcBareH.intConstField(1e-6,1e-6);
    tmp = op.getPopulations('ground');
    P(:,nn) = tmp(:,end);
end

%% Plot
figure(1);clf;
ax = gca;
grid on;
set(ax,'NextPlot','ReplaceChildren','LineStyleOrder',{'-','--','.-'})
plot(ax,f2*1e3,P,'linewidth',2);
legend(ax,op.getPopLegend('ground'));
xlabel('Frequency [kHz]');
ylabel('Population');