%% Example script
%
%   This demonstrates how to use the optical-bloch-solver classes to see
%   Raman spectroscopy line shapes for square pulses. This example shows
%   Rabi flopping between the clock states of Rb-87. This will fit the
%   resulting line shape to get the Rabi frequeny and the offset due to the
%   AC Stark shift
%
op = opticalSystem('Rb87','D2');
op.laser1.setGaussBeam(50e-3,10e-3)...
    .setPolarization([0,0,1],'spherical')...
    .setStates([2,0],[2,0],-2e3);
op.laser2.setGaussBeam(50e-3,10e-3)...
    .setPolarization([0,0,1],'spherical')...
    .setStates([1,0],[2,0],-2e3); %Note the AC Stark shift

th = 0*pi/180;
ph = 0;
op.setMagneticField(250e-3,[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]);
op.initPop(2) = 1;

tmp = 2*pi*1e2*ones(8,8)*0;
tmp = tmp - diag(diag(tmp));
op.decay(1:8,1:8) = tmp;


%% Solve for each detuning
f2 = (-1000:20:1000)*1e-3;
tau = 8e-6;
P = zeros(op.transition.ground.numStates,numel(f2));
for nn = 1:numel(f2)
    op.laser2.detuning = op.laser1.detuning + f2(nn);
    op.calcBareH.integrate(tau,tau);
    tmp = op.getPopulations('ground');
    P(:,nn) = tmp(:,end);
end

%% Plot
figure(1);clf;
ax = gca;
grid on;
set(ax,'NextPlot','ReplaceChildren','LineStyleOrder',{'-','--','.-'})
plot(ax,f2*1e3,P,'linewidth',0.5,'marker','.');
str = op.getPopLegend('ground');
legend(ax,str);
xlabel('Frequency [kHz]');
ylabel('Population');

%% Fit
nlf = nonlinfit(f2,P(2,:));
nlf.useErr = false;
nlf.setFitFunc(@(A,R,x0,x) A*(1 - 4*R.^2./(4*R.^2+(x-x0).^2).*sin(sqrt(4*R.^2+(x-x0).^2)*2*pi*tau*1e6/2).^2));
nlf.bounds([0,0,-5],[1,10,5],[0.95,0.025,0]);
nlf.fit;
hold(ax,'on');
plot(ax,nlf.x*1e3,nlf.f(nlf.x),'k-');
hold(ax,'off');
str{9} = 'Fit';
legend(ax,str);
fprintf(1,'Rabi Freq = %.3e kHz, Center = %.3e kHz\n',nlf.c(2,1)*1e3,nlf.c(3,1)*1e3);