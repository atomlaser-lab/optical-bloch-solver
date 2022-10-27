%% Example script
%
%   This demonstrates how to use the optical-bloch-solver classes to see
%   Raman spectroscopy line shapes for square pulses. This example shows
%   Rabi flopping between the clock states of Rb-87. This will fit the
%   resulting line shape to get the Rabi frequency and the offset due to the
%   AC Stark shift
%
op = opticalSystem('Rb87','D2');
op.laser1.setGaussBeam(15e-3,5e-3)...
    .setPolarization([0,0,1],'spherical')...
    .setStates([2,0],[2,0],15e9);
op.laser2.setGaussBeam(15e-3,5e-3)...
    .setPolarization([0,0,1],'spherical')...
    .setStates([1,0],[2,0],15e9); %Note the AC Stark shift

th = 0*pi/180;
ph = 0;
op.setMagneticField(5.5,[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]);
op.initPop(2) = 1;
%
% This bit of code attempts to calculate optical pumping rates and
% therefore decoherence associated with each laser addressing the other
% ground state manifold
%
R1 = op.getOffResonantPumping(op.laser1);
R2 = op.getOffResonantPumping(op.laser2);
R = R1 + R2;
op.makeTotalLindblad;
for gg = 1:op.transition.ground.numStates
    for ee = 1:op.transition.excited.numStates
        eShift = ee + op.transition.ground.numStates;
        op.decay(eShift,gg) = R(eShift,gg);
        op.lindblad = op.lindblad + op.makeLindblad(eShift,gg);
    end
end


%% Solve for each detuning
f2 = (-50:2:50)*1e3 - const.muBh*op.B*1e6*0;
tau = 50e-6;
P = zeros(op.transition.ground.numStates,numel(f2));
R = zeros(numel(f2),1);
for nn = 1:numel(f2)
    op.laser2.detuning = op.laser1.detuning + f2(nn);
    op.calcBareH.integrate(tau,tau);
    tmp = op.getPopulations('ground');
    P(:,nn) = tmp(:,end);
    tmp = sum(op.getScatteringRates,1);
    R(nn) = tau*tmp(end);
end

%% Plot
figure(3);clf;
ax = gca;
grid on;
set(ax,'NextPlot','ReplaceChildren','LineStyleOrder',{'-','--','.-'})
plot(ax,f2/1e3,P,'linewidth',0.5,'marker','.');
str = op.getPopLegend('ground');
legend(ax,str);
xlabel('Frequency [kHz]');
ylabel('Population');

%% Fit
nlf = nonlinfit(f2/1e6,P(2,:));
nlf.useErr = false;
nlf.setFitFunc(@(A,R,x0,x) A*(1 - 4*R.^2./(4*R.^2+(x-x0).^2).*sin(sqrt(4*R.^2+(x-x0).^2)*2*pi*tau*1e6/2).^2));
[m,idx] = min(nlf.y);
Rguess = asin(sqrt(1 - m))/(2*pi*tau*1e6);
nlf.bounds2('A',[0.9,1,0.97],'R',[0,2,Rguess],'x0',[min(nlf.x),max(nlf.x),nlf.x(idx)]);
nlf.fit;
hold(ax,'on');
plot(ax,nlf.x*1e3,nlf.f(nlf.x),'k-','linewidth',1);
hold(ax,'off');
str{9} = 'Fit';
legend(ax,str);
fprintf(1,'Rabi Freq = %.3e kHz, Center = %.3e kHz\n',nlf.c(2,1)*1e3,nlf.c(3,1)*1e3);