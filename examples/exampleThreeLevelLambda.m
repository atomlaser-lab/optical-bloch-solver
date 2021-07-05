%% Three-level lambda style examples
%
%   These examples show how to use the basic densityMatrix.m class to
%   simulate a basic three-level lambda-style system under different
%   conditions
%

%% Steady-state EIT
%
% This example solves the system under steady-state EIT conditions
%

%Define Rabi frequencies and detunings
rabip = 2*pi*2e5;
rabic = 2*pi*2e6;
deltac = 2*pi*0;
deltap = 2*pi*0;
% Create densityMatrix object with 3 states
d = densityMatrix(3);
% Bare Hamiltonian 
d.bare = [0     0                   0;
          0     deltac - deltap     0;
          0     0                   -deltap];
% Coupling to external fields
d.coupling = [0                0                  -rabip/2;
              0                0                  -rabic/2;
              -conj(rabip)/2   -conj(rabic)/2     0];
% Decay rates
g = 2*pi*6e6;
d.decay = [0 0 g;
           0 0 g;
           g g 0];

%Define probe frequency detuning from Raman resonance
f = 2*pi*linspace(-20e6,20e6,5e3);
P = zeros(3,numel(f));
susc = zeros(numel(f),1);
for nn = 1:numel(f)
    d.bare(2,2) = deltac - deltap - f(nn);
    d.bare(3,3) = -deltap - f(nn);
    d.solveSteadyState;
    P(:,nn) = d.getPopulations;
    susc(nn) = d.density(3,1);
end

%Plot the susceptibility
figure(1);clf;
plot(f/(2*pi*1e6),[real(susc),imag(susc)]/rabip,'.-');
xlabel('Frequency [MHz]');
ylabel('\rho_{31}/\Omega_p');
title('EIT Spectrum');
legend('Real(\rho_{31})','Imag(\rho_{31})');
grid on;

%% Steady-state Raman absorption
%
% This example solves the system under steady-state Raman absorption
% conditions, where the single photon detunings are very large compared t
% the natural linewidth
%

%
% Define Rabi frequencies and detunings - Raman resonance is at a two-photon
% detuning of 0.25*rabic^2/deltac
%
rabip = 2*pi*2e3;
rabic = 2*pi*2e7;
deltac = 2*pi*5e9;
deltap = 2*pi*5e9 + (0.5*rabic)^2/deltac;
% Create densityMatrix object with 3 states
d = densityMatrix(3);
% Bare Hamiltonian 
d.bare = [0     0                   0;
          0     deltac - deltap     0;
          0     0                   -deltap];
% Coupling to external fields
d.coupling = [0                0                  -rabip/2;
              0                0                  -rabic/2;
              -conj(rabip)/2   -conj(rabic)/2     0];
% Decay rates
g = 2*pi*6e6;
d.decay = [0 0 g;
           0 0 g;
           g g 0];

%Define probe frequency detuning from Raman resonance
f = 2*pi*linspace(-1e3,1e3,1e3);
P = zeros(3,numel(f));
susc = zeros(numel(f),1);
for nn = 1:numel(f)
    d.bare(2,2) = deltac - deltap - f(nn);
    d.bare(3,3) = -deltap - f(nn);
    d.solveSteadyState;
    P(:,nn) = d.getPopulations;
    susc(nn) = d.density(3,1);
end

%Plot the susceptibility
figure(2);clf;
plot(f/(2*pi*1e3),[real(susc),imag(susc)]/rabip,'.-');
xlabel('Frequency relative to Raman resonance [kHz]');
ylabel('\rho_{31}/\Omega_p');
title('Raman Spectrum');
legend('Real(\rho_{31})','Imag(\rho_{31})');
grid on;

%% Rabi flopping under Raman absorption
%
% This example solves the time-dependent system under conditions suitable
% for Raman absorption, showing that Raman absorption is a coherent process
%

%Define Rabi frequencies and detunings
rabip = 2*pi*2e5;
rabic = 2*pi*2e7;
deltac = 2*pi*5e9;
deltap = 2*pi*5e9 + (0.5*rabic)^2/deltac;
% Create densityMatrix object with 3 states
d = densityMatrix(3);
% Bare Hamiltonian 
d.bare = [0     0                   0;
          0     deltac - deltap     0;
          0     0                   -deltap];
% Coupling to external fields
d.coupling = [0                0                  -rabip/2;
              0                0                  -rabic/2;
              -conj(rabip)/2   -conj(rabic)/2     0];
% Decay rates
g = 2*pi*6e6;
d.decay = [0 0 g;
           0 0 g;
           g g 0];
% Set the initial state to be all atoms in state 1
d.initPop(1) = 1;
% Integrate with a time step of 100 ns to a maximum time of 100 us
d.intConstField(10e-6,20e-3);
% Plot populations - expected Rabi frequency is rabip*rabic/deltac^2
figure(3);clf;
d.plotPopulations;


