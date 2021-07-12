%% Example two-level Rabi flopping
%
%   This example shows how to use the basic densityMatrix.m class to
%   simulate very basic two-level Rabi flopping
%

% Create densityMatrix object with 2 states
d = densityMatrix(2);
% Bare Hamiltonian 
d.bare = [0 0;
          0 0];
% Coupling to external fields
rabi = 2*pi*10e3;
d.coupling = [0              rabi/2;
              conj(rabi)/2   0];
% Decay rates
d.decay = [0 0;
           0 0];
% Set the initial state to be all atoms in state 1
d.initPop(1) = 1;
% Integrate with a time step of 100 ns to a maximum time of 100 us
d.integrate(100e-9,100e-6);
% Plot populations
figure(1);clf;
d.plotPopulations;

