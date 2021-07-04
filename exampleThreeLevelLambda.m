%% Example three-level Rabi flopping
%
%   This example shows how to use the basic densityMatrix.m class to
%   simulate basic three-level system
%

%Define Rabi frequencies and detunings
rabip = 2*pi*1e6;
rabic = 2*pi*1e6;
deltac = 2*pi*0;
deltap = 2*pi*0;
% Create densityMatrix object with 3 states
d = densityMatrix(3);
% Bare Hamiltonian 
d.bare = [0     0                   0;
          0     deltac - deltap     0;
          0     0                   -deltap];
% Coupling to external fields
d.coupling = [0              0                  rabip/2;
              0              0                  rabic/2;
              conj(rabip)/2   conj(rabic)/2     0];
% Decay rates
g = 2*pi*6e6;
d.decay = [0 0 g;
           0 0 g;
           g g 0];
% Set the initial state to be all atoms in state 1
d.initPop(1) = 1;
% Integrate with a time step of 100 ns to a maximum time of 100 us
d.intConstField(1e-9,2e-6);
% Plot populations
figure(1);clf;
d.plotPopulations;

