# Optical Bloch Solver

This suite of MATLAB code solves the optical Bloch equations (OBEs) for either the steady state or for evolution for a given time.  It is mainly designed for solving problems involving optical transitions in alkali metal atoms, although it could be extended to other systems.  It is *not* well-designed for solving problems involving time-varying parameters such as frequency sweeps or pulsed lasers - if you need that, I suggest meta-programming as the solution for these kinds of problems.

## Density Matrix Class

The main class for solving the time-evolution of the density matrix for a given quantum system is the `densityMatrix` class.  The main idea behind this class is that the user provides a "bare" Hamiltonian which is the Hamiltonian in the absence of coupling between levels due to external fields; a "coupling" Hamiltonian which describes the interaction with external fields; and a matrix of "decay" rates which indicate the rate at which populations and/or coherences decay between states.  With these supplied, the class can then calculate the Lindblad term representing the loss of coherence and transfer of populations due to the decay rates, and then combine that with the rate of change of the density matrix due to the unitary evolution, to get a total time derivative of the density matrix.  This time derivative is represented as a super-operator that is applied to the density matrix.  So if we flatten the density matrix from NxN to N^2x1, where N is the number of eigenstates, then the super-operator is now a matrix of size N^2xN^2.  This means that for constant parameter problems, such as constant laser fields with constant detunings, the solution for any time is then the matrix exponential of the super-operator multiplied by the time step, left-multiplied with the initial density matrix (flattened).

A couple of example scripts are included which show how to use the simply `densityMatrix` class to solve master-equation style problems.  These scripts are `exampleTwoLevel.m`, `exampleThreeLevel.m`, and `exampleThreeLevelLambda.m`.  The relevant code from `exampleTwoLevel.m` is reproduced below.
```
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
d.intConstField(100e-9,100e-6);
% Plot populations
figure(1);clf;
d.plotPopulations;
```
As written above, this example demonstrates Rabi flopping between two quantum states at a frequency given by the Rabi frequency.  First, the `densityMatrix` object is created with dimension 2.  The bare and coupling parts of the Hamiltonian are then generated, and in this example we assume on-resonance coupling.  The decay terms are set to zero to see perfect Rabi flopping.  To see how decay from one level to the other affects Rabi flopping, set
```
d.decay = [0 [some value];
           [some value] 0];
```
where `[some value]` is a relevant decay rate.  After the decay rates are set, the initial population is defined as being all in the first quantum state, and then the equations of motion are integrated from 0 to 100 us in 100 ns steps.  As we are assuming that the Hamiltonian is time-independent, the length of the time step does not affect the accuracy as a matrix exponential is used for the calculation.  

The three level examples showcase pure three-level Rabi flopping (`exampleThreeLevel.m`), as might be seen when driving transitions between states in the same F manifold at low magnetic fields in alkali metal atoms, and lambda-style dynamics such as electromagnetically-induced transparency (EIT) and Raman absorption (`exampleThreeLevelLambda.m`).  The latter file also demonstrates using the `solveSteadyState` method, which computes the steady state solution for the system under the current parameters.  In the lambda examples, this method shows both EIT and Raman lineshapes in their respective regimes.

## Optical System Class

Constructing 

