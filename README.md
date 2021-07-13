# Optical Bloch Solver

This suite of MATLAB code solves the optical Bloch equations (OBEs) for either the steady state or for evolution for a given time.  It is mainly designed for solving problems involving optical transitions in alkali metal atoms, although it could be extended to other systems.  It is *not* well-designed for solving problems involving time-varying parameters such as frequency sweeps or pulsed lasers - if you need that, I suggest meta-programming as the solution for these kinds of problems.

### Dependencies

You should have the `const.m` field in your MATLAB path.  Obtain from https://github.com/ryan-james-thomas/matlab-common-functions.

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

When represented as a matrix equation, the time derivative of the density matrix is $\dot{\rho} = M\rho$ where we represent $\rho$ by a vector that is $N^2\times 1$ with $N$ the number of states or dimension of the system.  This means that the matrix $M$ is $N^2\times N^2$.  If we are only working with a two-level system, it is relatively straightforward to write down $M$ directly.  It's annoying but not too much trouble to write down $M$ for a three-level system, even though $M$ is now $9\times 9$.  If one wants to solve the optical Bloch equations for the full $^{87}$Rb D2 line, which has $N=24$, then one needs to write down all $N^2\times N^2 = 331776$ elements of $M$ by hand.  For alkali metal atoms with higher nuclear spin, such as $^{40}$K, this involves $N=54$ and thus $M$ is $2916\times 2916$.  This is clearly not something one wants to write down by hand.

The difficulty of writing down the solutions these equations is why this suite of code was developed.  While the `densityMatrix` class can be used by itself, most of the code is centered around the remaining classes which represent the relevant properties of optical transitions in certain alkali metal atoms.  The `opticalSystem` class, which is a subclass of `densityMatrix`, serves as the top-level class which contains instances of the remaining classes.  The `opticalSystem` class has the following properties
  - `laser1` and `laser2`: these are instances of the `laser` class and represent the laser fields incident on the system.  The assumption is that each laser interacts separately with the two ground state hyperfine manifolds.  This means that this suite of code cannot handle four-wave mixing problems.
  - 'transition': this is an instance of the `opticalTransition` class and contains the information that describes the optical transition, such as the type of atom and the ground and excited state hyperfine manifolds.
  - `B` and `Bdir`: a scalar and $3\times 1$ vector, respectively, representing the magnetic field in the system.  `B` is given in Gauss, as befits atomic physics, and `Bdir` is the direction of the field as a 3-component vector with components $(B_x,B_y,B_z)$.

One creates an `opticalSystem` object using
```
op = opticalSystem(atom,transition);
```
where `atom` is the atom type, given as a string.  Currently supported atom types are `'Rb87'`,`'K40'`, and `'K41'`, and the currently supported transitions are `'D1'` and `'D2'`.  The user then sets the laser parameters and the magnetic field and then solves for either the steady state solution of the constant $M$ solution.  As noted previously, varying parameters such as detuning or laser fields with time is not supported due to the complexity.  A number of examples are included in the `examples` folder.  A very simple example would be as follows
```
op = opticalSystem('Rb87','D2');                %Creates an optical system for the D2 line of Rb-87
op.laser1.setIntensity(1);                      %Sets the intensity in W/m^2
op.laser1.setPolarization([0,0,1],'spherical'); %Sets the polarization to by in x direction
op.laser1.setStates([2,2],[3,3],0);             %Tells the solver to assume that the field is always on resonance with the |2,2> -> |3,3> transition
op.setMagneticField(1);                         %Sets the magnetic field to 1 G in the z-direction
op.initPop(1:8) = 1;                            %Sets the initial population to all be equally distributed in all ground state levels
op.integrate(10e-9,20e-6);                      %Integrates the master equation with time steps of 10 ns up to a time of 20 us
op.plotPopulations('ground');                   %Plots all ground-state populations
```
The above example would calculate and plot the resulting population distribution of atoms after 20 us when a moderate strength laser field is applied with linear polarization in the $\hat{x}$ direction.  The initial density matrix has equal population in all ground states, and this shows optical pumping of atoms to the fully stretched state $|2,2\rangle$.  

The choices for the above options are rather opaque, so in the following sections we will go through the different options in detail.

### Setting laser parameters

The `laser` class describes the optical field incident on the atoms.  It has properties
  - `intensity`: the intensity in W/m^2 of the field
  - `field`: the electric field in V/m
  - `pol`: the normalized polarization in the linear basis as $(E_x,E_y,E_z)$.
  - `ground`: the $|F,m_F\rangle$ ground state to use as a reference for the detuning
  - `excited`: the $F,m_F\rangle$ excited state to use as a reference for the detuning
  - `detuning`: the detuning in Hz.

First, let's discuss the laser intensity and electric field.  The user should set the intensity using the `setIntensity(I)` method where the input argument `I` is the intensity.  When set this way, the electric field is automatically calculated.  Otherwise, the user can set the electric field directly.  A useful alternative is to use `setGaussbeam(P,w)` where `P` is the total optical power in W and `w` is the beam waist ($I = I_0\exp(-2r^2/w^2)$) in meters.  This calculates the peak intensity as $I = 2P/\pi w^2$ and sets the intensity and field based on the result.

The user should set the polarization using the method `setPolarization(pol,basis)`.  Here, `pol` is the 3 element vector specifying the polarization as seen along the optical axis, which is assumed to be along the $\hat{z}$ direction.  Note that the input `pol` is automatically normalized to unit length by the method.  If `basis` is omitted or set to `'linear'`, `pol` is assumed to be in the linear basis of $(\hat{x},\hat{y},\hat{z})$.  If `basis` is set to `'spherical'` then it is assumed to be in the spherical basis of $(\hat{e}_{-1},\hat{e}_0,\hat{e}_{+1})$ or left-handed circular, linear along $\hat{z}$, and right-handed circular polarization, respectively.  `pol` is stored only in the linear basis.

The final properties `ground`, `excited`, and `detuning` are more complicated.  Basically, in many calculations one wants to solve for time-dependent problems with a detuning relative to a particular transition.  However, one does not want to have to keep track of what the absolute frequency is if a magnetic field changes.  The properties `ground` and `excited` let the user fix what states the detuning should be calculated relative to.  Both `ground` and `excited` can be specified either as just the $F$ value or as a 2 element vector $|F,m_F\rangle$.  In `opticalSystem`, `laser` must have a ground state value set so that the solver knows how to calculate the rotating wave approximation; however, `ground` can be just an $F$ value.  If just an $F$ value is set, then the detuning is calculated with respect to the $B=0$ energy for the ground state and the fine structure frequency of the excited state.  If `ground` or `excited` is set to a 2 element vector, the detuning is calculated with respect to the relevant transition connecting the $|F,m_F\rangle\rightarrow |F\rangle$ or $|F,m_F\rangle\rightarrow|F',m_F'\rangle$ states.  Here are a couple of examples:
```
op.laser1.setStates(2,[],0);            %Laser 1 connects ground state with F = 2, and is detuned by 0 Hz from the F = 2 to F' transition
op.laser1.setStates([1,1],2,-1e9);      %Laser 1 connects the state |1,1> with the F' = 2 state with a detuning of -1 GHz. If the magnetic field is changed, the detuning relative to this transition is unchanged.
op.laser1.setStates([2,2],[3,3],0);     %Laser 1 connects the |2,2> -> |3,3> transition. Remains on resonance regardless of magnetic field
```

When using a second laser, simply set the `laser2` properties to have non-zero field strength and make sure that the ground state is a different $F$ value from the `laser` ground state property.  If `laser2` has no intensity or has no ground state property, then the solver assumes that `laser` couples to *both* ground states.

### Setting the magnetic field

The dynamics of the system change drastically depending on the relative alignment of the laser polarizations and the magnetic field direction.  This solver calculates states relative to the magnetic field direction.  This means that a right-hand circularly polarized beam with a magnetic field in the $\hat{x}$ direction will *not* pump atoms into the $|2,2\rangle$ ground state, since the decomposition of the light field in the basis along the magnetic field is not purely a right-hand circularly polarized field.  

The magnetic field should be set using the `setMagneticField(B,Bdir)` method of `opticalSystem`.  Here, `B` is the magnetic field in Gauss and `Bdir` is the direction as a 3-element vector.  `Bdir` is automatically normalized to unit length when using this method.  When this method is called, a rotation matrix is generated that is used to rotate the laser polarizations to the necessary frame.  It also automatically re-calculates the coupling of the laser fields to the states and the new detunings in the rotating wave approximation.

### Setting the initial density matrix

Generally speaking, one wants to solve these systems when the density matrix starts as a mixture of different state populations.  This can be done by setting the `initPop` value, which is an $N\times 1$ vector where $N$ is the total number of states.  This vector is organized such that ground states come first and then excited states.  Ground and excited states are ordered by increasing energy; this means that for Rb-87 the first state is $|F = 1,m_F = 1\rangle$, the second state is $|F = 1,m_F = 0\rangle$, the 6th state is $|F = 2,m_F = 0\rangle$, and the 20th state is $|F' = 3,m_F' = 0\rangle$.  Set the `initPop` vector to whatever initial mixture of populations you want; it is automatically normalized to unity.

### Solving the system of equations

One can solve the system of equations either as a function of time or for the steady state solution.  To solve as a function of time, use
```
op.integrate(dt,T);
```
which integrates the resulting system of equations up to time `T` using time steps `dt`.  Note that if we assume that the matrix `M` is constant then we can represent the solution of $\partial_t \rho = M\rho$ as $\rho(t) = \exp[M dt]\rho(0)$.  This involves no approximations, so `dt` can be set to any value less than or equal to `T`.  This can be useful if running simulations of Rabi spectroscopy, since `dt` can be set equal to `T` without loss of precision.

Alternatively, one can solve for the steady state solution using
```
op.solveSteadyState;
```
which solves for $\partial_t \rho = 0$.  This involves solving the system of equations with the additional constraint that all populations sum to 1.  Be warned that this can fail for very far detuned systems because of badly conditioned matrices, and even if it doesn't the results can be unexpected since the infinite time solution for many problems is rather different from intuition.

### Extracting information from the solution

Supposing that we have integrated the solutions as a function of time, how do we get the information out?  The solution is represented as the flattened density matrix $\mathbf{v}$ where $v_n = \rho_{ab}$ with the relationship between the indices being $n = a + (b-1)N$.  The flattened density matrix is stored in the property `densityVec`, which is an $N^2\times T$ matrix where $T$ is the number of time steps taken.  The density matrix can be recovered for each time step using `density = reshape(op.densityVec(:,idx),op.numStates,op.numStates);`.  

Populations can obviously be extracted from the re-constituted density matrix at each time, or one can use the method `getPopulation(type)` to retrieve them.  Here, `type` can be `'ground'`, `'excited'`, or `'all'` to retrieve the ground state, excited state, or all populations, respectively.  One can also use the `plotPopulations(type)` method to plot the populations as a function of time.

The population in the excited states gives the scattering rates when multiplied by the appropriate decay rate.  Use the method `getScatteringRates()` to return a matrix $N_{excited} \times T$ with the photon scattering rates for each excited state and for each time step.  Sum over the first dimension to get the total photon scattering rate.

Many times, though, we don't care just about the populations but also about the coherences, and especially how they impact the light fields.  From the Maxwell-Bloch equations, we know that the total polarization $\mathbf{P}$ is related to the sum over the upper-triangular part of the dipole matrix multiplied by the density matrix.  This functionality is conveniently provided as the method `getPolarization(basis,frame)`.  Here, `basis` is either `'linear'` or `'spherical'` and calculates the polarization in either the linear or spherical basis.  `frame` is the frame in which the polarization should be calculated.  If `frame` is `'lab'` (which is the default if `frame` isn't specified), then the polarization is calculated relative to the field propagation - this is what you would use as input to the Maxwell-Bloch equations.  If `frame` is `'field'` the polarization is calculated relative to the magnetic field direction.

The value returned by `getPolarization()` is a polarization per atom/m^3, so multiply by the number density to get a total polarization for a sample of atoms.


### Additional information

If you are interested in the various atomic parameters, several properties of the `opticalTransition` class may be of interest.  These are
  - `coupling`: the coupling matrix $\mathbf{d}/\langle J ||e\mathbf{r}|| J' \rangle$. Take a look at D. Steck's Rubidium reference papers for how $\langle J ||e\mathbf{r}|| J' \rangle$ is defined.
  - `qMatrix`: the matrix of $m_F' - m_F = q$ values for each transition where $|q| \leq 1$.  Note that in `qMatrix`, $q$ is shifted up by 1 to be compatible with MATLAB indexing.
  - `dipole`: the actual dipole matrix elements for each transition in Cm
  - `decay`: the decay rates for each transition.
Keep in mind that the coupling matrix and the `qMatrix` can change depending on the magnetic field due to mixing of hyperfine levels.

If you are looking for the matrix of Rabi frequencies, use `op.getLaserFieldMatrix` to return the matrix of $\mathbf{d}\cdot\mathbf{E}/\hbar$.

