classdef opticalSystem < densityMatrix
    %OPTICALSYSTEM Defines a class that describes the time-evolution of an
    %optical transition in an alkali metal atom.  It's a sub-class of
    %DENSITYMATRIX
    
    properties
        atom        %The atom to use, an instance of ALKALIATOM
        transition  %The transition to use, an instance of OPTICALTRANSITION
        ground      %The ground state to consider, an instance of FINESTRUCTURE
        excited     %The excited state to consider, an instance of FINESTRUCTURE
        
        laser1      %The primary laser field, an instance of LASER
        laser2      %The secondary laser field, an instance of LASER

        B           %The magnetic field to use in G
    end

    methods
        function self = opticalSystem(species,transition)
            %OPTICALSYSTEM Creates an instance of the class given the
            %species and transition to use
            %
            %   OP = OPTICALSYSTEM(SPECIES,TRANSITION) creates an instance
            %   for species SPECIES ('Rb87', 'K40', or 'K41') and
            %   transition 'D1' or 'D2'
            
            self = self@densityMatrix;  %This syntax is necessary for sub-classes
            self.atom = alkaliAtom(species);
            if strcmpi(transition,'D1')
                self.transition = self.atom.D1;
            elseif strcmpi(transition,'D2')
                self.transition = self.atom.D2;
            else
                error('Transition ''%s'' not supported!',transition);
            end
            %
            % Set the ground and excited states
            %
            self.ground = self.transition.ground;
            self.excited = self.transition.excited;
            self.setNumStates(self.transition.numStates,self.ground.numStates);
            %
            % Create blank instances of the primary and secondary lasers
            %
            self.laser1 = laser;
            self.laser2 = laser;
        end
        
        function str = getPopLegend(self,opt)
            %GETPOPLEGEND Returns a cell array that can be used as a legend
            %for a plot of populations
            %
            %   STR = D.GETPOPLEGEND() Returns a cell array of strings that
            %   is a legend for all populations
            %
            %   STR = D.GETPOPLEGEND(OPT) Returns a cell array of strings
            %   that is a legend for only those populations specified by
            %   OPT.  Valid values for OPT are the same as for
            %   GETPOPULATIONS
            %
            if isnumeric(opt)
                idx = opt;
            elseif all(ischar(opt))
                if strcmpi(opt,'ground')
                    idx = 1:self.numGroundStates;
                elseif strcmpi(opt,'excited')
                    idx = (self.numGroundStates+1):self.numStates;
                elseif strcmpi(opt,'all')
                    idx = 1:self.numStates;
                else
                    error('Population option not supported!');
                end
            end
            str = cell(length(idx),1);
            %
            % This creates labels in the |F,mF> basis
            %
            BV3 = [self.ground.BV3;self.excited.BV3];
            for nn=1:length(idx)
                str{nn} = sprintf('|%d,%d>',BV3(idx(nn),1),BV3(idx(nn),2));
            end
        end
        
        function plotScatteringRates(self)
            %PLOTSCATTERINGRATES Plots scattering rates as a function of
            %time
            %
            %   OP.PLOTSCATTERINGRATES() Plots the scattering rates as a
            %   function of time in the current axes
            R = self.getScatteringRates;
            set(gca,'NextPlot','ReplaceChildren','LineStyleOrder',{'-','--','.-'})
            plot(self.t,R,'linewidth',2);
            hold on;
            plot(self.t,sum(R,1),'linewidth',2,'color','k');
            hold off;
            s = self.getPopLegend('excited');
            s{end+1} = 'Total';
            legend(s);
        end
        
        %% Functions for setting atomic parameters        
        function self = refresh(self)
            %REFRESH Re-calculates the bare Hamiltonian and the coupling
            %terms
            %
            %   OP = OP.REFRESH() Re-calculates the bare and coupling
            %   Hamiltonians
            self.calcBareH;
            self.calcCoupling;
        end
        
        function self = setMagneticField(self,B)
            %SETMAGNETICFIELD Sets the magnetic field
            %
            %   OP = OP.SETMAGNETICFIELD(B) sets the magnetic field to B in
            %   Gauss and recalculates the bare and coupling Hamiltonians
            self.B = B;
            self.refresh;
        end

        function self = calcBareH(self,B)
            %CALCBAREH Calculates the "bare" Hamiltonian which is the
            %diagonal part of the Hamiltonian
            %
            %   OP = OP.CALCBAREH() Calculates the bare Hamiltonian using
            %   the internal magnetic field defined by property OP.B
            %
            %   OP = OP.CALCBAREH(B) Uses the new magnetic field B in Gauss
            
            if nargin==2
                self.B=B;
            end
            %
            % Diagonalize hyperfine Hamiltonians for the ground and excited
            % states in the optical transition under consideration
            %
            self.transition.setMagneticField(self.B);
            %
            % This shifts the bare Hamiltonian based on the laser detunings
            % and what reference states they are meant to track
            %
            if isempty(self.laser1.ground) || isempty(self.laser2.intensity) || self.laser2.intensity == 0
                %
                % Shift "bare" Hamiltonian based on the detuning of the
                % primary laser.  This only applies if laser2 is not set
                % and no ground state is set for the primary laser
                % (laser1). This also applies if the "ground" state of the
                % primary laser is empty, indicating that it is meant to
                % apply to both states
                %
                detuning1 = self.transition.calcNewDetuning(self.laser1);
                self.bare = 2*pi*1e6*blkdiag(self.ground.E,self.excited.E - detuning1*eye(self.excited.numStates));
            else                
                %
                % Shift the bare Hamiltonian based on the detunings of the
                % primary and secondary laser
                %
                idx = bsxfun(@eq,self.ground.BV3(:,1),self.laser2.ground(1)); %find ground states with same F as laser2 (the repump)
                g2 = zeros(self.ground.numStates,1);
                detuning1 = self.transition.calcNewDetuning(self.laser1);
                detuning2 = self.transition.calcNewDetuning(self.laser2);
                g2(idx) = detuning1 - detuning2;    %This is the two-photon detuning
                g2 = diag(g2);
                %
                % First part of matrix is the energies shifted by the
                % two-photon detuning, the second part is the excited state
                % energies shifted by the one-photon detuning
                %
                self.bare = 2*pi*1e6*blkdiag(self.ground.E - g2,self.excited.E - detuning1*eye(self.excited.numStates));
            end
        end
        
        
        %% Functions for getting coupling parameters
        function self = calcCoupling(self)
            %CALCCOUPLING Calculates the coupling matrices based on the
            %current transition and electric fields from the lasers. These
            %matrices are in the "internal" basis
            %
            %    D = D.CALCCOUPLING() calculate the coupling matrices and
            %    stores them as internal properties
            %
            self.transition.makeCoupling;
            groundU3int = self.ground.U31*self.ground.U1int;
            excitedU3int = self.excited.U31*self.excited.U1int;
            U3int = blkdiag(groundU3int,excitedU3int);
            self.coupling = U3int'*self.getLaserFieldMatrix*U3int;
            self.decay = self.transition.getDecayMatrix(U3int');
        end
       
        function omega = getLaserFieldMatrix(self)
            %GETLASERFIELDMATRIX Returns the matrix of d\cdot E/hbar in the
            %|F,mF> basis
            %
            %    OMEGA = D.GETLASERFIELDMATRIX() Returns the marix of d\cdot
            %    E/hbar in the |F,mF> basis
            %
            omega = zeros(self.numStates);    %omega=d.E/hbar
            %
            % Loop over all ground and excited states
            %
            for g = 1:self.numGroundStates
                for e = 1:self.numExcitedStates
                    eShift = e + self.numGroundStates;
                    q = self.transition.qMatrix(g,eShift);
                    if q ~= 0
                        if isempty(self.laser1.ground) || isempty(self.laser2.intensity) || (self.laser2.intensity == 0)
                            %
                            % This applies only if the primary laser is
                            % meant to couple both ground states or if the
                            % secondary laser is not set or has empty
                            % intensity
                            %
                            omega(g,eShift) = self.laser1.field*self.laser1.pol(q).*self.transition.dipole(g,eShift)/const.hbar;
                        else
                            %
                            % This applies if we are considering a
                            % two-laser system where the primary and
                            % secondary lasers couple the two ground F
                            % states independently
                            %
                            fStart = self.ground.BV3(g,1);
                            if fStart == self.laser1.ground(1)
                                omega(g,eShift) = self.laser1.field*self.laser1.pol(q).*self.transition.dipole(g,eShift)/hbar;
                            elseif fStart == self.laser2.ground(1)
                                omega(g,eShift) = self.laser2.field*self.laser2.pol(q).*self.transition.dipole(g,eShift)/hbar;
                            end
                        end
                    end
                end
            end
           omega = omega + omega';
        end
       
        function R = getScatteringRates(self)
            %GETSCATTERINGRATES Returns the scattering rates assuming that
            %they are R = \rho_{nn}*(decay rate for state n to all ground
            %states)
            %
            %   R = OP.GETSCATTERINGRATES() Returns the scattering rates R
            %   as a function of time
            %
            p = self.getPopulations('excited');
            d = sum(self.decay,1);
            d = d((self.ground.numStates+1):end);
            R = bsxfun(@times,p,d(:));
        end
    end
    

end