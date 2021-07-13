classdef opticalTransition < handle
    %OPTICALTRANSITION Defines a class that represents an optical
    %transition.
    %
    %   An optical transition consists of a ground and excited state, each
    %   labelled by L, J, and F.  The angular momentum is used to generate
    %   coupling matrices between the ground and excited states
    properties
        wavelength      %The wavelength of the transition in m
        totalDecay      %The total decay rate from the transition in rad s^{-1}
        crossSection    %The maximum absorption cross section
        coupling        %The matrix of couplings between ground and excited states
        dipole          %The dipole matrix, describing electric dipole transitions between ground and excited states
        qMatrix         %Labels the q value (q = mF' - mF) for each pair of ground and excited states
        decay           %Matrix labelling decay rates for each pair of ground and excited states
        transInfo       %Information about transitions between each allowed pair of sub-levels
        
        ground          %The ground state, instance of FineStructure
        excited         %The excited state, instance of FineStructure
        numStates       %The total number of states
    end
    
    methods
        function self = opticalTransition(ground,excited,wavelength,decay)
            %OPTICALTRANSITION Creates an instance of the OPTICALTRANSITION
            %class.  
            %
            %   T = OPTICALTRANSITION(GROUND,EXCITED,WAVELENGTH,DECAY)
            %   creates an instance T associated with the optical
            %   transition between ground state GROUND and excited state
            %   EXCITED with wavelength WAVELENGTH and total decay DECAY.
            %   GROUND and EXCITED are FINESTRUCTURE objects, WAVELENGTH is
            %   in m, and DECAY is in rad/s
            self.ground = ground;
            self.excited = excited;
            self.wavelength = wavelength;
            self.totalDecay = decay;
            self.numStates = self.ground.numStates+self.excited.numStates;
            self.crossSection = 3*self.wavelength^2./(2*pi);
        end
        
        function self = setMagneticField(self,B)
            %SETMAGNETICFIELD Sets the magnetic field and calculates the
            %new ground and excited state energies and eigenvectors
            %
            %   T = T.SETMAGNETICFIELD(B) sets the magnetic field to be B
            %   in G
            self.ground.solveHyperfine(B);
            self.excited.solveHyperfine(B);
        end
        
        function self = makeCoupling(self)
            %MAKECOUPLING Creates the coupling, dipole, q, and dipole
            %matrices
            %
            %   T = T.MAKECOUPLING Calculates the above matrices
            
            %
            % Pre-allocate arrays
            %
            self.coupling = zeros(self.numStates);
            self.qMatrix = zeros(self.numStates);
            I = self.ground.I;  %Shorten the name
            %
            % Loop over ground (g) and excited (e) states
            %
            for g = 1:self.ground.numStates
                for e = 1:self.excited.numStates
                    jStart = self.ground.J;
                    fStart = self.ground.BV3(g,1);
                    mStart = self.ground.BV3(g,2);
                    jEnd = self.excited.J;
                    fEnd = self.excited.BV3(e,1);
                    mEnd = self.excited.BV3(e,2);
                    q = mEnd-mStart;
                    if abs(q) <= 1 && abs(fEnd-fStart) <= 1
                        %
                        % Elements in the coupling and q matrices only
                        % exist when the momentum transferred is less than
                        % 1 unit
                        %
                        % These equations are straight out of Steck's
                        % Rubidium data
                        reducedF = (-1)^(fEnd+jStart+1+I)*sqrt((2*fEnd+1)*(2*jStart+1))*Wigner6jcoeff(jStart,jEnd,1,fEnd,fStart,I);
                        elementF = (-1)^(fEnd-1+mStart)*sqrt(2*fStart+1)*Wigner3j(fEnd,1,fStart,mEnd,-q,-mStart);
                        self.coupling(e+self.ground.numStates,g) = elementF*reducedF;
                        %
                        % The q matrix is offset by 2 so that -1 -> 1, 0 ->
                        % 2, and +1 -> 3, which means they can be used as
                        % indices when accessing other matrices
                        %
                        self.qMatrix(e+self.ground.numStates,g) = q+2;
                    end
                end
            end
            %
            % We've only calculated the upper halves of these matrices, so
            % now make them Hermitian
            %
            self.coupling = self.coupling+self.coupling';
            self.qMatrix = self.qMatrix+self.qMatrix';
            %
            % The dipole matrix is actually <a|er|b>, while the decay
            % matrix is |coupling|^2*TotalDecay
            %
            self.dipole = self.coupling*sqrt((3*pi*const.eps0*const.hbar*(self.wavelength/(2*pi))^3)*(2*self.excited.J+1)/(2*self.ground.J+1)*self.totalDecay);
            self.decay = abs(self.coupling).^2*self.totalDecay;
        end
        
        function decayOut = getDecayMatrix(self,U)
            %GETDECAYMATRIX Returns the decay matrix, possibly rotated into
            %a different basis
            %
            %   DECAY = T.GETDECAYMATRIX() returns the decay matrix in the
            %   |F,mF> basis
            %
            %   DECAY = T.GETDECAYMATRIX(U) returns the decay matrix in the
            %   basis defined the transformation unitary U that transforms
            %   vectors in the |F,mF> basis to the new basis
            %
            if nargin == 1
                decayOut = self.decay;
            else
                decayOut = abs(U*self.coupling*U').^2*self.totalDecay;
            end
        end
        
        function detuning = calcNewDetuning(self,laserIn)
            %CALCNEWDETUNING Calculates a new detuning given a laser
            %detuning
            %
            %   DETUNING = T.CALCNEWDETUNING(LASER) Takes an input LASER
            %   object and returns the new detuning given the current
            %   magnetic field and the states that are being used to
            %   reference the laser detuning
            %
            detuning = laserIn.detuning;
            if length(laserIn.ground) == 2
                %
                % Assuming that the ground state reference in the laser
                % object is |F,mF>, find the state corresponding to that
                % level and calculate the new detuning. Note the minus sign
                % since this is a ground state
                %
                idx = find(all(bsxfun(@eq,self.ground.BV3,laserIn.ground),2));
                if ~isempty(idx)
                    detuning = detuning-self.ground.E(idx,idx);
                end
            end
            if length(laserIn.excited) == 2
                %
                % Assuming that the excited state reference in the laser
                % object is |F,mF>, find the state corresponding to that
                % level and calculate the new detuning. Note the plus sign
                % since this is an excited state
                %
                idx = find(all(bsxfun(@eq,self.excited.BV3,laserIn.excited),2));
                if ~isempty(idx)
                    detuning = detuning+self.excited.E(idx,idx);
                end
            end
        end
        
        function self = getTransitionFrequencies(self)
            %GETTRANSITIONFREQUENCIES Calculates the transition frequencies
            %for every allowed transition and stores it in the internal
            %property transInfo
            self.makeCoupling;
            %
            % Makes a total transformation unitary from the "internal" basis
            % to the |F,mF> basis, and calculates couplings an decays from
            % that
            %
            U3int = blkdiag(self.ground.U3int,self.excited.U3int);
            D = U3int'*self.coupling*U3int;
            G = self.getDecayMatrix(U3int')/(2*pi);
            self.transInfo = transitionInformation;
            transCount = 1;
            %
            % Loop over every ground (g) and excited (e) state
            %
            for g = 1:self.ground.numStates
                for e=1:self.excited.numStates
                    eShift = e+self.ground.numStates;   %The shift of excited state indices since they start after the ground states
                    if G(g,eShift) > 1e-6
                        self.transInfo(transCount) = transitionInformation(-self.ground.E(g,g)+self.excited.E(e,e),...
                            G(g,eShift),D(g,eShift),g,e);
                        transCount = transCount+1;
                    end
                end
            end
        end
        
        function self = plotTransitionFreqs(self,groundF,excitedF,plotOpt)
            %PLOTTRANSITIONFREQS Plots transition frequencies as
            %Lorentzians with accurate widths and coupling strengths
            %
            %   T = T.PLOTTRANSITIONFREQS() Plots all possible transitions
            %   for this J -> J' transition
            %
            %   T = T.PLOTTRANSITIONFREQS(GROUND) Plots all possible
            %   transitions between the ground state labelled by F number
            %   GROUND. If GROUND is empty, plots all transitions
            %
            %   T = T.PLOTTRANSITIONFREQS(__,EXCITED) Plots all possible
            %   transitions to excited state labelled by F number
            %   EXCITED. If EXCITED is empty, plots all transitions
            %
            %   T = T.PLOTTRANSITIONFREQS(__,OPT) Plots transitions and
            %   includes a sum over all transition strengths if OPT = 'all'
            %
            
            %
            % Parse input arguments
            %
            if nargin == 1
                groundF = [];
                excitedF = [];
                plotOpt = '';
            elseif nargin == 2
                excitedF = [];
                plotOpt = '';
            elseif nargin==3
                plotOpt = '';
            end
            %Define Lorenztian profile
            L = @(x,A,x0,width) A./(1+4*((x-x0)/width).^2);
            %Get minimum and maximum frequencies
            fMin = min([self.transInfo.freq]);
            fMax = max([self.transInfo.freq]);
            %
            % Generate plot vectors for the total absorption spectrum
            %
            fTotal = linspace(fMin-20,fMax+20,1e4)';
            yTotal = zeros(size(fTotal));
            fMin = Inf;fMax=-Inf;
            %
            % Loop over all transitions
            %
            for nn = 1:numel(self.transInfo)
                %
                % Shorten names
                %
                g = self.transInfo(nn).g;
                e = self.transInfo(nn).e;
                %
                % Determine if the current ground and excited state should
                % be plotted based on input arguments
                %
                condGround = any(self.ground.BV3(g,1)==groundF) || isempty(groundF);
                condExcited = any(self.excited.BV3(e,1)==excitedF) || isempty(excitedF);
                %
                % Generates plot vectors
                %
                if condExcited && condGround
                    %
                    % Generate frequency vectors
                    %
                    fCenter = min(self.transInfo(nn).freq);
                    f = fCenter+linspace(-50e6,+50e6,1e3)';
                    fMin = min(fMin,min(f));
                    fMax = max(fMax,max(f));
                    %
                    % Generate absorption spectrum and plot
                    %
                    y = L(f,abs(self.transInfo(nn).coupling).^2,self.transInfo(nn).freq,self.transInfo(nn).decay);
                    yTotal = yTotal + L(fTotal,abs(self.transInfo(nn).coupling).^2,self.transInfo(nn).freq,self.transInfo(nn).decay);
                    plot(f,y);
                end
                hold on;
            end
            if strcmpi(plotOpt,'all')
                plot(fTotal,yTotal,'k-','linewidth',1.5);
            end
            hold off;
            xlim([fMin,fMax]);
        end
        
    end
    
end