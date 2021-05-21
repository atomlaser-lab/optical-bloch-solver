classdef transitionInformation
    %TRANSITIONINFORMATION Defines a class that contains information about
    %transitions between sublevels of a J -> J' transition
    
    properties
        freq        %The transition frequency in Mhz?
        decay       %The decay rate in rad/s
        coupling    %The coupling strength (dimensionless)
        g           %The ground state index
        e           %The excited state index
    end
    
    methods
        function self = transitionInformation(freq,decay,coupling,g,e)
            %TRANSITIONINFORMATION Creates an instance of the class
            %
            %   INFO = TRANSITIONINFORMATION Creates a blank instance
            %
            %   INFO = TRANSITIONINFORMATION(FREQ,DECAY,COUPLING,G,E)
            %   creates an instance with frequency FREQ, decay DECAY,
            %   coupling strength COUPLING, ground state index G, and
            %   excited state index E
            if nargin > 0
                self.freq = freq;
                self.decay = decay;
                self.coupling = coupling;
                self.g = g;
                self.e = e;
            end
        end
    end
    
end