classdef K41Atom < alkaliAtom
    %K41ATOM Sub-class of alkaliAtom containing necessary properties
    
    methods
        function self = K41Atom
            %K41ATOM Creates an instance of the K41ATOM object with
            %properties corresponding to K-41
            
            self.species = 'K41';
            I = 3/2;
            gI = -0.00007790600;
            self.ground = fineStructure(0,0.5,I,gI,127.0069352e6,0);
            self.excited1 = fineStructure(1,0.5,I,gI,15.245e6,0);
            self.excited2 = fineStructure(1,1.5,I,gI,3.363e6,3.351e6);
            
            self.D1 = opticalTransition(self.ground,self.excited1,770.108136507e-9,2*pi*5.956e6);
            self.D2 = opticalTransition(self.ground,self.excited2,766.700674872e-9,2*pi*6.035e6);
        end
    end
    
    methods(Static)
        function f = freq(transition,initState,finalState,B)
            %FREQ Computes the absolute frequency of a transition
            %
            %   F = FREQ(TRANSITION,INITSTATE,FINALSTATE) returns
            %   absolute frequency FREQ given [F,mF] states INITSTATE and
            %   FINALSTATE for TRANSITION (either 'D1' or 'D2')
            %
            %   F = FREQ(__,B) calculates the absolute
            %   frequency in magnetic field B
            
            a = K41Atom;
            if nargin < 4
                f = a.(transition).absoluteFreq(initState,finalState);
            else
                f = a.(transition).absoluteFreq(initState,finalState,B);
            end
        end
    end
end