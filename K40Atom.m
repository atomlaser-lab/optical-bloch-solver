classdef K40Atom < alkaliAtom
    %K40ATOM Sub-class of alkaliAtom containing necessary properties
    
    methods
        function self = K40Atom
            %K40ATOM Creates an instance of the K40ATOM object with
            %properties corresponding to K-40
            
            self.species = 'K40';
            I = 4;
            gI = 0.000176490;
            self.ground = fineStructure(0,0.5,I,gI,-285.7308e6,0);
            self.excited1 = fineStructure(1,0.5,I,gI,-34.523e6,0);
            self.excited2 = fineStructure(1,1.5,I,gI,-7.585e6,-3.445e6);
            
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
            
            a = K40Atom;
            if nargin < 4
                f = a.(transition).absoluteFreq(initState,finalState);
            else
                f = a.(transition).absoluteFreq(initState,finalState,B);
            end
        end
    end

end