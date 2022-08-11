classdef Rb85Atom < alkaliAtom
    %RB85ATOM Sub-class of alkaliAtom containing necessary properties
    
    methods
        function self = Rb85Atom
            %RB85ATOM Creates an instance of the RB85ATOM object with
            %properties corresponding to Rb-85
            
            self.species = 'Rb85';
            I = 5/2;
            gI = -0.000293640;
            self.ground = fineStructure(0,0.5,I,gI,1011.910831e6,0);
            self.excited1 = fineStructure(1,0.5,I,gI,120.527e6,0);
            self.excited2 = fineStructure(1,1.5,I,gI,25.002e6,25.790e6);
            
            self.D1 = opticalTransition(self.ground,self.excited1,794.7666414e-9,2*pi*5.75e6);
            self.D2 = opticalTransition(self.ground,self.excited2,780.241368271e-9,2*pi*6.0666e6);
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
            
            a = Rb85Atom;
            if nargin < 4
                f = a.(transition).absoluteFreq(initState,finalState);
            else
                f = a.(transition).absoluteFreq(initState,finalState,B);
            end
        end
    end
end