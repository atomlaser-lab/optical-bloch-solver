classdef Rb87Atom < alkaliAtom
    %RB87ATOM Sub-class of alkaliAtom containing necessary properties
    
    methods
        function self = Rb87Atom
            %RB87ATOM Creates an instance of the RB87ATOM object with
            %properties corresponding to Rb-87
            
            self.species = 'Rb87';
            I = 3/2;
            gI = -0.0009951414;
            self.ground = fineStructure(0,0.5,I,gI,6834.682610904e6/(3/2+1/2),0);
            self.excited1 = fineStructure(1,0.5,I,gI,408.328e6,0);
            self.excited2 = fineStructure(1,1.5,I,gI,84.7185e6,12.4965e6);
            
            self.D1 = opticalTransition(self.ground,self.excited1,794.9788509e-9,2*pi*5.746e6);
            self.D2 = opticalTransition(self.ground,self.excited2,780.241209686e-9,2*pi*6.065e6);
        end
    end
    
    
end