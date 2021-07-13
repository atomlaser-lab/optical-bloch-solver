classdef K39Atom < alkaliAtom
    %K39ATOM Sub-class of alkaliAtom containing necessary properties
    
    methods
        function self = K39Atom
            %K39ATOM Creates an instance of the K39ATOM object with
            %properties corresponding to K-39
            
            self.species = 'K39';
            I = 3/2;
            gI = -0.00014193489;
            self.ground = fineStructure(0,0.5,I,gI,230.8598601e6,0);
            self.excited1 = fineStructure(1,0.5,I,gI,27.775e6,0);
            self.excited2 = fineStructure(1,1.5,I,gI,6.093e6,2.786e6);
            
            self.D1 = opticalTransition(self.ground,self.excited1,770.108385049e-9,2*pi*5.956e6);
            self.D2 = opticalTransition(self.ground,self.excited2,766.700921822e-9,2*pi*6.035e6);
        end
    end
    
    
end