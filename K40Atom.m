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
    
    
end