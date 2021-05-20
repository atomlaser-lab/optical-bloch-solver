classdef alkaliAtom < handle
    properties
        species
        ground
        excited1
        excited2

        D1
        D2

    end

    methods
        function obj=alkaliAtom(species)
%             constants;
            obj.ground=fineStructure(species,0,0.5);
            obj.excited1=fineStructure(species,1,0.5);
            obj.excited2=fineStructure(species,1,1.5);
            if strcmpi(species,'Rb87'),
                obj.species='Rb87';
                obj.D1=opticalTransition(obj.ground,obj.excited1,795e-9,2*pi*6.065e6);
                obj.D2=opticalTransition(obj.ground,obj.excited2,780.241e-9,2*pi*6.065e6);
            elseif strcmpi(species,'K40'),
                obj.species='K40';
                obj.D1=opticalTransition(obj.ground,obj.excited1,770e-9,2*pi*6.035e6);
                obj.D2=opticalTransition(obj.ground,obj.excited2,766.7e-9,2*pi*6.035e6);
            elseif strcmpi(species,'K41'),
                obj.species='K41';
                obj.D1=opticalTransition(obj.ground,obj.excited1,770e-9,2*pi*6.035e6);
                obj.D2=opticalTransition(obj.ground,obj.excited2,766.7e-9,2*pi*6.035e6);
            else
                error('Unsupported species');
            end;
            
            
        end;    %End alkaliAtom
        
        function obj=setMagneticField(obj,B)
            obj.ground.solveHyperfine(B);
            obj.excited1.solveHyperfine(B);
            obj.excited2.solveHyperfine(B);
        end;    %end setMagneticField
        
        
        
    
    end %end methods
    
    methods(Static)
       
    end



end