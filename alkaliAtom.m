classdef alkaliAtom < handle
    %ALKALIATOM Defines a class that represents the optical transitions
    %within an alkali metal atom
    
    properties
        species     %The species of atoms as a character vector
        ground      %The ground state as an instance of FineStructure
        excited1    %The P1/2 excited state as an instance of FineStructure
        excited2    %The P3/2 excited state as an instance of FineStructure

        D1          %The D1 transition between ground -> excited1
        D2          %The D2 transition between ground -> excited2
    end

    methods
        function self = alkaliAtom(species)
            %ALKALIATOM Creates an instance of the object corresponding to
            %the species given at the input
            %
            %   ATOM = ALKALIATOM(SPECIES) creates an instance using
            %   SPECIES as the atom type.  Currently supports only 'Rb87',
            %   'K40', and 'K41'
            %
            self.ground = fineStructure(species,0,0.5);
            self.excited1 = fineStructure(species,1,0.5);
            self.excited2 = fineStructure(species,1,1.5);
            if strcmpi(species,'Rb87')
                self.species = 'Rb87';
                self.D1 = opticalTransition(self.ground,self.excited1,795e-9,2*pi*6.065e6);
                self.D2 = opticalTransition(self.ground,self.excited2,780.241e-9,2*pi*6.065e6);
            elseif strcmpi(species,'K40')
                self.species = 'K40';
                self.D1 = opticalTransition(self.ground,self.excited1,770e-9,2*pi*6.035e6);
                self.D2 = opticalTransition(self.ground,self.excited2,766.7e-9,2*pi*6.035e6);
            elseif strcmpi(species,'K41')
                self.species = 'K41';
                self.D1 = opticalTransition(self.ground,self.excited1,770e-9,2*pi*6.035e6);
                self.D2 = opticalTransition(self.ground,self.excited2,766.7e-9,2*pi*6.035e6);
            else
                error('Unsupported species');
            end
        end
        
        function self = setMagneticField(self,B)
            %SETMAGNETICFIELD Sets the magnetic field and solves for the
            %new hyperfine states given that field
            %
            %   ATOM = ATOM.SETMAGNETICFIELD(B) Uses the magnetic field B
            %   in Gauss to solve for the hyperfine states
            self.ground.solveHyperfine(B);
            self.excited1.solveHyperfine(B);
            self.excited2.solveHyperfine(B);
        end

    end

end