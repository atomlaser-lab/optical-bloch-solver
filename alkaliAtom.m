classdef (Abstract) alkaliAtom < handle
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
%         function self = alkaliAtom
%             %ALKALIATOM Creates an instance of the object corresponding to
%             %the species given at the input
%             %
%             %   ATOM = ALKALIATOM(SPECIES) creates an instance using
%             %   SPECIES as the atom type.  Currently supports only all Rb
%             %   and K isotopes
%             %
% 
%         end
        
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

    methods(Static,Abstract)
        freq = freq(transition,initState,finalState,B)
    end

end