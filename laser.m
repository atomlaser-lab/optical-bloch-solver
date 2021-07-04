classdef laser < handle
    %LASER Defines a class that represents a laser beam
    properties
        P           %Optical power in W
        w0          %Beam waist in m^2
        
        intensity   %Intensity in W/m^2
        field       %Electric field in V/m
        pol         %Polarization as 3-element vector in either spherical [sigma_-,pi,sigma_+] or linear [x,y,z] polarization
        
        ground      %|F,mF> ground state to use as a reference for the detuning
        excited     %|F',mF'> excited state to use a s a reference for the detuning
        detuning    %The detuning of the laser from the reference transition
    end
    
    properties(Constant)
        sphPolBasis = [1 1i 0;0 0 1;-1 1i 0]; %Converts (x,y,z) components to q=(-1,0,1) (spherical) components
    end
    
    methods
        function self = laser(P,w0)
            %LASER Creates a LASER object
            %
            %   L = LASER(I) Creates a LASER object with intensity I
            %
            %   L = LASER(P,W) Creates a LASER object corresponding to a
            %   Gaussian beam with total power P and beam waist W
            self.detuning = 0;
            if nargin == 1
                self.setIntensity(I);
            elseif nargin == 2
                self.setGaussBeam(P,w0);
            end
        end
        
        function self = setIntensity(self,I)
            %SETINTENSITY Sets the intensity of the laser field
            %
            %   L = L.SETINTENSITY(I) Sets the intensity of the laser field
            %   to I. Automatically calculates the new electric field
            self.intensity = I;
            self.field = self.calcField(self.intensity);
        end
        
        function self = setGaussBeam(self,P,w0)
            %SETGAUSSBEAM Sets intensity based on the assumption of a
            %Gaussian beam
            %
            %   L = L.SETGAUSSBEAM(P,W) sets the intensity and electric
            %   field assuming a Gaussian beam of total power P and beam
            %   waist W
            self.intensity = 2*P/(pi*w0^2);
            self.field = self.calcField(self.intensity);
        end
        
        function self = setPolarization(self,pol,polBasis)
            %SETPOLARIZATION Sets the polarization of the field. Function
            %automatically normalises the input polarization to have unit
            %length
            %
            %   L = L.SETPOLARIZATION(POL) Sets the polarization to the
            %   polarization POL assuming that POL is in the spherical
            %   basis q = [-1,0,1].
            %
            %   L = L.SETPOLARIZATION(__,BASIS) uses BASIS (either 'linear'
            %   or 'spherical') to set the polarization.  If 'linear',
            %   assumes polarization is in [x,y,z] order
            if nargin == 2 || strcmpi(polBasis,'spherical')
                self.pol = self.sphPolBasis'*pol(:);
            elseif strcmpi(polBasis,'linear')
                self.pol = pol(:);
            else
                error('Polarization basis not supported!');
            end
            self.pol = self.pol./sqrt(self.pol'*self.pol);
        end
        
        function self = setStates(self,ground,excited,detuning)
            %SETSTATES Sets the reference transition/states and the
            %detuning from that transition
            %
            %   L = L.SETSTATES(GROUND) Sets the reference ground state as
            %   a 2-element vector [F,mF].
            %
            %   L = L.SETSTATES(__,EXCITED) Sets the reference excited
            %   state as a 2-element vector [F,mF]
            %
            %   L = L.SETSTATUS(__,DETUNING) Sets the detuning from the
            %   reference states/transition to be DETUNING
           if nargin == 2
               self.ground = ground;
           elseif nargin == 3
               self.ground = ground;
               self.excited = excited;
           elseif nargin == 4
               self.ground = ground;
               self.excited = excited;
               self.detuning = detuning;
           end
        end
        
    end
    
    methods(Static)
        function E = calcField(I)
            %CALCFIELD Calculates the electric field for a given intensity
            %
            %   E = CALCFIELD(I) Calculates the electric field given
            %   intensity I
            E = sqrt(0.5*I./(const.eps0*const.c));
        end
    end
    
end