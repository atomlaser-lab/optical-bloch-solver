classdef laser < handle
    properties
        P
        w0
        
        intensity
        field
        pol
        
        %|F,mF> specifying state that the laser is detuned from
        ground  
        excited
        detuning
    end
    
    properties(Constant)
        sphPolBasis=[1 1i 0;0 0 1;-1 1i 0]; %Converts (x,y,z) components to q=(-1,0,1) (spherical) components
    end
    
    methods
        function obj=laser(P,w0)
            obj.detuning=0;
            if nargin~=0,
                obj.setGaussBeam(P,w0);
            end;
        end;
        
        function obj=setGaussBeam(obj,P,w0)
            if nargin==2,
                obj.P=P;
            elseif nargin==3,
                obj.P=P;
                obj.w0=w0;
            end;
            obj.field=sqrt(obj.P./(pi*const.eps0*const.c*obj.w0^2));
            obj.intensity=2*const.eps0*const.c*obj.field.^2;
        end;    %end setGaussBeam
        
        function obj=setPolarization(obj,pol,polBasis)
            if nargin==2,
                obj.pol=pol(:);
            elseif strcmpi(polBasis,'linear'),
                obj.pol=obj.sphPolBasis*pol(:);
            elseif strcmpi(polBasis,'spherical'),
                obj.pol=pol(:);
            else
                error('Polarization basis not supported!');
            end;
            obj.pol=obj.pol./sqrt(obj.pol'*obj.pol);
        end;    %end setPolarization
        
        function obj=setStates(obj,ground,excited,detuning)
           if nargin==2,
               obj.ground=ground;
           elseif nargin==3,
               obj.ground=ground;
               obj.excited=excited;
           elseif nargin == 4
               obj.ground = ground;
               obj.excited = excited;
               obj.detuning = detuning;
           end;
        end;    %end setDetuning
        
    end
    
end