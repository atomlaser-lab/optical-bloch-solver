classdef opticalSystem < densityMatrix
    properties
        atom
        transition
        ground
        excited
        
        laser1
        laser2

        B
        detuning

    end
    
    properties(Constant)
        sphPolBasis=[1 1i 0;0 0 1;-1 1i 0]; %Converts (x,y,z) components to q=(-1,0,1) (spherical) components
    end


    methods
        function obj=opticalSystem(species,transition)
            obj=obj@densityMatrix;
            if strcmpi(species,'Rb87'),
                species='Rb87';
            elseif strcmpi(species,'K40'),
                species='K40';
            elseif strcmpi(species,'K41'),
                species='K41';
            else
                error('Unsupported species');
            end;
            
            obj.atom=alkaliAtom(species);
            if strcmpi(transition,'D1'),
                obj.transition=obj.atom.D1;
            elseif strcmpi(transition,'D2'),
                obj.transition=obj.atom.D2;
            end;
            
            obj.ground=obj.transition.ground;
            obj.excited=obj.transition.excited;
            obj.setNumStates(obj.transition.numStates,obj.ground.numStates);
            obj.laser1=laser;
            obj.laser2=laser;
            
        end;    %end opticalSystem
        
        function a=getPopLegend(obj,opt)
            if isnumeric(opt),
                idx=opt;
            elseif all(ischar(opt)),
                if strcmpi(opt,'ground'),
                    idx=1:obj.numGroundStates;
                elseif strcmpi(opt,'excited'),
                    idx=(obj.numGroundStates+1):obj.numStates;
                elseif strcmpi(opt,'all'),
                    idx=1:obj.numStates;
                else
                    error('Population option not supported!');
                end;
            end;
            a=cell(length(idx),1);
            BV3=[obj.ground.BV3;obj.excited.BV3];
            for nn=1:length(idx),a{nn}=sprintf('|%d,%d>',BV3(idx(nn),1),BV3(idx(nn),2));end;
        end;
        
        function plotScatteringRates(obj,figNum)
            if nargin==2,
                figure(figNum);clf;cla;
            else
                clf;cla;
            end;
            R=obj.getScatteringRates;
            set(gca,'NextPlot','ReplaceChildren','LineStyleOrder',{'-','--','.-'})
            plot(obj.t,R,'linewidth',2);
            hold on;
            plot(obj.t,sum(R,1),'linewidth',2,'color','k');
            s=obj.getPopLegend('excited');
            s{end+1}='Total';
            legend(s);
        end;
        
        %% Functions for setting atomic parameters        
        function obj=refresh(obj)
            obj.calcBareH;
            obj.calcCoupling;
        end;
        
        %function for setting magnetic field
        function obj=setMagneticField(obj,B)
            obj.B=B;
            obj.calcBareH;
            obj.calcCoupling;
        end    %end setMagneticField
        
        %function for calculating diagonal elements of Hamiltonian
        function obj=calcBareH(obj,B)
            if nargin==2
                obj.B=B;
            end
            if isempty(obj.laser1.ground) || isempty(obj.laser2.intensity) || obj.laser2.intensity==0
                obj.transition.setMagneticField(obj.B);
                detuning1=obj.transition.calcNewDetuning(obj.laser1);
                obj.bare=2*pi*1e6*blkdiag(obj.ground.E,obj.excited.E-detuning1*eye(obj.excited.numStates));
            else                
                obj.transition.setMagneticField(obj.B); %update magnetic field
                idx=bsxfun(@eq,obj.ground.BV3(:,1),obj.laser2.ground(1)); %find ground states with same F as laser2 (the repump)
                g2=zeros(obj.ground.numStates,1);
                detuning1=obj.transition.calcNewDetuning(obj.laser1);
                detuning2=obj.transition.calcNewDetuning(obj.laser2);
                g2(idx)=detuning1-detuning2;
                g2=diag(g2);
                obj.bare=2*pi*1e6*blkdiag(obj.ground.E-g2,obj.excited.E-detuning1*eye(obj.excited.numStates));
            end;
        end;    %end calcBareH
        
        
       %% Functions for getting coupling parameters
       function obj=calcCoupling(obj)
           obj.transition.makeDipoleMatrix;
           groundU3int=obj.ground.U31*obj.ground.U1int;
           excitedU3int=obj.excited.U31*obj.excited.U1int;
           U3int=blkdiag(groundU3int,excitedU3int);
           obj.coupling=U3int'*obj.getLaserFieldMatrix*U3int;
           obj.decay=obj.transition.getDecayMatrix(U3int');
       end; %end calcCoupling
       
       %Calculates d.E/hbar in |F,mF> basis
       function omega=getLaserFieldMatrix(obj)
            hbar=6.62606957e-34/(2*pi);
            omega=zeros(obj.numStates);    %omega=d.E/hbar
            for g=1:obj.numGroundStates,
                for e=1:obj.numExcitedStates,
                    eShift=e+obj.numGroundStates;
                    q=obj.transition.qMatrix(g,eShift);
                    if q~=0,
                        if isempty(obj.laser1.ground) || isempty(obj.laser2.intensity) || obj.laser2.intensity==0
                            omega(g,eShift)=obj.laser1.field*obj.laser1.pol(q).*obj.transition.dipole(g,eShift)/hbar;
                        else
                            fStart=obj.ground.BV3(g,1);
                            if fStart==obj.laser1.ground(1),
                                omega(g,eShift)=obj.laser1.field*obj.laser1.pol(q).*obj.transition.dipole(g,eShift)/hbar;
                            elseif fStart==obj.laser2.ground(2),
                                omega(g,eShift)=obj.laser2.field*obj.laser2.pol(q).*obj.transition.dipole(g,eShift)/hbar;
                            end;
                        end;
                    end;
               end;
           end;
           omega=omega+omega';
       end; %end getLaserFieldMatrix
       
       function r=getScatteringRates(obj)
           p=obj.getPopulations('excited');
           d=sum(obj.decay,1);
           d=d((obj.ground.numStates+1):end);
           r=bsxfun(@times,p,d(:));
       end;

        
    end     %end methods
    

end