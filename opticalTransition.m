classdef opticalTransition < handle
    properties
        wavelength
        totalDecay
        crossSection
        coupling
        dipole
        qMatrix
        decay
        freqInfo
        
        ground
        excited
        numStates
    end
    
    methods
        function obj=opticalTransition(ground,excited,wavelength,decay)
            obj.ground=ground;
            obj.excited=excited;
            obj.wavelength=wavelength;
            obj.totalDecay=decay;
            obj.numStates=obj.ground.numStates+obj.excited.numStates;
            obj.crossSection=3*obj.wavelength^2./(2*pi);
        end;    %end opticalTransition
        
        function obj=setMagneticField(obj,B)
            obj.ground.solveHyperfine(B);
            obj.excited.solveHyperfine(B);
        end;    %end setMagneticField
        
        function obj=makeDipoleMatrix(obj)
            obj.coupling=zeros(obj.numStates);
            obj.qMatrix=zeros(obj.numStates);
            I=obj.ground.I;
            for g=1:obj.ground.numStates,
                for e=1:obj.excited.numStates,
                    jStart=obj.ground.J;
                    fStart=obj.ground.BV3(g,1);
                    mStart=obj.ground.BV3(g,2);
                    jEnd=obj.excited.J;
                    fEnd=obj.excited.BV3(e,1);
                    mEnd=obj.excited.BV3(e,2);
                    q=mEnd-mStart;
                    if abs(q)<=1 && abs(fEnd-fStart)<=1,
                        reducedF=(-1)^(fEnd+jStart+1+I)*sqrt((2*fEnd+1)*(2*jStart+1))*Wigner6jcoeff(jStart,jEnd,1,fEnd,fStart,I);
                        elementF=(-1)^(fEnd-1+mStart)*sqrt(2*fStart+1)*Wigner3j(fEnd,1,fStart,mEnd,-q,-mStart);
                        obj.coupling(e+obj.ground.numStates,g)=elementF*reducedF;
                        obj.qMatrix(e+obj.ground.numStates,g)=q+2;
                    end;
                end;
            end;
            obj.coupling=obj.coupling+obj.coupling';
            obj.qMatrix=obj.qMatrix+obj.qMatrix';
            obj.dipole=obj.coupling*sqrt((3*pi*const.eps0*const.hbar*(obj.wavelength/(2*pi))^3)*(2*obj.excited.J+1)/(2*obj.ground.J+1)*obj.totalDecay);
            obj.decay=abs(obj.coupling).^2*obj.totalDecay;
        end;    %end makeDipoleMatrix
        
        function decayOut=getDecayMatrix(obj,U)
            if nargin==1,
                decayOut=obj.decay;
            else
                decayOut=abs(U*obj.coupling*U').^2*obj.totalDecay;
            end;
        end;    %end makeDecayMatrix
        
        function detuning=calcNewDetuning(obj,laserIn)
            detuning=laserIn.detuning;
            if length(laserIn.ground)==2,
                idx=find(all(bsxfun(@eq,obj.ground.BV3,laserIn.ground),2));
                if ~isempty(idx),
                    detuning=detuning-obj.ground.E(idx,idx);
                end;
            end;
            if length(laserIn.excited)==2,
                idx=find(all(bsxfun(@eq,obj.excited.BV3,laserIn.excited),2));
                if ~isempty(idx),
                    detuning=detuning+obj.excited.E(idx,idx);
                end;
            end;
        end;    %end calcNewDetuning
        
        function obj=getTransitionFrequencies(obj)
            obj.makeDipoleMatrix;
            U3int=blkdiag(obj.ground.U3int,obj.excited.U3int);
            D=U3int'*obj.coupling*U3int;
            G=obj.getDecayMatrix(U3int')/(2*pi*1e6);
            obj.freqInfo=struct('F',0,'idx',[0,0],'decay',0,'coupling',0);
            freqCount=1;
            for g=1:obj.ground.numStates,
                for e=1:obj.excited.numStates,
                    eShift=e+obj.ground.numStates;
                    if G(g,eShift)>1e-6,
                        obj.freqInfo(freqCount).F=-obj.ground.E(g,g)+obj.excited.E(e,e);
                        obj.freqInfo(freqCount).idx=[g,e];
                        obj.freqInfo(freqCount).decay=G(g,eShift);
                        obj.freqInfo(freqCount).coupling=D(g,eShift);
                        freqCount=freqCount+1;
                    end;
                end;
            end;
        end;    %end getTransitionFrequencies
        
        function obj=plotTransitionFreqs(obj,groundF,excitedF,plotOpt,fNum)
            if nargin==1,
                groundF=[];
                excitedF=[];
                plotOpt='';
                figure;
            elseif nargin==2,
                excitedF=[];
                plotOpt='';
                figure;
            elseif nargin==3,
                plotOpt='';
                figure;
            elseif nargin==4,
                figure;
            elseif nargin==5,
                figure(fNum);clf;
            end;
            L=@(x,A,x0,width) A./(1+4*((x-x0)/width).^2);
            fMin=min([obj.freqInfo.F]);
            fMax=max([obj.freqInfo.F]);
            fTotal=linspace(fMin-20,fMax+20,1e4)';
            yTotal=zeros(size(fTotal));
            fMin=Inf;fMax=-Inf;
            for nn=1:numel(obj.freqInfo),
                g=obj.freqInfo(nn).idx(1);
                e=obj.freqInfo(nn).idx(2);
                condGround=any(obj.ground.BV3(g,1)==groundF) || isempty(groundF);
                condExcited=any(obj.excited.BV3(e,1)==excitedF) || isempty(excitedF);
                if condExcited && condGround,
                    fCenter=min(obj.freqInfo(nn).F);
                    f=fCenter+linspace(-50,+50,1e3)';
                    fMin=min(fMin,min(f));
                    fMax=max(fMax,max(f));
                    y=L(f,abs(obj.freqInfo(nn).coupling).^2,obj.freqInfo(nn).F,obj.freqInfo(nn).decay);
                    yTotal=yTotal+L(fTotal,abs(obj.freqInfo(nn).coupling).^2,obj.freqInfo(nn).F,obj.freqInfo(nn).decay);;
                    plot(f,y);
                end;
                hold on;
            end;
            if strcmpi(plotOpt,'all'),
                plot(fTotal,yTotal,'k-','linewidth',1.5);
            end;
            hold off;
            xlim([fMin,fMax]);
        end;    %end plotTransitionFreqs
        
    end
    
end