classdef densityMatrix < handle
    properties
        numStates
        numGroundStates
        numExcitedStates
        density
        densityVec
        initGroundPop
        initExcitedPop
        bare
        coupling
        decay
        lindblad
        
        dt
        nSteps
        t
    end
    
    methods
        function obj=densityMatrix(numStates,numGroundStates)
            if nargin==2,
                obj.setNumStates(numStates,numGroundStates);
            end;
        end;    %end densityMatrix
        
        function obj=setNumStates(obj,numStates,numGroundStates)
            obj.numStates=numStates;
            obj.numGroundStates=numGroundStates;
            obj.numExcitedStates=obj.numStates-obj.numGroundStates;
            obj.density=zeros(obj.numStates);
            obj.bare=zeros(obj.numStates);
            obj.coupling=zeros(obj.numStates);
            obj.decay=zeros(obj.numStates);
            obj.initGroundPop=zeros(obj.numGroundStates,1);
            obj.initExcitedPop=zeros(obj.numExcitedStates,1);
        end;    %end setNumStates
        
        function P=getPopulations(obj,opt)
            if isnumeric(opt),
                pTemp=obj.getPopFromVec;
                P=pTemp(opt(:),:);
            elseif all(ischar(opt)),
                pTemp=obj.getPopFromVec;
                if strcmpi(opt,'ground'),
                    P=pTemp(1:obj.numGroundStates,:);
                elseif strcmpi(opt,'excited'),
                    P=pTemp((obj.numGroundStates+1):end,:);
                elseif strcmpi(opt,'all'),
                    P=pTemp;
                else
                    error('Population option not supported!');
                end;
            end;    
            P=real(P);
        end;    %end getPopulations
        
        function plotPopulations(obj,opt,fnum)
            if nargin==3,
                figure(fnum);clf;cla;
            else
                clf;cla;
            end;
            P=obj.getPopulations(opt);
            set(gca,'NextPlot','ReplaceChildren','LineStyleOrder',{'-','--','.-'})
            plot(obj.t,P,'linewidth',2);
            legend(obj.getPopLegend(opt));
        end;    %end plotPopulations
        
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
            for nn=1:length(idx),a{nn}=sprintf('\rho_{%d%d}',idx(nn),idx(nn));end;
        end;
        
        
        function obj=intConstField(obj,dt,nSteps)
            obj.dt=dt;
            obj.nSteps=round(nSteps);
            obj.t=obj.dt*(0:(obj.nSteps-1));
            M=obj.makeM;
            obj.densityVec=zeros(obj.numStates^2,obj.nSteps);
            tmp=real([obj.initGroundPop(:);obj.initExcitedPop(:)]);
            obj.density=diag(tmp./sum(tmp));
            obj.densityVec(:,1)=obj.density(:);
            if any(isnan(obj.densityVec(:,1))),
                error('NaNs encountered in initial density state');
            end;
            %D_exp=(eye(size(M))-M*obj.dt/2)\(eye(size(M))+M*obj.dt/2);
            D_exp=expm(M*dt);
            for n=2:obj.nSteps,
                obj.densityVec(:,n)=D_exp*obj.densityVec(:,n-1);
            end;
        end;    %end intConstField
        
        function obj=solveSteadyState(obj)
            M=obj.makeM;
            a=diag(ones(obj.numStates,1));
            M(end+1,:)=a(:)';
            b=zeros(size(M,1),1);
            b(end)=1;
            obj.densityVec=M\b;
            obj.density=reshape(obj.densityVec,obj.numStates,obj.numStates);
        end;    %end solveSteadyState
        
        function M=makeM(obj)
            if numel(obj.lindblad)==0,
                obj.makeTotalLindblad;
            end;
            M=obj.flattenCoupling+obj.lindblad;
        end;    %end makeM
        
        function M=flattenCoupling(obj)
            M=zeros(obj.numStates^2);
            for n=1:obj.numStates,
                for m=1:obj.numStates,
                    row=m+(n-1)*obj.numStates;       
                    for j=1:obj.numStates,
                        col=j+(n-1)*obj.numStates;
                        M(row,col)=M(row,col)-1i*(obj.bare(j,m)+obj.coupling(j,m));
                        col=m+(j-1)*obj.numStates;
                        M(row,col)=M(row,col)+1i*(obj.bare(n,j)+obj.coupling(n,j));       
                    end;
                end;
            end;
        end;    %End flattenCoupling
        
        function L=makeTotalLindblad(obj)
            L=zeros(obj.numStates^2);
            for n=1:obj.numGroundStates,
                for m=(obj.numGroundStates+1):obj.numStates,
                    L=L+obj.makeLindblad(n,m);
                end;
            end;
            obj.lindblad=L;
        end;    %end makeTotalLindblad
        
        function L=makeLindblad(obj,g,e)
            L=zeros(obj.numStates^2);
            for n=1:obj.numStates,
                for m=1:obj.numStates,
                    row=m+(n-1)*obj.numStates;
                    if n==g && m==g,
                        col=e+(e-1)*obj.numStates;
                        L(row,col)=L(row,col)+obj.decay(g,e);
                    end;        
                    if n==e,
                        col=m+(e-1)*obj.numStates;
                        L(row,col)=L(row,col)-0.5*obj.decay(g,e);
                    end;        
                    if m==e,
                        col=e+(n-1)*obj.numStates;
                        L(row,col)=L(row,col)-0.5*obj.decay(g,e);
                    end;        
                end;
            end;
        end;    %end makeLindblad
        
        
    end;    %End methods
    
    methods(Access = private)
        function P=getPopFromVec(obj)
            P=obj.densityVec(1:(obj.numStates+1):obj.numStates^2,:);
        end;    %end getPopFromVec
    
    end;
    
    
end