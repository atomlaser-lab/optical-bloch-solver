classdef fineStructure < handle
    properties
        species
        numStates
        S
        L
        J
        I
        gI
        gJ
        A1
        A2
        
        E
        H0
        U1int
        U31
        U3int
        BV1
        BV3
    end

    methods
        function obj=fineStructure(species,L,J)
            obj.S=0.5;
            obj.L=L;
            obj.J=J;
            if strcmpi(species,'Rb87'),
                obj.species='Rb87';
                obj.I=3/2;
                obj.gI=-0.0009951414;
            elseif strcmpi(species,'K40'),
                obj.species='K40';
                obj.I=4;
                obj.gI=0.000176490;
            elseif strcmpi(species,'K41'),
                obj.species='K41';
                obj.I=3/2;
                obj.gI=0.000176490;
            end;
            
            if obj.L==0,
                obj.setGroundState;
            elseif obj.L==1 && obj.J==0.5,
                obj.setP1_2;
            elseif obj.L==1 && obj.J==1.5,
                obj.setP3_2;
            end;
            obj.gJ=obj.calcLandeJ(obj.S,obj.L,obj.J);
            obj.numStates=(2*obj.I+1)*(2*obj.J+1);
            obj.solveHyperfine(0);
        end;    %end fineStructure
        
        function obj=makeH0(obj)
            numI=2*obj.I+1;
            numJ=2*obj.J+1;
            nDim=obj.numStates;
            mJ=-obj.J:obj.J;
            mI=-obj.I:obj.I;
            %Uncoupled basis |mJ,mI>
            obj.BV1=[reshape(repmat(mJ,numI,1),nDim,1) repmat(mI(:),numJ,1)];
            
            %I-J coupled basis |F,mF>. Ordered by increasing energy for low magnetic fields
            obj.BV3=zeros(nDim,2);
            if obj.A1>0,
                F=abs(obj.I-obj.J):abs(obj.I+obj.J);
            else
                F=abs(obj.I+obj.J):-1:abs(obj.I-obj.J);
            end;
            mm=1;
            for nn=1:numel(F),
                gF=obj.calcLandeF(obj.I,obj.J,F(nn),obj.gI,obj.gJ);
                if gF>0,
                    mF=-F(nn):F(nn);
                else
                    mF=F(nn):-1:-F(nn);
                end;
                obj.BV3(mm:(mm+2*F(nn)),:)=[F(nn)*ones(2*F(nn)+1,1) mF(:)];
                mm=mm+2*F(nn)+1;
            end;
            
            %Transformation matrix from uncoupled |mJ,mI> basis to the coupled |F,mF> basis
            %Inverse is the hermitian conjugate
            obj.U31=zeros(nDim);
            for a=1:nDim,
                F=obj.BV3(a,1);
                mF=obj.BV3(a,2);
                for b=1:nDim,
                    mJ=obj.BV1(b,1);
                    mI=obj.BV1(b,2);
                    if abs(mI+mJ)>F || (mI+mJ)~=mF,
                        continue;
                    else
                        obj.U31(a,b)=ClebschGordan(obj.I,obj.J,F,mI,mJ,mF);
                    end;
                end;
            end;
            
            %% Bare Hamiltonian calculated in uncoupled basis
            H=zeros(nDim);
            F=abs(obj.I-obj.J):abs(obj.I+obj.J);
            for a=1:nDim,
                mJ1=obj.BV1(a,1);
                mI1=obj.BV1(a,2);
                for b=1:nDim,
                    mJ2=obj.BV1(b,1);
                    mI2=obj.BV1(b,2);
                    for c=1:numel(F),
                        if abs(mI1+mJ1)>F(c) || abs(mI2+mJ2)>F(c) || (mI1+mJ1)~=(mI2+mJ2),
                            continue;
                        end;
                            CG=ClebschGordan(obj.I,obj.J,F(c),mI1,mJ1,mI1+mJ1).*ClebschGordan(obj.I,obj.J,F(c),mI2,mJ2,mI2+mJ2);
                            K=F(c)*(F(c)+1)-obj.I*(obj.I+1)-obj.J*(obj.J+1);
                            if obj.J>=1.0 && obj.I>=1 && obj.L>0,
                                tmp=CG*(obj.A1/2*K+obj.A2*(1.5*K*(K+1)-2*obj.I*obj.J*(obj.I+1)*(obj.J+1))./(4*obj.I*obj.J*(2*obj.I-1)*(2*obj.J-1)));
                            else
                                tmp=CG*obj.A1/2*K;
                            end;
                            H(a,b)=H(a,b)+tmp;
                    end;
                end;
            end;
            obj.H0=H;
        end;    %%end makeH0
        
        function [E,U1int]=solveHyperfine(obj,B)
            mu_B=9.27400968e-24/(6.62606957e-34*1e4)/1e6;
            if numel(obj.H0)==0,
                obj.makeH0;
            end;
            %% Zeeman field
            %E is the eigen-energies as a matrix
            %U1int transforms a vector from the diagonal basis to the uncoupled basis
            if B==0,
                %No field case is diagonal in |F,mF> basis
                E=obj.U31*obj.H0*(obj.U31');
                U1int=obj.U31';
            else
                %With a field it is not diagonal in either basis
                HB=diag(mu_B*B*(obj.gJ*obj.BV1(:,1)+obj.gI*obj.BV1(:,2)));
                H2=obj.H0+HB;
                [U1int,E]=eig(H2);  
            end;
            
            obj.E=E;
            obj.U1int=U1int;
            obj.U3int=obj.U31*obj.U1int;
            if nargout<=1,
                E=diag(E);
            end;
        end;    %end solveHyperfine
        
        
    end

    methods(Access = protected)
        function obj=setGroundState(obj)
            obj.A2=0;
            if strcmpi(obj.species,'Rb87'),
                obj.A1=6834.682610904/(obj.I+obj.S);
            elseif strcmpi(obj.species,'K40'),
                obj.A1=-285.7308;
            elseif strcmpi(obj.species,'K41'),
                obj.A1=254.013872/(obj.I+obj.S);
            end;
        end;
        
        function obj=setP1_2(obj)
            obj.A2=0;
            if strcmpi(obj.species,'Rb87'),
                obj.A1=408.328;
            elseif strcmpi(obj.species,'K40'),
                obj.A1=-34.523;
            elseif strcmpi(obj.species,'K41'),
                obj.A1=15.245;
            end;
        end;    %end setP1_2
        
        function obj=setP3_2(obj)
            
            if strcmpi(obj.species,'Rb87'),
                obj.A1=84.7185;
                obj.A2=12.4965;
            elseif strcmpi(obj.species,'K40'),
                obj.A1=-7.585;
                obj.A2=-3.445;
            elseif strcmpi(obj.species,'K41'),
                obj.A1=3.363;
                obj.A2=3.351;
            end;
        end;    %end setP3_2
    end

    
    methods(Static)
        function gJ=calcLandeJ(S,L,J)
            gJ=(J.*(J+1)-S*(S+1)+L*(L+1))./(2*J.*(J+1))+2*(J.*(J+1)+S*(S+1)-L*(L+1))/(2*J*(J+1));
        end;    %end calcLandeJ
        
        function gF=calcLandeF(I,J,F,gI,gJ)
            gF=gJ*(F.*(F+1)-I*(I+1)+J*(J+1))./(2*F.*(F+1))+gI*(F.*(F+1)+I.*(I+1)-J.*(J+1))./(2*F.*(F+1));
        end;    %end calcLandeF
    end
    
    
    
end