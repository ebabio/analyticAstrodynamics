classdef flyBySolver
    % Solve lambert problem
    
    properties
        % Defining parameters
        localOrbit
        localStateIn % position is replaced with aimPoint
        localStateOut
        
        planetOrbit % of the assisting body        
        planetState
        
        outboundOrbit % of the assisted body after the flyby
        outboundState 
        
        aimPoint
    end
    
    methods
        function obj = flyBySolver(localBody, planetOrbit)
            % Class constructor, calling conventions:
            % 1) localBody, planetOrbit
            if(nargin == 2)
                obj.localOrbit = orbit();
                obj.localOrbit.gravP = localBody.gravP;
                if(isa(localBody,'orbit'))
                    obj.localOrbit.centralRadius = localBody.centralRadius;
                else
                    obj.localOrbit.centralRadius = localBody.radius;
                end
                obj.planetOrbit = planetOrbit;
            end
        end
        
        function obj = solveFlyBy(obj, stateIn, aimPoint)
            % aimPoint is a 3d Vector in local the following frame:
            % 1) V: velocity direction (irrelevant)
            % 2) N: maximum positive slope vector in the normal plane to V
            % 3) C: cross(V,N)
            
            % Define local states
            obj.planetState = obj.planetOrbit.computeState(stateIn.julian);
            
            obj.localStateIn = state(stateIn.rVec-obj.planetState.rVec, ...
                stateIn.vVec-obj.planetState.vVec);
            
            if(obj.localStateIn.r > 1e-2 * stateIn.r)
                warning('relative distance between local bodies may be high');
            end
            
            % Locate aimpoint in 3d
            vHat = obj.localStateIn.vVec/norm(obj.localStateIn.vVec);
            c = cross(vHat, [0; 0; 1]);
            cHat = c / norm(c);
            nHat = cross(cHat, vHat);
            
            CCM = [vHat nHat cHat];
            aimPoint(1) = 0;
            obj.localStateIn.rVec = CCM * aimPoint;
            obj.localStateIn.r = Inf;
            
            % Solve orbital parameters
            obj.localOrbit = orbit(obj.localStateIn, obj.localStateIn.julian, obj.localOrbit);
            
            % Compute escape velocity
            hHat = obj.localOrbit.hVec/obj.localOrbit.h;
            DCM = vrrotvec2mat([hHat; obj.localOrbit.delta]);
            obj.localStateOut.vVec = DCM * obj.localStateIn.vVec;
            
            % Compute outbound orbit
            obj.outboundState = state(obj.planetState.rVec, obj.planetState.vVec + obj.localStateOut.vVec);
            obj.outboundOrbit = orbit(obj.outboundState, obj.outboundState.julian, obj.planetOrbit);
        end
        
    end
end

