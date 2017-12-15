classdef state
    % State class dual with state class:
    % while orbit class defines orbit independent of time state class
    % defines body state in an instant in terms of cartesian coordinates
    % and angular coordinates in the plane of motion
    
    properties
        % Defining parameters
        rVec = nan;
        vVec = nan;
        julian = nan; % julian day
        
        % Auxiliary parameters
        r = nan;
        v = nan;
        trueA = nan;
        eccA = nan;
        meanA = nan;
        thetaA = nan;
        gamma = nan;
        
        initFlag = 0;
        % integer: 0 empty, 1 only position/vel, 2
    end
    
    methods
        function obj = state(varargin)
            % Two possible calling conventions:
            % 1) ~ : empty container
            % 2) rVec, vVec
            if(nargin == 0)
            elseif( nargin == 2)
                obj.rVec = varargin{1};
                obj.vVec = varargin{2};
                
                obj.r = norm(obj.rVec);
                obj.v = norm(obj.vVec);
                
                obj.initFlag = 1;
            end
        end
        
        % use cartesian components and a given orbit to find all angles of
        % interest
        function obj = updateAngles(obj, orbit)
            hHat = orbit.hVec/orbit.h;
            rHat = obj.rVec / obj.r;
            thetaHat = cross(hHat, rHat);
            
            obj.thetaA = atan2(rHat(3), thetaHat(3));
            obj.trueA = acos( (orbit.p / obj.r - 1)/orbit.e );
            if(orbit.e < 1)
                obj.eccA = acos( (1 - obj.r/orbit.a)/orbit.e );
                if(obj.rVec' * obj.vVec < 0)
                    obj.trueA = -obj.trueA;
                    obj.eccA = -obj.eccA;
                end
                obj.meanA = obj.eccA - orbit.e * sin(obj.eccA);
            end
            
            obj.initFlag = 2;
        end
        
        % visualization of the state in  given instant
        function plot(obj, quiverScale)
            hold on
            scatter3(obj.rVec(1), obj.rVec(2), obj.rVec(3),'x');
            quiver3(obj.rVec(1), obj.rVec(2), obj.rVec(3), ...
                obj.vVec(1), obj.vVec(2), obj.vVec(3), quiverScale);
            hold off
        end
    end
end

