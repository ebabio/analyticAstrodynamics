classdef lambertSolver
    % Solve lambert problem
    
    properties
        % Defining parameters
        state0
        state1
        
        transfer
    end
    
    methods
        function obj = lambertSolver(varargin)
            % Class constructor, calling conventions:
            % 1) orbit containing centralBody
            if(nargin ==1)
                obj.transfer = orbit();
                obj.transfer.gravP = varargin{1}.gravP;
                obj.transfer.centralRadius = varargin{1}.centralRadius;
            end
        end
        
        function obj = solveProblem(obj, state0, state1)
            obj.state0 = state0;
            obj.state1 = state1;
            
            obj = solve2dProblem(obj);
            obj = computeAuxiliary(obj);
        end
           
        % given initial and final position and time of flght use lambert's
        % theorem to solve for the conic arch
        function obj = solve2dProblem(obj)
            tof = 24*3600*(obj.state1.julian - obj.state0.julian);
            
            % Determine space triangle characteristics
            r0 = norm(obj.state0.r);
            r1 = norm(obj.state1.r);
            c = norm(obj.state0.rVec - obj.state1.rVec);
            s = 1/2 * (r0 + r1 + c);
            
            % Determine transfer
            % Type
            normal = cross(obj.state0.rVec, obj.state1.rVec);
            if(normal(3) > 0) % if 1 is less than pi rad CCW to 0 then ...
                % normal points upwards
                type = 1;
            else
                type = 2;
            end
            
            % Confirm transfer is elliptic
            if(type == 1)
                parabolicTof = 1/3 * sqrt(2/obj.transfer.gravP) * ...
                    (s^(3/2) - (s-c)^(3/2));
            else
                parabolicTof = 1/3 * sqrt(2/obj.transfer.gravP) * ...
                    (s^(3/2) + (s-c)^(3/2));
            end
            if(parabolicTof < tof)
                ellipseTransfer = 1;
            else
                ellipseTransfer = 0;
            end
            
            
            % Elliptic transfers
            if(ellipseTransfer)
                % Determine if A or B type
                aMin = s/2;
                alpha0 = @(a) (2 * asin( sqrt(s/(2*a)) ) );
                beta0 = @(a) (2 * asin( sqrt((s-c)/(2*a)) ) );
                
                if(type == 1)
                    beta = @(a) (beta0(a));
                else
                    beta = @(a) (-beta0(a));
                end
                tofMin = aMin^(3/2) / sqrt(obj.transfer.gravP) * ...
                    ( (alpha0(aMin) - beta(aMin)) - (sin(alpha0(aMin)) - sin(beta(aMin))) );
                
                if(tofMin > tof)
                    typeLetter = 'a';
                    alpha = @(a) (alpha0(a));
                else
                    typeLetter = 'b';
                    alpha = @(a) (2*pi-alpha0(a));
                end
                
                % Iterate to find a
                lambertsEq = @(a) (a^(3/2) / sqrt(obj.transfer.gravP) * ...
                    ( (alpha(a) - beta(a)) - (sin(alpha(a)) - sin(beta(a))) ) - tof);
                try
                    obj.transfer.a = fzero(lambertsEq, [aMin 3e2*aMin]);
                catch
                    obj.transfer.a = 1e3*s;
                    % disp('elliptic semi major axis too big')
                end
                
                % Find p
                p(1) = 4*obj.transfer.a*(s-r0)*(s-r1)/c^2 * ...
                    sin((alpha(obj.transfer.a)+beta(obj.transfer.a))/2)^2;
                p(2) = 4*obj.transfer.a*(s-r0)*(s-r1)/c^2 * ...
                    sin((alpha(obj.transfer.a)-beta(obj.transfer.a))/2)^2;
                if( (type==1 && typeLetter=='b') || (type==2 && typeLetter=='a') )
                    obj.transfer.p = min(p);
                else
                    obj.transfer.p = max(p);
                end
                
                % Find all other parameters
                obj.transfer.e = sqrt(1 - obj.transfer.p/obj.transfer.a);
                
                % Solve trigonometric ambiguity
                
                obj.state0.eccA = acos(1/obj.transfer.e - r0/(obj.transfer.a * obj.transfer.e) );
                obj.state1.eccA = acos(1/obj.transfer.e - r1/(obj.transfer.a * obj.transfer.e) );
                
                if( norm(...
                        wrapToPi(obj.state1.eccA - obj.state0.eccA) - ...
                        wrapToPi(alpha(obj.transfer.a) - beta(obj.transfer.a))) ...
                        < 1e-2)
                    %do nothing
                elseif ( norm(...
                        wrapToPi(obj.state1.eccA + obj.state0.eccA) - ...
                        wrapToPi(alpha(obj.transfer.a) - beta(obj.transfer.a))) ...
                        < 1e-2)
                    obj.state0.eccA = -obj.state0.eccA;
                elseif  ( norm(...
                        wrapToPi(-obj.state1.eccA - obj.state0.eccA) - ...
                        wrapToPi(alpha(obj.transfer.a) - beta(obj.transfer.a))) ...
                        < 1e-2)
                    obj.state1.eccA = -obj.state1.eccA;
                elseif  ( norm(...
                        wrapToPi(-obj.state1.eccA + obj.state0.eccA) - ...
                        wrapToPi(alpha(obj.transfer.a) - beta(obj.transfer.a))) ...
                        < 1e-2)
                    obj.state0.eccA = -obj.state0.eccA;
                    obj.state1.eccA = -obj.state1.eccA;
                else
                    warning('check was not correct')
                end
                
            else
                % Hyperbolic or Parabolic transfers
                obj.transfer.a = nan;
                obj.transfer.e = nan;
            end
            
        end
        
        function obj = computeAuxiliary(obj)
            perp = cross(obj.state0.rVec, obj.state1.rVec);
            perpHat = perp/norm(perp);
            if(perpHat(3) < 0) % so that it always has CCW rotation
                perpHat = -perpHat;
            end            
            rHat = obj.state0.rVec / obj.state0.r;
            thetaHat = cross(perpHat, rHat);
            
            % 3d parameters
            obj.state0.thetaA = atan2(rHat(3), thetaHat(3));
            obj.transfer.inclination = acos(perpHat(3));
            obj.transfer.lonAscNode = atan2(perpHat(1),-perpHat(2));
            obj.state0.trueA = 2*atan( sqrt((1+obj.transfer.e)/(1-obj.transfer.e)) * ...
                tan(obj.state0.eccA/2) );
            obj.transfer.lonPeriapsis = wrapToPi(obj.state0.thetaA - obj.state0.trueA);
            
            % Epoch mean anomaly
            obj.transfer = computeAuxiliaryParameters(obj.transfer);
            obj.state0.meanA = obj.state0.eccA - obj.transfer.e * sin(obj.state0.eccA);
            obj.transfer.meanA = wrapTo2Pi(obj.state0.meanA - obj.transfer.n*(24*3600)*obj.state0.julian);
            
            % Compute full states
            obj.state0 = obj.transfer.computeState(obj.state0.julian);
            obj.state1 = obj.transfer.computeState(obj.state1.julian);
        end
    end
end

