classdef orbit
    %Orbit implements common methods for orbits
    
    properties
        % Defining parameters
        lonAscNode
        inclination
        lonPeriapsis
        a
        e
        meanA % at J2000
        gravP
        
        % Auxiliary parameters
        centralRadius
        p
        hVec
        h
        energy
        n
        delta % only for hyperbolic orbits
    end
    
    methods
        function obj = orbit(varargin)
            % Two possible calling conventions:
            % 1) lonAscNode, inclination, lonPeriapsis, a, e, meanA@J200,
            % centralBody
            % 2) state object, julianDate, orbit w/ same centralBody
            % 3) orbitingBody, centralBody
            if(nargin == 7 || nargin == 2)
                if(nargin == 7)
                    % Initialize defining parameters
                    obj.lonAscNode = varargin{1};
                    obj.inclination = varargin{2};
                    obj.lonPeriapsis = varargin{3};
                    obj.a = varargin{4};
                    obj.e = varargin{5};
                    obj.meanA = wrapToPi(varargin{6});
                    obj.gravP = varargin{7}.gravP;
                    if(isa(varargin{7},'orbit'))
                        obj.centralRadius = varargin{7}.centralRadius;
                    else
                        obj.centralRadius = varargin{7}.radius;
                    end
                elseif(nargin == 2)
                    orbitingBody = varargin{1};
                    centralBody = varargin{2};
                    obj.lonAscNode = orbitingBody.lonAscNode;
                    obj.inclination = orbitingBody.inclination;
                    obj.lonPeriapsis = orbitingBody.lonPeriapsis;
                    obj.a = orbitingBody.a;
                    obj.e = orbitingBody.e;
                    obj.meanA = wrapToPi(orbitingBody.meanLon - orbitingBody.lonAscNode - orbitingBody.lonPeriapsis);
                    
                    obj.gravP = centralBody.gravP;
                    if(isa(varargin{2},'orbit'))
                        obj.centralRadius = varargin{2}.centralRadius;
                    else
                        obj.centralRadius = varargin{2}.radius;
                    end
                end
                
                obj = obj.computeAuxiliaryParameters();
                
            elseif(nargin == 3)
                % handle provided arguments
                state0 = varargin{1};
                state0.julian = varargin{2};
                obj.gravP = varargin{3}.gravP;
                if(isa(varargin{3},'orbit'))
                    obj.centralRadius = varargin{3}.centralRadius;
                else
                    obj.centralRadius = varargin{3}.radius;
                end
                
                % Initialize parameters
                obj.energy = state0.v^2/2 - obj.gravP/state0.r;
                obj.a = -obj.gravP / (2*obj.energy);
                obj.n = sqrt(obj.gravP / obj.a^3);
                obj.hVec = cross(state0.rVec, state0.vVec);
                obj.h = norm(obj.hVec);
                obj.p = obj.h^2/obj.gravP;
                obj.e = sqrt(1 - obj.p/obj.a);
                
                hHat = obj.hVec/obj.h;
                obj.inclination = acos(hHat(3));
                obj.lonAscNode = atan2(hHat(1),-hHat(2));
                
                if(state0.initFlag < 2)
                    state0 = state0.updateAngles(obj);
                end
                
                obj.lonPeriapsis = state0.thetaA - state0.trueA;
                
                if(obj.e < 1)
                    obj.meanA = wrapTo2Pi(state0.meanA - obj.n*(24*3600)*state0.julian);
                elseif(obj.e > 1)
                    obj.delta = pi- 2 * acos(1/obj.e);
                end
            end
        end
        
        function obj = computeAuxiliaryParameters(obj)
            % Initialize auxiliary parameters
            obj.p = obj.a * (1 - obj.e^2);
            obj.energy = -obj.gravP / (2*obj.a);
            obj.n = sqrt(obj.gravP / obj.a^3);
            obj.h = sqrt(obj.gravP * obj.p);
            hHat = [sin(obj.lonAscNode)*sin(obj.inclination);
                -cos(obj.lonAscNode)*sin(obj.inclination);
                cos(obj.inclination)];
            obj.hVec = obj.h * hHat;
            
            if(obj.e > 1)
                obj.delta = pi- 2 * acos(1/obj.e);
            end
        end
        
        function state0 = computeState(obj, julian)
            state0 = state();
            
            % Compute elapsed days since J2000
            state0.julian = julian;
            
            % Compute instantaneuous state
            state0.meanA = wrapToPi(obj.n * (24*3600) * state0.julian + obj.meanA);
            
            keplersEq = @(E) (E - obj.e * sin(E) - state0.meanA);
            state0.eccA = fzero(keplersEq, 0);
            
            state0.trueA = 2 * atan( sqrt( (1+obj.e)/(1-obj.e) ) * tan(state0.eccA/2) );
            state0.r = obj.a * (1 - obj.e * cos(state0.eccA));
            state0.v = sqrt( 2 *(obj.energy + obj.gravP / state0.r) );
            state0.gamma = acos( (obj.h/state0.r) / state0.v );
            
            DCM = rotZXZ(obj.lonAscNode, obj.inclination, obj.lonPeriapsis+state0.trueA);
            state0.rVec = DCM * [state0.r; 0; 0];
            
            if(state0.trueA <0)
                state0.gamma = -1 * state0.gamma;
            end
            DCM = rotZXZ(obj.lonAscNode, obj.inclination, obj.lonPeriapsis+state0.trueA+pi/2-state0.gamma);
            state0.vVec = DCM * [state0.v; 0,;0];
        end
        
        function plot(obj)
            % Compute Orbit radius
            if(obj.e > 1)
                thetaMax = .9* acos(-1/obj.e);
            else
                thetaMax = pi;
            end
            theta = linspace(-thetaMax, thetaMax, 1e3);
            r = obj.p ./ ( 1 + obj.e * cos(theta) );
            
            % Generate in-plane projection
            x = r .* cos(theta);
            y = r .* sin(theta);
            
            % Generate rotations
            DCM = rotZXZ(obj.lonAscNode, obj.inclination, obj.lonPeriapsis);
            
            % Planet dependencies
            [xPlanet, yPlanet, zPlanet] =  sphere(40);
            xPlanet =  obj.centralRadius .*xPlanet;
            yPlanet =  obj.centralRadius .*yPlanet;
            zPlanet =  obj.centralRadius .*zPlanet;
            
            % Plot
            hold on
            
            % Check if figure has been initialized
            persistent initializedFigures;
            figHangle = gcf;
            if(isempty(find(initializedFigures == figHangle.Number, 1)))
                rMax = max(r);
                axisScale = rMax / 2;
                % Plot planet
                surf(xPlanet , yPlanet , zPlanet ,'EdgeColor', 'none')
                colormap(winter)
                % Plot axes
                plot3([0, axisScale, nan, 0, 0, nan, 0, 0], [0, 0, nan, 0, axisScale, nan, 0, 0], [0, 0, nan, 0, 0, nan, 0, axisScale] );
                text([axisScale, 0, 0], [0, axisScale, 0], [0, 0, axisScale], ['X';'Y';'Z']);
                
                initializedFigures = [initializedFigures, figHangle.Number];
            end
            
            % Plot orbit
            rVec = DCM * [x; y; 0 * r];
            plot3(rVec(1,:), rVec(2,:), rVec(3,:));
            axis equal
            grid on
            hold off
        end
        
    end
end

