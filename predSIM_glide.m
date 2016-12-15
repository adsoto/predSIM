function dy = predSIM_glide(t,y,s)

% predSIM_turn: computes the derivatives involved in solving a system of ODEs
% which describe the position and orientation of a rigid body. This model
% simulates a predator fish as it glides. Here we are only modeling the forces
% as the predator glides through the water so there is no turning moment. 
%
% The input variable 's' is a structure with parameter values

% Unpack state variables
Vbod_x = y(2);
Vbod_y = y(4);
theta  = y(5);
dTheta = y(6);
pitch  = y(7);
heave  = y(8);

% Components of drag on body
drag_x      = - 0.5 * s.cDrag * s.rho * s.SA * abs(Vbod_x) * Vbod_x;
drag_y      = - 0.5 * s.cDrag * s.rho * s.SA * abs(Vbod_y) * Vbod_y;

% Viscous drag that resists rotation
drag_theta  = - s.cDrag_rot * (s.bodyL/2)^3 * s.visc * dTheta^2;
% drag_theta  = - s.SA*s.rho*y(6)*abs(y(6));
% drag_theta  = - 0.5*s.cDrag_rot*y(6)*abs(y(6));

% Preallocate derivative vector for system of equations
dy = zeros(8,1);

% x-velocity of body
dy(1) = Vbod_x;

% x-acceleration of body
dy(2) = (drag_x) ./ s.mass;

% y-velocity of body
dy(3) = y(4);

% y-acceleration of body
dy(4) = (drag_y) ./ s.mass;

% Rate of body rotation (assume COM of body is anterior to its midpoint)
dy(5) = dTheta;

% Body rotational acceleration
dy(6) = (drag_theta) ./ (s.bodyI);

% Rate of change in fin pitch 
dy(7) = 0;

% Rate of change in fin heave
dy(8) = 0;

end


