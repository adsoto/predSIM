function dy = predSIM_glide(t,y)

% predSIM_turn: computes the derivatives involved in solving a system of ODEs
% which describe the position and orientation of a rigid body. This model
% simulates a predator fish as it glides. Here we are only modeling the forces
% as the predator glides through the water so there is no turning moment. 
%
% The glabal variale 's' is a structure with parameter values
%

global s


% Thrust perpendicular to long axis of body
F_norm_val = 0;

% Thrust parallel to long axis of body
F_parl_val = 0;

% Components of drag
drag_x      = - 0.5*s.cDrag*s.rho*s.SA*(sqrt(y(2)^2 + y(4)^2))*y(2);
drag_y      = - 0.5*s.cDrag*s.rho*s.SA*(sqrt(y(2)^2 + y(4)^2))*y(4);
drag_theta  = - s.SA*s.rho*y(6)*abs(y(6));

% Preallocate derivative vector for system of equations
dy = zeros(6,1);

% Equations for x-coordinate (F(2) is along x-direction)
dy(1) = y(2);
dy(2) = (F_parl_val*cos(y(5)) + drag_x) ./ s.mass;

% Equations for y-coordinate (F(1) is along y-direction)
dy(3) = y(4);
dy(4) = (F_parl_val*sin(y(5)) + drag_y) ./ s.mass;

% Equations for theta (assume COM of body is anterior to its midpoint)
dy(5) = y(6);
dy(6) = (0.7 * s.bodyL * F_norm_val + drag_theta) ./ (s.bodyI);

end


