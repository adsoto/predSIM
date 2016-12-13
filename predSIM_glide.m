function dy = predSIM_glide(t,y,s)

% predSIM_turn: computes the derivatives involved in solving a system of ODEs
% which describe the position and orientation of a rigid body. This model
% simulates a predator fish as it glides. Here we are only modeling the forces
% as the predator glides through the water so there is no turning moment. 
%
% The input variable 's' is a structure with parameter values

% Modified tail kinematics during glide (tail doesn't move)
f = s;
f.pitch0 = 0;
f.h0 = 0;

% Get speed of fin
[f,lift,torque,finVel] = fin_kine(f,y,t);

% Update fin position
s.finPos = f.finPos;
s.finPos_body = f.finPos_body;

% Components of drag
drag_x      = - 0.5*s.cDrag*s.rho*s.SA*(sqrt(y(2)^2 + y(4)^2))*y(2);
drag_y      = - 0.5*s.cDrag*s.rho*s.SA*(sqrt(y(2)^2 + y(4)^2))*y(4);
drag_theta  = - s.SA*s.rho*y(6)*abs(y(6));
% drag_theta  = - 0.5*s.cDrag_rot*y(6)*abs(y(6));

% Preallocate derivative vector for system of equations
dy = zeros(8,1);

% Equations for x-coordinate 
dy(1) = y(2);
dy(2) = (drag_x) ./ s.mass;

% Equations for y-coordinate 
dy(3) = y(4);
dy(4) = (drag_y) ./ s.mass;

% Equations for theta (assume COM of body is anterior to its midpoint)
dy(5) = y(6);
dy(6) = (drag_theta) ./ (s.bodyI);

% Fin velocity (so that fin position is returned by the solution)
dy(7) = finVel(:,1);
dy(8) = finVel(:,2);

end


