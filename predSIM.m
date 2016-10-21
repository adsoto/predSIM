function dy = predSIM(t,y,s)

% predSIM: computes the derivatives involved in solving a system of ODEs
% which describe the position and orientation of a rigid body. This
% simulates a fish swimming through water. 
%
% The glabal variale 's' is a structure with parameter values

% global s

% Components of drag on body
drag_x  = - 0.5*s.cDrag*s.rho*s.SA*(sqrt(y(2)^2 + y(4)^2))*y(2);
drag_y  = - 0.5*s.cDrag*s.rho*s.SA*(sqrt(y(2)^2 + y(4)^2))*y(4);
drag_theta  = - s.SA*s.rho*y(6)*abs(y(6));

% Lift and Drag on fin
[beta, fin_L, fin_D] = fin_kine(s,y,t);

% Forward thrust
thrust_fwd = fin_L * sin(beta) - fin_D * cos(beta);

% Lateral thrust
thrust_lat = -fin_L * cos(beta) + fin_D * sin(beta);

% Preallocate derivative vector for system of equations
dy = zeros(6,1);

% Equations for x-coordinate 
dy(1) = y(2);
dy(2) = (thrust_fwd*cos(y(5)) + drag_x) ./ s.mass;

% Equations for y-coordinate 
dy(3) = y(4);
dy(4) = (thrust_fwd*sin(y(5)) + drag_y) ./ s.mass;

% Equations for theta (assume COM of body is anterior to its midpoint)
dy(5) = y(6);
dy(6) = (0.7 * s.bodyL * thrust_lat + drag_theta) ./ (s.bodyI);

end


