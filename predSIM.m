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
% drag_theta  = - s.SA*s.rho*y(6)*abs(y(6));
drag_theta  = - 0.5*s.cDrag_rot*y(6)*abs(y(6));

% Lift and Drag on fin
[s, fin_L, ~] = fin_kine(s,y,t);

% Forward thrust (aligned with long axis of fish)
thrust_fwd = fin_L(1) * cos(y(5)) + fin_L(2) * sin(y(5));
% thrust_fwd = fin_L * sin(beta) - fin_D * cos(beta);

% Lateral thrust (perpendicular to long axis of fish)
thrust_lat = -fin_L(1) * sin(y(5)) + fin_L(2) * cos(y(5));

% x-component of lift on body
lift_x = thrust_fwd*cos(y(5));

% y-component of lift on body
lift_y = thrust_fwd*sin(y(5));

% Torque due to lift force (cross product of fin position and lift vector)
fin_pos = [s.finPos_body(:,1),s.finPos_body(:,2),zeros(length(thrust_fwd),1)];
lift_vec = [fin_L(:,1), fin_L(:,2), zeros(length(thrust_fwd),1)];
fin_pos_cross_lift = cross(fin_pos,lift_vec,2);

torque_lift = fin_pos_cross_lift(:,3); 

% Preallocate derivative vector for system of equations
dy = zeros(6,1);

% Equations for x-coordinate 
dy(1) = y(2);
dy(2) = (lift_x + drag_x) ./ s.mass;

% Equations for y-coordinate 
dy(3) = y(4);
dy(4) = (lift_y + drag_y) ./ s.mass;

% Equations for theta (assume COM of body is anterior to its midpoint)
dy(5) = y(6);
% dy(6) = (0.7 * s.bodyL * thrust_lat + drag_theta) ./ (s.bodyI);
% dy(6) = (s.d_bodyfin * thrust_lat + drag_theta) ./ (s.bodyI);
dy(6) = (torque_lift + drag_theta) ./ (s.bodyI);
end


