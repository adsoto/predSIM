function dy = predSIM_turn(t,y,tF,F_norm,F_parl)

% predSIM_turn: computes the derivatives involved in solving a system of ODEs
% which describe the position and orientation of a rigid body. This
% simulates the predator as it initiates a turn. Here we include the full
% dynamic model for the predator 
%
% The glabal variale 's' is a structure with parameter values
%
% The direction of the turn is dictated by the sign of the bearing angle
% and is encoded in the parameter: s.turnDirec

global s

% Interpolate the data set (fT,F_norm) at time t (relative to heading)
    % Thrust perpendicular to long axis of body
    F_norm_val = interp1(tF,F_norm,t);


% Interpolate the data set (fT,F_parl) at time t (relative to heading)
    % Thrust parallel to long axis of body
    F_parl_val = interp1(tF,F_parl,t);
    
% NOTE: Consider including an added mass term to drag equations        

% Components of drag
drag_x  = - 0.5*s.cDrag*s.rho*s.SA*(sqrt(y(2)^2 + y(4)^2))*y(2);
drag_y  = - 0.5*s.cDrag*s.rho*s.SA*(sqrt(y(2)^2 + y(4)^2))*y(4);

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
dy(6) = (0.7 * s.bodyL * F_norm_val - s.SA*s.rho*y(6)*abs(y(6))) ./ (s.bodyI);

end


