function dy = predSIM(t,y,tF,F_norm,F_parl)

% predSIM: computes the derivatives involved in solving a system of ODEs
% which describe the position and orientation of a rigid body. This
% simulates a fish swimming through water. 
%
% The glabal variale 's' is a structure with parameter values

global s

% Interpolate the data set (fT,Fy) at time t (relative to heading)
if t < (s.simDur / 1)
    % Thrust perpendicular to long axis of body
    F_norm_val = interp1(tF,F_norm,t);
else
    F_norm_val = 0;
%     F_norm_val = 1*interp1(tF,F_norm,t);
end

% Interpolate the data set (fT,Fx) at time t (relative to heading)
if t < (s.simDur / 1)
    % Thrust parallel to long axis of body
    F_parl_val = interp1(tF,F_parl,t);
else
    F_parl_val = 0;
end
    
% Rotation matrix (rotates force vector CW by given angle)
% R = [cos(y(5)) -sin(y(5)); sin(y(5)) cos(y(5))]';

% Thrust vector (inertial FOR)
% F = R * [F_norm_val; F_parl_val];

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


