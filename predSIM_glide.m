function dy = predSIM_glide(t,y,s)

% predSIM_turn: computes the derivatives involved in solving a system of ODEs
% which describe the position and orientation of a rigid body. This model
% simulates a predator fish as it glides. Here we are only modeling the forces
% as the predator glides through the water so there is no turning moment. 
%
% The input variable 's' is a structure with parameter values

% Unpack state variables
Vbod_x  = y(2);
Vbod_y  = y(4);
% theta   = y(5);
theta   = atan2(Vbod_y,Vbod_x);
dTheta  = y(6);
pitch   = y(7);
heave   = y(8);
d_prime = y(9);
h_prime = y(10);

% Compute drag on body
[~,~,drag,drag_theta] = fin_kine(s,theta,dTheta,pitch,heave,...
                                         d_prime,h_prime,Vbod_x,Vbod_y);

% Preallocate derivative vector for system of equations
dy = zeros(10,1);

% Equations for x-coordinate 
dy(1) = Vbod_x;
dy(2) = drag(:,1) ./ s.mass;

% Equations for y-coordinate 
dy(3) = Vbod_y;
dy(4) = drag(:,2) ./ s.mass;

% Equations for theta (assume COM of body is anterior to its midpoint)
dy(5) = dTheta;
dy(6) = (drag_theta) ./ (s.bodyI);

% Rate of change in fin pitch 
dy(7) = d_prime;

% Rate of change in fin heave
dy(8) = h_prime;

% Acceleration of fin pitch
dy(9) = 0;

% Acceleration of heaving
dy(10) = 0;

end


