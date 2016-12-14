function animate_sol(sol,p)



num_time = 500;

bod_clr = .8.*[1 1 1];

play_interval = 0.01;

% Local coordinates for body
theta   = linspace(0,2*pi,300)';
xBodL   = p.bodyL/2 .* cos(theta);
yBodL   = p.bodyW/2 .* sin(theta);
%TODO:Fin



% Time vector
t = linspace(sol.t(1),sol.t(end),num_time)';

% Interpolate all necessary parameters
x       = interp1(sol.t,sol.x,t);
y       = interp1(sol.t,sol.y,t);
theta   = interp1(sol.t,sol.theta,t);
xFin    = interp1(sol.t,sol.finX,t);
yFin    = interp1(sol.t,sol.finY,t);

% Range and limits for axes
rng = max([range(x) range(y)]);
lims = [min([min(x) min(y)])-p.bodyL min([min(x) min(y)])+rng+p.bodyL];

f = figure('DoubleBuffer','on');

% Loop thru time
for i = 1:length(t)
    
    % Current rotation matrix
    R = local_system([0 0],[cos(theta(i)) sin(theta(i))]);
    
    [xBodG,yBodG] = local_to_global([x(i) y(i)],R,xBodL,yBodL);
    
    h = fill(xBodG,yBodG,bod_clr);
    set(h,'EdgeColor','none')
    axis square
    xlim(lims);ylim(lims)
    
    pause(play_interval);
    ttt=3;
end









function R = local_system(origin,rost)
% Defines rotation matrix for a coordinate system

% Check dimensions
if size(origin,1)~=1 || size(origin,2)~=2 || size(rost,1)~=1 || size(rost,2)~=2 
    error('inputs have incorrect dimensions')
end

% Retrieve local x axis to determine coordinate system
xaxis(1,1) = rost(1) - origin(1);
xaxis(1,2) = rost(2) - origin(2);
xaxis(1,3) = 0;

% Normalize to create a unit vector
xaxis = xaxis./norm(xaxis);

%Determine local y axis
%Short hand of cross product of inertial z axis and local x axis
yaxis = [-xaxis(2) xaxis(1) 0];

% Normalize to create a unit vector
yaxis = yaxis./norm(yaxis);

%Determine local z axis
zaxis = cross(xaxis,yaxis);

% Normalize to create a unit vector
zaxis = zaxis./norm(zaxis);

%Create rotation matrix (from inertial axes to local axes)
R = [xaxis(1:2); yaxis(1:2)];
%end

function [xT,yT] = global_to_local(origin,R,x,y)
% Assumes columns vectors for coordinates

pts = [x y];

% Translate
pts(:,1) = pts(:,1) - origin(1);
pts(:,2) = pts(:,2) - origin(2);

% Rotate points
ptsT = [R * pts']';

% Extract columns of points
xT = ptsT(:,1);
yT = ptsT(:,2);
%end

function [xT,yT] = local_to_global(origin,R,x,y)
% Assumes columns vectors for coordinates

pts = [x y];

% Rotate points
ptsT = [inv(R) * pts']';

% Translate global coordinates wrt origin
ptsT(:,1) = ptsT(:,1) + origin(1);
ptsT(:,2) = ptsT(:,2) + origin(2);

% Extract columns of points
xT = ptsT(:,1);
yT = ptsT(:,2);
%end
