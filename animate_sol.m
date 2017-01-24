function animate_sol(sol,p,save_movie)


if nargin < 3
    save_movie = 0;
    
end

% Number of time points (frames)
num_time = 300;

% Rate to play exported video
playrate = 30;

%v.bod_clr = .8.*[1 1 1];

play_interval = 0.001;

% Local coordinates for body
theta   = linspace(0,2*pi,300)';
v.xBodL   = p.bodyL/2 .* cos(theta);
v.yBodL   = p.bodyW/2 .* sin(theta);
%TODO:Fin


% Time vector
v.t = linspace(sol.t(1),sol.t(end),num_time)';

% Interpolate all necessary parameters
v.x       = interp1(sol.t,sol.x,v.t);
v.y       = interp1(sol.t,sol.y,v.t);
v.theta   = interp1(sol.t,sol.theta,v.t);
v.h       = interp1(sol.t,sol.heave,v.t);
v.p       = interp1(sol.t,sol.pitch,v.t);
v.force   = sqrt(interp1(sol.t,sol.lift(:,1),v.t).^2 + ...
                 interp1(sol.t,sol.lift(:,2),v.t).^2 );
v.drag   = sqrt(interp1(sol.t,sol.drag(:,1),v.t).^2 + ...
                interp1(sol.t,sol.drag(:,2),v.t).^2 );
         
%xFin    = interp1(sol.t,sol.finX,t);
%yFin    = interp1(sol.t,sol.finY,t);

% Range and limits for axes
%rng = max([range(x) range(y)]);
%lims = [min([min(x) min(y)])-p.bodyL min([min(x) min(y)])+rng+p.bodyL];

f = figure('DoubleBuffer','on','Visible','off','Color','w');

% Plot series of frames
plot_body(v,1,p,5); hold on
plot_body(v,round(length(v.t)/4),p,5);hold on
plot_body(v,round(length(v.t)/2),p,5);hold on
plot_body(v,round(length(v.t)*.75),p,5);hold on
plot_body(v,length(v.t),p,5);hold on
axis equal

% Plot prey position
plot(p.preyX,p.preyY,'r.','MarkerSize',18)
hold off


% Save axis limits
xL = xlim; yL = ylim;

% Make some space
xL(1) = xL(1) - range(xL)/10;
xL(2) = xL(2) + range(xL)/10;
yL(1) = yL(1) - range(xL)/10;
yL(2) = yL(2) + range(xL)/10;

set(f,'Visible','on')

% Prompt for where to save movie            
if save_movie
    [FileName,PathName] = uiputfile;
    if isempty(FileName)
        return
    end
    outputVideo = VideoWriter(fullfile(PathName,[FileName '.m4v']));
    outputVideo.FrameRate = playrate;
    open(outputVideo)
end
   
% Get size of tail and prey
tailThick = ceil(0.1*range(xlim)*100);
preySize = round(3*range(xlim)*100);

% Loop thru time
for i = 1:length(v.t)
    
    % Make figure visible
    set(f,'Visible','on')

    % Plot prey position
    plot(sol.preyPos(1),sol.preyPos(2),'r.','MarkerSize',...
        preySize)
    hold on
    
    % Plot current frame
    plot_body(v,i,p,tailThick)
    hold off
    
    % Set axes
    xlim(xL);ylim(yL)
    
    if save_movie
        % Grab frame
        img = getframe(f);
        
        writeVideo(outputVideo,img)
    end
    
    % Wait before next frame
    pause(play_interval);
end

if save_movie
    close(outputVideo)
    close(f)
end




function plot_body(v,frame,p,tailThick)

% Distance from COM to posterior margin of trunk
tr_len = 0.7*p.bodyL;

% Length of peduncle
pd_len = p.pedL;

% Distance of COP along chord length
tl_len = p.finL;

% Current rotation matrix
R = local_system([0 0],[cos(v.theta(frame)) sin(v.theta(frame))]);
    
x     = v.x(frame);
y     = v.y(frame);
theta = v.theta(frame);
heave = v.h(frame);
pitch = v.p(frame);
force = v.force(frame);
drag  = v.drag(frame);

% Colormap for tail
cmap = colormap('parula');

% Color position for tail force
cpos = force./max([max(v.force) max(v.drag)]) * size(cmap,1);

tailclr = [interp1(1:size(cmap,1),cmap(:,1),cpos) ...
           interp1(1:size(cmap,1),cmap(:,2),cpos) ...
           interp1(1:size(cmap,1),cmap(:,3),cpos)];
if max(isnan(tailclr))
    tailclr = cmap(1,:);
end

% Color position for body
cpos = drag./max([max(v.force) max(v.drag)]) * size(cmap,1);

bodclr = [interp1(1:size(cmap,1),cmap(:,1),cpos) ...
           interp1(1:size(cmap,1),cmap(:,2),cpos) ...
           interp1(1:size(cmap,1),cmap(:,3),cpos)];
if max(isnan(bodclr))
    bodclr = cmap(1,:);
end
       
% Transform body in global FOR
[xBodG,yBodG] = local_to_global([v.x(frame) v.y(frame)],R,v.xBodL,v.yBodL);
    

% Coordinates of trailing edge of trunk
tr_pos(:,1) = x - tr_len.*cos(theta);
tr_pos(:,2) = y - tr_len.*sin(theta);

% Coordinates of peduncle 
pd_pos(:,1) = [tr_pos(:,1) - pd_len.*cos(theta+heave)];
pd_pos(:,2) = [tr_pos(:,2) - pd_len.*sin(theta+heave)];

% Coordinates of qtr-chord point           
fin_pos(:,1) = pd_pos(:,1) -  tl_len.*cos(theta+heave+pitch);
fin_pos(:,2) = pd_pos(:,2) -  tl_len.*sin(theta+heave+pitch);

% Coordinate of leading edge of fin




% Current position of fin quarter-chord point
% finPos(1,1) = x - 0.7*p.bodyL*cos(theta) ...
%                  - p.pedL*cos(theta+heave) ...
%                  - 0.25*p.finL*cos(theta+heave+pitch);
% finPos(1,2) = y - 0.7*p.bodyL*sin(theta) ...
%                  - p.pedL*sin(theta+heave) ...
%                  - 0.25*p.finL*sin(theta+heave+pitch);

% % Fin leading edge coordinates             
% leadEdge(1,1) = x - 0.7*p.bodyL*cos(theta) ...
%                  - p.pedL*cos(theta+heave); 
% leadEdge(1,2) = y - 0.7*p.bodyL*sin(theta) ...
%                   - p.pedL*sin(theta+heave); 
% 
% % Fin trailing edge coordinates               
% trailEdge(1,1) = x - 0.7*p.bodyL*cos(theta) ...
%                    - p.pedL*cos(theta+heave) ...
%                    - p.finL*cos(theta+heave+pitch);
% trailEdge(1,2) = y - 0.7*p.bodyL*sin(theta) ...
%                    - p.pedL*sin(theta+heave) ...
%                    - p.finL*sin(theta+heave+pitch);              

% Draw body    
h = fill(xBodG,yBodG,bodclr);
set(h,'EdgeColor','none')
hold on


h = plot([pd_pos(1) fin_pos(1)],[pd_pos(2) fin_pos(2)],'k-');
set(h,'LineWidth',tailThick,'Color',tailclr);
hold off

% Colorbar
%h = colorbar;
%set(h,'TickLabels','','Box','off')

% Axes properties
set(gca,'Color','w','XColor','w','YColor','w')
set(gca,'Position',[0 0 1 1])
%set(gca,'Units','Pixels')
%set(gcf,'Units','Pixels')

title(['t = ' num2str(v.t(frame),'%10.2f\n')],'Color',.5.*[1 1 1])

%     axis e
%     xlim(lims);ylim(lims)


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
