function [s, fin_L,fin_D] = fin_kine(s,y,t)
% FIN_KINE sets up the motion of the caudal fin used for propulsion.
% Both heave and pitch are modeled as sine functions.
% The angle of attack, alpha, is a function of heave and pitch.
%
% INPUTS:   -p: structure with fin motion parameters
%           -t: time point (ideally, should be a vector)

%% Default parameters (if no input)
if nargin<3
    
    t = 1.0 * linspace(0,1)';
    
    % Set fish body velocity (m/s)
    u_fish = 0.1 * ones(length(t),1);
    
    if nargin<1
        s.finL      = 0.3;
        s.bodyL     = 1;
        s.h0        = 0.0 * s.finL;
        s.pitch0    = 30*pi/180;
        s.tailFreq  = 1;
        s.psi       = 70*pi/180;
        s.cD_parl   = 0.3;
        s.cD_perp   = 0.1;
        s.pedL      = 0.5 * s.finL;
        s.rho       = 1000;
        s.finA      = s.finL*(0.25*s.finL); 
        
    end
else
    % Check that input y, is size n x 6.
    sz_in = size(y);
    if sz_in(2)~=6
        y = y';
    end
        
    % Linear speed of fish
    u_fish = sqrt(y(:,2).^2 + y(:,4).^2);
%     u_fish = u_fish';
end



%% Motion Parameters

% Max heave amplitude (m)
h0 = s.h0;

% Max pitch amplitude (rad)
pitch0 = s.pitch0;

% Tail-beat frequency of fin (s^-1)
freq = s.tailFreq;

% Angular frequency of fin (rad/s);
omega = 2*pi*freq;

% Phase lag (between pitch and heave, pitch leads)
psi = s.psi;

% % Check for zero fish velocity
% if u_fish==0
%     % Set speed of fish to max. heave velocity
%     u_fish = h0*omega;
% end

%% Fin motion Equations

% Fish velocities and heading angle
if nargin < 2
    hd_ang = 0*pi/180 .* ones(length(t),1);
    x_vel = u_fish .* cos(hd_ang);
    y_vel = u_fish .* sin(hd_ang);
    hd_vel = zeros(length(t),1);
    x_fish = x_vel .* t; 
    y_fish = y_vel .* t;
else
    x_vel = y(:,2);
    y_vel = y(:,4);
    hd_ang = y(:,5);
    hd_vel = y(:,6);
    x_fish = y(:,1);
    y_fish = y(:,3);
    % Heading (derived from velocity)
%     hd_ang = atan2(y_vel,x_vel);
end

% Heave equation (peduncle angle relative to fish body midline)
heave = h0 * sin(omega * t);

% Pitch equation (tail angle relative to peduncle midline)
pitch = pitch0 * sin(omega * t + psi);

% Derivative of heave
h_prime = h0 * omega * cos(omega * t);

% Derivative of pitch
p_prime = pitch0 * cos(omega * t + psi) * omega;

% Angle of attack (needs velocity of fish, u_fish or max heave velocity)
alpha = -atan(h_prime ./ u_fish) + (pitch + hd_ang);

% Sum of alpha and pitch
beta = heave + pitch + hd_ang; %alpha + pitch;

%% Compute forces acting on fin
   
% Current position of fin quarter-chord point (w.r.t. inertial FOR)
s.finPos(:,1) = x_fish-0.7*s.bodyL*cos(hd_ang) - s.pedL*cos(hd_ang+heave) ...
    -0.25*s.finL*cos(hd_ang+heave+pitch);

s.finPos(:,2) = y_fish-0.7*s.bodyL*sin(hd_ang) - s.pedL*sin(hd_ang+heave) ...
    -0.25*s.finL*sin(hd_ang+heave+pitch);

% Current position of fin 1/4 chord point (inertial FOR)
fp = s.finPos;

% Current position of fin quarter-chord point (w.r.t. body COM)
s.finPos_body(:,1) = -0.7*s.bodyL*cos(hd_ang) - s.pedL*cos(hd_ang+heave) ...
    -0.25*s.finL*cos(hd_ang+heave+pitch);

s.finPos_body(:,2) = -0.7*s.bodyL*sin(hd_ang) - s.pedL*sin(hd_ang+heave) ...
    -0.25*s.finL*sin(hd_ang+heave+pitch);

% Speed of fin in x-direction (interial frame)
u_parl = x_vel + 0.7*s.bodyL*sin(hd_ang).*hd_vel + ...
        s.pedL*sin(hd_ang+heave).*(hd_vel+h_prime) + ...
        0.25*s.finL*sin(hd_ang+heave+pitch).*(hd_vel+h_prime+p_prime);
    
% Speed of fin in y-direction (interial frame)
u_perp = y_vel - 0.7*s.bodyL*cos(hd_ang).*hd_vel - ...
        s.pedL*cos(hd_ang+heave).*(hd_vel+h_prime) - ...
        0.25*s.finL*cos(hd_ang+heave+pitch).*(hd_vel+h_prime+p_prime);
    
% Unit vector pointed in direction of fin leading edge (w.r.t. fixed axes)
unit_f = [cos(hd_ang+heave+pitch), sin(hd_ang+heave+pitch),...
            zeros(length(heave),1)];

% Construct velocity vector for cross product calculation
vel_vec = [u_parl, u_perp, zeros(length(u_parl),1)];  

% Cross products of velocity with unit vectors (may be size Nx3)
x_prod = cross(vel_vec,unit_f,2);

% Lift
% fin_L = (0.5*s.rho*s.cD_parl).*abs((u_perp)).*(u_perp);
fin_L = pi*s.rho*s.finA .* (cross(x_prod,vel_vec,2));

% if fin_L(1) > 1
%     pause
% end

% Drag
fin_D = (0.5*s.rho*s.cD_perp).*abs((u_parl)).*(u_parl);

if 0
    % Plots to verify velocity (of fin) and lift are orthogonal
    figure
    
    % Plot the orientation of the fin starting at the quarter chord point
    quiver(fp(:,1),fp(:,2),unit_f(:,1),unit_f(:,2),0.5)
    hold on
    % Plot the lift force vectors
    quiver(fp(:,1),fp(:,2),fin_L(:,1),fin_L(:,2))
    
    % Plot the quarter chord point velocity vector
    quiver(fp(:,1),fp(:,2),u_parl(:),u_perp(:),0.5)
    axis equal
    
    % Norm of unit vectors (should be 1)
    unit_norm = sqrt(sum(unit_f.^2,2));
    
    % Norm of velocity vector
    vel_norm = sqrt(sum(vel_vec.^2,2));
    
    % Dot product between orientation of fin and velocity
    orient_dot_vel = dot(unit_f,vel_vec,2);
    
    % Cross product between orientation of fin and velocity
    orient_x_vel = cross(unit_f,vel_vec,2);
    
    % Norm of cross products
    cross_norm = sqrt(sum(orient_x_vel.^2,2));
    
    % Angle between fin orientation and velocity vector
    atck_ang = atan2(cross_norm,orient_dot_vel);
    
    figure, plot(t,atck_ang*180/pi)
    
end

% function [u_parl,u_perp] = getSpeed(u_fish,pitch,heave,h_prime,fp,t)
% 
% % Rigid Body transformation matrix for current time
% rot_mat = [cos(pitch) sin(pitch) -u_fish*t;...
%     -sin(pitch) cos(pitch) heave;...
%     0 0 1];
% 
% warning off
% % Inverse of transformation
% % rot_inv = inv(rot_mat);
% 
% warning on
% 
% % Time derivative of transformation
% rot_prime = [-sin(pitch) cos(pitch) -u_fish;...
%              -cos(pitch) -sin(pitch) h_prime;...
%              0 0 0];
% 
% % Angular velocity matrix
% omega_mat = rot_prime / rot_mat;
% 
% % Velocity of center of mass of fin
% u_parl = omega_mat(1,3) + fp(2);
% u_perp = omega_mat(2,3) + fp(1);



%% For test simulations

if 0
    % Define airfoil geometry
    x_foil = linspace(0,s.finL,500) - s.fp(1,1);
    y_foil = 0.6*(0.2969*sqrt(x_foil) - 0.1260*x_foil ...
        - 0.3516*x_foil.^2 + 0.2843*x_foil.^3 - 0.1036*x_foil.^4);
    
    x_foil = -x_foil;
    
    % Close al figures
    close all
    
    % Create figure
    h = figure;
    
    % Set figure size
    h.Position = [100, 100, 1080, 900];
    
    % Set figure size for saving
    h.PaperUnits = 'inches';
    h.PaperPosition = [0 0 6 6];
    
    % Plot foil (before motion)
    h1 = subplot(2,2,[1,2]);
    plot(x_foil,y_foil,'k-','LineWidth',2)
    hold on
    plot(x_foil,-y_foil,'k-','LineWidth',2)
    h1.XLim = ([-1.05 1.05]);
    h1.YLim = ([-0.75 0.75]);
    hold off
    xlabel('x-position (cm)','FontSize',16)
    ylabel('y-position (cm)','FontSize',16)
    
    % Plot heave curve
    subplot(2,2,3);
    plot(t,heave,'LineWidth',2)
    hold on
    ax = gca;
    h_line = line([t(1),t(1)],[ax.YLim(1),ax.YLim(2)],'Color','k');
    hold off
    xlabel('time (s)','FontSize',16)
    ylabel('Heave Amplitude (cm)','FontSize',16)
    
    % Plot pitch curve
    subplot(2,2,4);
    plot(t,pitch*180/pi,'LineWidth',2)
    hold on
    ax = gca;
    p_line = line([t(1),t(1)],[ax.YLim(1),ax.YLim(2)],'Color','k');
    hold off
    xlabel('time (s)','FontSize',16)
    ylabel('Pitch Angle (deg)','FontSize',16)
    
    drawnow
    % Set path for saving movie stills
    %     [~, pathMovie] = uiputfile;
    
    % Save first frame
    %     print(h,[pathMovie filesep 'airfoilMovie' num2str(0)],'-djpeg')
    
    for j=1:length(t)
        
        % current heave value
        heave_curr = heave(j);
        
        % current pitch angle
        pitch_curr = pitch(j);
        
        % rotation + translation matrix (homogeneous coordinates)
        rot_mat = [cos(pitch_curr) sin(pitch_curr) u_fish(j)*t(j);...
            -sin(pitch_curr) cos(pitch_curr) heave_curr;...
            0 0 1];
        
        % Transformed x & y values
        top_foil = rot_mat * [x_foil ; y_foil; ones(size(x_foil))];
        bot_foil = rot_mat * [x_foil ; -y_foil; ones(size(x_foil))];
        
        x_new = top_foil(1,:);
        y_top = top_foil(2,:);
        y_bot = bot_foil(2,:);

        fin_l = sqrt((x_new(1) - x_new(end))^2 + (y_top(1) - y_top(end))^2);
        
        % Update foil position
        plot(h1,x_new,y_top,'k','LineWidth',2)
        hold(h1,'on')
        plot(h1,x_new,y_bot,'k','LineWidth',2)
        h1.XLim = ([-1.05 1.05]);
        h1.YLim = ([-0.75 0.75]);
        hold(h1,'off')
        xlabel(h1,'x-position (cm)','FontSize',16)
        ylabel(h1,'y-position (cm)','FontSize',16)
        
        % Update vertical line position for heave
        h_line.XData = [t(j),t(j)];

        % Update vertical line position for pitch
        p_line.XData = [t(j),t(j)];

        drawnow
%         pause
        
        if 0%makeMovie
            
            % Save figure window as jpeg
            print(h,[pathMovie filesep 'airfoilMovie' num2str(j)],'-djpeg')
        end
    end
    
end


