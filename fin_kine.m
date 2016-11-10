function [beta,fin_L,fin_D] = fin_kine(s,y,t)
% FIN_KINE sets up the motion of the caudal fin used for propulsion.
% Both heave and pitch are modeled as sine functions.
% The angle of attack, alpha, is a function of heave and pitch.
%
% INPUTS:   -p: structure with fin motion parameters
%           -t: time point (ideally, should be a vector)

%% Default parameters (if no input)
if nargin<2
    
    t = linspace(0,1)';
    
    % Set fish body velocity (m/s)
    u_fish = 1 * zeros(length(t),1);
    
    if nargin<1
        s.h0        = 0.1;
        s.pitch0    = 30*pi/180;
        s.tailFreq  = 1;
        s.psi       = 90*pi/180;
        s.cD_parl   = 0.3;
        s.cD_perp   = 0.1;
        s.finL      = 2*s.h0;
        s.rho       = 1000;
        
    end
else
    % Linear speed of fish
    u_fish = sqrt(y(2)^2 + y(4)^2);
end



%% Motion Parameters

% Max heave amplitude (Length)
h0 = s.h0;

% Max pitch amplitude (rads)
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

% Heave equation (peduncle angle relative to body midline)
heave = h0 * sin(omega * t);

% Pitch equation (tail angle relative to inertial horizontal axis)
pitch = pitch0 * sin(omega * t + psi);

% Derivative of heave
h_prime = h0 * omega * cos(omega * t);

% Angle of attack (needs velocity of fish, u_fish or max heave velocity)
alpha = -atan(h_prime ./ u_fish) + pitch;

% Sum of alpha and pitch
beta = alpha + pitch;

%% Compute forces acting on fin

% Current position of fin 1/4 chord point
fp = s.fp;

[u_parl,u_perp] = getSpeed(u_fish,pitch,heave,h_prime,fp,t);

% % Fin speed parallel to longitudinal axis of fin
% uFin_parl = u_fish .* cos(pitch);
% 
% % Fin speed perpendicular to longitudinal axis of fin
% uFin_perp = u_fish .* sin(pitch);

% Lift
fin_L = (0.5*s.rho*s.cD_parl).*abs((u_perp)).*(u_perp);

% Drag
fin_D = (0.5*s.rho*s.cD_perp).*abs((u_parl)).*(u_parl);

function [u_parl,u_perp] = getSpeed(u_fish,pitch,heave,h_prime,fp,t)

% Rigid Body transformation matrix for current time
rot_mat = [cos(pitch) sin(pitch) -u_fish*t;...
    -sin(pitch) cos(pitch) heave;...
    0 0 1];

warning off
% Inverse of transformation
% rot_inv = inv(rot_mat);

warning on

% Time derivative of transformation
rot_prime = [-sin(pitch) cos(pitch) -u_fish;...
             -cos(pitch) -sin(pitch) h_prime;...
             0 0 0];

% Angular velocity matrix
omega_mat = rot_prime / rot_mat;

% Velocity of center of mass of fin
u_parl = omega_mat(1,3) + fp(2);
u_perp = omega_mat(2,3) + fp(1);


%% For test simulations

if 0
    % Define airfoil geometry
    x_foil = linspace(0,1,500);
    y_foil = 0.6*(0.2969*sqrt(x_foil) - 0.1260*x_foil ...
        - 0.3516*x_foil.^2 + 0.2843*x_foil.^3 - 0.1036*x_foil.^4);
    
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
    h1.XLim = ([-0.5 1.05]);
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
        rot_mat = [cos(pitch_curr) sin(pitch_curr) -u_fish(j)*t(j);...
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
        h1.XLim = ([-0.5 1.05]);
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


