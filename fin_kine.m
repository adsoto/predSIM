function [gama,v_n,a_n] = fin_kine(p,t)
% FIN_KINE sets up the motion of the caudal fin used for propulsion.
% Both heave and pitch are modeled as sine functions.
% The angle of attack, alpha, is a function of heave and pitch. 
%
% INPUTS:   -p: structure with fin motion parameters
%           -t: time point (ideally, should be a vector) 

%% Default parameters (if no input)
if nargin<1
    p.h0        = 0.5;
    p.theta0    = 30*pi/180;
    p.tailFreq  = 2;
    p.psi       = 60*pi/180;
    
    t = linspace(0,1);
end

% Set fish body velocity (m/s)
u_fish = 0.02; 

%% Motion Parameters
% Max heave amplitude (Length)
h0 = p.h0;

% Max pitch amplitude (rads)
theta0 = p.theta0;

% Tail-beat frequency of fin (s^-1)
freq = p.tailFreq;

% Circular frequency of fin (rad/s);
omega = 2*pi*freq;

% Phase lag (between pitch and heave, pitch leads)
psi = p.psi;

%% Equations

% Heave equation
heave = h0 * sin(omega * t);

% Pitch equation
theta = theta0 * sin(omega * t + psi);

% Derivative of heave
h_prime = h0 * omega * cos(omega * t);

% Angle of attack (needs velocity of fish, u_fish or max heave velocity)
alpha = -atan(h_prime ./ u_fish) + theta;

%% For test simulations

if 1
    % Define airfoil geometry
    x_foil = linspace(0,1,500);
    y_foil = 0.6*(0.2969*sqrt(x_foil) - 0.1260*x_foil ...
                 - 0.3516*x_foil.^2 + 0.2843*x_foil.^3 - 0.1036*x_foil.^4);

     % Create figure
     h = figure;

     % Set figure size
     h.Position = [100, 100, 1080, 900];

     % Set figure size for saving
     h.PaperUnits = 'inches';
     h.PaperPosition = [0 0 7.2 6];
    
    % Plot foil (before motion)
    plot(x_foil,y_foil,'k-','LineWidth',2)
    hold on
    plot(x_foil,-y_foil,'k-','LineWidth',2)
   	xlim([-0.5 1.5])
    ylim([-1 1])
    hold off
    drawnow
    
    % Set path for saving movie stills
    [~, pathMovie] = uiputfile;
    
    % Save first frame
    print(h,[pathMovie filesep 'airfoilMovie' num2str(0)],'-djpeg')
    
    for j=1:length(t)
        
        % current heave value
        heave_curr = heave(j);
        
        % current pitch angle
        pitch = theta(j);
        
        % rotation + translation matrix (homogeneous coordinates)
        rot_mat = [cos(pitch) sin(pitch) 0;...
                   -sin(pitch) cos(pitch) heave_curr;...
                    0 0 1];
        
        % Transformed x & y values
        top_foil = rot_mat * [x_foil; y_foil; ones(size(x_foil))];
        bot_foil = rot_mat * [x_foil; -y_foil; ones(size(x_foil))];

        x_new = top_foil(1,:);
        y_top = top_foil(2,:);
        y_bot = bot_foil(2,:);
        
        plot(x_new,y_top,'k','LineWidth',2)
        hold on
        plot(x_new,y_bot,'k','LineWidth',2)
        xlim([-0.5 1.5])
        ylim([-1 1])
        hold off
        drawnow 
        
        if 1%makeMovie 
            
            % Save figure window as jpeg
            print(h,[pathMovie filesep 'airfoilMovie' num2str(j)],'-djpeg')
        end
    end
    
end
              

