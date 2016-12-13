function sol = run_predSIM(xPrey,yPrey,kp)
% Runs a numerical simulation of a swimming fish predator.
% INPUTS:   - xPrey, yPrey: coordinates of initial prey position
%           - kp: gain parameter for proportional control
%
% OUTPUT:   - sol: structure with solution output


%%% TO DO: Think of a way to modulate the duration of a glide based on...
% distance and/or bearing angle

%% Simulation Parameters 

% Turn on figures
plotOn = 1;

% Prey initial position (from input)
if nargin < 2
    p.preyX = 1;                       % (m)
    p.preyY = 1;                      % (m)
else
    p.preyX = xPrey;
    p.preyY = yPrey;
end

% Gain Parameter
if nargin < 3
    kp = 3.5e-3;
end

% Time span (sec)
p.simDur = 20;

% Maximum step size of simulation (s)
p.maxStep   = 1e-2;

% Relative tolerence of the simulation
p.rel_tol = 1e-7;

%% Morphological and mechanical parameters
% Scaling relations come from McHenry & Lauder (2006)

% Density of fluid (kg/m^3)
p.rho = 1000;

% Body length for a small adult (mm)
bodyL   = 10^1.4;
p.bodyL = bodyL * 10^-3;          % (m)

% Body width (mm)
bodyW   = (6.22e-2) * bodyL^(1.56);
p.bodyW = bodyW * 10^-3;          % (m)

% Body mass (g)
mass    = (4.14E-6) * bodyL^(3.17);
p.mass  = mass * 10^-3;            % (kg)

% Wetted surface area (mm^2)
surfA   = 3.06E-1 * bodyL^(2.16);
p.SA    = surfA * 10^-6;               % (m^2)

% Body moment of inertia---for a solid ellipsoid about z-axis---(kg m^2)
p.bodyI   = (p.mass/5) * (p.bodyL^2 + p.bodyW^2) + p.mass*(0.2*p.bodyL)^2;

% Drag coefficent for coasting zebrafish (dimensionless)
cDrag   = 1.44E2 * bodyL^(-2.34);
p.cDrag = cDrag * 10^-3;

% Rotational drag (dimensionless)
p.cDrag_rot = 0.02;

% Pred initial position
p.predX = -0.1;                       % (m)
p.predY = -0.2;                       % (m)

% Pred initial heading
p.theta = pi/10;                        % (rad)

% Distance threshold
p.dThresh = 2 * p.bodyL;       % (m)

%% Global variables declared
% These variables are passed to the governing function during the
% simulation

global s

%% Scale input parameter values for numerical stability
% All parameters used by the model are rescaled, made dimensionless, and
% stored in the 's' structure.

% Scaling factors
sL = 1 / p.bodyL;
sM = 1 / p.mass;
sT = 10^0;

% Dimensionless parameters
s.cDrag     = p.cDrag;
s.cDrag_rot = p.cDrag_rot;
s.rel_tol   = p.rel_tol;
s.theta     = p.theta;

% Linear/Area dimensions
s.bodyL     = p.bodyL   * sL;
s.bodyW     = p.bodyW   * sL;
s.SA        = p.SA      * sL^2;
s.preyX     = p.preyX   * sL;
s.preyY     = p.preyY   * sL;
s.predX     = p.predX   * sL;
s.predY     = p.predY   * sL;
s.dThresh   = p.dThresh * sL;

% Mechanical properties
s.mass      = p.mass    * sM;
s.bodyI     = p.bodyI   * sM    * sL^2;
s.rho       = p.rho     * sM    / sL^3;

% Time
s.simDur    = p.simDur  * sT;
s.maxStep   = p.maxStep * sT;

% Indicator variable for capture
s.capture = 0;

capInd = 0;

%% Controller parameters

% Vector from pred to prey (range vector)
rangeX = s.preyX - (s.predX + (0.3*s.bodyL)*cos(s.theta));
rangeY = s.preyY - (s.predY + (0.3*s.bodyL)*sin(s.theta));

% Angle of range angle (intertial FOR)
alpha = atan2(rangeY,rangeX);

% Bearing angle (positive is to left of pred, negative to the right)
phi = atan2(sin(alpha - s.theta), cos(alpha - s.theta));

% Direction of the turn (1=CCW, -1=CW)
s.turnDirec = sign(phi);


%% Run ODE solver in a loop

refine = 4;

% Solver options
opts = odeset('Events',@turnEvents,'Refine',refine,'RelTol', s.rel_tol);

% Time span for simulation
tspan = [0 s.simDur];

% Initial conditions in the form: [x, x', y, y', theta, theta']
init = [s.predX, 0, s.predY, 0, s.theta, 0]';

% Initial distance
[~,~,distInit] = controlParams(init);

% Initial thrust parameters
[tF,F_parl,F_norm] = thrustFnc(tspan(1),s.turnDirec,phi,kp);

% Create empty output vectors and counters
iter    = 0;
tout    = 0;
yout    = init';
teout   = [];
yeout   = [];
ieout   = [];
phiPre  = phi;
phiPost = [];
 
while ~s.capture
    
    % Iteration counter, keep track of beat-glide events
    iter = iter + 1; 
    
    % Solve ODE during turn (time dependent terms passed as params)
    [t,y,te,ye,ie] = ode15s(@(t,y) predSIM_turn(t,y,tF,F_norm,F_parl),...
        tspan, init, opts);
    
    % Accumulate output.  This could be passed out as output arguments.
    nt      = length(t);
    tout    = [tout; t(2:nt)];
    yout    = [yout; y(2:nt,:)];
    teout   = [teout; te];          % Events at tstart are never reported.
    yeout   = [yeout; ye];
    ieout   = [ieout; ie];

    % Set the new initial conditions.
    init = y(nt,:);
    
    % Set the new start time
    tspan(1) = t(nt);
    
    % Bearing angle after a turn
    [~,phiTurn,~] = controlParams(init); 
    phiPost = [phiPost; phiTurn];
    
    % check for a distance threshold event (ieout will contain a 1)
    distEvnt = ieout<2;
    if ~isempty(ieout(distEvnt))
        capInd = 1;
        disp('   Target captured')
        break
    else
    end
  
    % A good guess of a valid first timestep is the length of the last valid
    % timestep, so use it for faster computation.  'refine' is 4 by default.
    %     opts = odeset(opts,'InitialStep',t(nt)-t(nt-refine),...
    %         'MaxStep',t(nt)-t(1));
    
%     % Set glide duration based on current distance 
%     [~,~,distCurr] = controlParams(init);
%     glideDur = distCurr 
    
    % Solve ODE during glide for 0.4 sec
    [t,y,te,ye,ie] = ode45(@(t,y) predSIM_turn(t,y),...
        [t(nt),t(nt)+0.4], init, opts);
    
    % Accumulate output.
    nt      = length(t);
    tout    = [tout; t(2:nt)];
    yout    = [yout; y(2:nt,:)];
    teout   = [teout; te];          % Events at tstart are never reported.
    yeout   = [yeout; ye];
    ieout   = [ieout; ie];
    
    % Set the new initial conditions.
    init = y(nt,:);
    
    % Set the new start time
    tspan(1) = t(nt);
    
    % Controller parameters, computed with current state variable values
    [turnDirec,phi,dist] = controlParams(init); 
    
    % Store bearing angle after a glide
    phiPre = [phiPre; phi];
    
    % check for a distance threshold event (ieout will contain a 1)
    distEvnt = ieout<2;
    if ~isempty(ieout(distEvnt))
        capInd = 1;
        disp('   Target captured')
        break
    else
    end
    
    % Check time interval
    if t(nt)>=tspan(2)
        break
    else
        % Thrust parameters for next iteration (begins with turn)
        [tF,F_parl,F_norm] = thrustFnc(t(nt),turnDirec,phi,kp);
    end
    
end

% Store results
sol.t       = tout              ./ sT;
sol.x       = yout(:,1) ./ sL;
sol.y       = yout(:,3) ./ sL;
sol.theta   = yout(:,5);
sol.dx      = yout(:,2) ./ sL   .* sT;
sol.dy      = yout(:,4) ./ sL   .* sT;
sol.dtheta  = yout(:,6)         .* sT;
sol.phiPre  = phiPre;
sol.phiPost = phiPost;
sol.distInit= distInit  ./ sL;
sol.turns   = iter;
sol.capture = capInd;

% Clear others
% clear t y tspan init s sT sL sM

%% Plot solutions

close all

if plotOn
    
    % Plot heading angle (from solution)
    figure,
    plot(sol.t, sol.theta*180/pi,'LineWidth', 2)
    ylabel('Heading (deg)')
    xlabel('time (s)')
    
    set(findall(gcf,'-property','FontSize'),'FontSize',14)
    
    figure,
    
    % Plot position
    subplot(2,2,[1;3])
    plot(sol.x,sol.y,'LineWidth', 2)
    hold on, plot(p.preyX,p.preyY,'or'), hold off
    ylabel('y position')
    xlabel('x position')
    
    % Plot heading (derived from velocity)
    heading = atan2(sol.dy,sol.dx);
    subplot(2,2,2)
    plot(sol.t,unwrap(heading)*180/pi,'LineWidth', 2)
    ylabel('Velocity Direction (deg)')
    xlabel('time (s)')
    
    % Plot speed
    speed = sqrt(sum([sol.dx,sol.dy].^2,2));
    subplot(2,2,4)
    plot(sol.t,speed,'LineWidth', 2)
    ylabel('Speed (m/s)')
    xlabel('time (s)')
    
    set(findall(gcf,'-property','FontSize'),'FontSize',14)

end


% -----------------------------------------------------------------------
% Nested functions -- problem parameters provided by the outer function.
%
    function [value,isterminal,direction] = turnEvents(t,y)
        % Locate the time when a turn is completed or when the distance
        % threshold is satisfied
        
        % Get current distance to prey
        [~,~,dist]  = controlParams(y);
        
        % Detect distance threshold 
        dThresh     = (dist - s.dThresh) - 1e-5;
        
        % Detect rotational velocity = 0; (turn completed) 
        % look at absolute value so that crossings are from negative direc.
        rotVel      = abs(y(6)) - 1e-2;  
        
        % Value contains both events that are checked for zero crossings
        value       = [dThresh; rotVel];
        
        % stop the integration if either event is detected
        isterminal  = [1; 1]; 
        
        % zero can be approached from either direction for distance
        % threshold and negative direction (decreasing) for rot. velocity
        direction   = [0; -1];      
    end

% -----------------------------------------------------------------------

    function [turnDirec,phi,dist] = controlParams(y)
        % controlParams computes the bearing angle which is then used for 
        % the control input that computes the required thrust parameters
        %
        % INPUT: y contains the current value of all state variables
        
%         global s

        % Heading angle (velocity direction)
%         heading = y(5);
        heading = atan2(y(4),y(2));
        
        % Vector from pred to prey (range vector)
        rangeX = s.preyX - (y(1) + (0.3*s.bodyL)*cos(heading));
        rangeY = s.preyY - (y(3) + (0.3*s.bodyL)*sin(heading));
        
        % Distance to prey (scaled units)
        dist = norm([rangeX, rangeY]);
        
        alpha = atan2(rangeY,rangeX);
        
        % Bearing angle
        phi = atan2(sin(alpha - heading), cos(alpha - heading));
        
        % Check that bearing is between -pi and pi
%         if phi >= pi
%             phi = 2*pi - phi;
%         elseif phi < -pi
%             phi = 2*pi + phi;
%         end
        
        % Direction of the turn (1=CCW, -1=CW)
        turnDirec = sign(phi);
        
    end

% -----------------------------------------------------------------------

    function [tF,F_parl,F_norm] = thrustFnc(tStart,turnDirec,phi,kp)
        % thrustFnc defines the thrust force (time dependent & relative to heading)
        % from the control parameters
        
%         global s
        
        % Gain constant for controller (N/rad)
        if nargin<4
            kp = 2.0e-3;
        end
        
        % Magnitude of thrust (N) (Estimate from Weihs 1972)
        % NOTE: This is also the (P) control variable in the model
        % p.fMag = 20e-3;
        p.fMag = kp * abs(phi);
        
        % Thrust (rescaled)
        s.fMag = p.fMag * sM * sL / sT^2;
        
        % Time vector for thrust pulse
        tF = (tStart:1/500:s.simDur)';
        
        % Ratio between ON duration and OFF duration (<= 1)
        % This parameter is related to the duration between tail beats
        r = (1/3) / s.simDur;
        
        % Number of ON intervals (i.e., number of thrust pulses)
        nInt = 1;
        
        % Define a rectangular pulse (alternating 'on' and 'off' intervals)
        w   = s.simDur / (1+1/r) / nInt;        % ON pulse width
        d   = w/2 + tF(1):w*(1+1/r):s.simDur;   % delay vector, defines ON periods
        f_pulse  = pulstran(tF,d,'rectpuls',w); % pulse train
        
        % Force perpendicular to long axis (generates turning moment)
        F_norm = (3*s.fMag/5) * f_pulse * turnDirec;
        
        % Forward thrust parallel to long axis
        % F_parl  = s.fMag * cos(phi);
        F_parl = (2*s.fMag/5) * f_pulse;
        
    end

% -----------------------------------------------------------------------

end


% function [transVel,isterminal,direction] = glideEvents(t,y)
% % Locate the time when height passes through zero in a decreasing direction
% % and stop integration.
% % Locate the time when a turn is completed
%
% transVel      = y(6);  % detect rotational velocity = 0
% isterminal  = 1;     % stop the integration
% direction   = 0;     % zero can be approached from either direction


