function [p,sol] = run_predSIM_fin(xPrey,yPrey,kp)
% Runs a numerical simulation of a swimming fish predator.
% INPUTS:   - xPrey, yPrey: coordinates of initial prey position
%           - kp: gain parameter for proportional control
%
% OUTPUT:   - sol: structure with solution output


%%% TO DO: Think of a way to modulate the duration of a glide based on...
% distance and/or bearing angle

%% Simulation Parameters 

% Turn on figures
plotOn = 0;

% Plot force data
plotForce = 1;

% Prey initial position (from input)
if nargin < 2
    p.preyX = 0.1;                       % (m)
    p.preyY = 0.1;                      % (m)
else
    p.preyX = xPrey;
    p.preyY = yPrey;
end


% Gain Parameter
if nargin < 3
    kp = 4.5e-3;
end

% Time span (sec)
p.simDur = 2;

% Maximum step size of simulation (s)
p.maxStep   = 1e-2;

% Relative tolerence of the simulation
p.rel_tol = 1e-4;


%% Morphological and mechanical parameters
% Scaling relations come from McHenry & Lauder (2006)

% Density of fluid (kg/m^3)
p.rho = 1000;

% Dynamic viscosity of water (Pa s)
p.visc = 8.9e-4;

% Body length for a small adult (mm)
bodyL   = 10^1.5;
p.bodyL = bodyL * 10^-3;                % (m)

% Body width (mm)
bodyW   = (6.22e-2) * bodyL^(1.56);
p.bodyW = (bodyW * 10^-3)/1.2;                % (m)

% Body mass (g)
mass    = (4.14E-6) * bodyL^(3.17);
p.mass  = mass * 10^-3;                 % (kg)

% Wetted surface area (mm^2)
surfA   = 3.06E-1 * bodyL^(2.16);
p.SA    = surfA * 10^-6;                % (m^2)

% Body moment of inertia---for a solid ellipsoid about z-axis---(kg m^2)
p.bodyI = (p.mass/5) * (p.bodyL^2 + p.bodyW^2) + p.mass*(0.2*p.bodyL)^2;

%p.bodyI = p.bodyI/10;

% Drag coefficent for coasting zebrafish (dimensionless)
%cDrag   = 1.44E2 * bodyL^(-2.34);
%p.cDrag = cDrag * 10^-3;
p.cDrag = 0.07;

% Lift coefficient of the tail
p.cLift  = 2*pi;

% Rotational drag (dimensionless)
%p.cDrag_rot = 0.02;
p.cDrag_rot = 8*pi;

% Pred initial position
p.predX = 0;                         % (m)
p.predY = 0;                         % (m)

% Pred initial heading
p.theta0 = 0*pi/180;                  % (rad)

% Distance threshold
p.dThresh = 0.5 * p.bodyL;           % (m)

% Initial speed                      % (m/s) 
p.U0 = 0.01;

% Speed of tail beat (m/s)
%p.beatSpd = 0.02;

% Max amplitude of tail heaving (rad)
p.maxHeave = 60 * pi/180;


%% Caudal fin parameters

% Fin length (m); estimate based on literature (Plaut, 2000)
p.finL      = p.bodyL * 0.19 ;

% Peduncle length (m); estimate based on observation & anatomy
p.pedL      = p.bodyL * 0.12;

% Fin height (m); estimate based on literature (Plaut, 2000)
p.finH      = p.bodyL * 0.18;

% Fin span (m^2)
p.finSpan   = p.finH^2;

% Fin surface area (m^2), estimate based on literature (Plaut, 2000)
p.finA      = p.finSpan / 2.05; 

% Heave amplitude (rad)
p.h0        = 0*pi/180;

% Pitch amplitude (rad)
p.pitch0    = 15*pi/180;

% Tail-beat frequency (Hz)
p.tailFreq  = 6;

% Glide duration (s)
p.glideDur = 0.4;

% Phase lag (pitch leads heave) (rad)
p.psi       = 0*pi/180;

% Drag on fin
p.cD_parl   = 0.3;
p.cD_perp   = 0.1;


%% Global variables declared
% These variables are passed to the governing function during the
% simulation

% global s


%% Scale input parameter values for numerical stability
% All parameters used by the model are rescaled, made dimensionless, and
% stored in the 's' structure.

% Scaling factors
sL = 1 / p.bodyL;
sM = 1 / p.mass;
sT = 10^0;

% Store scaling factors in 's' structure
s.SL = sL;
s.sM = sM;
s.sT = sT;

% Scaling factor for force & pressure
sF = sM * sL / sT^2;
sP = sF / sL^2;

% Dimensionless parameters
s.cDrag     = p.cDrag;
s.cDrag_rot = p.cDrag_rot;
s.cLift     = p.cLift;
s.rel_tol   = p.rel_tol;
s.theta0    = p.theta0;
s.psi       = p.psi; 
s.pitch0    = p.pitch0;
s.h0        = p.h0;
s.cD_parl   = p.cD_parl;
s.cD_perp   = p.cD_perp;

% Linear/Area dimensions
s.bodyL     = p.bodyL   * sL;
s.bodyW     = p.bodyW   * sL;
s.SA        = p.SA      * sL^2;
s.preyX     = p.preyX   * sL;
s.preyY     = p.preyY   * sL;
s.predX     = p.predX   * sL;
s.predY     = p.predY   * sL;
s.dThresh   = p.dThresh * sL;
s.finL      = p.finL    * sL;
s.pedL      = p.pedL    * sL;
s.finA      = p.finA    * sL^2;

% Mechanical properties
s.mass      = p.mass    * sM;
s.bodyI     = p.bodyI   * sM    * sL^2;
s.rho       = p.rho     * sM    / sL^3;
s.visc      = p.visc    * sP * sT;

% Time
s.simDur    = p.simDur  * sT;
s.maxStep   = p.maxStep * sT;
s.tailFreq  = p.tailFreq / sT;

% Kinematics
s.U0        = p.U0 * sL / sT;
%s.beatSpd   = p.beatSpd * sL / sT;
s.maxHeave  = p.maxHeave;
s.glideDur  = p.glideDur * sT;

% Indicator variable for capture
s.capture = 0;
s.preyX   = p.preyX *sL;
s.preyY   = p.preyY *sL;

capInd    = 0;

%clear p xPrey yPrey 
 

%% Controller parameters

% Vector from pred to prey (range vector)
rangeX = s.preyX - (s.predX + (0.3*s.bodyL)*cos(s.theta0));
rangeY = s.preyY - (s.predY + (0.3*s.bodyL)*sin(s.theta0));

% Angle of range angle (intertial FOR)
alpha = atan2(rangeY,rangeX);

% Bearing angle (positive is to left of pred, negative to the right)
phi = alpha - s.theta0;

% Deal with zero bearing with a slight offset
if phi==0
    phi = pi/10000;
end

% Direction of the turn (1=CCW, -1=CW)
turnDirec = sign(phi);

% Set turn direction parameter
s.turnDirec = turnDirec;


%% Run ODE solver in a loop

refine = 4;

% Solver options for turning phase
opts = odeset('Events',@turnEvents,'Refine',refine,'RelTol', s.rel_tol);

% Solver options for glide phase
opts2 = odeset('Events',@turnEvents2,'Refine',refine,'RelTol', s.rel_tol);

% Time span for simulation
%tspan = [0 s.simDur];

% Initial conditions in the form: [x, x', y, y', theta, theta']
%init = [s.predX, 1.1*sL, s.predY, 0.0*sL, s.theta, 0];
init = [s.predX, s.U0*cos(s.theta0), s.predY, s.U0*sin(s.theta0), s.theta0, 0, 0, 0];

% Get initial position of fin (saved in 's' structure)
%[s,~] = fin_kine(s,init,tspan(1));

% Distance from body COM to fin quarter-chord point
s.d_bodyfin = 0.7*s.bodyL+s.pedL+0.25*s.finL;

% Initial conditions in the form: [x, x', y, y', theta, theta',xFin,yFin]
%init = [init, s.finPos(1), s.finPos(2)]';
init = [init, -s.d_bodyfin*cos(s.theta0), -s.d_bodyfin*sin(s.theta0)]';

% Initial distance to prey
[~,~,distInit] = controlParams(init);

% Create empty output vectors and counters
iter    = 0;
tout    = 0;
yout    = init';
teout   = [];
yeout   = [];
ieout   = [];
phiPre  = phi;
phiPost = [];

tspan(1) = 0;
 
while ~s.capture

    % Iteration counter, keeps track of beat-glide events
    iter = iter + 1;
    
    % Current time 
    s.tCurr = tspan(1);
    
    % Set turn direction parameter
    s.turnDirec = turnDirec;
    
    % Speed of tail beat
    %TODO: Make this a control parameter
    s.beatSpd = (s.bodyL * s.tailFreq)*6;
    
    % Generate fin kinematics for tail beat 
    s = gen_kinematics(s);
    
    % Simulation period for beat
    tspan(1,2) = tspan(1) + 1/s.tailFreq; 
    
    % Solve ODE (during fin oscillation, 1 sec)
    [t,y,te,ye,ie] = ode15s(@(t,y) predSIM(t,y,s),...
        tspan, init, opts);%[tspan(1),tspan(1)+1], init, opts);
      
    % Accumulate output.  
    nt      = length(t);
    tout    = [tout; t(2:nt)];
    yout    = [yout; y(2:nt,:)];
    teout   = [teout; te];          % Events at tstart are never reported.
    yeout   = [yeout; ye];
    ieout   = [ieout; ie];

    % Set the new initial conditions.
    init = y(nt,:);
    
    % Set the new simulation period for glide
    tspan(1) = t(nt);
    tspan(2) = t(nt) + s.glideDur;
    
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
    
    % Solve ODE (during glide for 0.5 sec)
    [t,y,te,ye,ie] = ode45(@(t,y) predSIM_glide(t,y,s),...
        tspan, init, opts2);
    
    % Accumulate output.
    nt      = length(t);
    tout    = [tout; t(2:nt)];
    yout    = [yout; y(2:nt,:)];
    teout   = [teout; te];          % Events at tstart are never reported.
    yeout   = [yeout; ye];
    ieout   = [ieout; ie];
    
    % Set the new initial conditions.
    init = y(nt,:)';
    
    % Set the new start time
    tspan(1) = t(nt);
    
    % Controller parameters, computed with current state variable values
    [turnDirec,phi,dist] = controlParams(init); 
    
    % Store bearing angle after a glide
%     phiPre = [phiPre; phi];
    
    % check for a distance threshold event (ieout will contain a 1)
    distEvnt = ieout<2;
    if ~isempty(ieout(distEvnt))
        capInd = 1;
        disp('   Target captured')
        break
    else
    end
    
    % Check time interval
    if t(nt)>=s.simDur
        break
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
sol.pitch   = yout(:,7);
sol.heave   = yout(:,8);
sol.dpitch  = yout(:,9)         .* sT;
sol.dheave  = yout(:,10)        .* sT;
sol.phiPre  = phiPre;
sol.phiPost = phiPost;
sol.distInit= distInit  ./ sL;
sol.turns   = iter;
sol.capture = capInd;
sol.params  = s;

% Calculate forces
[sol.lift,sol.torque,sol.drag,sol.drag_theta] = ...
    fin_kine(p,sol.theta,sol.dtheta,sol.pitch,sol.heave,...
             sol.dpitch,sol.dheave,sol.dx,sol.dy);

% Clear others
% clear t y tspan init s sT sL sM


%% Plot force data

if plotForce
   
    subplot(5,1,1)
    plot(sol.t,sol.heave.*180/pi,'-',sol.t,sol.pitch.*180/pi,'-')
    xlabel('t (s)')
    ylabel('Tail angle (deg)')
    legend('h','p')
    grid on
       
    subplot(5,1,2)
    plot(sol.t,sol.lift(:,1).*1000,'-',sol.t,sol.lift(:,2).*1000,'-')                    
    xlabel('t (s)')
    ylabel('Thrust (mN)')
    legend('x','y')
    grid on
    
    subplot(5,1,3)
    plot(sol.t,sol.theta.*180/pi,'-')                    
    xlabel('t (s)')
    ylabel('Heading (deg)')
    grid on

    subplot(5,1,4)
    plot(sol.t,sol.x.*100,'-',sol.t,sol.y.*100,'-')                    
    xlabel('t (s)')
    ylabel('Position (cm)')
    legend('x','y')
    grid on
    
    subplot(5,1,5)
    plot(sol.t,sol.dx.*100,'-',sol.t,sol.dy.*100,'-')                    
    xlabel('t (s)')
    ylabel('Speed (cm/s)')
    legend('x','y')
    grid on
    
    % Trajectory
    figure;
    plot(sol.x.*100,sol.y.*100,'-',sol.x(1).*100,sol.y(1).*100,'o')
    axis equal
    xlabel('x (cm)'); ylabel('y (cm)')
end


%% Plot solutions

%close all

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
    axis equal 
    
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
%% Nested functions -- problem parameters provided by the outer function.
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
        
        %TODO: Fix this.  It's shutting off the solver in the middle of a beat
        
        % Value contains both events that are checked for zero crossings
        value       = [dThresh; rotVel];
        
        % stop the integration if either event is detected (set both to 1)
        %isterminal  = [1; 1]; 
        isterminal = [0; 0];
           
        
        % zero can be approached from either direction for distance
        % threshold and negative direction (decreasing) for rot. velocity
        direction   = [0; -1];      
    end

% -----------------------------------------------------------------------
%
    function [value,isterminal,direction] = turnEvents2(t,y)
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
        
        % stop the integration if either event is detected (set both to 1)
        %isterminal  = [1; 0]; 
        isterminal = [0; 0];
        
        % zero can be approached from either direction for distance
        % threshold and negative direction (decreasing) for rot. velocity
        direction   = [0; 0];      
    end

% -----------------------------------------------------------------------


    function [turnDirec,phi,dist] = controlParams(y)
        % controlParams computes the bearing angle which is then used for 
        % the control input that computes the required thrust parameters
        %
        % INPUT: y contains the current value of all state variables
        
%         global s

        % Heading angle (velocity direction)
        heading = y(5);
%         heading = atan2(y(4),y(2));
        
        % Vector from pred to prey (range vector)
        rangeX = s.preyX - (y(1) + (0.3*s.bodyL)*cos(heading));
        rangeY = s.preyY - (y(3) + (0.3*s.bodyL)*sin(heading));
        
        % Distance to prey (scaled units)
        dist = norm([rangeX, rangeY]);
        
        % Angle of range vector
        alpha = atan2(rangeY,rangeX);
        
        % Bearing angle
        phi = atan2(sin(alpha - heading), cos(alpha - heading));
        
        % Direction of the turn (1=CCW, -1=CW)
        turnDirec = sign(phi);
        
    end


end


function s = gen_kinematics(s)

% Number of end and center points
numendpts  = 1000;
numcntrpts = 300;

% Default duty cycle for heaving of tail
dCycle = 0.2;

% Phase lag of pitch after heave (fraction of tail beat)
lag = 0.05;

% Tail-beat period (includes zero padding)
tau = 1./s.tailFreq;

% Peak tail speed
Vmax = s.beatSpd;

% Duration of zero-padding at start
zStartDur = tau/10;

% Duration of zero-padding at end
zEndDur = tau/10;

% Duration of beat outward
pwr_dur = dCycle * (tau-zStartDur-zEndDur);

% If the tail cannot move far enough . . .
if pwr_dur*Vmax > s.maxHeave
    % Set amplitude at the max
    h_amp = s.maxHeave;
    
    % Calculate a new powerstroke duration
    pwr_dur = s.maxHeave / Vmax;
    
    % And a new duty cycle (Later, we could put a limit on this too)
    dCycle = pwr_dur / (tau-zStartDur-zEndDur);
else
    % Set heave amplitude
    h_amp = pwr_dur*Vmax;
end

% Duration of recovery stroke
rcvr_dur = tau-zStartDur-zEndDur-pwr_dur;

% Time values
t = [linspace(0,zStartDur,numendpts)';...
     linspace(zStartDur,tau-zEndDur,numcntrpts)';...
     linspace(tau-zEndDur,tau,numendpts)'];
 
% Initialize heave and pitch as series of zeros
p = t.*0;
h = t.*0;

% Index of when tail is moving
idx = (t > zStartDur) & (t<(tau-zEndDur));

% Discrete heave values
h(idx) = s.turnDirec * h_amp .* ...
    (sawtooth(2*pi*(t(idx)-zStartDur)./(tau-zStartDur-zEndDur),dCycle)./2+0.5);

% Delay of pitching after heaving
delay = lag * (pwr_dur + rcvr_dur);

% Check that zero-padding is long enough
if zEndDur < delay
    error(['The duration of zero padding must exceed the delay '...
           'between heaving and pitching'])
end

% Index of power stroke
iPwr = (t >= (zStartDur + delay)) & (t<(zStartDur + delay + pwr_dur));

% Index of recovery stroke
iRecov = (t>=(zStartDur + delay + pwr_dur)) & (t<(tau-zEndDur+delay));

% Discrete pitch values during power stroke
p(iPwr) = s.turnDirec * h_amp/2 .* ...
    (sawtooth(2*pi*(t(iPwr)-zStartDur-delay)./(2*pwr_dur)-pi/2,0.5));

% Discrete pitch values during recovery stroke
p(iRecov) = -s.turnDirec * h_amp/2 .* ...
    (sawtooth(2*pi*(t(iRecov)-zStartDur-delay-pwr_dur)./(2*rcvr_dur)-pi/2,0.5));

% Fourier fit to heave data
s.fHeave = fit(t,h,'fourier8');

% Fourier fit to pitch data
s.fPitch = fit(t,p,'fourier8');

% Plot fits
if 0
   figure
   subplot(2,1,1)
   plot(s.fHeave,t,h)
   grid on
   xlabel('t');ylabel('h (rad)')
   
   subplot(2,1,2)
   plot(s.fPitch,t,p)
   xlabel('t');ylabel('p (rad)')
   grid on
end

end

