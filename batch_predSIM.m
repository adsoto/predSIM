function batch_predSIM

tic 

% batch_predSIM runs batch simulations of predSIM. The initial position of
% the prey item is defined using ANG and RAD. 

% angles to define initial prey position
ang = linspace(0,pi/2 - pi/12,20);                       % (rad)

% radial distances to define initial prey position
rad = linspace(0.01, 0.3, 20);                    % (m) 

% Initial gain parameter
% kp = 3.0e-3;

% Get size of initial position vectors
N1 = length(rad);
N2 = length(ang);

% Set up data containers
distInit    = zeros(N1,N2);
phiInit     = zeros(N1,N2);
phiFinal    = zeros(N1,N2);
captureInd  = zeros(N1,N2);
totTurns    = zeros(N1,N2);
gainParam   = zeros(N1,N2);
preyAng     = zeros(N1,N2);

% sol         = cell(N1,N2);

% Outer loop
for j=1:N1
    
    % Set current prey radial distance
    currDist = rad(j);
    
    % Inner loop
    parfor k=1:N2
        
        % Initial gain parameter
        kp = 3e-3;
        
        % Set current angle
        currAng = ang(k);

        % Prey initial position
        preyX = currDist * cos(currAng);       % (m)
        preyY = currDist * sin(currAng);       % (m)
        
        % Run the simulation
        sol = run_predSIM(preyX,preyY,kp);
        
        % Store initial gain
        gain_kp = kp;
        
        % Store capture indicator
        capture_ind = sol.capture;
        
        while ~sol.capture && kp < 11e-3 %
            % Adjust gain parameter
            kp = kp * 1.05;
            
            % Collect gain parameter
            gain_kp = [gain_kp; kp];
            
            % rerun simulation with new gain 
            sol = run_predSIM(preyX,preyY,kp);
            
            % Collect capture indicator
            capture_ind = [capture_ind; sol.capture];
        end
        
        % Store simulation results (of last run)
        
        % Initial Distance
        distInit(j,k)   = sol.distInit;
        
        % Initial bearing angle
        phiInit(j,k)    = sol.phiPre(1);
        
        % Final bearing angle
        phiFinal(j,k)   = sol.phiPost(end);
        
        % Capture inidicator
        captureInd(j,k) = sol.capture;
        
        % Total # of turns
        totTurns(j,k)   = sol.turns;
        
        % Gain parameter
        gainParam(j,k)  = kp;
        
        % Angle to define initial prey position
        preyAng(j,k)   = currAng;
        
        % Full solution variables
        solFull{j,k}    = sol;
        
        % All gain parameters
        gainParamAll{j,k}  = gain_kp;
        
%         % Reset gain parameter
%         kp = 3.0e-3;
        
                % Clear results for next iteration
%         clear sol preyX preyY
    end
%     disp(' Running next distance')
end

dOut.dist = distInit;
dOut.phiI = phiInit;
dOut.phiF = phiFinal;
dOut.cap  = captureInd;
dOut.turns= totTurns; 
dOut.gain = gainParam;
dOut.pAng = preyAng;

save('batch_predSim.mat','dOut','solFull','gainParamAll')

toc

end
% 