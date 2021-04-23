clc
close all
clear all
%% QUESTION: Target is now moving along x-axis with velocity of 20m/s Target starting position is 10 m along x-axis with
%%%% sampling time T=0.01 second, total time=5 seconds. Measurement noise have zero mean and 1.2 standarddeviation.
%%************************************************************************%%%
%%%                                  Steps                                %%%
%%% 1.Global variable delaration and common variables                     %%%
%%% 2. Kalman Filter                                                      %%%
%%%   2.1. The prediction part                                            %%%
%%%   2.2. The estimation Part                                            %%%                                                           %%%
%%% 3. Averaging the results to eleminate uncertainity (MonteCarlo Runs)  %%%
%%% 4. Statistics (Root mean Sauare Error)                                %%%
%%% 5. Plotting                                                           %%%
%%%***********************************************************************%%%

%% Global variable delaration and common variables

Runs = 5;
T = 0.01;             % sampling Time
total_time = 5;
t =0: T: total_time;        % Time vector
TotalScans = length(t);

mean = 0;            % measurement noise mean
meas_noise_sd = 1.2; % measurement noise standard deviation

nstates = 2;         % numbers of States
nmeas = 1;           % number of measurements
init_pos = [110]; % True inital Position (m) of the target on x and y axis [x; y]
V = [20];          % True Velocity (m/s) of the target on x and y axis [vx; vy]

%% Kalman Filter
F=[eye(1) T*eye(1);zeros(1) eye(1)]; % Transition matrix for constant velocity model

H=[eye(nmeas) zeros(nstates-nmeas)];   %Output Coeffecient Matrix 

R=(meas_noise_sd^2)*eye(nmeas); % measurement covariance matrix


for run = 1 : Runs
    
    targetX(1) = init_pos(1);
   
    target_velocityX = V(1);
   
    
    for scan = 2 : TotalScans
               
        targetX(scan) = targetX(scan - 1) + T*target_velocityX;    % True/actual position: all scans
       

    end
    
    Ztrue = [targetX ];    % True/Actual target position
    
    meas_noise = mean + R*randn(nmeas, TotalScans);    % generating measurement noise (gaussian noise) 
    Znoisy = Ztrue + meas_noise;                       % The noisy measurment vector
    
    X0 = [Znoisy(:,1) ; V+diag(R).*randn( size(V) )];                   % inital State for first scan only  
    P0 = [R*eye(nmeas)  zeros(nmeas) ; zeros(nmeas)  (R^2)*eye(nmeas)/T];  % inital State Covariacne for first scan only 

    for scan = 1 : TotalScans

        if scan > 1
            
            Xpred=F*Xest;       % Predicted state vector
Ppred=F*Pest*F';    % Predicted state covariance matrix

S = H*Ppred*H' + R;     % Innovation covariance
K = Ppred*H'*inv(S);    % optimal kalman gain

Xest = Xpred + K*( Znoisy(:,scan) - H*Xpred );    % a-posterori state estimate
Pest = ( eye(nstates) - K*H )*Ppred;        % a-posterori estimate Covariance
            
           
        else
            
            Xest = X0;      % inital State declaration (first scan only)
            Pest = P0;      % inital State covariance declaration (first scan only)
            
        end
        
        EstimatedX(run,scan) = Xest(1);
     
        
        Estimated_Vx(run,scan) = Xest(2);
      
        
    end
    
    Zx(run,:)=Znoisy(1,:);
   
end

% Averaging the results to eleminate uncertainity (MonteCarlo Runs)

Monte_estimateX = sum(EstimatedX)/Runs;

Monte_estimateVx = sum(Estimated_Vx)/Runs;


Monte_Znoisyx = sum(Zx)/Runs;


%% Statistics (Root mean Sauare Error)
 TruevelocityX = repmat(V(1),Runs, TotalScans);  % generating veloxity vector for x-axis
  
    
    RMSEX = sqrt( sum( (repmat(Ztrue(1,:),Runs,1) - EstimatedX).^2 )/Runs ); % Computing Root Mean Square error in x-position
    
    RMSEVx = sqrt( sum( (TruevelocityX - Estimated_Vx).^2 )/Runs ); % Computing Root Mean Square error in x-velocity
   
    PositionRMSE = sqrt( RMSEX.^2); % computing the total or resultant position RMSE (both in x and y position)
    VelocityRMSE = sqrt( RMSEVx.^2); % computing the total or resultant velocity RMSE (Both in Vx and Vy)
    
    RMSE = [ PositionRMSE; VelocityRMSE];
    

%% Plotting 
Estimate=[Monte_estimateX;  Monte_estimateVx]; 
Monte_Znoisy=[Monte_Znoisyx];
 figure
% subplot(211)
  plot(t,Ztrue,'r','LineWidth',1.5)
 hold on
 plot(t,Monte_Znoisyx,'g','LineWidth',1.5)
  hold on
  plot(t,Monte_estimateX,'b','LineWidth',1.25)
  legend('True Position','Measured Position','Estimated Position')
   xlabel('Time (s)')
    ylabel('Position (m)')
    title('Target position')
    hold off
%      plot(t,Ztrue,'--b',t,Monte_estimateX,'r*',t,Monte_estimateX,'g+')
grid on
% subplot(212)
figure
      plot(t,V,'r','LineWidth',1.5)
      hold on
  plot(t,Monte_estimateVx,'b','LineWidth',1.5)
 legend('Target Velcoity' , 'Estimated Velocity')
   xlabel('Time (s)')
    ylabel('Velocity (m/s)')
    title('Target Velocity')
    grid on
    hold off
    figure
    subplot(211)
    plot(t, RMSE(1,:),'LineWidth',1.5)
    
    xlabel('Time (s)')
    ylabel('RMSE (m)')
    title('RMSE in Position')
    grid on
    
    subplot(212)
    plot(t, RMSE(2,:),'LineWidth',1.5)
    xlabel('Time (s)')
    ylabel('RMSE (m/s)')
    title('RMSE in Velocity')
    grid on