function [Spots, statSpots] = OrnsteinUhlenbeck(nDays, nSims, mu, s0, vol, theta)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$
    %
    % Ornstein Uhlenbeck process:
    %
    % $dX_t = \theta (\mu - X_t)dt + \sigma dW_t$
    %
    % Discretise:
    %
    % $x_t-x_(t-1) = \theta (\mu - X_t-1)\Delta t + \sigma \epsilon_i\sqrt(t)$
    % 
    % valid for small t where 
    %
    % $\epsilon_i \epsilon N(0,1)$ .
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inputs
    %        nDays - number of days or steps to simulate
    %        nSims - number of simulations
    %        Seed - seed for RNG
    %        mu - parameter in OU
    %        s0 - starting point
    %        vol - parameter in OU
    %        theta - parameter in OU
    %
    % Outputs:
    %        Spots - simulated values
    %        statSpots - stats - mean and std
    %
    % Ahmos Sansom - August 2014
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Seed for RNG
    Seed    = 12345;
    % set dt to be unity - could add as input...
    dt = 1.0;
 
    % initialise random number generator
    s = RandStream('twister','Seed',Seed);
    RandStream.setGlobalStream(s)
    
    % Initialise spot sims and stats
    Spots = zeros(nDays,nSims);
    statSpots = zeros(nDays,4);
    
    % pre set vol*sqrt(dt) vol to reduce run time
    volsqrdt = vol*sqrt(dt);
    
    % Loop round days
    for j=1:nDays
        
        % Store some randns:
        myRandn = randn(nSims);
        
        for k=1:nSims % Each days has nSims
            
            if j == 1 % Initial Step
               Spots(j,k) = s0 + theta*(mu - s0)*dt + myRandn(k)*volsqrdt;
            else
               Spots(j,k) = Spots(j-1,k) + theta*(mu - Spots(j-1,k))*dt + myRandn(k)*volsqrdt;
            end
        end
        
        % Collect stats:
        statSpots(j,1) = mean(Spots(j,:));
        statSpots(j,2) = var(Spots(j,:));
        
        statSpots(j,3) = s0 * exp(-theta*j*dt) + mu *(1-exp(-theta*j*dt));
        statSpots(j,4) = 0.5 * vol * vol * (1.0 - exp(-2.0*theta*j*dt))/ theta;
           
    end
        
end
