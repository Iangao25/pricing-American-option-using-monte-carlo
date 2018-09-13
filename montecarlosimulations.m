dt = 1/50; %Lenght of the interval of time
N = T/dt;  %Number of periods to simulate the price
NSim = 100000;  %Number of simulations

dBt = sqrt(dt)*randn(NSim,N); % Brownian motion
St = zeros(NSim,N); % Initialize matrix
St(:,1) = S0*ones(NSim,1); % vector of initial stock price per simulation

for t = 2:N;
    St(:,t) = St(:,t-1).*exp((r-D-.5*sigma^2)*dt+sigma*dBt(:,t)); %simulation of prices
end

SSit = St; % just change the name
NSim = size(SSit,1); % Number of simulations