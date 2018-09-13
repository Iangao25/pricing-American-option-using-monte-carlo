%stock
sigma = .2;
S0 = 100;

%interest rate and dividend yield
r = .04;
D = 0;

%option
T = 2; %time to maturity
KP = 100; % Strike Price





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






% Work Backwards
% Initialize CashFlow Matrix
MM = NaN*ones(NSim,N);
MM(:,N) = max(KP-SSit(:,N),0);

for tt = N:-1:3;
    % disp('Time to Maturity')
    % disp(1-tt/N)
    % Step 1: Select the path in the money at time tt-1
    I = find(KP-SSit(:,tt-1)>0);
    ISize = length(I);
    % Step 2: Project CashFlow at time tt onto basis function at time tt-1
    if tt == N
        YY = (ones(ISize,1)*exp(-r*[1:N-tt+1]*dt)).*MM(I,tt:N);
    else
        YY = sum(((ones(ISize,1)*exp(-r*[1:N-tt+1]*dt)).*MM(I,tt:N))')';
    end
    
    SSb3 = SSit(I,tt-1);
    SSb = SSb3/S0;
    XX  = [ones(ISize,1),exp(-SSb./2),exp(-SSb./2).*(1-SSb),...
        exp(-SSb./2).*(1-SSb.*2+SSb.^2./2)];
    BB = XX'*XX\XX'*YY;
    

    % Find when the option is exercised:
    IStop = max(KP - SSb3, 0) > XX*BB;
    % Find when the option is not exercised:
    ICon = setdiff([1:NSim],I(IStop));
    % Replace the payoff function with the value of the option (zeros when
    % not exercised and values when exercised):
    MM(I(IStop),tt-1) = KP-SSit(I(IStop),tt-1);
    MM(I(IStop),tt:N) = zeros(length(I(IStop)),N-tt+1);
    MM(ICon,tt-1) = zeros(length(ICon),1);
    
end

YY = sum(((ones(NSim,1)*exp(-r*[1:N-1]*dt)).*MM(:,2:N))')';

StopTime = zeros(NSim,1);
Euro = zeros(NSim,1);
EuroPayoff = zeros(NSim,1);
MM(:,1)=[];

for i=1:NSim
    if sum(MM(i,:)) == 0
        StopTime(i) = 0;
        Euro(i) = 0;
        EuroPayoff(i) = 0;
    else
        StopTime(i) = find(MM(i,:)>0);
        StopTime(i) = StopTime(i)/N;
        d1 = (log((KP-sum(MM(i,:)))/KP)+(r+sigma^2/2)*(T-StopTime(i)))/(sigma*sqrt(T-StopTime(i)));
        d2 = (log((KP-sum(MM(i,:)))/KP)+(r-sigma^2/2)*(T-StopTime(i)))/(sigma*sqrt(T-StopTime(i)));
        Euro(i) = KP*exp(-r*(T-StopTime(i)))*normcdf(-d2)-(KP-sum(MM(i,:)))*normcdf(-d1);
        EuroPayoff(i) = exp(-r*(T-StopTime(i)))*max((KP - SSit(i,N)), 0);
    end
    
    
    
end

MatCov = cov(EuroPayoff, YY);
c = -MatCov(1,2)/var(EuroPayoff);


ValueWithControl = YY + c * (EuroPayoff - Euro);

d1 = (log(S0/KP)+(r+sigma^2/2)*T)/(sigma*sqrt(T));
d2 = (log(S0/KP)+(r-sigma^2/2)*T)/(sigma*sqrt(T));
P_bseu = KP*exp(-r*T)*normcdf(-d2)-S0*normcdf(-d1)

valuewithcontrol = mean(ValueWithControl)
Sterr = std(ValueWithControl)/sqrt(NSim)

