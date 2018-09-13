clear all
delta=0;% Dividend Yield
T=1;% Time to maturity
sigma=0.2;% Volatility
t=[0;1/3;2/3;1]; % Exercise opportunities for the Bermudan option
r=0.06;% Risk-free rate
reiteration=25;% Number of estimates
b=1000;% Branching parameter
S0=[70;80;90;100;110;120]; % Initial stock price
K0=[100;100;100;100;100;100];% Strike Price
K=K0./K0;
S=S0./K0;
truevalue=[0.121;0.670;2.303;5.731;11.341;20.000];
o=length(S0);
% variable generation
option=zeros(reiteration,o);
estimator=zeros(1,o);
serr=zeros(1,o);
realerror=zeros(1,o);
RMSE=zeros(1,o);
time=zeros(1,o);
percitm=zeros(1,o);
inthemoneyperc=zeros(1,o);
for j=1:o
    tic
    for f=1:reiteration
        %Generate the stock paths
        S1=zeros(b,1);
        S2=zeros(b,1);
        for i=1:b/2
            x=normrnd(0,1,1,1);
            S1(i)=S(j)*exp((r-delta-sigma^2/2)*(1-2/3)+sigma*sqrt(1-2/3)*x);
            S1(i+b/2)=S(j)*exp((r-delta-sigma^2/2)*(1-2/3)+sigma...
                *sqrt(1-2/3)*(-x));
            x=normrnd(0,1,1,1);
            S2(i)=S1(i)*exp((r-delta-sigma^2/2)*(1-2/3)+sigma*sqrt(1-2/3)*x);
            S2(i+b/2)=S1(i+b/2)*exp((r-delta-sigma^2/2)*(1-2/3)+sigma...
                *sqrt(1-2/3)*(-x));
        end
        % Calculate the option price
        % At T-1 the calculations are straightforward thanks to pruning
        Yt2=blsprice(S2,K(j),r,1/3,sigma,delta);
        Lt2=max(max(S2-K(j),0),blsprice(S2,K(j),r,1/3,sigma,delta));
        %At time T-2 we need first to see where the option is in the money
        indicator=zeros(b,1);
        inthemoney=max(S1-K(j),0);
        for i=1:b
            if inthemoney(i)>0
                indicator(i)=1;
            end
        end
        %Generate the basic functions
        X2=zeros(b,4);
        Y=Lt2*exp(-r*1/3);%F
        Y1=Yt2*exp(-r*1/3); %Y
        Y12=(Yt2.^2)*exp(-r*1/3);%Y^2
        Y13=(Lt2.*Yt2)*exp(-r*1/3);%FY
        S22=S1;
        for i=1:b
            if indicator(i)==1
                X2(i,1)=K(j);
                X2(i,2)=S22(i);
                X2(i,3)=blsprice(S22(i),K(j),r,2/3,sigma,delta);
                X2(i,4)=S22(i)*blsprice(S22(i),K(j),r,2/3,sigma,delta);
            end
        end
        
        %These lines delete the zero rows from the regression matrices
        X3 = X2(any(X2,2),:);
        Y2= Y(any(X2,2),:);%F
        Y21= Y1(any(X2,2),:);%Y
        Y22= Y12(any(X2,2),:);%Y^2
        Y23= Y13(any(X2,2),:);%FY
        beta=X3'*X3\X3'*Y2; %F
        beta1=X3'*X3\X3'*Y21;%Y
        beta2=X3'*X3\X3'*Y22;%Y^2
        beta3=X3'*X3\X3'*Y23;%FY
        Continuation=X3*beta;%F
        Continuation1=X3*beta1;%Y
        Continuation2=X3*beta2;%Y^2
        Continuation3=X3*beta3;%FY
        coeff=-(Continuation3-Continuation.*Continuation1)./...
            (Continuation2-(Continuation1.^2));
        % control variate adjusted continuation value
        Contval=Continuation+coeff.*(Continuation1-X3(:,3));
        continuation1=zeros(b,1);
        S21=S1;
        S21(indicator==0)=[];
        for i=1:length(Y2)
            n=find(S1==S21(i));
            continuation1(n)=Contval(i);
        end
        Lt1=Lt2*exp(-r*1/3);
        Yt1=Yt2*exp(-r*1/3);
        for i=1:b
            if max(S1(i)-K(j),0)>0
                if max(S1(i)-K(j),0)>continuation1(i)
                    Lt1(i,1)=max(S1(i)-K(j),0);
                    Yt1(i,1)=blsprice(S1(i),K(j),r,2/3,sigma,delta);
                end
            end
        end
        Lt0=Lt1*exp(-r*1/3);
        Yt0=Yt1*exp(-r*1/3);
        est=mean(Lt0);
        Y0=mean(Yt0);
        LY=mean(Lt0.*Yt0);
        Ysquared=mean(Yt0.^2);
        coefficient=-(LY-est*Y0)/(Ysquared-(Y0^2));
        % control variate adjusted continuation value
        estimator1=est+coefficient*(Y0-blsprice(S(j),K(j),r,1,sigma,delta));
        % Point estimate
        option(f,j)=max(max(S(j)-K(j),0),estimator1)*100;
        % n. of In-the-money paths
        inthemoneyperc(f,j)=sum(indicator);
    end
    % Final Estimator
    estimator(1,j)= mean(option(:,j));
    % SE Estimator
    serr(1,j)=std(option(:,j))/sqrt(reiteration);
    realerror(1,j)=abs(estimator(1,j)-truevalue(j))/truevalue(j);
    % RMSE
    RMSE(1,j)=sqrt(mean((option(:,j)-truevalue(j)).^2));
    % In-the-money paths as a percentage of the branching parameter
    percitm(j)=mean(inthemoneyperc(:,j))/b;
    % cpu time
    time(j)=toc;
end
%Building the table
dataLength = length(S0);
%Building the table
%These are lines used to space between an output and the preceding one
disp(' ');
disp(' ');
disp(' ');
name = 'Branching Parameter';
str1 = [name,' ', num2str(b)];
disp(str1);
name = 'Number of Estimates';
str = [name,' ', num2str(reiteration)];
disp(str);
disp(' ');
% This code is used to generate the table of results.
label = char('70','80','90','100','110','120');
fprintf(['Stock Price Estimator SE True Value',...
    ' Rel Error RMSE Average time per Estimate(in sec)\n']);
for i=1:dataLength
    fprintf(['%8s %-8.4g %-8.3g %-8.4g',...
        ' %-8.3g %-8.3g %-8.4g\n'],...
        label(i,:),estimator(i),serr(i),truevalue(i),realerror(i),RMSE(i),time(i))
end
disp(' ');
% Generate the output for the total cpu time computation
name = 'Time Elapsed for constructing the entire table (in sec)';
str1 = [name,' ', num2str(sum(time))];
disp(str1)