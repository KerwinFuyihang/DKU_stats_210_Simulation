
clear,clc,close all,format compact
%% Question1
Ntrials = [100,1000,10000,10000];
PIs = zeros(length(Ntrials),1);
for i = 1:length(Ntrials)
    N = Ntrials(i);%No.of trials
    x = rand(N,1)*pi/2;
    y = rand(N,1)*1/2;
    p = mean(y<(0.5*sin(x)));
end

figure()
semilogx(Ntrials,PIs,'LineWidth',2)
hold on
plot(Ntrials,Ntrials*0+pi,'--','LineWidth',2)
xlabel('Trials')
title('\pi estimation by buffon needle simulation')

%% Question 2
% Equation to integrate over. Integral of this from a to b.
f = @(x) 1./x+1;

% Parameters
N = 10000;  % number of points to generate randomly
a = 2;      % integral lower x limit - manual
b = 4;      % integral upper x limit - manual
M = 1.4*max(f(linspace(a,b)));     % upper y bound = c*max y_value from equation

% Generate Dots
for i = 1:N
   
    % generate random pt
    x_val = rand(1)*(b-a) + a;
    y_val = rand(1)*M;
    
    % compare random against the curve
    fx = f(x_val);
    
    % logic statement
    if y_val < fx           % under the curve - blue
        under(i,1) = x_val;
        under(i,2) = y_val;
    else                    % above the curve - red
        above(i,1) = x_val;
        above(i,2) = y_val;
    end
end

% Filter out zeros
under2(:,1) = nonzeros(under(:,1));
under2(:,2) = nonzeros(under(:,2));
above2(:,1) = nonzeros(above(:,1));
above2(:,2) = nonzeros(above(:,2));

% Plotting
plot(above2(:,1),above2(:,2),'ro','MarkerFaceColor','r')
hold on
plot(under2(:,1),under2(:,2),'bo','MarkerFaceColor','b')
title('Monte Carlo Integration'), xlabel('x'), ylabel('y')
legend('above','under')

% Integral Calcs
MonteCarlo_Integral = length(under2) / N * (M*(b-a));
MATLAB_Integral = integral(f,a,b);
PercentError = abs(MATLAB_Integral - MonteCarlo_Integral)/MATLAB_Integral * 100;

%% Question 3


