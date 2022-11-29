
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
f = @(x) 1./(x+1);

% Parameters
N = 10000;  % number of points to generate randomly
a = 0;      % integral lower x limit - manual
b = 1;      % integral upper x limit - manual
M = max(f(linspace(a,b)));     % upper y bound
% Generate Dots
for i = 1:N
   
    % generate random pt
    x_val = rand(1);
    y_val = rand(1);
    
    % compare random against the curve
    fx = f(x_val);
    
    % logiccl statement
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
plot(above2(:,1),above2(:,1),'ro','MarkerFaceColor','r')
hold on
plot(under2(:,1),under2(:,1),'bo','MarkerFaceColor','b')
title('Monte Carlo Integration'), xlabel('x'), ylabel('y')
legend('above','under')

% Integral Calcs
MonteCarlo_Integral = length(under2) / N * (b-a);
MATLAB_Integral = integral(f,a,b);
PercentError = abs(MATLAB_Integral - MonteCarlo_Integral)/MATLAB_Integral * 100;

%% Question 3
% Real calculation
p =0.25;
X =4;
p_theory = (1-p)^(X-1);

% Simulation

N = 100000;
number = rand(N,X-1) <p;
count = sum(max(number,[],2));
p_sim = 1-count/N;
fprintf ('Calcluated P(X>=4) is %.4f, estimated P(x>=4) is %.4f\n', p_theory,p_sim )

%% Question 4

p = 0.5;
N = 10000;
nMax = 100;
[~, X] = max(rand(N, nMax) < p, [], 2); % find first success position in each row
[~, Y] = max(rand(N, nMax) < p, [], 2); % find first success position in each row
Z = X + Y;

k = 2:max(Z);
pmf = p^2*(k-1).*(1-p).^(k-2);

counts = pmf*0;
for i = 1:length(k)
    counts(i) = sum(Z == k(i));
end
pmf_est = counts/N;

figure()
bar(k, pmf_est)
hold on
plot(k, pmf, 'r-', 'Linewidth', 2)
xlabel('Z')
legend(' PMF of Z', 'True PMF')
%% Question 5
close all
r = normrnd(1,1,[1,500]);
histogram(r);
histfit(r);
mean_r = mean(r);





