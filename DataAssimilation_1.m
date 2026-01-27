%obtaining the data
rng(1);
A = readtable('meanlogprice.csv');
y1 = A{:,2}; %extracting y1 for n=1,..,30
y2 = A{:,3}; %extracting y2 for n=1,..,30
Nsamples = numel(y1); %number of observations

R = [1 0.5; 0.5 1]; %correlation matrix
Rinverse = (1/0.75)*[1 -0.5; -0.5 1]; %inverse of R

%gamma(0.01,0.01) prior on γ
alpha = 0.01;
beta = 0.01;

%initializations for Gibbs sampling
B = 5000; % burn-in period
MHsamples = 10000; %number of Metropolis-Hastings samples
thinning = 1; %thinning
B1 = B + MHsamples*thinning;

gammavec = zeros(1, B1);
muvec    = zeros(1, B1);

%initializing the given priors
gammavec(1) = gamrnd(alpha, 1/beta);
muvec(1) = normrnd(0,sqrt(100/gammavec(1)));

for i = 2:B1
    %updating γ|mu,L
    mu_prev = muvec(i-1);
    %quadratic form of the sum of matrices
    z = [y1 - mu_prev,  y2 - 2];
    quad_sum = sum(sum((z*Rinverse) .* z, 2));
    gamma_shape = Nsamples + 0.51; %shape parameter of conditional posterior Q5 
    gamma_rate  = 0.5 * quad_sum + (mu_prev^2)/200 + 0.01;  %rate parameter of conditional posterior
    gammavec(i) = gamrnd(gamma_shape, 1/gamma_rate); 
    %updating mu|γ,L
    sum_x = sum(y1-0.5*(y2-2)); %where x = y1-0.5*(y2-2)
    post_mean = (sum_x/0.75)/((Nsamples/0.75)+(1/100)); %posterior mean derived in Q4
    post_var = 1/(gammavec(i)*(Nsamples/0.75)+(1/100)); %posterior variance
    muvec(i) = normrnd(post_mean,sqrt(post_var)); %draws mu(i)
end

keep_index = (B+1):thinning:B1;
mu = muvec(keep_index).'; %column vectors
gamma = gammavec(keep_index).';

%numerical summaries
summary = @(x) [mean(x), std(x), quantile(x,0.025), quantile(x,0.975)];
mu_stats = summary(mu);
gamma_stats = summary(gamma);

num_summary = table(["mu"; "gamma"], [mu_stats(1); gamma_stats(1)], [mu_stats(2); gamma_stats(2)], ...
    [mu_stats(3);gamma_stats(3)], [mu_stats(4);    gamma_stats(4)], 'VariableNames', {'Parameter','PosteriorMean','PosteriorSD','CI-2p5','CI-97p5'});

disp('Numerical summaries:');
disp(num_summary);

%graphical evidence that the Gibbs chain has mixed satisfactorily
%mu trace plot
figure;
plot(mu); grid on
xlabel('Number of iterations'); ylabel('\mu'); title('Trace Plot for \mu');

%gamma trace plot
figure;
plot(gamma); grid on
xlabel('Number of iterations'); ylabel('\gamma'); title('Trace Plot for \gamma');

%posterior histogram of mu
figure;
histogram(mu); grid on
xlabel('\mu'); ylabel('frequency'); title('Posterior Histogram of \mu');

%posterior histogram of gamma
figure;
histogram(gamma); grid on
xlabel('\gamma'); ylabel('frequency'); title('Posterior Histogram of \gamma');
