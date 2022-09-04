% % K Plot
% N = 100;
% samples = 80;
% y0 = unifrnd(0,2*pi,[N,samples]);
% 
% 
% A = binornd(1,0.3,N,N);
%     for i = 1:N
%         A(i,i) = 0;
%     end
% 
%     for i = 1:N
%         for j = 1:N
%             A(j,i) = A(i,j);
%         end
%     end
% 
% T = 5;
% n = 51;
% tspan = linspace(0, T, n);
% dt = T/(n-1);
% omega = normrnd(1,0.5,[N,1]);
% % K = 0.1;
% 
% K01 = Autoinf(A,N,T,n,omega,y0,0.1,samples,0.86);
% K5 = Autoinf(A,N,T,n,omega,y0,5,samples,0.86);
% K10 = Autoinf(A,N,T,n,omega,y0,10,samples,0.86);
% 
% % Plot
% x = 1:samples;
% plot(x,K01,'b-',x,K5,'g-', x,K10,'m-');
% xlabel('Number of Perturbations');
% ylabel('Accuracy');
% legend({'k = 0.1','k = 5','k = 10'},'Location','northeast','Orientation','vertical');






% K plot with error bar
N = 15;
samples = 30;
experiments = 10;
y0 = unifrnd(0,2*pi,[N,samples]);
omega = normrnd(0,0.5,[N,experiments]);

A = binornd(1,0.2,N,N);
for i = 1:N
    A(i,i) = 0;
end

for i = 1:N
    for j = 1:N
        A(j,i) = A(i,j);
    end
end

T = 5;
n = 51;
tspan = linspace(0, T, n);
dt = T/(n-1);
K = 80;
threshold = 0.86;

% for recording accuracy
global y_cel_global
y_cel_global = cell(samples,experiments);
V = zeros([experiments,samples]);
for exp = 1:experiments
    V(exp,:) = Autoinf(exp,A,N,T,n,omega(:,exp),y0,K,samples,threshold);
end

% get mean and se
mean = sum(V)/experiments;
se = std(V)/sqrt(experiments);


% for estimating r
s = 0;
for exp = 1:experiments
    for sample = 1:samples
        s = s + r(y_cel_global{sample,exp},n,N);
    end
end
r_ave = s/(samples*experiments);


% % plot the generated data
% plot(0:0.1:5,y_cel_global{1,1});
% xlabel('Time');
% ylabel('Angular displacement \theta');

mean_K01 = mean;
se_K01 = se;
r_K01 = r_ave;

mean_K5 = mean;
se_K5 = se;
r_K5 = r_ave;

mean_K10 = mean;
se_K10 = se;
r_K10 = r_ave;

mean_K20 = mean;
se_K20 = se;
r_K20 = r_ave;

mean_K40 = mean;
se_K40 = se;
r_K40 = r_ave;

mean_K80 = mean;
se_K80 = se;
r_K80 = r_ave;





% plot mean with error bars
errorbar(mean_K01, se_K01); hold on;
errorbar(mean_K5, se_K5); hold on;
errorbar(mean_K10, se_K10); hold on;
errorbar(mean_K20, se_K20); hold on;
errorbar(mean_K40, se_K40); hold on;
errorbar(mean_K80, se_K80);

xlabel('Number of Perturbations');
ylabel('Accuracy');
legend({'k = 0.1, r = 0.24','k = 5, r = 0.44','k = 10, r = 0.58','k = 20, r = 0.74','k = 40, r = 0.86','k = 80, r = 0.92'},'Location','southeast','Orientation','vertical');


%% plot 2nd order transitions between r and k
rrr = zeros([1,100]);
for kkk = 1:1:100
    global y_cel_global
    y_cel_global = cell(samples,experiments);
    V = zeros([experiments,samples]);
    for exp = 1:experiments
        V(exp,:) = Autoinf(exp,A,N,T,n,omega(:,exp),y0,kkk,samples,threshold);
    end
    s = 0;
    for exp = 1:experiments
        for sample = 1:samples
            s = s + r(y_cel_global{sample,exp},n,N);
        end
    end
    rrr(kkk) = s/(samples*experiments);
end
ngrrr = zeros([1,21]);
for kkk = -20:1:0
    global y_cel_global
    y_cel_global = cell(samples,experiments);
    V = zeros([experiments,samples]);
    for exp = 1:experiments
        V(exp,:) = Autoinf(exp,A,N,T,n,omega(:,exp),y0,kkk,samples,threshold);
    end
    s = 0;
    for exp = 1:experiments
        for sample = 1:samples
            s = s + r(y_cel_global{sample,exp},n,N);
        end
    end
    ngrrr(kkk+21) = s/(samples*experiments);
end

plot(-20:1:100,[ngrrr,rrr])
xlabel('k');
ylabel('r');

%% plot noise for a single reconstruction
K = 0.1;
lambda = 0.01;

global y_cel_global_noise
y_cel_global_noise = cell(samples,experiments);
V_noise = zeros([experiments,samples]);
for exp = 1:experiments
    V_noise(exp,:) = Autoinf_noise(lambda,exp,A,N,T,n,omega(:,exp),y0,K,samples,threshold);
end

% get mean and se
mean = sum(V_noise)/experiments;
se = std(V_noise)/sqrt(experiments);

% for estimating r
s = 0;
for exp = 1:experiments
    for sample = 1:samples
        s = s + r(y_cel_global_noise{sample,exp},n,N);
    end
end
r_ave = s/(samples*experiments);


% lambda = 0.01
mean_N001_K01=mean;
se_N001_K01=se;
r_N001_K01=r_ave;
% plot a difference between N001K01 and K01
errorbar(mean_K01, se_K01);  hold on;
errorbar(mean_N001_K01, se_N001_K01);
xlabel('Number of Perturbations');
ylabel('Accuracy');
legend({'k = 0.1, \lambda = 0','k = 0.1, \lambda = 0.01'},'Location','southeast','Orientation','vertical');

%% plot multiple noise level
% now we need 100 samples instead of 30
samples = 100;
y0 = unifrnd(0,2*pi,[N,samples]);
lambda = 1;

global y_cel_global_noise
y_cel_global_noise = cell(samples,experiments);
V_noise = zeros([experiments,samples]);
for exp = 1:experiments
    V_noise(exp,:) = Autoinf_noise(lambda,exp,A,N,T,n,omega(:,exp),y0,K,samples,threshold);
end

% get mean and se
mean = sum(V_noise)/experiments;
se = std(V_noise)/sqrt(experiments);

% for estimating r
s = 0;
for exp = 1:experiments
    for sample = 1:samples
        s = s + r(y_cel_global_noise{sample,exp},n,N);
    end
end
r_ave = s/(samples*experiments);


% lambda = 0.01
mean_N001_K01=mean;
se_N001_K01=se;
r_N001_K01=r_ave;

% lambda = 0.05
mean_N005_K01=mean;
se_N005_K01=se;
r_N005_K01=r_ave;

% lambda = 0.1
mean_N01_K01=mean;
se_N01_K01=se;
r_N01_K01=r_ave;

% lambda = 0.5
mean_N05_K01=mean;
se_N05_K01=se;
r_N05_K01=r_ave;

% lambda = 1
mean_N1_K01=mean;
se_N1_K01=se;
r_N1_K01=r_ave;

% plot out
errorbar(mean_N001_K01, se_N001_K01); hold on;
errorbar(mean_N005_K01, se_N005_K01); hold on;
errorbar(mean_N01_K01, se_N01_K01); hold on;
errorbar(mean_N05_K01, se_N05_K01); hold on;
errorbar(mean_N1_K01, se_N1_K01); 
xlabel('Number of Perturbations');
ylabel('Accuracy');
legend({'\lambda = 0.01','\lambda = 0.05','\lambda = 0.1','\lambda = 0.5','\lambda = 1'},'Location','east','Orientation','vertical');


%% 100-network case
% K plot
N = 50;
samples = 30;
experiments = 10;
T = 5;
n = 51;
tspan = linspace(0, T, n);
dt = T/(n-1);
threshold = 0.86;

K = 1;

% for recording accuracy
global y_cel_global
y_cel_global = cell(samples,experiments);
V = zeros([experiments,samples]);
for exp = 1:experiments
    y0 = unifrnd(0,2*pi,[N,samples]);
    omega = normrnd(0,0.5,[N,experiments]);
    A = binornd(1,0.05,N,N);
    for i = 1:N
        A(i,i) = 0;
    end

    for i = 1:N
        for j = 1:N
            A(j,i) = A(i,j);
        end
    end
    V(exp,:) = Autoinf(exp,A,N,T,n,omega(:,exp),y0,K,samples,threshold);
end

% get mean and se
mean = sum(V)/experiments;
se = std(V)/sqrt(experiments);


% for estimating r
s = 0;
for exp = 1:experiments
    for sample = 1:samples
        s = s + r(y_cel_global{sample,exp},n,N);
    end
end
r_ave = s/(samples*experiments);

% store the value

% 
% mean_K50 = mean;
% se_K50 = se;
% r_K50 = r_ave;

% mean_K100 = mean;
% se_K100 = se;
% r_K100 = r_ave;


% mean_K150 = mean;
% se_K150 = se;
% r_K150 = r_ave;

% mean_K200 = mean;
% se_K200 = se;
% r_K200 = r_ave;

% mean_K400 = mean;
% se_K400 = se;
% r_K400 = r_ave;

% mean_K600 = mean;
% se_K600 = se;
% r_K600 = r_ave;

% mean_K1000 = mean;
% se_K1000 = se;
% r_K1000 = r_ave;

% now plot the error bars & marking their K values
errorbar(mean_K1, se_K1); hold on;
% errorbar(mean_K10, se_K10); hold on;
% errorbar(mean_K20, se_K20); hold on;
errorbar(mean_K50, se_K50); hold on;
errorbar(mean_K100, se_K100); hold on;
% errorbar(mean_K300, se_K300); hold on;
% errorbar(mean_K150, se_K150); hold on;
% errorbar(mean_K400, se_K400); hold on;
% errorbar(mean_K600, se_K600); hold on;
% errorbar(mean_K1000, se_K1000); hold on;
errorbar(mean_K2000, se_K2000); hold on;


xlabel('Number of Perturbations');
ylabel('Accuracy');
legend({'k = 1, r = 0.1315','k = 50, r = 0.4059','k = 100, r = 0.5126','k = 2000, r = 0.8558'},'Location','east','Orientation','vertical');
%% 100-network noise plot

K=10;

lambda = 1;

% global y_cel_global_noise
% y_cel_global_noise = cell(samples,experiments);
V_noise = zeros([experiments,samples]);
for exp = 1:experiments
    y0 = unifrnd(0,2*pi,[N,samples]);
    omega = normrnd(0,0.5,[N,experiments]);
    A = binornd(1,0.05,N,N);
    V_noise(exp,:) = Autoinf_noise(lambda,exp,A,N,T,n,omega(:,exp),y0,K,samples,threshold);
end

% get mean and se
mean = sum(V_noise)/experiments;
se = std(V_noise)/sqrt(experiments);

% %for estimating r
% s = 0;
% for exp = 1:experiments
%     for sample = 1:samples
%         s = s + r(y_cel_global_noise{sample,exp},n,N);
%     end
% end
% r_ave = s/(samples*experiments);


% lambda = 0.01
mean_N1=mean;
se_N1=se;
%r_N001_K20=r_ave;

% % lambda = 0.05
% mean_N005_K20=mean;
% se_N005_K20=se;
% %r_N005_K20=r_ave;

% % lambda = 0.1
% mean_N01_K20=mean;
% se_N01_K20=se;
% %r_N01_K20=r_ave;

% % lambda = 0.5
% mean_N05_K20=mean;
% se_N05_K20=se;
% %r_N05_K20=r_ave;

% % lambda = 1
% mean_N1_K20=mean;
% se_N1_K20=se;
% %r_N1_K20=r_ave;

% % lambda = 2
% mean_N2_K20=mean;
% se_N2_K20=se;
% %r_N2_K20=r_ave;

% % lambda = 5
% mean_N5_K20=mean;
% se_N5_K20=se;
%r_N5_K20=r_ave;

% plot out
errorbar(mean_N001, se_N001); hold on;
errorbar(mean_N005, se_N005); hold on;
errorbar(mean_N01, se_N01); hold on;
% errorbar(mean_N05_K20, se_N05_K20); hold on;
errorbar(mean_N1, se_N1); hold on;
% errorbar(mean_N2_K20, se_N2_K20); hold on;
% errorbar(mean_N5_K20, se_N5_K20);
xlabel('Number of Perturbations');
ylabel('Accuracy');
legend({'\lambda = 0.01','\lambda = 0.05','\lambda = 0.1','\lambda = 1'},'Location','northwest','Orientation','vertical');

