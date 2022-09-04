%% general parameter
N = 50;
T = 5;
n = 51;
tspan = linspace(0, T, n);
samples = 30;


%% scale free network generation
N_seed = 5;
seed = binornd(1,0.5,N_seed,N_seed);
for i = 1:N_seed
    seed(i,i) = 0;
end

for i = 1:N_seed
    for j = 1:N_seed
        seed(j,i) = seed(i,j);
    end
end
SF = double(SFNG(N, 1, seed));
% issymmetric(SF)
% sum(sum(SF))

omega_SF = normrnd(0,0.5,[N,1]);
y0_SF = unifrnd(0,2*pi,[N,samples]);

%% random graph generation
RG = binornd(1,0.2,N,N);
for i = 1:N
    RG(i,i) = 0;
end

for i = 1:N
    for j = 1:N
        RG(j,i) = RG(i,j);
    end
end
% issymmetric(RG)
% sum(sum(RG))

omega_RG = omega_SF;
y0_RG = unifrnd(0,2*pi,[N,samples]);

%% K values
K_RG = zeros([samples,1]);
K_SF = zeros([samples,1]);
for i = 1:samples
    K_RG(i) = 200;
    K_SF(i) = 200;
end
r_RG = zeros([samples,1]);
r_SF = zeros([samples,1]);
for i = 1:samples
[t,y_RG] = ode45(@(t,theta)ode(t,theta,omega_RG,K_RG(i),N,RG),tspan,y0_RG(:,i));
[t,y_SF] = ode45(@(t,theta)ode(t,theta,omega_SF,K_SF(i),N,SF),tspan,y0_SF(:,i));
r_RG(i) = r(y_RG,n,N);
r_SF(i) = r(y_SF,n,N);
end

% when K_SF(1)=200, r_SF(1)=xx; then wanna find K_RG(1) s.t. r_RG(1)=r_SF(1)
% do this for every i

% set same k=200
% if r_RG(i)<r_SF(i), while <: K_RG=K_RG+1 until r_RG(i)>r_SF(i)
% else: K_RG=K_RG-1 until r_RG(i)<r_SF(i)

% Then since we exceed it by less than 1
% while abs(det_r)~=0: if det_r>0, K_RG=K_RG-0.1;if det_r<0, K_RG=K_RG+0.1

% by observation r_RG always seems lower than r_SF by our setting
% so we just need to increase K_RG (no need for decrease for now)
for i = 1:samples

        while r_RG(i)>r_SF(i)
            K_RG(i)=K_RG(i)-5;
            [t,y_RG] = ode45(@(t,theta)ode(t,theta,omega_RG,K_RG(i),N,RG),tspan,y0_RG(:,i));
            r_RG(i)=r(y_RG,n,N);
        end

        while r_RG(i)<r_SF(i)
            K_RG(i)=K_RG(i)+0.1;
            [t,y_RG] = ode45(@(t,theta)ode(t,theta,omega_RG,K_RG(i),N,RG),tspan,y0_RG(:,i));
            r_RG(i)=r(y_RG,n,N);
        end
end

% max difference on a scale of 10^-3
max(abs(r_RG-r_SF))

%% recon
w_RG = Autoinf_2(RG,N,T,n,omega_RG,y0_RG,K_RG,samples,0.86);
w_SF = Autoinf_2(SF,N,T,n,omega_SF,y0_SF,K_SF,samples,0.86);



%% multiple exp
N = 100;
N_seed = 5;
T = 5;
n = 51;
tspan = linspace(0, T, n);
samples = 50;
experiments = 5;
omega = normrnd(0,0.5,[N,experiments]);
% y0_RG = unifrnd(0,2*pi,[N,samples,experiments]);
y0_RG = cell([1,experiments]); % inside each cell it's an N*samples matrix
% y0_SF = unifrnd(0,2*pi,[N,samples,experiments]);
y0_SF = cell([1,experiments]); % inside each cell it's an N*samples matrix
for exp = 1:experiments
    y0_RG{exp} = unifrnd(0,2*pi,[N,samples]);
    y0_SF{exp} = unifrnd(0,2*pi,[N,samples]);
end
K_RG = zeros([samples,experiments]);
K_SF = zeros([samples,experiments]);
for i = 1:samples
    for j = 1:experiments
        K_RG(i,j) = 500;
        K_SF(i,j) = 500;
    end
end
r_RG = zeros([samples,experiments]);
r_SF = zeros([samples,experiments]);


s = 0;
for exp = 1:200
    seed = binornd(1,0.5,N_seed,N_seed);
    for i = 1:N_seed
        seed(i,i) = 0;
    end
    for i = 1:N_seed
        for j = 1:N_seed
            seed(j,i) = seed(i,j);
        end
    end
    SF = double(SFNG(N, 3, seed));
    s = s + sum(sum(SF));
end
ave = s/(N*200);


for exp=1:experiments
    % generate SF
    seed = binornd(1,0.5,N_seed,N_seed);
    for i = 1:N_seed
        seed(i,i) = 0;
    end
    for i = 1:N_seed
        for j = 1:N_seed
            seed(j,i) = seed(i,j);
        end
    end
    SF = double(SFNG(N, 3, seed));

    % generate RG
    RG = binornd(1,ave/N,N,N);
    for i = 1:N
        RG(i,i) = 0;
    end
    for i = 1:N
        for j = 1:N
            RG(j,i) = RG(i,j);
        end
    end

    for sample = 1:samples
        [t,y_RG] = ode45(@(t,theta)ode(t,theta,omega(:,exp),K_RG(sample,exp),N,RG),tspan,y0_RG{exp}(:,sample));
        [t,y_SF] = ode45(@(t,theta)ode(t,theta,omega(:,exp),K_SF(sample,exp),N,SF),tspan,y0_SF{exp}(:,sample));
        r_RG(sample,exp) = r(y_RG,n,N);
        r_SF(sample,exp) = r(y_SF,n,N);
    end
end



tic
for exp = 1:experiments
    % matching r
    for i = 1:samples
        if r_RG(i,exp)>r_SF(i,exp)
            while r_RG(i,exp)>r_SF(i,exp)
                K_RG(i,exp)=K_RG(i,exp)-10;
                [t,y_RG] = ode45(@(t,theta)ode(t,theta,omega(:,exp),K_RG(i,exp),N,RG),tspan,y0_RG{exp}(:,i));
                r_RG(i,exp)=r(y_RG,n,N);
            end
        else
            while r_RG(i,exp)<r_SF(i,exp)
                K_RG(i,exp)=K_RG(i,exp)+10;
                [t,y_RG] = ode45(@(t,theta)ode(t,theta,omega(:,exp),K_RG(i,exp),N,RG),tspan,y0_RG{exp}(:,i));
                r_RG(i,exp)=r(y_RG,n,N);
            end
        end


        if r_RG(i,exp)<r_SF(i,exp)
            while r_RG(i,exp)<r_SF(i,exp)
                K_RG(i,exp)=K_RG(i,exp)+0.1;
                [t,y_RG] = ode45(@(t,theta)ode(t,theta,omega(:,exp),K_RG(i,exp),N,RG),tspan,y0_RG{exp}(:,i));
                r_RG(i,exp)=r(y_RG,n,N);
            end
        else
            while r_RG(i,exp)>r_SF(i,exp)
                K_RG(i,exp)=K_RG(i,exp)-0.1;
                [t,y_RG] = ode45(@(t,theta)ode(t,theta,omega(:,exp),K_RG(i,exp),N,RG),tspan,y0_RG{exp}(:,i));
                r_RG(i,exp)=r(y_RG,n,N);
            end

        end
    end
end
toc

% report the max difference of r as a result of the matching
r_dif = abs(r_RG-r_SF); % it is supposed to stay below 0.1
max(max(r_dif))
mean(mean(r_SF))
mean(mean(r_RG))


w_RG = zeros([experiments,samples]);
w_SF = zeros([experiments,samples]);
tic
for exp = 1:experiments
    w_RG(exp,:) = Autoinf_2(RG,N,T,n,omega(:,exp),y0_RG{exp}(:,:),K_RG(:,exp),samples,0.86);
    w_SF(exp,:) = Autoinf_2(SF,N,T,n,omega(:,exp),y0_SF{exp}(:,:),K_SF(:,exp),samples,0.86);
end
toc

mean_RG = sum(w_RG)/experiments;
se_RG = std(w_RG)/sqrt(experiments);
mean_SF = sum(w_SF)/experiments;
se_SF = std(w_SF)/sqrt(experiments);



meanRG_r095 = mean_RG;
seRG_r095 = se_RG;
meanSF_r095 = mean_SF;
seSF_r095 = se_SF;

layout = tiledlayout(3,1);
nexttile
errorbar(meanSF_r027,seSF_r027); hold on;
errorbar(meanRG_r027,seRG_r027); hold on;
xlabel('Perturbations','FontSize', 12);
ylabel("Accuracy",'FontSize', 12);
legend('SF,r=0.27','RG,r=0.27',Orientation='vertical',Location='east')
nexttile
errorbar(meanSF_r062,seSF_r062); hold on;
errorbar(meanRG_r062,seRG_r062); hold on;
xlabel('Perturbations','FontSize', 12);
ylabel("Accuracy",'FontSize', 12);
legend('SF,r=0.62','RG,r=0.62',Orientation='vertical',Location='east')
% nexttile
% errorbar(meanSF_r089,seSF_r089); hold on;
% errorbar(meanRG_r089,seRG_r089); hold on;
nexttile
errorbar(meanSF_r095,seSF_r095); hold on;
errorbar(meanRG_r095,seRG_r095); hold on;
xlabel('Perturbations','FontSize', 12);
ylabel("Accuracy",'FontSize', 12);
legend('SF,r=0.95','RG,r=0.95',Orientation='vertical',Location='east')
% nexttile
% errorbar(meanSF_r078,seSF_r078); hold on;
% errorbar(meanRG_r078,seRG_r078); hold on;
% 
% nexttile
% errorbar(meanSF_r096,seSF_r096); hold on;
% errorbar(meanRG_r096,seRG_r096); hold on;

% spacing
layout.TileSpacing = 'compact';
layout.Padding = 'compact';

% errorbar(meanSF_r031,seSF_r031); hold on;
% errorbar(meanRG_r031,seRG_r031); hold on;
% errorbar(meanSF_r045,seSF_r045); hold on;
% errorbar(meanRG_r045,seRG_r045); hold on;
% errorbar(mean(w_SF),se_SF); hold on;
% errorbar(mean(w_SF),se_SF); hold on;
% errorbar(mean(w_SF),se_SF); hold on;
% errorbar(mean(w_SF),se_SF); 
xlabel('Perturbations','FontSize', 12);
ylabel("Accuracy",'FontSize', 12);
legend('SF,r\approx0.31','RG,r\approx0.31',Orientation='vertical',Location='east')



