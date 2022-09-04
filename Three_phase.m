% generated data for different K/three phase

N = 15;
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
omega = normrnd(0,0.5,[N,1]);
y0 = unifrnd(0,2*pi,[N,1]);

%% Plot unsync
K = 0.1;

% this is our data generated
[t,y] = ode45(@(t,theta)ode(t,theta,omega,K,N,A),tspan,y0);

% unsynch
% y_k01 = mod(y,2*pi);
y_k01 = y;
r_k01 = r(y_k01,n,N);


%% Plot partly sync
K = 5;

% this is our data generated
[t,y] = ode45(@(t,theta)ode(t,theta,omega,K,N,A),tspan,y0);

% unsynch
%y_k5 = mod(y,2*pi);
y_k5 = y;
r_k5 = r(y_k5,n,N);

%% Plot full sync
K = 50;

% this is our data generated
[t,y] = ode45(@(t,theta)ode(t,theta,omega,K,N,A),tspan,y0);

% unsynch
%y_k50 = mod(y,2*pi);
y_k50 = y;
r_k50 = r(y_k50,n,N);

%% Plot out
layout = tiledlayout(3,1);
% 1
nexttile
plot(0:1:200,y_k01,'LineWidth',0.2)
xlabel('Time','FontSize', 12);
ylabel("\theta for k = 0.1, r = 0.2432",'FontSize', 12);

% 2
nexttile
plot(0:1:200,y_k5,'LineWidth',0.2)
xlabel('Time','FontSize', 12);
ylabel("\theta for k = 5, r = 0.5732",'FontSize', 12);


%3
nexttile
plot(0:1:200,y_k50,'LineWidth',0.2)
xlabel('Time','FontSize', 12);
ylabel("\theta for k = 30, r = 0.9711",'FontSize', 12);


% spacing
layout.TileSpacing = 'compact';
layout.Padding = 'compact';

%% r-k plot for 10 experiments with different networks of the same size
% N=15
yy = cell([10,100]);
N=15;

for i = 1:10
    A = binornd(1,0.2,N,N);

    for a = 1:N
        A(a,a) = 0;
    end

    for b = 1:N
        for c = 1:N
            A(b,c) = A(c,b);
        end
    end

    omega = normrnd(0,0.5,[N,1]);
    y0 = unifrnd(0,2*pi,[N,1]);


    for k = 1:100
        [t,y] = ode45(@(t,theta)ode(t,theta,omega,k,N,A),tspan,y0);
        yy{i,k} = y;
    end
end
 
yy15 = yy;
rr15 = zeros([10,100]);
for i = 1:10
    for j = 1:100
        rr15(i,j) = r(yy15{i,j},n,N);
    end
end
mean15 = sum(rr15)/10;
se15 = std(rr15)/sqrt(10);


% N = 30
N = 30;
yy = cell([10,100]);
for i = 1:10
    A = binornd(1,0.2,N,N);

    for a = 1:N
        A(a,a) = 0;
    end

    for b = 1:N
        for c = 1:N
            A(b,c) = A(c,b);
        end
    end

    omega = normrnd(0,0.5,[N,1]);
    y0 = unifrnd(0,2*pi,[N,1]);


    for k = 1:100
        [t,y] = ode45(@(t,theta)ode(t,theta,omega,k,N,A),tspan,y0);
        yy{i,k} = y;
    end
end
 
yy30 = yy;
rr30 = zeros([10,100]);
for i = 1:10
    for j = 1:100
        rr30(i,j) = r(yy30{i,j},n,N);
    end
end
mean30 = sum(rr30)/10;
se30 = std(rr30)/sqrt(10);


% N = 50
N = 50;
yy = cell([10,100]);
for i = 1:10
    A = binornd(1,0.2,N,N);

    for a = 1:N
        A(a,a) = 0;
    end

    for b = 1:N
        for c = 1:N
            A(b,c) = A(c,b);
        end
    end

    omega = normrnd(0,0.5,[N,1]);
    y0 = unifrnd(0,2*pi,[N,1]);


    for k = 1:100
        [t,y] = ode45(@(t,theta)ode(t,theta,omega,k,N,A),tspan,y0);
        yy{i,k} = y;
    end
end
 
yy50 = yy;
rr50 = zeros([10,100]);
for i = 1:10
    for j = 1:100
        rr50(i,j) = r(yy50{i,j},n,N);
    end
end
mean50 = sum(rr50)/10;
se50 = std(rr50)/sqrt(10);

% N = 100
N = 100;
yy = cell([10,100]);
for i = 1:10
    A = binornd(1,0.2,N,N);

    for a = 1:N
        A(a,a) = 0;
    end

    for b = 1:N
        for c = 1:N
            A(b,c) = A(c,b);
        end
    end

    omega = normrnd(0,0.5,[N,1]);
    y0 = unifrnd(0,2*pi,[N,1]);


    for k = 1:100
        [t,y] = ode45(@(t,theta)ode(t,theta,omega,k,N,A),tspan,y0);
        yy{i,k} = y;
    end
end
 
yy100 = yy;
rr100 = zeros([10,100]);
for i = 1:10
    for j = 1:100
        rr100(i,j) = r(yy100{i,j},n,N);
    end
end
mean100 = sum(rr100)/10;
se100 = std(rr100)/sqrt(10);



% N = 200
N = 200;
yy = cell([10,100]);
for i = 1:10
    A = binornd(1,0.2,N,N);

    for a = 1:N
        A(a,a) = 0;
    end

    for b = 1:N
        for c = 1:N
            A(b,c) = A(c,b);
        end
    end

    omega = normrnd(0,0.5,[N,1]);
    y0 = unifrnd(0,2*pi,[N,1]);


    for k = 1:100
        [t,y] = ode45(@(t,theta)ode(t,theta,omega,k,N,A),tspan,y0);
        yy{i,k} = y;
    end
end
 
yy200 = yy;
rr200 = zeros([10,100]);
for i = 1:10
    for j = 1:100
        rr200(i,j) = r(yy200{i,j},n,N);
    end
end
mean200 = sum(rr200)/10;
se200 = std(rr200)/sqrt(10);


% N = 500
N = 500;
yy = cell([10,100]);
for i = 1:10
    A = binornd(1,0.2,N,N);

    for a = 1:N
        A(a,a) = 0;
    end

    for b = 1:N
        for c = 1:N
            A(b,c) = A(c,b);
        end
    end

    omega = normrnd(0,0.5,[N,1]);
    y0 = unifrnd(0,2*pi,[N,1]);


    for k = 1:100
        [t,y] = ode45(@(t,theta)ode(t,theta,omega,k,N,A),tspan,y0);
        yy{i,k} = y;
    end
end
 
yy500 = yy;
rr500 = zeros([10,100]);
for i = 1:10
    for j = 1:100
        rr500(i,j) = r(yy500{i,j},n,N);
    end
end
mean500 = sum(rr500)/10;
se500 = std(rr500)/sqrt(10);



errorbar(mean15,0.25*se15); hold on;
%errorbar(mean30,0.25*se30); hold on;
%errorbar(mean50,0.25*se50); hold on;
errorbar(mean100,0.25*se100); hold on;
errorbar(mean200,0.25*se200); hold on;
errorbar(mean500,0.25*se500);
% all p=0.2; 10 exp; 0.25*se; time 5 dt 0.1; each exp a different network;
xlabel('Coupling Strength K','FontSize', 14);
ylabel("Order Parameter r",'FontSize', 14);
legend({'N=15','N=100','N=200','N=500'},Location="southeast",Orientation="vertical")




