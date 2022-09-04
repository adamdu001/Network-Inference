% tuning omega, multi network
N = 5;
samples = 1;
T = 5;
n = 51;
tspan = linspace(0, T, n);
dt = T/(n-1);

K = 0.1;
config_exp = 1000;
config_run = 50;
compare_step = 1;
y0cell = cell([config_exp,1]);
% book = []; % record those changed position
om = cell([config_exp,config_run+1]); % record all of config

for i=1:config_exp
    om{i,1} = normrnd(0,0.5,[N,1]);
    y0cell{i} = unifrnd(0,2*pi,[N,samples]);
end

w = cell([config_exp,config_run+1]); %config_exp experiments

global y_samples_vs_run
y_samples_vs_run = cell(samples,config_run,config_exp);


tic
for exp = 1:config_exp
    A = binornd(1,1,N,N);
    for i = 1:N
        A(i,i) = 0;
    end

    for i = 1:N
        for j = 1:N
            A(j,i) = A(i,j);
        end
    end
    for i = 1:config_run

        w{exp,i} = Autoinf_config_r(exp,i,A,N,5,51,om{exp,i},y0cell{exp},K,samples,0.86);

        swap = randi([1 N],[2,1]);
        omega_mirror = om{exp,i}; % create a copy of current config to operate on
        omega_mirror([swap(1),swap(2)]) = omega_mirror([swap(2),swap(1)]);

        w{exp,i+1} = Autoinf_1(A,N,5,51,omega_mirror,y0cell{exp},K,samples,0.86);

        if w{exp,i+1}(compare_step) > w{exp,i}(compare_step)
            %             book = [book,swap];
            om{exp,i+1} = omega_mirror;
            % update the omega config to the better performing one

        else
            %             book = [book,zeros([2,1])]; % otherwise let it be 0
            om{exp,i+1} = om{exp,i};
            % and w{i} stay as w{i}
        end
    end
end
toc
%% how large is the accuracy increase

colsum_acc = zeros([config_exp,config_run]);
for i = 1:config_run
    for j = 1:config_exp
        colsum_acc(j,i) = w{j,i}(compare_step);
    end
end
acc_rise_mean = sum(colsum_acc)/config_exp; %the mean
acc_rise_se = std(colsum_acc)/sqrt(config_exp); % the se


%% see how r has been changed
r_config = zeros([config_exp,config_run]);
for j = 1:config_run
    for exp = 1:config_exp
        for i = 1:samples
            r_config(exp,j) = r_config(exp,j) + r(y_samples_vs_run{i,j,exp},n,N);
        end
    end
end
r_config = r_config/samples; % take average all samples
r_config_mean = sum(r_config)/config_exp;
r_config_se = std(r_config)/sqrt(config_exp);



%% plot
layout = tiledlayout(1,3);

nexttile
errorbar(acc_rise_mean,acc_rise_se,'LineWidth',0.2)
xlabel('Number of Swaps','FontSize', 12);
ylabel("Accuracy rate at pertubation = "+compare_step,'FontSize', 12);
%title("Accuracy versus swaps in the \omega configuration at pertubation = "+compare_step);

nexttile
errorbar(r_config_mean,r_config_se)
xlabel('Number of Swaps','FontSize', 12);
ylabel("Average r",'FontSize', 12);

nexttile
scatter(r_config(:,1),r_config(:,50),10,'filled');hold on;
plot([0,1],[0,1]); 
xlim([0,1])
ylim([0,1])
xlabel('Starting r','FontSize', 12);
ylabel("Final r",'FontSize', 12);

% spacing
layout.TileSpacing = 'compact';
layout.Padding = 'compact';


