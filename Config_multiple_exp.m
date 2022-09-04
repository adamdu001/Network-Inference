% config test with multiple experiments
N = 50;
samples = 1;
y0 = unifrnd(0,2*pi,[N,samples]);
omega = normrnd(0,0.5,[N,1]);

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

%% compare acc and r part1
K = 20;
config_exp = 100;
config_run = 200;
compare_step = 1;
book = []; % record those changed position
om = cell([config_exp,config_run+1]); % record all of config
for i=1:config_exp
    om{i,1} = omega;
end

w = cell([config_exp,config_run+1]); %config_exp experiments

% global y_samples_vs_run
% y_samples_vs_run = cell(samples,config_run);

global y_samples_vs_run
% y_samples_vs_run = cell(samples,config_run);
y_samples_vs_run = cell(samples,config_run,config_exp);

for exp = 1:config_exp
    for i = 1:config_run

        w{exp,i} = Autoinf_config_r(exp,i,A,N,5,51,om{exp,i},y0,K,samples,0.86);

        swap = randi([1 N],[2,1]);
        omega_mirror = om{exp,i}; % create a copy of current config to operate on
        omega_mirror([swap(1),swap(2)]) = omega_mirror([swap(2),swap(1)]);

        w{exp,i+1} = Autoinf_1(A,N,5,51,omega_mirror,y0,K,samples,0.86);

        if w{exp,i+1}(compare_step) > w{exp,i}(compare_step)
            book = [book,swap];
            om{exp,i+1} = omega_mirror;
            % update the omega config to the better performing one

        else
            book = [book,zeros([2,1])]; % otherwise let it be 0
            om{exp,i+1} = om{exp,i};
            % and w{i} stay as w{i}
        end
    end
end
%% how large is the accuracy increase

colsum_acc = zeros([config_exp,config_run]);
for i = 1:config_run
    for j = 1:config_exp
        colsum_acc(j,i) = w{j,i}(compare_step);
    end
end
acc_rise_mean = sum(colsum_acc)/config_exp; %the mean
acc_rise_se = std(colsum_acc)/sqrt(experiments); % the se

subplot(211)
errorbar(acc_rise_mean,acc_rise_se,'LineWidth',0.2)
xlabel('Number of Swaps','FontSize', 12);
ylabel("Accuracy rate at pertubation = "+compare_step,'FontSize', 12);
%title("Accuracy versus swaps in the \omega configuration at pertubation = "+compare_step);


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


subplot(212)
errorbar(r_config_mean,r_config_se,'LineWidth',0.1)
xlabel('Number of Swaps','FontSize', 12);
ylabel("Oder parameter r at pertubation = "+compare_step,'FontSize', 12);
%title("Order parameter versus swaps in the \omega configuration at pertubation = "+compare_step);

%% study confi details (not needed in tuning freq as we consider multiple networks)
bookcell = cell([1,config_exp]);
for i=0:(config_exp-1)
    bookcell{i+1} = book(:,config_run*i+(1:config_run));
end

conf_final = zeros([N,config_exp]);
for exp = 1:config_exp
    conf_final(:,exp) = (1:N)';
end

for exp = 1:config_exp
    bookp1 = bookcell{exp};
    bookp1_refined = [];

    % effective permutations bookp1_refined
    for i = 1:length(bookp1)
        if bookp1(:,i) ~= [0;0]
            bookp1_refined = [bookp1_refined,bookp1(:,i)];
        end
    end

    a = conf_final(:,exp);

    for i = 1:length(bookp1_refined)
        a([bookp1_refined(1,i),bookp1_refined(2,i)]) = a([bookp1_refined(2,i),bookp1_refined(1,i)]);
    end
    conf_final(:,exp) = a;

end

%% draw the network before and after

% before
names = cell([10,1]);
for i = 1:10
    names{i} = strcat('\omega_{',num2str(i),'}=',num2str(round(omega(i),3)));
end
G = graph(A,names);
p = plot(G);


% afterwards, no.9, 99% accuracy
names_after = cell([10,1]);
for i = 1:10
    names_after{i} = strcat('\omega_{',num2str(conf_final(i,9)),'}=',num2str(round(omega(conf_final(i,9)),3)));
end
G_after = graph(A,names_after);
p_after = plot(G_after);



layout = tiledlayout(2,1);
% 1
nexttile
plot(G,NodeFontSize=10);

% 2
nexttile
plot(G_after,NodeFontSize=10);

% spacing
layout.TileSpacing = 'compact';
layout.Padding = 'compact';
