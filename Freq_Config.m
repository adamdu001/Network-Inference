% Now consider different config of natural frequencies

% fix the variables
N = 20;
samples = 30;
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

        % if acc higher then save it to book
        % update omega according to latest col of book
        % keep climbing hills until no longer higher for a while

%% compare acc and r part1
K = 0.1;
config_run = 500;
compare_step = 1;
book = []; % record those changed position
om = cell([1,config_run+1]); % record all of config
om{1} = omega;


w = cell([1,config_run+1]);

global y_samples_vs_run
y_samples_vs_run = cell(samples,config_run);

for i = 1:config_run
    
    w{i} = Autoinf_config_r(i,A,N,5,51,om{i},y0,K,samples,0.86);

    swap = randi([1 N],[2,1]);
    omega_mirror = om{i}; % create a copy of current config to operate on
    omega_mirror([swap(1),swap(2)]) = omega_mirror([swap(2),swap(1)]);
    
    w{i+1} = Autoinf_1(A,N,5,51,omega_mirror,y0,K,samples,0.86);

    if w{i+1}(compare_step) > w{i}(compare_step)
        book = [book,swap];
        om{i+1} = omega_mirror;
        % update the omega config to the better performing one
        
    else 
        book = [book,zeros([2,1])]; % otherwise let it be 0
        om{i+1} = om{i};
        % and w{i} stay as w{i}
    end
end

%% how large is the accuracy increase

acc_rise = zeros([1,config_run]);
for i = 1:config_run
    acc_rise(i) = w{i}(compare_step);
end
subplot(121)
plot(1:config_run,acc_rise,'LineWidth',1);
xlabel('Number of Swaps','FontSize', 12);
ylabel("Accuracy rate at pertubation = "+compare_step,'FontSize', 12);
%title("Accuracy versus swaps in the \omega configuration at pertubation = "+compare_step);


%% see how r has been changed
r_config = zeros([1,config_run]);
for j = 1:config_run
    for i = 1:samples
        r_config(j) = r_config(j) + r(y_samples_vs_run{i,j},n,N);
    end
end
r_config = r_config/samples;
subplot(122)
plot(r_config, 'LineWidth',1)
xlabel('Number of Swaps','FontSize', 12);
ylabel("Oder parameter r at pertubation = "+compare_step,'FontSize', 12);
%title("Order parameter versus swaps in the \omega configuration at pertubation = "+compare_step);



%%
bookp1 = book;
omp1 = om;
conf_p1 = (1:N)';
bookp1_refined = [];
for i = 1:length(bookp1)
    if bookp1(:,i) ~= [0;0]
        bookp1_refined = [bookp1_refined,bookp1(:,i)];
    end
end

for i = 1:length(bookp1_refined)
    conf_p1([bookp1_refined(1,i),bookp1_refined(2,i)]) = conf_p1([bookp1_refined(2,i),bookp1_refined(1,i)]);
end
% 
% 
% bookp2 = book;
% omp2 = om;
% conf_p2 = (1:N)';
% bookp2_refined = [];
% for i = 1:length(bookp2)
%     if bookp2(:,i) ~= [0;0]
%         bookp2_refined = [bookp2_refined,bookp2(:,i)];
%     end
% end
% 
% for i = 1:length(bookp2_refined)
%     conf_p2([bookp2_refined(1,i),bookp2_refined(2,i)]) = conf_p2([bookp2_refined(2,i),bookp2_refined(1,i)]);
% end
% 
% 
% 
% 
% 
