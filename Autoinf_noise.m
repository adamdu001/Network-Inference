% Autoinf with noise
function f = Autoinf_noise(lambda,exp,A,N,T,n,omega,y0,K,samples,threshold)

    tspan = linspace(0, T, n);
    dt = T/(n-1);
    y_cel = cell(samples,1);
    for i = 1:samples
        [t,y] = ode45(@(t,theta)ode(t,theta,omega,K,N,A),tspan,y0(:,i));
        y_cel{i} = y + lambda*normrnd(0,1,[n,N]);
%         global y_cel_global_noise
%         y_cel_global_noise{i,exp} = y_cel{i}; % for recording y
    end


    y_dot_cel = cell(samples,1);
    for i = 1:samples
        y_dot_cel{i} = zeros([n-1,N]);
        for j = 2:n
            y_dot_cel{i}(j-1,:) = (y_cel{i}(j,:)-y_cel{i}(j-1,:))/dt;
        end
    end


    Y_cel = cell(samples,1);
    for i = 1:samples
        Y_cel{i} = y_dot_cel{i}' - omega;
    end

    G = cell(samples,n-1);
    % for sample a, time b, row i, col j
    for a = 1:samples
        for b = 1:n-1
            for i = 1:N % for each row
                for j = 1:N
                    G{a,b}(i,j) = (K/N)*sin(y_cel{a}(b,j)-y_cel{a}(b,i));
                end
            end
        end
    end

    % now add small value to diagonal
    for a = 1:samples
        for b = 1:n-1
            for i = 1:N
                G{a,b}(i,i) = 0.0001;
            end
        end
    end

    f = CFS(A,G,Y_cel,threshold,N,n,samples);
end
