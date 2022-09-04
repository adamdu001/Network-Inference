function f = CFS(A,G,Y_cel,threshold,N,n,samples)
    AA = zeros([N,samples]);
    for oscillator = 1:N
        GG = []; YY = [];
        for sample = 1:samples
            for time = 1:n-1
                GG = [GG;G{sample,time}(oscillator,:)];%row concat for all time
                % then keep row concat GG for each sample i.e. 50, 100, 150, ...
            end
            YY = [YY,Y_cel{sample}(oscillator,:)];
            % col concat YY for each sample, i.e. 50, 100, 150, ...
            estimate = YY/GG'; % this is a row vector: 1 by N
            for element = 1:N
                if estimate(element) > threshold
                    estimate(element) = 1;
                else 
                    estimate(element) = 0;
                end
            end
            AA(oscillator,sample) = sum(estimate==A(oscillator,:));
        end
    end
    f = sum(AA)/N^2;
end