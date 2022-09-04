function accuracy = acc(A,A_hat,N)
         m = A == A_hat;
         s = sum(sum(m));
         accuracy = s/(N^2);
end