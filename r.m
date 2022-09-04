% order parameter r
function orderpara = r(y,n,N)
         sinvalue = sum(sin(y'));
         cosvalue = sum(cos(y'));
         s = 0;
         for t = 1:n
             s = s + (1/N)*sqrt(sinvalue(t)^2+cosvalue(t)^2);
         end
         orderpara = s/n;
end

