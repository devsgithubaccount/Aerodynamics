%Simpson's rule Function
function [s] = simp(funk,a,b,N)
    
    h = (b-a)/N;
    sum1 = 0;
    sum2 = 0;
        for i = 1:2:(N-1)
            x = a + h*i;
            sum1 = sum1 + funk(x);
        end
        for i = 0:2:N
            x = a + h*i;
            sum2 = sum2 + funk(x);
        end
    s = (h/3) * (funk(a) + 4*sum1 + 2*sum2 +funk(b));

end