%Trapazoidal Rule Function
function [s] = trap(funk,a,b,N)

    h = (b-a)/N;
    sumstart = funk(a) + funk(b);
    sum = 0;
    for i = 1:N-1
        x = a + i*h;
        sum = sum + funk(x);
    end
    s = (h/2) * (sumstart +2*sum);

end