%my_pdf=@(x) (exp(-(((x-1)/1).^2)/2))/1/sqrt(2*pi);
%my_pdf=@(x) (x*x);

n = 50000;
x = zeros (n,1);
x(1) = rand;
for i = 1:n-1
    x_c = normrnd (x(i), 0.1);
    if rand < min (1, mypdf(x_c)/mypdf(x(i))) x(i+1) = x_c;
    else x(i+1) = x(i);
    end
end

histogram (x)