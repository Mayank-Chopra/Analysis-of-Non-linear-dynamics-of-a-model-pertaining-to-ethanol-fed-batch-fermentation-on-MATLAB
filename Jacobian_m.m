function gradf = Jacobian_m(f,x_in)
len = length(x_in);
for i=1:len
    eps = x_in(i)/10000;
    if eps == 0
        eps = 1/10000;
    end
    xi = x_in;
    xf = x_in;
    xf(i) = xf(i) + eps;
    fi = f(xi(1),xi(2),xi(3),xi(4));
    ff = f(xf(1),xf(2),xf(3),xf(4));
    gradf(:,i) = (ff-fi)/eps;
end
end