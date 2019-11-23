clc
syms x y lamda;
b = [x;y];
f = 2*x^2 + y^2;
c= [1;1];
max = 10;
epsilon = 0.1;
grad = gradient(f,b);
k = 0;
while (k < max)
    g = eval(subs(grad,b,c));
    d = -g;
    nor = norm(d);
    if (nor > epsilon)
        s = c + lamda*d;
    else
        break
    end
    dif = eval(gradient(subs(f,b,s)));
    lamda_new = solve(dif == 0);
    c = c + lamda_new*d
    k = k + 1;
end