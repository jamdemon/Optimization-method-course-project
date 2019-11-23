clc
syms a0 a1 a2 stepsize
Omiga = [0.5236, 1.0472, 1.5708, 2.0944, 2.6180];
P = [0.6821, 4.3232, 3.7540, 0.4368, 0.1988];
%P_hat = 1 ./ abs(a0 + a1 * exp(-i * Omiga) + a2 * exp(-2 * i * Omiga)) .^ 2;
P_hat = 1 ./ ((a0 + a1 * cos(-Omiga) + a2 * cos(-2 * Omiga)) .^ 2 + ...
    (a1 * sin(-Omiga) + a2 * sin(-2 * Omiga)) .^ 2);
e = P ./ P_hat - log(P ./ P_hat) - 1;
E = sum(e);
 
X = [a0; a1; a2];
point = [1; 0; 0];
optimal = [1; -0.5161; 0.9940];
% grad_a0 = vpa(diff(E, a0), 10);
% grad_a1 = vpa(diff(E, a1), 10);
% grad_a2 = vpa(diff(E, a2), 10);
% grad = [grad_a0; grad_a1; grad_a2];
grad = vpa(gradient(E),10);
value = [];
distance = [];
while true
    g = eval(subs(grad, X, point));
    d = -g;
    X_new = point + stepsize * d;
    f = subs(E, X, X_new);
    stepsize_new = golden_section_search(f);
    point = point + stepsize_new * d;
    
    temp_value = eval(subs(E, X, point));
    value = [value; temp_value];
    distance = [distance; norm(point - optimal) ^ 2];
    if norm(stepsize_new * d) < 0.0001
        break
    end
end
Iteration_number = length(value)
Point = point
E_value = temp_value
Distance = norm(point- optimal) ^ 2


figure(1)
semilogy((1: length(value)), value);
title('E vs Iteration Number');
xlabel('Iteration Number');
ylabel('Value of E');

figure(2)
semilogy((1: length(distance)), distance);
title('Distance vs Iteration Number');
xlabel('Iteration Number');
ylabel('Distance');
 
function y = golden_section_search(f)
syms stepsize
a = -1;
b = 1;
x1 = a + 0.382 * (b - a);
x2 = a + 0.618 * (b - a);
while true
if b - a < 0.001
    break
end
f1 = eval(subs(f, stepsize, x1));
f2 = eval(subs(f, stepsize, x2));
if f2 < f1
    a = x1;
    b = b;
    x1 = x2;
    x2 = a + 0.618 * (b - a);
else
    a = a;
    b = x2;
    x2 = x1;
    x1 = a + 0.382 * (b - a);
end
end
y = (x1 + x2) / 2;
end

