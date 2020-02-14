function [newbracketL, newbracketU] = bisection(f, a, b, varargin)
% Bisection Method for Root Finding
% Kochenderfer & Wheeler, Algorithms for Optimization, Algorithm 3.6

if isempty(varargin)
    eps = 1e-10;
else
    eps = varargin{1};
end

if a > b % ensure a < b
    a = b;
    b = a;
end
ya = f(a);
yb = f(b);
if ya == 0
    b = a;
end
if yb == 0
    a = b;
end
while b - a > eps
    x = (a+b)/2;
    y = f(x);
    if y == 0
        a = x;
        b = x;
    elseif sign(y) == sign(ya)
        a=x;
    else
        b=x;
    end
end
newbracketL = a;
newbracketU = b;
end