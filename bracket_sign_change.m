function newbracket = bracket_sign_change(f, a, b)
% Bracket Sign Change for Root Finding
% Kochenderfer & Wheeler, Algorithms for Optimization, Algorithm 3.7

k=2;
if a > b %ensure a < b
    a = b;
    b = a;
end
center = (b+a)/2;
half_width = (b-a)/2;

while f(a)*f(b) > 0
    half_width = half_width*k;
    a = center - half_width;
    b = center + half_width;
end
newbracket = [a,b];
end