function bracket_sign_change(f, a, b)
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
return [a,b]
end