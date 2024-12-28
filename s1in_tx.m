function F = s1in_tx(t, x)
% s_1^in(t,x)
F = fsolve(@(s) x+s-t+(cos(2*pi*t)-cos(2*pi*s))/(4*pi), t);

end