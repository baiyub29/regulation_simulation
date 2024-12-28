function F = s2in_tx(t, x)
% s_2^in(t,x)
F = fsolve(@(s) x-1.5*(s-t)+(sin(2*pi*t)-sin(2*pi*s))/(2*pi)-1, t);

end