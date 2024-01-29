function F = s2(t)
% s_2^out(t,1)
F = fsolve(@(x) 1.5*x+sin(2*pi*x)/(2*pi)-1.5*t-sin(2*pi*t)/(2*pi)-1, t);

end