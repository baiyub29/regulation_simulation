function F = u(t,s2_t)
% control
F = 2*(sin(0.5*t)+cos(0.5*t)-sin(0.5*s2_t))+(cos(0.8*s2_t+0.8)+sin(0.5*s2_t)-3*sin(0.5*s2_t+0.5)-2*sin(0.5*s2_t+0.5))./(1+0.5*cos(s2_t));

end