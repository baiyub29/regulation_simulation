function F = u_s(t,s2,w,v1,v2)
% control
F = -w+2*(1-sin(0.5*(s2-t))).*v1+2*(1-cos(0.5*(s2-t))).*v2 +...
    (cos(0.8*(s2+1)) + (sin(0.5*(s2-t))-5*sin(0.5*(s2-t+1))).*v1 + ...
    (cos(0.5*(s2-t))-5*cos(0.5*(s2-t+1))).*v2) ./ (1+0.5*cos(s2));

end