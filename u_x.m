function [res] = u_x(phi, u)
res = [0; u * cos(phi); u * sin(phi)];
end
