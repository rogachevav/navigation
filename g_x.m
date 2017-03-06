function [res] = g_x(phi, h, a, e2)
res = [0; 0; -g_phi(phi, h, a, e2)];
end
