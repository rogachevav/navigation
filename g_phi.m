function [g] = g_phi(phi, h, a, e2)
g = 9.78030 * (1 - 2 * h / a + (4/5) * e2 * sin(phi) ^ 2);
end
