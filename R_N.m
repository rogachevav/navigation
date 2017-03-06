function [res] = R_N(phi, a, e2)
res = a * (1 - e2) / ((1 - e2 * (sin(phi) ^ 2)) ^ 1.5);
end
