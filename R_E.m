function [res] = R_E(phi, a, e2)
res = a / ((1 - e2 * (sin(phi) ^ 2)) ^ 0.5);
end