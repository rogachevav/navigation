function [x, y, z] = convertation(coord, R_earth)    
x = (R_earth + coord(:, 3)).*cos(coord(:, 1)).*cos(coord(:, 2));    
y = (R_earth + coord(:, 3)).*cos(coord(:, 1)).*sin(coord(:, 2));    
z = (R_earth + coord(:, 3)).*sin(coord(:, 1)); 
end
