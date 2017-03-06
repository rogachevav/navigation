A = load('./imu.dat');
B = load('./trj.dat');
%загрузили показания
A(:, 2:4) = A(:, 2:4).*(pi / 180);
B(:, 2:3) =B(:, 2:3).* (pi / 180);
B(:, 8:10) = B(:, 8:10).* (pi / 180);
%перевели в радианы
Init = A(:, 2).*A(:, 5) +  A(:, 3).*A(:, 6) +  A(:, 4).*A(:, 7);
A_init = A(A(:, 1) < 60, :);
%отрезали кусок для выставки
R_earth = 6371000;
a = 6378137.0; 
b = 6356752.0; 
e2 = 6.6943799901413 * 10 ^ (-3);     %> lim n-> inf (ksi + ksi - n*a)/std sqrt n
delta_t = 0.01;
phi_0 = B(1, 2); 
lambda_0 = B(1, 3); 
h_0 = B(1, 4);
u = 2 * pi / 86164.090530833;
%константы+н.у.
%L_0 = initial_alignment(A_init, phi_0, h_0, a, e2, u);
disp(L_0), disp(det(L_0))
%Z = 
%L_0 =initial_alignment(A_init, phi_0, h_0, a, e2, u);
A_x_0= eye(3);
A_z_0= L_0;
%phi_0= B(1, 2);
%lambda_0 = B(1, 3); h_0 = B(1, 4);
V_0 =[0; 0; 0];
alpha = [0, 0, 0];
coordinates = [phi_0, lambda_0, h_0, V_0'];
for row = A'
	t = row(1)';
	if mod(t, 250.0) == 0
		disp(t)
    end
	w_z         = row(2:4)';
	w_z_matrix  = create_matrix(w_z);
	A_z         = matrix_multiplier(w_z, w_z_matrix, delta_t) * A_z_0;
	w_x         = u_x(phi_0, u) + Omega(V_0, phi_0, h_0, a, e2);
	w_x_matrix  = create_matrix(w_x);
	A_x         = matrix_multiplier(w_x, w_x_matrix, delta_t) * A_x_0;
	L           = A_z * A_x';
    if mod(t, 10.0) == 0
        temp = L_0 * L' - eye(size(L));
		alpha(end+1, :) = [temp(2, 3), -temp(1, 3), temp(1, 2)];
    end
	f_z         = row(5:7);
	V           = V_0 + (create_matrix(Omega(V_0, phi_0, h_0, a, e2) + 2*u_x(phi_0, u))*V_0 + g_x(phi_0, h_0, a, e2) + L_0'*f_z)*delta_t;
	phi         = phi_0 + V_0(2) / (R_N(phi_0, a, e2) + h_0) * delta_t;
	lambda      = lambda_0 + V_0(1) / ((R_E(phi_0, a, e2) + h_0) * cos(phi_0)) * delta_t;
	h           = h_0 + V_0(3) * delta_t;
	coordinates(end + 1, :) = [phi, lambda, h, V'];
    
	A_z_0 = A_z; A_x_0 = A_x; L_0 = L; V_0 = V; phi_0 = phi; lambda_0 = lambda; h_0 = h;
end;
hold on;
[X11, X12, X13] = convertation(B(:, 2:4), R_earth);
[X21, X22, X23] = convertation(coordinates, R_earth);
hold on;
%new_coordinates = coordinates;
%new_coordinates(:,1) = (coordinates(:,2) - lambda_0)*R_E(phi_0, a, e2)*cos(phi_0);
%new_coordinates(:,2) = (coordinates(:,1) - phi_0)*R_N(phi_0, a, e2);
%new_coordinates(:,3) = coordinates(:,3) - h_0;
grid on;
daspect([1 1 1]);
hold on;
%true_coord = B(:,2:4);
%true_coord(:,1) = (true_coord(:,1) - lambda_0)*R_E(phi_0, a, e2)*cos(phi_0);
%true_coord(:,2) = (true_coord(:,2) - phi_0)*R_N(phi_0, a, e2);
%true_coord(:,3) = true_coord(:,3) - h_0;
%plot3(true_coord(:,1), true_coord(:,2), true_coord(:,3),'red')
%hold on;
%plot3(new_coordinates(:,1), new_coordinates(:,2), new_coordinates(:,3),'blue')
hold on;
%hold off;
%new_coordinates(2:end, 1)
plot3(X11, X12, X13, X21, X22, X23)
grid on
daspect([1 1 1])