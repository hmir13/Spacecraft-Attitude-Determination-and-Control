close all ; clc ; clear

%% Problem Constants
% Earth constants
mu_e = 398600.4418; % km^3/s^2

% initial state
r_0 = [1029.7743
       6699.3469
       3.7896]; % km
v_0 = [-6.2119
        0.9524
        4.3946]; % km/s

z_0 = [0.685
       0.695
       0.153
       0.153];
q_0 = z_0 / norm(z_0);

w_0 = [0.53
       0.53
       0.053]; % deg/s
w_0 = deg2rad(w_0);

h_w_0 = [0
         0
         0]; % kg-m^2/s

% desired final state
q_c = [0
       0
       0
       1];

w_c = [0
       0
       0]; % rad/s

% inertial properties
I = [6400   -76.4   -25.6
     -76.4   4730    -40
     -25.6   -40     8160]; % kg-m^2

% time span
t_0 = 0; % s
dt = 5; % s
t_f = 2 * 3600; % s
t_span = t_0:dt:t_f; % s

%% Attitude reorientation maneuver using reaction wheels

x_0 = [r_0;
       v_0];
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
[~, out] = ode45(@(t,x) keplers_2body(t,x,mu_e), t_span, x_0, opts);

X = out(:,1);
Y = out(:,2);
Z = out(:,3);
X_dot = out(:,4);
Y_dot = out(:,5);
Z_dot = out(:,6);

r_eci = [X, Y, Z];

figure();
sgtitle("Position vs. time");

subplot(3,1,1);
plot(t_span, X, DisplayName="X(t)");
ylabel("X [km]");

subplot(3,1,2);
plot(t_span, Y, DisplayName="Y(t)");
ylabel("Y [km]");

subplot(3,1,3);
plot(t_span, Z, DisplayName="Z(t)");
xlabel("time [s]");
ylabel("Z [km]");
snapnow

figure();
plot3(X, Y, Z);
title("Orbit in the inertial frame");
xlabel("X [km]");
ylabel("Y [km]");
zlabel("Z [km]");
grid;
snapnow

k_p = 0;
k_d = 0;
x_0 = [q_0 ; w_0 ; h_w_0];

[~, out] = ode45(@(t,x) RW_controller(t,x,I,q_c,w_c,k_p,k_d), t_span, x_0);

q = out(:, 1:4);
w = out(:, 5:7);

figure();
sgtitle("Quaternion vs. time")

subplot(4,1,1)
plot(t_span/60, q(:, 1), DisplayName="q_1")
ylabel("q_1")

subplot(4,1,2)
plot(t_span/60, q(:, 2), DisplayName="q_2")
ylabel("q_2")

subplot(4,1,3)
plot(t_span/60, q(:, 3), DisplayName="q_3")
ylabel("q_3")

subplot(4,1,4)
plot(t_span/60, q(:, 4), DisplayName="q_4")
ylabel("q_4")
xlabel("time [min]")
    
figure();
sgtitle("Angular velocity vs. time")

subplot(3,1,1)
plot(t_span/60, w(:, 1), DisplayName="w_1")
ylabel("\omega_1 (rad/s)")

subplot(3,1,2)
plot(t_span/60, w(:, 2), DisplayName="w_2")
ylabel("\omega_2 (rad/s)")

subplot(3,1,3)
plot(t_span/60, w(:, 3), DisplayName="w_3")
ylabel("\omega_3 (rad/s)")
xlabel("time [min]")

k_p = 10;
k_d = 150;
x_0 = [q_0 ; w_0 ; h_w_0];

[~, out] = ode45(@(t,x) RW_controller(t,x,I,q_c,w_c,k_p,k_d), t_span, x_0, opts);

q = out(:, 1:4);
w = out(:, 5:7);
h_w = out(:, 8:10);

figure();
sgtitle("Quaternion vs. time")

subplot(4,1,1)
plot(t_span/60, q(:, 1), DisplayName="q_1")
ylabel("q_1")

subplot(4,1,2)
plot(t_span/60, q(:, 2), DisplayName="q_2")
ylabel("q_2")

subplot(4,1,3)
plot(t_span/60, q(:, 3), DisplayName="q_3")
ylabel("q_3")

subplot(4,1,4)
plot(t_span/60, q(:, 4), DisplayName="q_4")
ylabel("q_4")
xlabel("time [min]")

figure();
sgtitle("Angular velocity vs. time")

subplot(3,1,1)
plot(t_span/60, w(:, 1), DisplayName="w_1")
ylabel("\omega_1 (rad/s)")

subplot(3,1,2)
plot(t_span/60, w(:, 2), DisplayName="w_2")
ylabel("\omega_2 (rad/s)")

subplot(3,1,3)
plot(t_span/60, w(:, 3), DisplayName="w_3")
ylabel("\omega_3 (rad/s)")
xlabel("time [min]")

figure();
sgtitle("RW Angular momentum vs. time")

subplot(3,1,1)
plot(t_span/60, h_w(:, 1), DisplayName="h_1")
ylabel("h_1 (kg-m^2/s)")

subplot(3,1,2)
plot(t_span/60, h_w(:, 2), DisplayName="h_2")
ylabel("h_2 (kg-m^2/s)")

subplot(3,1,3)
plot(t_span/60, h_w(:, 3), DisplayName="h_3")
ylabel("h_3 (kg-m^2/s)")
xlabel("time [min]")

t_star = 10 * 60;
t_star_idx = find(t_span == t_star);
fprintf("q at t=10 min: [%f, %f, %f, %f]\n", q(t_star_idx, :));
fprintf("w at t=10 min: [%f, %f, %f]\n", w(t_star_idx, :));
fprintf("h_w at t=10 min: [%f, %f, %f]\n\n", h_w(t_star_idx, :));

%% Momentum dumping using magnetic torquers

k_p = 10;
k_d = 150;
k_m = 0.001;
x_0 = [q_0 ; w_0 ; h_w_0];
b_eci = readmatrix("bvec_ECI_rec.xlsx");

t_span1 = t_0:dt:t_star;
[~, out1] = ode45(@(t,x) RW_controller(t,x,I,q_c,w_c,k_p,k_d), t_span1, x_0, opts);

x_0 = out1(end, :);
t_span2 = t_star:dt:t_f;
[~, out2] = ode45(@(t,x) RW_Mag_controller(t,x,I,q_c,w_c,b_eci,k_p,k_d,k_m,t_span), t_span2, x_0, opts);

out = [out1; out2(2:end, :)];

q = out(:, 1:4);
w = out(:, 5:7);
h_w = out(:, 8:10);

figure;
sgtitle("RW Angular momentum vs. time")

subplot(3,1,1)
plot(t_span/60, h_w(:, 1), DisplayName="h_1")
ylabel("h_1 (kg-m^2/s)")

subplot(3,1,2)
plot(t_span/60, h_w(:, 2), DisplayName="h_2")
ylabel("h_2 (kg-m^2/s)")

subplot(3,1,3)
plot(t_span/60, h_w(:, 3), DisplayName="h_3")
ylabel("h_3 (kg-m^2/s)")
xlabel("time [min]")

m_array = zeros(length(t_span), 3);
L_m_array = zeros(length(t_span), 3);
k_m = 0;
for i=1:length(t_span)
    A_q = quat_att_matrix(q(i, :));
    b_body = A_q * b_eci(:, i);
    b_body_norm = norm(b_body);
    b = b_body / b_body_norm;

    if (t_span(i) == t_star)
        k_m = 0.001;
    end

    m_array(i,:) = (k_m / b_body_norm) * cross(h_w(i, :), b);

    L_m_array(i,:) = -k_m * (eye(3) - b * b') * h_w(i, :)';
end

figure;
sgtitle("Magnetic Torque vs. time")

subplot(3,1,1)
plot(t_span/60, L_m_array(:, 1))
ylabel("\tau_1 (Nm)")

subplot(3,1,2)
plot(t_span/60, L_m_array(:, 2))
ylabel("\tau_2 (Nm)")

subplot(3,1,3)
plot(t_span/60, L_m_array(:, 3))
ylabel("\tau_3 (Nm)")

figure;
sgtitle("Commanded dipole vs. time")

subplot(3,1,1)
plot(t_span/60, m_array(:, 1))
ylabel("m_1 (Am^2)")

subplot(3,1,2)
plot(t_span/60, m_array(:, 2))
ylabel("m_2 (Am^2)")

subplot(3,1,3)
plot(t_span/60, m_array(:, 3))
ylabel("m_3 (Am^2)")

h_tot = I * w' + h_w';
h_tot_norm = zeros(length(t_span),1);
for i=1:length(t_span)
    h_tot_norm(i) = norm(h_tot(:,i));
end

figure;
plot(t_span'/60, h_tot_norm)
ylim([0, 100])
title("Total Angular momentum norm vs. time")
ylabel("||h_{tot}|| (kg-m^2/s)")
xlabel("time [min]")

figure;
ylim([0, 100])
title("Total Angular momentum norm vs. time")
ylabel("||h_{tot}|| (kg-m^2/s)")
xlabel("time [min]")
hold on;
k_m_values = [0.0001, 0.001, 0.01, 0.02];
for i=1:length(k_m_values)
    k_m = k_m_values(i);

    t_span1 = t_0:dt:t_star;
    [~, out1] = ode45(@(t,x) RW_controller(t,x,I,q_c,w_c,k_p,k_d), t_span1, x_0, opts);
    
    x_0 = out1(end, :);
    t_span2 = t_star:dt:t_f;
    [~, out2] = ode45(@(t,x) RW_Mag_controller(t,x,I,q_c,w_c,b_eci,k_p,k_d,k_m,t_span), t_span2, x_0, opts);
    
    out = [out1; out2(2:end, :)];
    
    w = out(:, 5:7);
    h_w = out(:, 8:10);

    h_tot = I * w' + h_w';
    h_tot_norm = zeros(length(t_span),1);
    for j=1:length(t_span)
        h_tot_norm(j) = norm(h_tot(:,j));
    end

    plot(t_span'/60, h_tot_norm, 'DisplayName', ['k_m = ', num2str(k_m)])
end
hold off;
legend;

%% Spacecraft detumbling using magnetic torquers

k_D = 10;
q_0 = (1/sqrt(2)) * [1 0 0 1]';
w_0 = [0.01 0.01 0.01]';
x_0 = [q_0 ; w_0];

[~, out] = ode45(@(t,x) Mag_controller(t,x,I,b_eci,k_D,t_span), t_span, x_0, opts);

q = out(:, 1:4);
w = out(:, 5:7);

L_D_array = zeros(length(t_span), 3);
for i=1:length(t_span)
    A_q = quat_att_matrix(q(i, :));
    b_body = A_q * b_eci(:, i);
    b_body_norm = norm(b_body);
    b = b_body / b_body_norm;

    L_D_array(i,:) = -k_D * (eye(3) - b * b') * w(i, :)';
end

figure;
sgtitle("Angular velocities vs. time")

subplot(3,1,1)
plot(t_span/60, w(:, 1), DisplayName="w_1")
ylabel("\omega_1 (rad/s)")

subplot(3,1,2)
plot(t_span/60, w(:, 2), DisplayName="w_2")
ylabel("\omega_2 (rad/s)")

subplot(3,1,3)
plot(t_span/60, w(:, 3), DisplayName="w_3")
ylabel("\omega_3 (rad/s)")

figure;
sgtitle("Detumbling torque vs. time")

subplot(3,1,1)
plot(t_span/60, L_D_array(:, 1))
ylabel("\tau_1 (Nm)")

subplot(3,1,2)
plot(t_span/60, L_D_array(:, 2))
ylabel("\tau_2 (Nm)")

subplot(3,1,3)
plot(t_span/60, L_D_array(:, 3))
ylabel("\tau_3 (Nm)")

figure;
sgtitle("Angular velocities vs. time")
subplot(3,1,1)
ylabel("\omega_1 (rad/s)")
subplot(3,1,2)
ylabel("\omega_2 (rad/s)")
subplot(3,1,3)
ylabel("\omega_3 (rad/s)")
k_D_values = [5, 10, 50, 100];
for i=1:length(k_D_values)
    k_D = k_D_values(i);
    
    [~, out] = ode45(@(t,x) Mag_controller(t,x,I,b_eci,k_D,t_span), t_span, x_0, opts);
    
    w = out(:, 5:7);

    subplot(3,1,1)
    hold on;
    plot(t_span/60, w(:, 1), 'DisplayName', ['k_D = ', num2str(k_D)])
    
    subplot(3,1,2)
    hold on;
    plot(t_span/60, w(:, 2), 'DisplayName', ['k_D = ', num2str(k_D)])
    
    subplot(3,1,3)
    hold on;
    plot(t_span/60, w(:, 3), 'DisplayName', ['k_D = ', num2str(k_D)])
end
hold off;
legend;

%% Functions

% Quaternion functions
function d_q = quat_error(q,q_c)
    sigma_q_c = quat_tilde_matrix(q_c);

    dq13 = sigma_q_c' * q;
    dq4 = q_c' * q;

    d_q = [dq13; dq4];
end

function sigma_q = quat_tilde_matrix(q)
    q1 = q(1); 
    q2 = q(2); 
    q3 = q(3); 
    q4 = q(4);

    sigma_q = [ q4, -q3,  q2
                q3,  q4, -q1
               -q2,  q1,  q4
               -q1, -q2, -q3];
end

function A_q = quat_att_matrix(q)
    q1 = q(1); 
    q2 = q(2); 
    q3 = q(3); 
    q4 = q(4);

    A_q = [q1^2 - q2^2 - q3^2 + q4^2, 2 * (q1 * q2 + q3 * q4)   , 2 * (q1 * q3 - q2 * q4)   
           2 * (q2 * q1 - q3 * q4)  , -q1^2 + q2^2 - q3^2 + q4^2, 2 * (q2 * q3 + q1 * q4)   
           2 * (q3 * q1 + q2 * q4)  , 2 * (q3 * q2 - q1 * q4)   , -q1^2 - q2^2 + q3^2 + q4^2];
end

% RW controller
function dx = RW_controller(~,x,I,q_c,w_c,k_p,k_d)
    q = x(1:4);
    w = x(5:7);
    h_w = x(8:10);

    dq = quat_error(q,q_c);
    dw = w - w_c;

    L_s = -k_p * sign(dq(4)) * dq(1:3) - k_d * dw;
    L_w = -L_s;

    sigma_q = quat_tilde_matrix(q);
    
    q_dot = (1/2) * sigma_q * w;
    w_dot = inv(I) *(L_s - cross(w,I*w));
    h_dot = L_w - cross(w,h_w);

    dx = [q_dot ; w_dot ; h_dot];

end

function dx = keplers_2body(~,x,mu)
    r = x(1:3);
    v = x(4:6);
 
    r_norm = sqrt(r(1)^2 + r(2)^2 + r(3)^2);

    a = -mu / r_norm^3 * r;

    dx = [v;a];
end

% RW and MagTorquer controller
function b = local_body_geomagnetic_field(q, b_eci, t_idx)
    A_q = quat_att_matrix(q);

    b_body = A_q * b_eci(:, t_idx);
    
    b = b_body / norm(b_body);
end

function dx = RW_Mag_controller(t,x,I,q_c,w_c,b_eci,k_p,k_d,k_m,t_span)
    q = x(1:4);
    w = x(5:7);
    h_w = x(8:10);

    dq = quat_error(q,q_c);
    dw = w - w_c;

    [~, t_idx] = min(abs(t_span - t));
    A_q = quat_att_matrix(q);
    b_body = A_q * b_eci(:, t_idx);
    b = b_body / norm(b_body);

    L_s = -k_p * sign(dq(4)) * dq(1:3) - k_d * dw;
    L_m = -k_m * (eye(3) - b * b') * h_w;
    L_w = -(L_s);

    sigma_q = quat_tilde_matrix(q);
    
    q_dot = (1/2) * sigma_q * w;
    w_dot = inv(I) *(L_s + L_m - cross(w,I*w));
    h_dot = L_w - cross(w,h_w);

    dx = [q_dot ; w_dot ; h_dot];

end

% MagTorquer controller
function dx = Mag_controller(t,x,I,b_eci,k_D,t_span)
    q = x(1:4);
    w = x(5:7);

    [~, t_idx] = min(abs(t_span - t));
    A_q = quat_att_matrix(q);
    b_body = A_q * b_eci(:, t_idx);
    b = b_body / norm(b_body);

    L_D = -k_D * (eye(3) - b * b') * w;

    sigma_q = quat_tilde_matrix(q);
    
    q_dot = (1/2) * sigma_q * w;
    w_dot = inv(I) *(L_D - cross(w,I*w));

    dx = [q_dot ; w_dot];

end