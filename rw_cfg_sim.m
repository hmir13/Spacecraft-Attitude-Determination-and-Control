close all ; clc ; clear

% Reaction Wheel Configurations Simulation
%% Quaternion-based Attitude Control Using Reaction Wheels
I = [6400 -76.4 -25.6;
     -76.4 4730 -40;
     -25.6 -40 8160];

q_0 = (1/sqrt(2)) * [1 0 0 1]';
w_0 = [0.01 0.01 0.01]';
h_w_0 = [0 0 0]';

q_c = [0 0 0 1]';
w_d = [0 0 0]';

k_p = 10;
k_d = 150;

x_0 = [q_0 ; w_0 ; h_w_0];

t = 20 * 60;
tspan = 0:1:t;
[t, out] = ode45(@(t,x) controller(t,x,I,q_c,k_p,k_d), tspan, x_0);

q = out(:,1:4);
w = out(:, 5:7);
h = out(:, 8:10);

% Quaternion
figure;
plot(t, q);
title('Quaternion Time History q(t)');
xlabel('Time (s)');
ylabel('Quaternion components');
ylim([-0.2,1.2]);
legend('q_1', 'q_2', 'q_3', 'q_4');

% Angular velocity
figure;
plot(t, w);
title('Angular Velocity Time History \omega(t)');
xlabel('Time (s)');
ylabel('Angular velocity (rad/s)');
legend('\omega_1', '\omega_2', '\omega_3');

% Wheel momentum
figure;
plot(t, h);
xlabel('Time (s)');
ylabel('Wheel Momentum (kg-m^2/s)');
legend('h_1', 'h_2', 'h_3');
title('Time History of Wheel Momentum h_w(t)');

%% Problem 2 Redundant Reaction Wheel Configurations 
% Pyramid config
w_p_1 = (1/sqrt(2)) * [1 1 0]';
w_p_2 = (1/sqrt(2)) * [-1 1 0]';
w_p_3 = (1/sqrt(2)) * [0 1 1]';
w_p_4 = (1/sqrt(2)) * [0 1 -1]';

w_p = [w_p_1 w_p_2 w_p_3 w_p_4];
w_p_inv = pinv(w_p);
h_p = w_p_inv * h';

figure;
plot(t, h_p);
title('Angular Momentum for Pyramid Configuration');
xlabel('Time (s)');
ylabel('Angular Momentum (kgm/s)');
legend('h_w1', 'h_w2', 'h_w3', 'h_w4');

% NASA config
w_n_1 = [1 0 0]';
w_n_2 = [0 1 0]';
w_n_3 = [0 0 1]';
w_n_4 = (1/sqrt(3)) * [1 1 1]';

w_n = [w_n_1 w_n_2 w_n_3 w_n_4];
w_n_inv = pinv(w_n);
h_n = w_n_inv * h';

figure;
plot(t, h_n);
title('Angular Momentum for NASA Standard Configuration');
xlabel('Time (s)');
ylabel('Angular Momentum (kg-m/s)');
legend('h_w1', 'h_w2', 'h_w3', 'h_w4');

figure;
norm_p = sqrt(sum(h_p.^2, 1));
norm_n = sqrt(sum(h_n.^2, 1));
plot(t, norm_p, t, norm_n);
title('Norm of the Wheel Momentum');
xlabel('Time (s)');
ylabel('Norm of Momentum (kg-m/s)');
legend('Pyramid Config', 'NASA Standard Config'); 

steady_p = norm_p(end);  
steady_n = norm_n(end);        

disp(['Steady-state wheel momentum norm for Pyramid configuration: ', num2str(steady_p)]);
disp(['Steady-state wheel momentum norm for NASA configuration: ', num2str(steady_n)]);

function dq = quat_error(q,qc)
    q1 = qc(1); 
    q2 = qc(2); 
    q3 = qc(3); 
    q4 = qc(4);

    A = [q4, -q3, q2;
          q3, q4, -q1;
          -q2, q1, q4;
          -q1, -q2, -q3];

    dq13 = A' * q;
    dq4 = qc' * q;

    dq = [dq13; dq4];
end

function dxdt = controller(~,x,I,q_c,k_p,k_d)
    q = x(1:4);
    w = x(5:7);
    h = x(8:10);

    dq = quat_error(q,q_c);

    Ls = -k_p * sign(dq(4)) * dq(1:3) - k_d * w;
    Lw = -Ls;

    q1 = q(1); 
    q2 = q(2); 
    q3 = q(3); 
    q4 = q(4);

    Aq = [q4, -q3, q2;
          q3, q4, -q1;
          -q2, q1, q4;
          -q1, -q2, -q3];
    
    q_dot = (1/2) * Aq * w;
    w_dot = inv(I) *(Ls - cross(w,I*w));
    h_dot = Lw - cross(w,h);

    dxdt = [q_dot ; w_dot ; h_dot];

end
