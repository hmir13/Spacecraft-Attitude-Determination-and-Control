close all ; clc ; clear

% Attitude Determination Algorithms
%% TRIAD Algorithm

% body frame measurements
b1 = [0.2500
      0.5177
      0.8311];
b2 = [0.8479
      0.5040
      0.2018];

% reference frame measurements
r1 = [0.5637
      0.3054
      0.7674];
r2 = [0.2569
      0.9337
      0.2495];

% determine attitude matrix
w1 = b1;
c1 = cross(b1, b2);
w2 = c1 / norm(c1);
w3 = cross(w1, w2);

v1 = r1;
c2 = cross(r1, r2);
v2 = c2 / norm(c2);
v3 = cross(v1, v2);

A = [w1 w2 w3] * [v1 v2 v3]';
disp('Attitude Matrix/DCM:');
disp(A);

%% Quaternion Solution of Wahbas Problem: Davenportas q Method
close all ; clc ; clear

bvecs = readmatrix("bvec_Meas.xls");
rvecs = readmatrix("rvec_Refs.xls");
a = 0.01 * ones(1, size(bvecs,2));

B = zeros(3);
for i=1:size(bvecs,2)
    B = B + a(i) * bvecs(:, i) * rvecs(:, i)';
end
z = [B(2,3) - B(3,2)
     B(3,1) - B(1,3)
     B(1,2) - B(2,1)];
K = [B + B' - trace(B) * eye(3) z
     z'                         trace(B)];

% find maximum eigenvalue and respective eigenvector
[eigvecs, eigvals] = eig(K);
eigvals = diag(eigvals);

[max_eigval, idx] = max(eigvals);
max_eigvec = eigvecs(:, idx);

% find optimal quaternion and corresponding attitude matrix
q_opt = max_eigvec;
disp('Estimated Optimal Quaternion:');
disp(q_opt);

A = [q_opt(1)^2 - q_opt(2)^2 - q_opt(3)^2 + q_opt(4)^2, 2 * (q_opt(1) * q_opt(2) + q_opt(3) * q_opt(4))   , 2 * (q_opt(1) * q_opt(3) - q_opt(2) * q_opt(4))   ;
     2 * (q_opt(2) * q_opt(1) - q_opt(3) * q_opt(4))  , -q_opt(1)^2 + q_opt(2)^2 - q_opt(3)^2 + q_opt(4)^2, 2 * (q_opt(2) * q_opt(3) + q_opt(1) * q_opt(4))   ;
     2 * (q_opt(3) * q_opt(1) + q_opt(2) * q_opt(4))  , 2 * (q_opt(3) * q_opt(2) - q_opt(1) * q_opt(4))   , -q_opt(1)^2 - q_opt(2)^2 + q_opt(3)^2 + q_opt(4)^2];
disp('Corresponding Attitude Matrix/DCM:');
disp(A);

% check if attitude matrix is correct
valid = true;
tol = 1e-1;
for i=1:size(bvecs,2)
    if (abs(norm(A*rvecs(:,i)) - norm(bvecs(:,i))) > tol)
        valid = false;
    end
end
if (valid)
    disp("Attitude Matrix/DCM is valid (A*r_i ~= b_i)")
else
    disp("Attitude Matrix/DCM is not valid (A*r_i != b_i)")
end

