%% clean workspace
clc; clear; close all;

%% Filter
% time
dt = 0.01;
tf = 20; 
t = [0:dt:tf]';
m = length(t);

% system parameters
sigw = 0.5624;%3.51e-2; % process noise
sigv = 4.224e-4; % measurement noise
m_b = 410; % sprung mass [kg]
m_w = 45; % unsprung mass [kg]
k_t = 1.83e5; % tire stiffness [N/m]
k_s = 2e4; % suspension spring stiffness [N/m]
c = 2e3; % suspension damping [Ns/m]

% model simulation and measurements
A = [0 0 1 -1; 0 0 0 1;...
    -k_s/m_b 0 -c/m_b c/m_b; k_s/m_w -k_t/m_w c/m_w -c/m_w];
B = [0; -1; 0; 0];
C = [-k_s/m_b 0 -c/m_b c/m_b; k_s/m_w -k_t/m_w c/m_w -c/m_w];
sys = ss(A,B,C,0);
wk = sqrt(sigw)*randn(m,1);
[y,tout,xt] = lsim(sys,wk,t);
ym = y + sqrt(sigv)*randn(m,2);

% discrete-time system
Q = sigw^2/10*eye(4); % process noise covariance 
R = sigv^2*eye(2); % measurement noise covariance 
sysd = c2d(sys,dt,'zoh');
phi = sysd.A;
gam = sysd.B;
H = sysd.C;

% Initial Covariance
p = eye(4); %sigv^2*
pcov = zeros(m,4);
pcov(1,:) = diag(p);

% Initial Condition and H Matrix
x0 = [0; 0; 0; 0];
xe = zeros(m,4);
Ks1 = zeros(m,4);
Ks2 = zeros(m,4);
xe(1,:) = x0';
x = x0;

% Main Loop
for i = 1:m-1
    % Kalman gain
    K = p*H'*inv(H*p*H' + R);

    % update
    x = x + K*(ym(i) - H*x);
    p = [eye(4) - K*H]*p;

    % propagate
    x = phi*x + gam*wk(i);
    p = phi*p*phi' + Q;

    % store state and covariance 
    xe(i+1,:) = x';
    pcov(i+1,:) = diag(p)';
    Ks1(i,:) = K(:,1)';
    Ks2(i,:) = K(:,2)';
end
% 3-sigma bound
sig3 = pcov.^(0.5)*3;

% plot 
% plot output vs measurement
figure;
subplot(2,1,1);
plot(t,y(:,1));
title('Actual')
grid;
xlabel('Time (s)');
ylabel('Sprung Mass Acceleration (m/s/s)')
subplot(2,1,2); 
plot(t,ym(:,1));
title('Measurement')
grid;
xlabel('Time (s)');
ylabel('Sprung Mass Acceleration (m/s/s)')

figure;
subplot(2,1,1);
plot(t,y(:,2));
title('Actual')
grid;
xlabel('Time (s)');
ylabel('Unsprung Mass Acceleration (m/s/s)')
subplot(2,1,2); 
plot(t,ym(:,2));
title('Measurement')
grid;
xlabel('Time (s)');
ylabel('Unsprung Mass Acceleration (m/s/s)')

% plot state errors with sigma bounds
figure;
plot(t,[sig3(:,1) xe(:,1)-xt(:,1) -sig3(:,1)])
xlabel('Time (s)');
ylabel('Rattle Space Error');
figure;
plot(t,[sig3(:,2) xe(:,2)-xt(:,2) -sig3(:,2)])
xlabel('Time (s)');
ylabel('Tire Deflection Error');
figure;
plot(t,[sig3(:,3) xe(:,3)-xt(:,3) -sig3(:,3)])
xlabel('Time (s)');
ylabel('Sprung Mass Velocity Error');
figure;
plot(t,[sig3(:,4) xe(:,4)-xt(:,4) -sig3(:,4)])
xlabel('Time (s)');
ylabel('Unsprung Mass Velocity Error');

% plot states
figure;
subplot(2,1,1);
plot(t,xe(:,1));
title('Estimate')
grid;
xlabel('Time (s)');
ylabel('Rattle Space (m)')
subplot(2,1,2); 
plot(t,xt(:,1));
title('Actual')
grid;
xlabel('Time (s)');
ylabel('Rattle Space (m)')

figure;
subplot(2,1,1);
plot(t,xe(:,2));
title('Estimate')
grid;
xlabel('Time (s)');
ylabel('Tire Deflection (m)')
subplot(2,1,2); 
plot(t,xt(:,2));
title('Actual')
grid;
xlabel('Time (s)');
ylabel('Tire Deflection (m)')

figure;
subplot(2,1,1);
plot(t,xe(:,3));
title('Estimate')
grid;
xlabel('Time (s)');
ylabel('Sprung Mass Velocity (m/s)')
subplot(2,1,2); 
plot(t,xt(:,3));
title('Actual')
grid;
xlabel('Time (s)');
ylabel('Sprung Mass Velocity (m/s)')

figure;
subplot(2,1,1);
plot(t,xe(:,4));
title('Estimate')
grid;
xlabel('Time (s)');
ylabel('Unsprung Mass Velocity (m/s)')
subplot(2,1,2); 
plot(t,xt(:,4));
title('Actual')
grid;
xlabel('Time (s)');
ylabel('Unsprung Mass Velocity (m/s)')

%plot Kalman gains
figure;
plot(t,[Ks1]);grid
xlabel('Time (s)');
ylabel('Kalman Gain Column 1')
figure;
plot(t,[Ks2]);grid
xlabel('Time (s)');
ylabel('Kalman Gain Column 2')

%% Monte Carlo for DKF
n_sim = 100;
for n = 1:n_sim
    % time
    dt = 0.01;
    tf = 20;
    t = [0:dt:tf]';
    m = length(t);

    % system parameters
    sigw = 0.5624;%3.51e-2; % process noise
    sigv = 4.224e-4; % measurement noise
    m_b = 410; % sprung mass [kg]
    m_w = 45; % unsprung mass [kg]
    k_t = 1.83e5; % tire stiffness [N/m]
    k_s = 2e4; % suspension spring stiffness [N/m]
    c = 2e3; % suspension damping [Ns/m]

    % model simulation and measurements
    A = [0 0 1 -1; 0 0 0 1;...
        -k_s/m_b 0 -c/m_b c/m_b; k_s/m_w -k_t/m_w c/m_w -c/m_w];
    B = [0; -1; 0; 0];
    C = [-k_s/m_b 0 -c/m_b c/m_b; k_s/m_w -k_t/m_w c/m_w -c/m_w];
    sys = ss(A,B,C,0);
    wk = sqrt(sigw)*randn(m,1);
    [y,tout,xt] = lsim(sys,wk,t);
    ym = y + sqrt(sigv)*randn(m,2);

    % discrete-time system
    Q = sigw^2/10*eye(4); % process noise covariance
    R = sigv^2*eye(2); % measurement noise covariance
    sysd = c2d(sys,dt,'zoh');
    phi = sysd.A;
    gam = sysd.B;
    H = sysd.C;

    % Initial Covariance
    p = eye(4); %sigv^2*
    pcov = zeros(m,4);
    pcov(1,:) = diag(p);

    % Initial Condition and H Matrix
    x0 = [0; 0; 0; 0];
    xe = zeros(m,4);
    Ks1 = zeros(m,4);
    Ks2 = zeros(m,4);
    xe(1,:) = x0';
    x = x0;

    % Main Loop
    for i = 1:m-1
        % Kalman gain
        K = p*H'*inv(H*p*H' + R);

        % update
        x = x + K*(ym(i,1) - H*x);
        p = [eye(4) - K*H]*p;

        % propagate
        x = phi*x + gam*wk(i);
        p = phi*p*phi' + Q;

        % store state and covariance
        xe(i+1,:) = x';
        pcov(i+1,:) = diag(p)';
        Ks1(i,:) = K(:,1)';
        Ks2(i,:) = K(:,2)';
    end
    % 3-sigma bound
    sig3 = pcov.^(0.5)*3;

    % calculate number of times estimate exceeds 3 sig boundary
    est_exceed1 = [];
    est_exceed2 = [];
    est_exceed3 = [];
    est_exceed4 = [];

    for in = 1:m
        if xe(in,1)-xt(in,1) > sig3(in,1) || xe(in,1)-xt(in,1) < -sig3(in,1)
            est_exceed1(in) = 1;
        else
            est_exceed1(in) = 0;
        end

        if xe(in,2)-xt(in,2) > sig3(in,2) || xe(in,2)-xt(in,2) < -sig3(in,2)
            est_exceed2(in) = 1;
        else
            est_exceed2(in) = 0;
        end

        if xe(in,3)-xt(in,3) > sig3(in,3) || xe(in,3)-xt(in,3) < -sig3(in,3)
            est_exceed3(in) = 1;
        else
            est_exceed3(in) = 0;
        end

        if xe(in,4)-xt(in,4) > sig3(in,4) || xe(in,4)-xt(in,4) < -sig3(in,4)
            est_exceed4(in) = 1;
        else
            est_exceed4(in) = 0;
        end
    end
    exceed_1(n) = sum(est_exceed1);
    exceed_2(n) = sum(est_exceed2);
    exceed_3(n) = sum(est_exceed3);
    exceed_4(n) = sum(est_exceed4);

    % plot state errors with sigma bounds
    figure(101);
    plot(t,[sig3(:,1) xe(:,1)-xt(:,1) -sig3(:,1)])
    xlabel('Time (s)');
    ylabel('Rattle Space Error');
    hold on
    figure(102);
    plot(t,[sig3(:,2) xe(:,2)-xt(:,2) -sig3(:,2)])
    xlabel('Time (s)');
    ylabel('Tire Deflection Error');
    hold on
    figure(103);
    plot(t,[sig3(:,3) xe(:,3)-xt(:,3) -sig3(:,3)])
    xlabel('Time (s)');
    ylabel('Sprung Mass Velocity Error');
    hold on
    figure(104);
    plot(t,[sig3(:,4) xe(:,4)-xt(:,4) -sig3(:,4)])
    xlabel('Time (s)');
    ylabel('Unsprung Mass Velocity Error');
    hold on
end

avg_exceed1 = mean(exceed_1);
avg_exceed2 = mean(exceed_2);
avg_exceed3 = mean(exceed_3);
avg_exceed4 = mean(exceed_4);
disp(avg_exceed1)
disp(avg_exceed2)
disp(avg_exceed3)
disp(avg_exceed4)
[row1,~] = find(exceed_1>0);
[row2,~] = find(exceed_2>0);
[row3,~] = find(exceed_3>0);
[row4,~] = find(exceed_4>0);
disp(length(row1)/100)
disp(length(row2)/100)
disp(length(row3)/100)
disp(length(row4)/100)

%% parameter estimation with EKF (just c)
% time
dt = 0.01;
tf = 20; 
t = [0:dt:tf]';
m = length(t);

% system parameters
sigw = 0.5624;%3.51e-2; % process noise
sigv = 4.224e-4; % measurement noise
m_b = 410; % sprung mass [kg]
m_w = 45; % unsprung mass [kg]
k_t = 1.83e5; % tire stiffness [N/m]
k_s = 2e4; % suspension spring stiffness [N/m]
c = 2e3; % suspension damping [Ns/m]

% State and Initialize
wk = sqrt(sigw)*randn(m,1);
noise = sqrt(sigv)*randn(m,2);
C = [-k_s/m_b 0 -c/m_b c/m_b; k_s/m_w -k_t/m_w c/m_w -c/m_w];
h=[-k_s/m_b 0 -c/m_b c/m_b 0; k_s/m_w -k_t/m_w c/m_w -c/m_w 0]; 
R = sigv^2*eye(2); % measurement noise covariance
x=zeros(m,4); ym=zeros(m,2);
x0=[0; 0; 0; 0]; x(1,:)=x0'; xe(1,:)=[0;0;0;0;2.5e3]'; 
p0=100*eye(5); p=p0; p_cov = zeros(m,5); p_cov(1,:)=diag(p0)'; 

% Truth and Measurements
for i=1:m-1
    f1=dt*susp(x(i,:),c,k_t,k_s,m_b,m_w,wk(i));
    f2=dt*susp(x(i,:)+0.5*f1',c,k_t,k_s,m_b,m_w,wk(i));
    f3=dt*susp(x(i,:)+0.5*f2',c,k_t,k_s,m_b,m_w,wk(i));
    f4=dt*susp(x(i,:)+f3',c,k_t,k_s,m_b,m_w,wk(i));
    x(i+1,:)=x(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');
    ymi = C*x(i,:)' + noise(i,:)';
    ym(i,:) = ymi';
end

% Extended Kalman Filter
for i=1:m-1
    % Kalman Update
    gain=p*h'*inv(h*p*h'+R);
    p=(eye(5)-gain*h)*p;
    xe(i,:)=xe(i,:)+(gain*(ym(i,1)-h*xe(i,:)'))';

    % Propagation
    f1=dt*susp_id(xe(i,:),k_t,k_s,m_b,m_w,wk(i));
    f2=dt*susp_id(xe(i,:)+0.5*f1',k_t,k_s,m_b,m_w,wk(i));
    f3=dt*susp_id(xe(i,:)+0.5*f2',k_t,k_s,m_b,m_w,wk(i));
    f4=dt*susp_id(xe(i,:)+f3',k_t,k_s,m_b,m_w,wk(i));
    xe(i+1,:)=xe(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');
    F=[0 0 1 -1 0
       0 0 0  1 0 
       -k_s/m_b 0 -xe(i,5)/m_b xe(i,5)/m_b -xe(i,3)/m_b+xe(i,4)/m_b
       k_s/m_w -k_t/m_w xe(i,5)/m_w -xe(i,5)/m_w xe(i,3)/m_w-xe(i,4)/m_w
        0 0 0 0 0];
    G=[0;0;0;0;1];
    phi=c2d(F,G,dt);
    q=0.001*[0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 1 0;0 0 0 0 1];%eye(5)
    p=phi*p*phi'+q*dt;
    p_cov(i+1,:)=diag(p)';
end

% 3-Sigma Outlier
sig3_id=(p_cov.^(0.5))*3;

% plot states
figure;
subplot(2,1,1);
plot(t,xe(:,1));
title('Estimate')
grid;
xlabel('Time (s)');
ylabel('Rattle Space (m)')
subplot(2,1,2); 
plot(t,x(:,1));
title('Actual')
grid;
xlabel('Time (s)');
ylabel('Rattle Space (m)')

figure;
subplot(2,1,1);
plot(t,xe(:,2));
title('Estimate')
grid;
xlabel('Time (s)');
ylabel('Tire Deflection (m)')
subplot(2,1,2); 
plot(t,x(:,2));
title('Actual')
grid;
xlabel('Time (s)');
ylabel('Tire Deflection (m)')

figure;
subplot(2,1,1);
plot(t,xe(:,3));
title('Estimate')
grid;
xlabel('Time (s)');
ylabel('Sprung Mass Velocity (m/s)')
subplot(2,1,2); 
plot(t,x(:,3));
title('Actual')
grid;
xlabel('Time (s)');
ylabel('Sprung Mass Velocity (m/s)')

figure;
subplot(2,1,1);
plot(t,xe(:,4));
title('Estimate')
grid;
xlabel('Time (s)');
ylabel('Unsprung Mass Velocity (m/s)')
subplot(2,1,2); 
plot(t,x(:,4));
title('Actual')
grid;
xlabel('Time (s)');
ylabel('Unsprung Mass Velocity (m/s)')

figure;
plot(t,xe(:,5),t,c*ones(m));
title('Estimated Damping Coefficient vs Actual')
grid;
xlabel('Time (s)');
ylabel('Damping Coefficient [Ns/m]')

% plot errors
figure;
plot(t,xe(:,1)-x(:,1),t,sig3_id(:,1),t,-sig3_id(:,1))
xlabel('Time (s)')
ylabel('Rattle Space Error')
figure;
plot(t,xe(:,2)-x(:,2),t,sig3_id(:,2),t,-sig3_id(:,2))
xlabel('Time (s)')
ylabel('Tire Deflection Error')
figure;
plot(t,xe(:,3)-x(:,3),t,sig3_id(:,3),t,-sig3_id(:,3))
xlabel('Time (s)')
ylabel('Sprung Mass Velocity Error')
figure;
plot(t,xe(:,4)-x(:,4),t,sig3_id(:,4),t,-sig3_id(:,4))
xlabel('Time (s)')
ylabel('Unsprung Mass Velocity Error')
figure;
plot(t,xe(:,5)-c,t,sig3_id(:,5),t,-sig3_id(:,5))
xlabel('Time (s)')
ylabel('Damping Coefficient Error')

%% parameter estimation with EKF (c, k_t, k_s)
% time
dt = 0.01;
tf = 20; 
t = [0:dt:tf]';
m = length(t);

% system parameters
sigw = 0.5624;%3.51e-2; % process noise
sigv = 4.224e-4; % measurement noise
m_b = 410; % sprung mass [kg]
m_w = 45; % unsprung mass [kg]
k_t = 1.83e5; % tire stiffness [N/m]
k_s = 2e4; % suspension spring stiffness [N/m]
c = 2e3; % suspension damping [Ns/m]

% State and Initialize
wk = sqrt(sigw)*randn(m,1);
noise = sqrt(sigv)*randn(m,2);
C = [-k_s/m_b 0 -c/m_b c/m_b; k_s/m_w -k_t/m_w c/m_w -c/m_w];
h=[-k_s/m_b 0 -c/m_b c/m_b 0 0 0; k_s/m_w -k_t/m_w c/m_w -c/m_w 0 0 0]; 
R = sigv^2*eye(2); % measurement noise covariance
x=zeros(m,4); ym=zeros(m,2);
x0=[0; 0; 0; 0]; x(1,:)=x0'; xe(1,:)=[0;0;0;0;2.5e3;2.05e4;1.835e5]'; 
p0=100*eye(7); p=p0; p_cov = zeros(m,7); p_cov(1,:)=diag(p0)'; 
% best k_t results: p0 = 100*eye(7), q = 0.0001
% best k_s, c results: p0 = 100*eye(7), q = 0.001
% Truth and Measurements
for i=1:m-1
    f1=dt*susp(x(i,:),c,k_t,k_s,m_b,m_w,wk(i));
    f2=dt*susp(x(i,:)+0.5*f1',c,k_t,k_s,m_b,m_w,wk(i));
    f3=dt*susp(x(i,:)+0.5*f2',c,k_t,k_s,m_b,m_w,wk(i));
    f4=dt*susp(x(i,:)+f3',c,k_t,k_s,m_b,m_w,wk(i));
    x(i+1,:)=x(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');
    ymi = C*x(i,:)' + noise(i,:)';
    ym(i,:) = ymi';
end

% Extended Kalman Filter
for i=1:m-1
    % Kalman Update
    gain=p*h'*inv(h*p*h'+R);
    p=(eye(7)-gain*h)*p;
    xe(i,:)=xe(i,:)+(gain*(ym(i,1)-h*xe(i,:)'))';

    % Propagation
    f1=dt*susp_id2(xe(i,:),m_b,m_w,wk(i));
    f2=dt*susp_id2(xe(i,:)+0.5*f1',m_b,m_w,wk(i));
    f3=dt*susp_id2(xe(i,:)+0.5*f2',m_b,m_w,wk(i));
    f4=dt*susp_id2(xe(i,:)+f3',m_b,m_w,wk(i));
    xe(i+1,:)=xe(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');
    F=[0 0 1 -1 0 0 0
       0 0 0  1 0 0 0
       -xe(i,6)/m_b 0 -xe(i,5)/m_b xe(i,5)/m_b -xe(i,3)/m_b+xe(i,4)/m_b -xe(i,1)/m_b 0
       xe(i,6)/m_w -xe(i,7)/m_w xe(i,5)/m_w -xe(i,5)/m_w xe(i,3)/m_w-xe(i,4)/m_w xe(i,1)/m_w -xe(i,2)/m_w
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0];
    G=[0;0;0;0;1;1;1];
    phi=c2d(F,G,dt);
    q=0.001*[0 0 0 0 0 0 0;0 0 0 0 0 0 0;0 0 0 0 0 0 0;
        0 0 0 1 0 0 0;0 0 0 0 1 0 0;0 0 0 0 0 1 0;0 0 0 0 0 0 1];
    p=phi*p*phi'+q*dt;
    p_cov(i+1,:)=diag(p)';
end

% 3-Sigma Outlier
sig3_id=(p_cov.^(0.5))*3;

% plot states
figure;
subplot(2,1,1);
plot(t,xe(:,1));
title('Estimate')
grid;
xlabel('Time (s)');
ylabel('Rattle Space (m)')
subplot(2,1,2); 
plot(t,x(:,1));
title('Actual')
grid;
xlabel('Time (s)');
ylabel('Rattle Space (m)')

figure;
subplot(2,1,1);
plot(t,xe(:,2));
title('Estimate')
grid;
xlabel('Time (s)');
ylabel('Tire Deflection (m)')
subplot(2,1,2); 
plot(t,x(:,2));
title('Actual')
grid;
xlabel('Time (s)');
ylabel('Tire Deflection (m)')

figure;
subplot(2,1,1);
plot(t,xe(:,3));
title('Estimate')
grid;
xlabel('Time (s)');
ylabel('Sprung Mass Velocity (m/s)')
subplot(2,1,2); 
plot(t,x(:,3));
title('Actual')
grid;
xlabel('Time (s)');
ylabel('Sprung Mass Velocity (m/s)')

figure;
subplot(2,1,1);
plot(t,xe(:,4));
title('Estimate')
grid;
xlabel('Time (s)');
ylabel('Unsprung Mass Velocity (m/s)')
subplot(2,1,2); 
plot(t,x(:,4));
title('Actual')
grid;
xlabel('Time (s)');
ylabel('Unsprung Mass Velocity (m/s)')

figure;
plot(t,xe(:,5),t,c*ones(m));
title('Estimated Damping Coefficient vs Actual')
grid;
xlabel('Time (s)');
ylabel('Damping Coefficient')

figure;
plot(t,xe(:,6),t,k_s*ones(m));
title('Estimated Suspension Spring Stiffness vs Actual')
grid;
xlabel('Time (s)');
ylabel('Suspension Spring Stiffness [N/m]')

figure;
plot(t,xe(:,7),t,k_s*ones(m));
title('Estimated Tire Stiffness vs Actual')
grid;
xlabel('Time (s)');
ylabel('Tire Stiffness [N/m]')

% plot error
figure;
plot(t,xe(:,1)-x(:,1),t,sig3_id(:,1),t,-sig3_id(:,1))
xlabel('Time (s)')
ylabel('Rattle Space Error')
figure;
plot(t,xe(:,2)-x(:,2),t,sig3_id(:,2),t,-sig3_id(:,2))
xlabel('Time (s)')
ylabel('Tire Deflection Error')
figure;
plot(t,xe(:,3)-x(:,3),t,sig3_id(:,3),t,-sig3_id(:,3))
xlabel('Time (s)')
ylabel('Sprung Mass Velocity Error')
figure;
plot(t,xe(:,4)-x(:,4),t,sig3_id(:,4),t,-sig3_id(:,4))
xlabel('Time (s)')
ylabel('Unsprung Mass Velocity Error')
figure;
plot(t,xe(:,5)-c,t,sig3_id(:,5),t,-sig3_id(:,5))
xlabel('Time (s)')
ylabel('Damping Coefficient Error')
figure;
plot(t,xe(:,6)-k_s,t,sig3_id(:,6),t,-sig3_id(:,6))
xlabel('Time (s)')
ylabel('Suspension Spring Stiffness Error')
figure;
plot(t,xe(:,7)-k_t,t,sig3_id(:,7),t,-sig3_id(:,7))
xlabel('Time (s)')
ylabel('Tire Stiffness Error')
