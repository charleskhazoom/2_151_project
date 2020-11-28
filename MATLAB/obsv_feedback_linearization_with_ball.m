function dz_hat = obsv_feedback_linearization_with_ball(y,x_hat,u,p,p_obsv)

A  = p_obsv.A;
B  = p_obsv.B;
C  = p_obsv.C;
L  = p_obsv.L;
zf = p_obsv.zf;

% undo feedback linearization
M = A_arm(x_hat, p); % Mass matrix
V = Corr_arm(x_hat, p); % Coriolis matrix
G = Grav_arm(x_hat, p); % Gravity matrix
w = (M^-1)*(u -V - G);

% compute observer dynamics to be integrated
dz_hat = (A - L*C)*(x_hat-zf) + B*w + L*(y-C*zf);