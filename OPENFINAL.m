%% Complete Cruise-Phase Optimization of Commercial Aircraft
clear all; clc;
global W_x_func W_y_func g_func dg_dx_func dg_dy_func ...
    dWx_dx_func dWx_dy_func dWy_dx_func dWy_dy_func ...
    C_D_func dC_D_dv_func C_s_func dC_s_dv_func
%% User Inputs for Domain Dimensions
x_f = input('Enter the x-dimension of the flight box (x_f[m]); e.g., 3*10^6: ');
y_f = input('Enter the y-dimension of the flight box (y_f[m]); e.g., 3*10^6: ');
font_size = 12;
line_width = 1;
%% Choose Input Method for Scattered Points
disp('Select input method for flight-sensitive area points:');
disp('1: Interactive clicking (ginput)');
disp('2: Load from file (CSV, TXT, or XLSX)');
method = input('Enter 1 or 2: ');
switch method
    case 1
        figure (1)
        axis([0 x_f 0 y_f]);
        grid on
        title('Click to select points; press ENTER when finished');
        xlabel('x'); ylabel('y');
        hold on;
        Xs = [];
        while true
            [x, y, button] = ginput(1);  % Get one point
            if isempty(button)
                break;  % Stop if ENTER is pressed
            end
            Xs = [Xs; x, y];  % Append the new point
            plot(x, y, 'k.', 'MarkerSize', 25);  % Display the clicked point
        end
        hold off;       
    case 2
        [fileName, pathName] = uigetfile({'*.csv;*.txt;*.xlsx', 'Data Files (*.csv, *.txt, *.xlsx)'}, ...
            'Select the data file containing scattered points');
        if isequal(fileName, 0)
            error('No file selected. Exiting.');
        else
            fullFileName = fullfile(pathName, fileName);
            Xs = readmatrix(fullFileName);
        end
        
    otherwise
        error('Invalid input method. Please restart and choose either 1 or 2.');
end
Xs(:,1) = min(max(Xs(:,1), 0), x_f);
Xs(:,2) = min(max(Xs(:,2), 0), y_f);
%% Ask for the Number of Clusters (Ellipses)
K = input('Enter the number of clusters (ellipses): ');
% Ensure we have enough points for clustering: N ≥ 4K
N = size(Xs,1);
if N < 4*K
    error('Insufficient points! You provided %d points, but at least %d points are required for %d clusters. Please add more points or reduce the number of clusters.', N, 4*K, K);
end

%% Automatically Determine the Scaling Factor 'k'
% For a 85% confidence ellipse in 2D:
k = sqrt(chi2inv(0.85, 2));
%% Clustering Using k-Means
[idx, centroids] = kmeans(Xs, K);
%% Initialize Arrays for Ellipse Parameters
a     = zeros(K, 1);  % Semi-axis lengths
b     = zeros(K, 1);
alpha = zeros(K, 1);  % Orientation angles (in radians)
%% Plot the Scattered Points and Cluster Centroids
figure (2)
hold on;
plot(Xs(:,1), Xs(:,2), 'k.', 'MarkerSize', 15);
set(gca, 'FontSize', font_size, 'TickLabelInterpreter', 'latex', 'GridColor', 'r');
plot(centroids(:,1), centroids(:,2), 'r.', 'MarkerSize', 20, 'LineWidth', 2);
set(gca, 'FontSize', font_size, 'TickLabelInterpreter', 'latex', 'GridColor', 'r');
%% Compute Ellipse Parameters for Each Cluster
for i = 1:K
    % Extract points in the i-th cluster
    cluster_points = Xs(idx == i, :);    
    % Compute the covariance matrix
    Sigma = cov(cluster_points);   
    % Eigenvalue decomposition of the covariance matrix
    [V, D] = eig(Sigma);
    % Eigenvalues (variances along principal directions)
    lambda1 = D(1,1);
    lambda2 = D(2,2);    
    % Compute semi-axes using the scaling factor 'k'
    a(i) = k * sqrt(lambda1);
    b(i) = k * sqrt(lambda2);   
    % Determine orientation angle from the first eigenvector
    alpha(i) = atan2(V(2,1), V(1,1));   
    % --- Plot the Ellipse ---
    theta = linspace(0, 2*pi, 100);
    ellipse_x = a(i) * cos(theta);
    ellipse_y = b(i) * sin(theta);    
    % Rotate the ellipse using the computed angle
    R = [cos(alpha(i)) -sin(alpha(i)); sin(alpha(i)) cos(alpha(i))];
    ellipse_rotated = (R * [ellipse_x; ellipse_y])';    
    % Translate the ellipse to the cluster centroid
    ellipse_rotated(:,1) = ellipse_rotated(:,1) + centroids(i,1);
    ellipse_rotated(:,2) = ellipse_rotated(:,2) + centroids(i,2);    
    % Plot the ellipse with a thicker line
    figure (2)
    plot(ellipse_rotated(:,1), ellipse_rotated(:,2), 'LineWidth', 2);
    hold on
end
xlabel('$x(t)$','Interpreter', 'latex', 'FontSize', font_size);
ylabel('$y(t)$', 'Interpreter', 'latex', 'FontSize', font_size);
title('Flight-Sensitive Areas and the Optimal Flight Path', 'Interpreter', 'latex', 'FontSize', font_size);
xlim([-0.2*x_f, 1.2*x_f]);  % Limit x-axis
ylim([-0.2*y_f, 1.2*y_f]);  % Limit y-axis
grid on;
set(gca, 'FontSize', font_size, 'TickLabelInterpreter', 'latex', 'GridColor', 'r');
%% User Inputs for Wind Functions W_x(x,y) and W_y(x,y)
disp('Define the wind functions W_x(x,y) and W_y(x,y).');
wx_str = input('Enter W_x(x,y) as a function of x and y; e.g., 20: ', 's');
wy_str = input('Enter W_y(x,y) as a function of x and y; e.g., -20: ', 's');
% wx_str=0; wy_str=0;
% Convert user input strings into function handles
W_x = str2func(['@(x,y) ' wx_str]);
W_y = str2func(['@(x,y) ' wy_str]);


%% --- Additional Calculations: g(x,y) and its Partial Derivatives ---
% Prompt user for the coefficients c_{s,i} for each ellipse (cluster)
disp('Enter the weights c_{s,i} for each ellipse (as a vector of length K).');
c_s = input(sprintf('Enter weights for each of the %d ellipses (e.g., [0.1 0.2 ...]): ', K));
if length(c_s) ~= K
    error('The number of coefficients must equal the number of clusters (ellipses)!');
end

% Define symbolic variables
syms x y

% Initialize symbolic expression for g(x,y)
g_sym = sym(0);

% Loop over each ellipse (cluster) to define G_{s,i}(x,y) and sum its weighted inverse
for i = 1:K
    % Get the centroid and ellipse parameters for the i-th cluster
    cxi = centroids(i,1);
    cyi = centroids(i,2);
    ai = a(i);
    bi = b(i);
    alphai = alpha(i);
    
    % Define rotated coordinates:
    % (x-cxi) and (y-cyi) are the differences from the centroid.
    G_si = (( (x - cxi)*cos(alphai) + (y - cyi)*sin(alphai) )^2 / (ai^2)) + ...
        (( (x - cxi)*sin(alphai) - (y - cyi)*cos(alphai) )^2 / (bi^2));
    
    syms kappa
    kappa = 0.001;
    % Add the term (c_{s,i} / G_si) to g(x,y)
    g_sym = g_sym + c_s(i)/(G_si+kappa);
end

% Compute the partial derivatives of g(x,y)
dg_dx_sym = simplify(diff(g_sym, x));
dg_dy_sym = simplify(diff(g_sym, y));

%% --- Additional Calculations: Partial Derivatives of Wind Functions ---
% Convert the user-defined wind functions to symbolic expressions
W_x_sym = str2sym(wx_str);
W_y_sym = str2sym(wy_str);

% Compute the partial derivatives
dWx_dx = simplify(diff(W_x_sym, x));
dWx_dy = simplify(diff(W_x_sym, y));
dWy_dx = simplify(diff(W_y_sym, x));
dWy_dy = simplify(diff(W_y_sym, y));

%% --- Converting symbolic expressions to MATLAB function handles ---

% Convert g(x,y) and its derivatives into functions
g_func    = matlabFunction(g_sym, 'vars', [x, y]);
dg_dx_func = matlabFunction(dg_dx_sym, 'vars', [x, y]);
dg_dy_func = matlabFunction(dg_dy_sym, 'vars', [x, y]);

% Convert wind functions and their derivatives into functions
W_x_func   = matlabFunction(W_x_sym, 'vars', [x, y]);
W_y_func   = matlabFunction(W_y_sym, 'vars', [x, y]);
dWx_dx_func = matlabFunction(dWx_dx, 'vars', [x, y]);
dWx_dy_func = matlabFunction(dWx_dy, 'vars', [x, y]);
dWy_dx_func = matlabFunction(dWy_dx, 'vars', [x, y]);
dWy_dy_func = matlabFunction(dWy_dy, 'vars', [x, y]);

%% =================== Begin Model Definitions =========================
%% Atmospheric Model and Cruise Altitude
% Default constants (Table 2)
P0      = 101325;         % Sea-level pressure (Pa)
c0      = 1.4;            % Specific heat ratio
s       = 283.3;          % Wing area (m^2)
T0      = 5e5;            % Reference thrust (N)
C_s0    = 9e-6;           % Fuel consumption constant (kg/N/s)
R_air   = 287.04;         % Specific gas constant (J/kg/K)
g       = 9.81;           % Gravity (m/s^2)
Theta0  = 288.15;         % Sea-level temperature (K)
beta    = 0.0065;         % Temperature lapse rate (K/m)

% Ask user for cruise altitude (in meters)
h_bar = input('Enter the cruise altitude (in meters); e.g., 9000: ');

% Atmospheric Constants
Theta= Theta0 - beta * h_bar;
P= P0 * (Theta/Theta0)^(g/(beta*R_air));
rho= P / (R_air * Theta);

%% Aerodynamic Drag Coefficient Model
% Default constants
% For c(M,v): use k0p
k01 = 0.0067;   k02 = -0.1861;  k03 = 2.2420;   k04 = -6.4350;  k05 = 6.3428;
% For b(M,v): use k1p
k11 = 0.0962;   k12 = -0.7602;  k13 = -1.2870;  k14 = 3.7925;   k15 = -2.7672;
% For a(M,v): use k2p
k21 = -0.1317;  k22 = 1.3427;   k23 = -1.2839;  k24 = 5.0164;   k25 = 0.0000;
% Constant offsets for c, b, a
C_D0i = 0.01322;  C_D1i = -0.00610;  C_D2i = 0.06000;

%% --- Aerodynamic and Thrust Models Setup ---
% Ask the user if they want to use the default models or provide custom formulas
useDefault = input('Do you want to use the default formulas for C_D, C_s, and T_max? (Y/N): ', 's');

% Define symbolic variables (common to both cases)
syms v_sym M_sym m_sym  % v: speed, M: Mach number, m: mass variable for C_D

% dM/dv is common (from the atmospheric model)
dM_sym_dv = 1 / sqrt(c0 * R_air * Theta);

if strcmpi(useDefault, 'Y')
    %% DEFAULT MODELS
    % --- Aerodynamic Drag Coefficient Model: Default Constants
    % Define barK as a function of Mach number (symbolically)
    barK_sym = ((M_sym - 0.4)^2) / sqrt(1 - M_sym^2);
    
    % Default symbolic expressions for a, b, and c (which depend on M only)
    a_sym = C_D2i + k21 * barK_sym + k22 * barK_sym^2 + k23 * barK_sym^3 + k24 * barK_sym^4 + k25 * barK_sym^5;
    b_sym = C_D1i + k11 * barK_sym + k12 * barK_sym^2 + k13 * barK_sym^3 + k14 * barK_sym^4 + k15 * barK_sym^5;
    c_sym = C_D0i + k01 * barK_sym + k02 * barK_sym^2 + k03 * barK_sym^3 + k04 * barK_sym^4 + k05 * barK_sym^5;
    
    % Define A and B (which incorporate a and b) symbolically;
    % Note that these expressions depend on the velocity v_sym explicitly.
    A_sym = (4 * g^2 * a_sym) / (rho^2 * s^2 * v_sym^4);
    B_sym = (2 * g * b_sym) / (rho * s * v_sym^2);
    
    % Define the symbolic drag coefficient:
    % Here, m_sym is an additional variable (e.g., mass) that appears in the quadratic model
    C_D_sym = A_sym * m_sym^2 + B_sym * m_sym + c_sym;
    
    % Compute the derivative of C_D with respect to v using the chain rule.
    % Note: Since a_sym, b_sym, and c_sym depend only on M_sym in the default case,
    % the v-dependence comes only through the explicit v_sym factors and via M_sym (through dM_sym_dv)
    dC_D_dv = (diff(A_sym, v_sym) + diff(A_sym, M_sym)*dM_sym_dv)*m_sym^2 + ...
        (diff(B_sym, v_sym) + diff(B_sym, M_sym)*dM_sym_dv)*m_sym + ...
        (diff(c_sym, v_sym) + diff(c_sym, M_sym)*dM_sym_dv);
    
    % --- Fuel Consumption Model: Default ---
    % Default formula for C_s (assumed to depend on M_sym only in the default)
    C_s_sym = C_s0 * (Theta/Theta0)^(0.5) * (1 + 1.2 * M_sym);
    dC_s_dv = (diff(C_s_sym, v_sym) + diff(C_s_sym, M_sym)*dM_sym_dv);
    
    % --- Thrust Model: Default ---
    T_max_func = @(M) (P * Theta0 / (P0 * Theta)) * T0 * ...
        (1 + ((c0 - 1)/2) * M.^2).^(c0/(c0-1)) .* (1 - 0.49 * sqrt(M));
    
else
    %% CUSTOM MODELS
    % For the drag coefficient, ask the user to provide formulas for a(v,M), b(v,M), and c(v,M)
    disp('Enter your custom formulas for the aerodynamic drag model:');
    a_str = input('Enter a(v,M) as a function of v and M (e.g., ''v.^0 + M''): ', 's');
    b_str = input('Enter b(v,M) as a function of v and M: ', 's');
    c_str = input('Enter c(v,M) as a function of v and M: ', 's');
    
    % Convert the strings to symbolic expressions in v_sym and M_sym
    a_sym = str2sym(a_str);
    b_sym = str2sym(b_str);
    c_sym = str2sym(c_str);
    
    % Define A and B using the user-provided functions (note the explicit v_sym dependence)
    A_sym = (4 * g^2 * a_sym) / (rho^2 * s^2 * v_sym^4);
    B_sym = (2 * g * b_sym) / (rho * s * v_sym^2);
    
    % Construct the symbolic drag coefficient:
    C_D_sym = A_sym * m_sym^2 + B_sym * m_sym + c_sym;
    
    % Compute its derivative with respect to v (again, applying the chain rule for any M_sym dependence)
    dC_D_dv = (diff(A_sym, v_sym) + diff(A_sym, M_sym)*dM_sym_dv)*m_sym^2 + ...
        (diff(B_sym, v_sym) + diff(B_sym, M_sym)*dM_sym_dv)*m_sym + ...
        (diff(c_sym, v_sym) + diff(c_sym, M_sym)*dM_sym_dv);
    
    % For fuel consumption, ask for a custom function C_s(v,M)
    C_s_str = input('Enter your custom formula for C_s(v,M) (e.g., ''1e-6*(v+M)''): ', 's');
    C_s_sym = str2sym(C_s_str);
    dC_s_dv = (diff(C_s_sym, v_sym) + diff(C_s_sym, M_sym)*dM_sym_dv);
    
    % For the thrust model, ask for a custom function T_max(M)
    T_max_str = input('Enter your custom formula for T_max(M) (e.g., ''1e5*(1 - M)''): ', 's');
    T_max_sym = str2sym(T_max_str);
    % Convert the T_max symbolic expression to a function later.
end

%% --- Convert Symbolic Expressions to Function Handles ---
% Convert C_D and its derivative: note that these functions take (v, M, m) as inputs.
C_D_func    = matlabFunction(C_D_sym, 'Vars', [v_sym, M_sym, m_sym]);
dC_D_dv_func = matlabFunction(dC_D_dv, 'Vars', [v_sym, M_sym, m_sym]);

% Convert fuel consumption functions (C_s and its derivative) that depend on v and M.
C_s_func    = matlabFunction(C_s_sym, 'Vars', [v_sym, M_sym]);
dC_s_dv_func = matlabFunction(dC_s_dv, 'Vars', [v_sym, M_sym]);

% If using custom formulas, T_max_sym is defined symbolically; otherwise, T_max_func is already defined.
if strcmpi(useDefault, 'Y')
    % Default T_max_func is already set above.
else
    T_max_func = matlabFunction(T_max_sym, 'Vars', M_sym);
end

%% Optimal Solution
ct = input('Enter the final time weight (ct); e.g., 0.1:');
% cm = input('Enter the final mass weight (cm); e.g., -1:');
cm=-1;
m_0 = input('Enter the initil mass; e.g., 140000:');
% Initial shooting parameters:
x0 = 0;  y0 = 0;
bestD = inf;
q0    = 0;
for i = 1:K
    % ellipse params
    cx = centroids(i,1);  cy = centroids(i,2);
    ai = a(i);            bi = b(i);
    alph = alpha(i);
    
    % 1) transform to ellipse frame
    ui =  (x0-cx)*cos(alph) + (y0-cy)*sin(alph);
    vi = -(x0-cx)*sin(alph) + (y0-cy)*cos(alph);
    
    % 2) project radially to boundary
    ti = atan2( vi/bi, ui/ai );
    
    % boundary point in world coords
    ub = ai*cos(ti);  vb = bi*sin(ti);
    xb = cx + ub*cos(alph) - vb*sin(alph);
    yb = cy + ub*sin(alph) + vb*cos(alph);
    
    % 3) check distance
    d = hypot(xb - x0, yb - y0);
    if d < bestD
        bestD = d;
        % 4) compute raw tangent derivative dy/dx
        dxdt = -ai*sin(ti)*cos(alph) - bi*cos(ti)*sin(alph);
        dydt = -ai*sin(ti)*sin(alph) + bi*cos(ti)*cos(alph);
        
        q0 = abs(dxdt / dydt);
    end
end
Init = zeros(3,1);
Init(2) = q0;        % initial q(0)
Init(3) = (sqrt(x_f^2+y_f^2))/300;     % final time, t_f
Init(1)=-(1+ct)*(1/sqrt(1+Init(2)^2))/200;  % initial lambda_x(0)
Init(1)=-0.1;

N_p = 150;

tic
% Define the constraint and objective functions
Fconst = @(x_opt) constraints(x_opt, c0, R_air, Theta, x_f, y_f, rho, s, ct,cm,m_0,N_p);
Fmin   = @(x_opt) minsolves(x_opt);

typicalX = [Init(1);Init(2);Init(3)];
options = optimoptions(@fmincon, ...
    'Algorithm','interior-point', ...
    'MaxFunctionEvaluations', 1e6, ...
    'MaxIterations',100, ...
    'ScaleProblem','obj-and-constr', ...
    'TypicalX',typicalX, ...
    'ConstraintTolerance', 1e-3);
[x_opt, fval, exitflag, output, lambda] = fmincon(Fmin, Init, [], [], [], [], [], [], Fconst, options);
toc
%% Visualization and Post-Processing
% State vector: [x, y, m, z, lambda_x, q]
X = zeros(N_p, 6);
Time=zeros(N_p,1);
MaxThrust=zeros(N_p,1);
Drag=zeros(N_p,1);
Pi=zeros(N_p,1);

v_opt = zeros(N_p, 1);
X(1,3) = m_0;
X(1,5) = x_opt(1);
X(1,6) = x_opt(2);

t_f = x_opt(3);
dt = t_f / (N_p - 1);
Time(N_p,1)=t_f;

for i = 1:N_p-1
    % Extract current state variables:
    Time(i,1)=(i-1)*dt;
    x_var  = X(i,1);
    y_var  = X(i,2);
    m_var  = X(i,3);
    z_var = X(i,4);
    lx_var = X(i,5);
    q_var  = X(i,6);
    
    % Define function handles (all depend on v):
    Qx_v = @(v) W_x_func(x_var, y_var) + v/sqrt(1 + q_var^2);
    Qy_v = @(v) W_y_func(x_var, y_var) + v*q_var/sqrt(1 + q_var^2);
    Qz_v = g_func(x_var, y_var);  % g_func is independent of v
    
    % Mach number as a function of v:
    M = @(v) v / sqrt(c0 * R_air * Theta);
    
    % Define T_1 and T_2 with evaluation at v:
    T_1 = @(v) (1 / C_s_func(v, M(v))) * dC_s_dv_func(v, M(v));
    T_2 = @(v) (1 / C_D_func(v, M(v), m_var)) * dC_D_dv_func(v, M(v), m_var);
    
    % Algebraic equation for optimal speed v:
    Qv = @(v) lx_var * sqrt(1 + q_var^2) - ( ct + Qz_v + lx_var * Qx_v(v) + lx_var*q_var*Qy_v(v) ) ...
        * ( T_1(v) + 2/v + T_2(v) );
    
    %     vinit = 190;  % initial guess for v
    %     v = fzero(Qv, vinit);
    %     v_opt(i) = v;
    
    residual = @(v) abs(Qv(v));
    v_lo = 150;
    v_hi = 280;
    v = fminbnd(residual, v_lo, v_hi);
    v_opt(i) = v;
    
    MaxThrust(i,1)=T_max_func(M(v));
    Drag(i,1)=0.5 * rho * s * v^2 * C_D_func(v, M(v), m_var);
    Pi(i,1)=Drag(i,1)/MaxThrust(i,1);
    %%RK-3rd-Order Block
    x_var_ref=x_var;
    y_var_ref=y_var;
    m_var_ref=m_var;
    z_var_ref=z_var;
    lx_var_ref=lx_var;
    q_var_ref=q_var;
    
    Qx=zeros(3,1);
    Qy=zeros(3,1);
    Qm=zeros(3,1);
    Qz=zeros(3,1);
    Qlx=zeros(3,1);
    Qq=zeros(3,1);
    for RK_Iteration=1:3
        Qx(RK_Iteration) = W_x_func(x_var, y_var) + v/sqrt(1 + q_var^2);
        Qy(RK_Iteration) = W_y_func(x_var, y_var) + v*q_var/sqrt(1 + q_var^2);
        Qm(RK_Iteration) = -0.5 * rho * s * v^2 * C_D_func(v, M(v), m_var) * C_s_func(v, M(v));
        Qz(RK_Iteration) = g_func(x_var, y_var);
        Qlx(RK_Iteration) = - dg_dx_func(x_var, y_var) - lx_var*(dWx_dx_func(x_var, y_var) + q_var*dWy_dx_func(x_var, y_var));
        Qq(RK_Iteration)  = - dWx_dy_func(x_var, y_var) + (dWx_dx_func(x_var, y_var) - dWy_dy_func(x_var, y_var))*q_var + ...
            dWy_dx_func(x_var, y_var)*q_var^2 + (1/lx_var)*(q_var*dg_dx_func(x_var, y_var) - dg_dy_func(x_var, y_var));
        
        if RK_Iteration==1
            x_var=x_var_ref+Qx(RK_Iteration)*dt/2;
            y_var=y_var_ref+Qy(RK_Iteration)*dt/2;
            m_var=m_var_ref+Qm(RK_Iteration)*dt/2;
            z_var=z_var_ref+Qz(RK_Iteration)*dt/2;
            lx_var=lx_var_ref+Qlx(RK_Iteration)*dt/2;
            q_var=q_var_ref+Qq(RK_Iteration)*dt/2;
        end
        
        if RK_Iteration==2
            x_var=x_var_ref-Qx(RK_Iteration-1)*dt+2*Qx(RK_Iteration)*dt;
            y_var=y_var_ref-Qy(RK_Iteration-1)*dt+2*Qy(RK_Iteration)*dt;
            m_var=m_var_ref-Qm(RK_Iteration-1)*dt+2*Qm(RK_Iteration)*dt;
            z_var=z_var_ref-Qz(RK_Iteration-1)*dt+2*Qz(RK_Iteration)*dt;
            lx_var=lx_var_ref-Qlx(RK_Iteration-1)*dt+2*Qlx(RK_Iteration)*dt;
            q_var=q_var_ref-Qq(RK_Iteration-1)*dt+2*Qq(RK_Iteration)*dt;
        end
        
    end
    
    % Update states:
    X(i+1,1) = X(i,1) + (Qx(1)+4*Qx(2)+Qx(3)) * dt/6;
    X(i+1,2) = X(i,2) + (Qy(1)+4*Qy(2)+Qy(3)) * dt/6;
    X(i+1,3) = X(i,3) + (Qm(1)+4*Qm(2)+Qm(3)) * dt/6;
    X(i+1,4) = X(i,4) + (Qz(1)+4*Qz(2)+Qz(3)) * dt/6;
    X(i+1,5) = X(i,5) + (Qlx(1)+4*Qlx(2)+Qlx(3)) * dt/6;
    X(i+1,6) = X(i,6) + (Qq(1)+4*Qq(2)+Qq(3)) * dt/6;
end
v_opt(N_p,1)=v_opt(N_p-1,1);
Pi(N_p,1)=Pi(N_p-1,1);

figure (2)
plot(X(:,1), X(:,2), 'k', 'LineWidth', line_width)
hold on
bullet_1 = [0, x_f];
bullet_2 = [0, y_f];
plot(bullet_1, bullet_2, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); 
hold on

figure (3)
plot(Time, X(:,3), 'k', 'LineWidth', line_width)
hold on
xlabel('$t[s]$', 'Interpreter', 'latex', 'FontSize', font_size);
ylabel('$m(t)[kg]$', 'Interpreter', 'latex', 'FontSize', font_size);
title('Optimal Mass-Time Plot', 'Interpreter', 'latex', 'FontSize', font_size);
grid on;
set(gca, 'FontSize', font_size, 'TickLabelInterpreter', 'latex', 'GridColor', 'r');

figure (4)
plot(Time, X(:,5), 'k', 'LineWidth', line_width)
hold on
xlabel('$t[s]$', 'Interpreter', 'latex', 'FontSize', font_size);
ylabel('$\lambda_x(t)$', 'Interpreter', 'latex', 'FontSize', font_size);
title('Optimal $\lambda_x$-Time Plot', 'Interpreter', 'latex', 'FontSize', font_size);
grid on;
set(gca, 'FontSize', font_size, 'TickLabelInterpreter', 'latex', 'GridColor', 'r');

figure (5)
plot(Time, atan(X(:,6)), 'k', 'LineWidth', line_width)
hold on
xlabel('$t[s]$', 'Interpreter', 'latex', 'FontSize', font_size);
ylabel('$\chi(t)[rad]$', 'Interpreter', 'latex', 'FontSize', font_size);
title('Optimal Heading Angle-Time Plot', 'Interpreter', 'latex', 'FontSize', font_size);
grid on;
set(gca, 'FontSize', font_size, 'TickLabelInterpreter', 'latex', 'GridColor', 'r');

figure (6)
plot(Time, v_opt, 'k', 'LineWidth', line_width)
hold on
xlabel('$t[s]$', 'Interpreter', 'latex', 'FontSize', font_size);
ylabel('$v(t)[m/s]$', 'Interpreter', 'latex', 'FontSize', font_size);
title('Optimal Speed-Time Plot', 'Interpreter', 'latex', 'FontSize', font_size);
grid on;
set(gca, 'FontSize', font_size, 'TickLabelInterpreter', 'latex', 'GridColor', 'r');

figure (7)
plot(Time, Pi, 'k', 'LineWidth', line_width)
hold on
xlabel('$t[s]$', 'Interpreter', 'latex', 'FontSize', font_size);
ylabel('$\Pi(t)$', 'Interpreter', 'latex', 'FontSize', font_size);
title('Optimal Throttle Setting-Time Plot', 'Interpreter', 'latex', 'FontSize', font_size);
grid on;
set(gca, 'FontSize', font_size, 'TickLabelInterpreter', 'latex', 'GridColor', 'r');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constraint Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cin, ceq] = constraints(x_opt, c0, R_air, Theta, x_f, y_f, rho, s,ct,cm,m_0,N_p)
global W_x_func W_y_func g_func dg_dx_func dg_dy_func ...
    dWx_dx_func dWx_dy_func dWy_dx_func dWy_dy_func ...
    C_D_func dC_D_dv_func C_s_func dC_s_dv_func

% Preallocate state trajectory matrix
% State vector: [x, y, m, z, lambda_x, q]
X = zeros(N_p, 6);
SaveQq=zeros(N_p-1,1);
X(1,3) = m_0;
X(1,5) = x_opt(1); % initial lambda_x(0)
X(1,6) = x_opt(2); % initial q(0)
V_u=zeros(N_p,1);

t_f = x_opt(3);
dt = t_f / (N_p - 1);

for i = 1:N_p-1
    % Extract current state variables:
    x_var  = X(i,1);
    y_var  = X(i,2);
    m_var  = X(i,3);
    z_var = X(i,4);
    lx_var = X(i,5);
    q_var  = X(i,6);
    
    % Define function handles (all depend on v):
    Qx_v = @(v) W_x_func(x_var, y_var) + v/sqrt(1 + q_var^2);
    Qy_v = @(v) W_y_func(x_var, y_var) + v*q_var/sqrt(1 + q_var^2);
    Qz_v = g_func(x_var, y_var);  % g_func is independent of v
    
    % Mach number as a function of v:
    M = @(v) v / sqrt(c0 * R_air * Theta);
    
    % Define T_1 and T_2 with evaluation at v:
    T_1 = @(v) (1 / C_s_func(v, M(v))) * dC_s_dv_func(v, M(v));
    T_2 = @(v) (1 / C_D_func(v, M(v), m_var)) * dC_D_dv_func(v, M(v), m_var);
    
    % Algebraic equation for optimal speed v:
    Qv = @(v) lx_var * sqrt(1 + q_var^2) - ( ct + Qz_v + lx_var * Qx_v(v) + lx_var*q_var*Qy_v(v) ) ...
        * ( T_1(v) + 2/v + T_2(v) );
    
    residual = @(v) abs(Qv(v));
    v_lo = 150;
    v_hi = 280;
    v = fminbnd(residual, v_lo, v_hi);
    V_u(i,1) = v;
    
    %%RK-3rd-Order Block
    x_var_ref=x_var;
    y_var_ref=y_var;
    m_var_ref=m_var;
    z_var_ref=z_var;
    lx_var_ref=lx_var;
    q_var_ref=q_var;
    
    Qx=zeros(3,1);
    Qy=zeros(3,1);
    Qm=zeros(3,1);
    Qz=zeros(3,1);
    Qlx=zeros(3,1);
    Qq=zeros(3,1);
    for RK_Iteration=1:3
        Qx(RK_Iteration) = W_x_func(x_var, y_var) + v/sqrt(1 + q_var^2);
        Qy(RK_Iteration) = W_y_func(x_var, y_var) + v*q_var/sqrt(1 + q_var^2);
        Qm(RK_Iteration) = -0.5 * rho * s * v^2 * C_D_func(v, M(v), m_var) * C_s_func(v, M(v));
        Qz(RK_Iteration) = g_func(x_var, y_var);
        Qlx(RK_Iteration) = - dg_dx_func(x_var, y_var) - lx_var*(dWx_dx_func(x_var, y_var) + q_var*dWy_dx_func(x_var, y_var));
        Qq(RK_Iteration)  = - dWx_dy_func(x_var, y_var) + (dWx_dx_func(x_var, y_var) - dWy_dy_func(x_var, y_var))*q_var + ...
            dWy_dx_func(x_var, y_var)*q_var^2 + (1/lx_var)*(q_var*dg_dx_func(x_var, y_var) - dg_dy_func(x_var, y_var));
        
        if RK_Iteration==1
            x_var=x_var_ref+Qx(RK_Iteration)*dt/2;
            y_var=y_var_ref+Qy(RK_Iteration)*dt/2;
            m_var=m_var_ref+Qm(RK_Iteration)*dt/2;
            z_var=z_var_ref+Qz(RK_Iteration)*dt/2;
            lx_var=lx_var_ref+Qlx(RK_Iteration)*dt/2;
            q_var=q_var_ref+Qq(RK_Iteration)*dt/2;
        end
        
        if RK_Iteration==2
            x_var=x_var_ref-Qx(RK_Iteration-1)*dt+2*Qx(RK_Iteration)*dt;
            y_var=y_var_ref-Qy(RK_Iteration-1)*dt+2*Qy(RK_Iteration)*dt;
            m_var=m_var_ref-Qm(RK_Iteration-1)*dt+2*Qm(RK_Iteration)*dt;
            z_var=z_var_ref-Qz(RK_Iteration-1)*dt+2*Qz(RK_Iteration)*dt;
            lx_var=lx_var_ref-Qlx(RK_Iteration-1)*dt+2*Qlx(RK_Iteration)*dt;
            q_var=q_var_ref-Qq(RK_Iteration-1)*dt+2*Qq(RK_Iteration)*dt;
        end
        
    end
    
    % Update states:
    X(i+1,1) = X(i,1) + (Qx(1)+4*Qx(2)+Qx(3)) * dt/6;
    X(i+1,2) = X(i,2) + (Qy(1)+4*Qy(2)+Qy(3)) * dt/6;
    X(i+1,3) = X(i,3) + (Qm(1)+4*Qm(2)+Qm(3)) * dt/6;
    X(i+1,4) = X(i,4) + (Qz(1)+4*Qz(2)+Qz(3)) * dt/6;
    X(i+1,5) = X(i,5) + (Qlx(1)+4*Qlx(2)+Qlx(3)) * dt/6;
    X(i+1,6) = X(i,6) + (Qq(1)+4*Qq(2)+Qq(3)) * dt/6;
    SaveQq(i,1)=(Qq(1)+4*Qq(2)+Qq(3))/6;
end
V_u(N_p,1)=V_u(N_p-1,1);
Qx_end=(Qx(1)+4*Qx(2)+Qx(3))/6;
Qy_end=(Qy(1)+4*Qy(2)+Qy(3))/6;
Qz_end=(Qz(1)+4*Qz(2)+Qz(3))/6;
Qm_end=(Qm(1)+4*Qm(2)+Qm(3))/6;
cin (1,1) = -t_f;

ceq(1,1)=X(N_p,1)-x_f;
ceq(2,1)=X(N_p,2)-y_f;
ceq(3,1) = -(ct + Qz_end + lx_var*(Qx_end + q_var*Qy_end)) / Qm_end - cm;
end
%% Objective Function (Dummy)
function Y = minsolves(x_opt)
Y = 1;
end