clc; clear;

%% ────────── 已有几何部分 ──────────
syms x y L L1 L2 L5 theta real     % 符号量

% 方程
eq1 = (x - L*sin(theta))^2 + (y + L1 - L*cos(theta))^2 == L5^2;
eq2 = x^2 + y^2 == (x - L2*sin(theta))^2 + (y + L1 - L2*cos(theta))^2;

% 解析解
sol = solve([eq1, eq2], [x, y], 'Real', true);

% 替换数值
L_val  = 12;
L1_val = 23;
L2_val = 24;
L5_val = 10;
theta_val = deg2rad(83.66);

vars = [L, L1, L2, L5, theta];
vals = [L_val, L1_val, L2_val, L5_val, theta_val];

x_all = double(subs(sol.x, vars, vals));
y_all = double(subs(sol.y, vars, vals));

% 取 x>0 的那一支
idx     = find(x_all > 0, 1);
x_val   = x_all(idx);
y_val   = y_all(idx);

% 坐标
A = [0, 0];
B = [0, -L1_val];
C = [L2_val * sin(theta_val), -L1_val + L2_val * cos(theta_val)];
D = [L_val * sin(theta_val),  -L1_val + L_val  * cos(theta_val)];
E = [x_val, y_val];

AE = norm(E - A);
EC = norm(E - C);

fprintf('x  = %.4f, y  = %.4f\n', x_val, y_val);
fprintf('AE = %.4f\n', AE);
fprintf('EC = %.4f\n\n', EC);

%% ────────── 追加：求解 Fec, Fde, Fae ──────────
syms Fec Fde Fae real    % 未知杆力
W = 10;                  % 外荷载

% ── 力臂中需要的 sin(∠BCE) 与 sin(∠BDE) ──
vecCB  = B - C;                       % C→B
vecCE  = E - C;                       % C→E
sinBCE = abs(det([vecCB; vecCE])) / (norm(vecCB)*norm(vecCE));

vecDB  = B - D;                       % D→B
vecDE  = E - D;                       % D→E
sinBDE = abs(det([vecDB; vecDE])) / (norm(vecDB)*norm(vecDE));

% ── 三个独立方程 ──
eqF(1) = Fec*L2_val*sinBCE + Fde*L_val*sinBDE == W*L2_val*sin(theta_val);
eqF(2) = Fae*(E(1)-A(1))/AE  + Fde*(E(1)-D(1))/L5_val == Fec*(C(1)-E(1))/EC;
eqF(3) = Fae*(A(2)-E(2))/AE  == Fde*(E(2)-D(2))/L5_val + Fec*(E(2)-C(2))/EC;

solF   = solve(eqF, [Fec, Fde, Fae]);

Fec_val = double(solF.Fec);
Fde_val = double(solF.Fde);
Fae_val = double(solF.Fae);

fprintf('Fec = %.4f\n', Fec_val);
fprintf('Fde = %.4f\n', Fde_val);
fprintf('Fae = %.4f\n\n', Fae_val);

%% ────────── 绘图（保持不变）──────────
figure; hold on; axis equal
plot([A(1), B(1)], [A(2), B(2)], 'k-', 'LineWidth', 2)
plot([B(1), C(1)], [B(2), C(2)], 'k-', 'LineWidth', 2)
plot([A(1), E(1)], [A(2), E(2)], 'r--','LineWidth',1.5)   % AE
plot([E(1), C(1)], [E(2), C(2)], 'b--','LineWidth',1.5)   % EC
plot([D(1), E(1)], [D(2), E(2)], 'g--','LineWidth',1.5)   % DE
scatter([A(1),B(1),C(1),D(1),E(1)], [A(2),B(2),C(2),D(2),E(2)],60,'filled')
text(A(1),A(2),' A'); text(B(1),B(2),' B');
text(C(1),C(2),' C'); text(D(1),D(2),' D');
text(E(1),E(2),' E');
title('Structure'); grid on;
