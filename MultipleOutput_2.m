clc; clear;

%% ────────── 定义符号量与几何方程 ──────────
syms x y L L1 L2 L5 theta real

% 几何方程
eq1 = (x - L*sin(theta))^2 + (y + L1 - L*cos(theta))^2 == L5^2;
eq2 = x^2 + y^2 == (x - L2*sin(theta))^2 + (y + L1 - L2*cos(theta))^2;

% 参数赋值
L_val  = 18;
L1_val = 25.2;
L2_val = 26.8;
L5_val = 0.35 * L2_val;
theta_val = deg2rad(180-145);

%% ────────── 判断几何约束是否满足 ──────────
exprLHS = (2*L*L1*cos(theta) - 2*L*L2 - L1^2 + L2^2)^2;
exprRHS = 4*L5^2 * (L1^2 - 2*L1*L2*cos(theta) + L2^2);

valid_geom = double(subs(exprLHS <= exprRHS, ...
                  [L, L1, L2, L5, theta], ...
                  [L_val, L1_val, L2_val, L5_val, theta_val]));

if ~valid_geom
    fprintf('[×] 几何约束不满足：无解存在。\n');
    return
end

%% ────────── 解方程组 ──────────
sol = solve([eq1, eq2], [x, y], 'Real', true);
x_all = double(subs(sol.x, [L, L1, L2, L5, theta], [L_val, L1_val, L2_val, L5_val, theta_val]));
y_all = double(subs(sol.y, [L, L1, L2, L5, theta], [L_val, L1_val, L2_val, L5_val, theta_val]));

%% ────────── 判断合法解 ──────────
% 三点坐标
A = [0, 0];
B = [0, -L1_val];
C = [L2_val * sin(theta_val), -L1_val + L2_val * cos(theta_val)];

xC = C(1); yC = C(2);
xB = B(1); yB = B(2);
xA = A(1); yA = A(2);

% 直线BC 和 AC 的表达式： y = kx + b
kBC = (yC - yB) / (xC - xB); bBC = yB;
kAC = (yC - yA) / (xC - xA); bAC = yA;

y_BC = kBC * x_all + bBC;
y_AC = kAC * x_all + bAC;

valid_idx = find((x_all > 0) & (y_all > y_BC) & (y_all < y_AC));

if isempty(valid_idx)
    fprintf('[×] 无合法几何解（不满足 Ex>0，Ey > BC 且 Ey < AC）。\n');
    return
end

% 多解中取 Ex 最大
[~, best_i] = max(x_all(valid_idx));
idx = valid_idx(best_i);

x_val = x_all(idx);
y_val = y_all(idx);
E = [x_val, y_val];

%% ────────── 输出长度信息 ──────────
D = [L_val * sin(theta_val), -L1_val + L_val * cos(theta_val)];

AE = norm(E - A);
EC = norm(E - C);

fprintf('✔ 合法解选定：\n');
fprintf('x  = %.4f, y  = %.4f\n', x_val, y_val);
fprintf('AE = %.4f\n', AE);
fprintf('EC = %.4f\n\n', EC);

%% ────────── 力平衡求解 Fec, Fde, Fae ──────────
syms Fec Fde Fae real
W = 10;

vecCB  = B - C;
vecCE  = E - C;
sinBCE = abs(det([vecCB; vecCE])) / (norm(vecCB)*norm(vecCE));

vecDB  = B - D;
vecDE  = E - D;
sinBDE = abs(det([vecDB; vecDE])) / (norm(vecDB)*norm(vecDE));

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

%% ────────── 绘图 ──────────
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
