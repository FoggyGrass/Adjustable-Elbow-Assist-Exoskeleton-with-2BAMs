clc; clear; close all;

%% --- 参数定义 ---
L1 = 23; L2 = 24;              % AB、BC 长度
l3_0 = 16; l4_0 = 16; l5_0 = 39;
r = 0.5;                       % 最大收缩率
l3_vals = linspace(l3_0*(1 - r), l3_0, 6);
l4_vals = linspace(l4_0*(1 - r), l4_0, 6);
l5_range = [l5_0*(1 - r), l5_0];

A = [0; 0];
B = A + [0; -L1];

% 初始化数据
record = [];  % 每行：[angleABC, l3, l4, l5]
C_all = []; D_all = [];

%% --- 主循环 ---
for theta = linspace(0.01, pi - 0.01, 200)
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    C = B + R * [0; -L2];

    for l3 = l3_vals
        for l4 = l4_vals
            [D1, D2, valid] = circle_intersection(B, l4, C, l3);
            if ~valid, continue; end

            for D = [D1, D2]
                D = D(:);
                BA = A - B; BC = C - B; BD = D - B;

                if ~is_in_angle_region(BA, BC, BD), continue; end

                l5 = norm(D - A);
                if l5 < l5_range(1) || l5 > l5_range(2), continue; end

                angleABC = acos(dot(BA, BC)/(norm(BA)*norm(BC)));
                record(end+1,:) = [rad2deg(angleABC), l3, l4, l5];
                C_all(:,end+1) = C;
                D_all(:,end+1) = D;
            end
        end
    end
end

%% --- 图像绘制：C点工作空间 + 最小/最大角构型 ---
figure; hold on; axis equal;
title('C点工作空间 + 最小/最大角构型');
xlabel('X (cm)'); ylabel('Y (cm)');

fill(C_all(1,:), C_all(2,:), [0.6 0.8 1], 'FaceAlpha', 0.4, 'EdgeColor','none');
plot(A(1), A(2), 'ko', 'MarkerFaceColor', 'k'); text(A(1)-1, A(2)+1, 'A');
plot(B(1), B(2), 'ko', 'MarkerFaceColor', 'g'); text(B(1)-1, B(2), 'B');

% === 最小角构型绘制 ===
[min_angle, idx_min] = min(record(:,1));
C_min = C_all(:,idx_min); D_min = D_all(:,idx_min);
l3_min = record(idx_min,2); l4_min = record(idx_min,3); l5_min = record(idx_min,4);

plot(C_min(1), C_min(2), 'ro', 'MarkerFaceColor', 'r');
plot(D_min(1), D_min(2), 'mo', 'MarkerFaceColor', 'm');
plot([A(1) B(1)], [A(2) B(2)], 'k-', 'LineWidth', 2);
plot([B(1) C_min(1)], [B(2) C_min(2)], 'b-', 'LineWidth', 2);
plot([C_min(1) D_min(1)], [C_min(2) D_min(2)], 'c--');
plot([B(1) D_min(1)], [B(2) D_min(2)], 'm--');
plot([A(1) D_min(1)], [A(2) D_min(2)], 'g--');
text(D_min(1)+1, D_min(2), sprintf('\\theta_{min} = %.1f°', min_angle));

% === 最大角构型绘制 ===
[max_angle, idx_max] = max(record(:,1));
C_max = C_all(:,idx_max); D_max = D_all(:,idx_max);
l3_max = record(idx_max,2); l4_max = record(idx_max,3); l5_max = record(idx_max,4);

plot(C_max(1), C_max(2), 'rs', 'MarkerFaceColor', 'r');
plot(D_max(1), D_max(2), 'ms', 'MarkerFaceColor', 'm');
plot([B(1) C_max(1)], [B(2) C_max(2)], 'b-', 'LineWidth', 2);
plot([C_max(1) D_max(1)], [C_max(2) D_max(2)], 'c--');
plot([B(1) D_max(1)], [B(2) D_max(2)], 'm--');
plot([A(1) D_max(1)], [A(2) D_max(2)], 'g--');
text(D_max(1)+1, D_max(2), sprintf('\\theta_{max} = %.1f°', max_angle));

legend('C 工作空间', 'A', 'B', ...
       'C_{min}', 'D_{min}', ...
       'C_{max}', 'D_{max}', ...
       'AB','BC','CD','BD','AD');

%% --- 收缩率标注 ---
% 最小角构型
mid_CD = (C_min + D_min)/2;
mid_BD = (B + D_min)/2;
mid_AD = (A + D_min)/2;
text(mid_CD(1), mid_CD(2), sprintf('%.1f%%', (1 - l3_min/l3_0)*100), 'Color','c', 'FontSize', 9);
text(mid_BD(1), mid_BD(2), sprintf('%.1f%%', (1 - l4_min/l4_0)*100), 'Color','m', 'FontSize', 9);
text(mid_AD(1), mid_AD(2), sprintf('%.1f%%', (1 - l5_min/l5_0)*100), 'Color','g', 'FontSize', 9);

% 最大角构型
mid_CD = (C_max + D_max)/2;
mid_BD = (B + D_max)/2;
mid_AD = (A + D_max)/2;
text(mid_CD(1), mid_CD(2), sprintf('%.1f%%', (1 - l3_max/l3_0)*100), 'Color','c', 'FontSize', 9);
text(mid_BD(1), mid_BD(2), sprintf('%.1f%%', (1 - l4_max/l4_0)*100), 'Color','m', 'FontSize', 9);
text(mid_AD(1), mid_AD(2), sprintf('%.1f%%', (1 - l5_max/l5_0)*100), 'Color','g', 'FontSize', 9);

%% --- 连续性检查 ---
angles = sort(record(:,1));
angle_diff = diff(angles);
threshold = 5;  % 容差
discontinuities = find(angle_diff > threshold);
if isempty(discontinuities)
    disp('✅ 角ABC从最大角到最小角连续可达，无中断区间。');
else
    disp('⚠️ 检测到角ABC存在不可达区域：');
    for i = discontinuities'
        fprintf('缺口：%.2f° -> %.2f°\n', angles(i), angles(i+1));
    end
end

%% --- 控制台输出 ---
fprintf('\n------ 最小角度构型数据 ------\n');
fprintf('∠ABC = %.2f°\n', min_angle);
fprintf('l3 (CD) = %.2f cm, 收缩率 = %.1f%%\n', l3_min, (1 - l3_min / l3_0)*100);
fprintf('l4 (BD) = %.2f cm, 收缩率 = %.1f%%\n', l4_min, (1 - l4_min / l4_0)*100);
fprintf('l5 (AD) = %.2f cm, 收缩率 = %.1f%%\n', l5_min, (1 - l5_min / l5_0)*100);

fprintf('\n------ 最大角度构型数据 ------\n');
fprintf('∠ABC = %.2f°\n', max_angle);
fprintf('l3 (CD) = %.2f cm, 收缩率 = %.1f%%\n', l3_max, (1 - l3_max / l3_0)*100);
fprintf('l4 (BD) = %.2f cm, 收缩率 = %.1f%%\n', l4_max, (1 - l4_max / l4_0)*100);
fprintf('l5 (AD) = %.2f cm, 收缩率 = %.1f%%\n', l5_max, (1 - l5_max / l5_0)*100);

%% --- 几何函数 ---
function inside = is_in_angle_region(BA, BC, BD)
    BA = BA / norm(BA);
    BC = BC / norm(BC);
    BD = BD / norm(BD);
    angle_ABC = acos(dot(BA, BC));
    angle_ABD = acos(dot(BA, BD));
    angle_DBC = acos(dot(BD, BC));
    inside = abs(angle_ABD + angle_DBC - angle_ABC) < 1e-3;
end

function [P1, P2, valid] = circle_intersection(c1, r1, c2, r2)
    d = norm(c2 - c1);
    if d > r1 + r2 || d < abs(r1 - r2)
        P1 = []; P2 = []; valid = false; return;
    end
    a = (r1^2 - r2^2 + d^2)/(2*d);
    h = sqrt(max(r1^2 - a^2, 0));  % 防止负数
    p2 = c1 + a*(c2 - c1)/d;
    offset = h * [0 -1; 1 0] * (c2 - c1)/d;
    P1 = p2 + offset;
    P2 = p2 - offset;
    valid = true;
end
