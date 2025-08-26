clc; clear;

% 固定参数
L1 = 25.2;
L2 = 26.8;
theta_list = deg2rad(1:1:179);      % 弧度
theta_deg_list = rad2deg(theta_list);
L_list = (1: 0.25: L2-1);            % L 值范围
L5_ratios = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7,0.75];

% 初始化面积结果
area_result = zeros(length(L5_ratios), 1);

for r = 1:length(L5_ratios)
    L5 = L5_ratios(r) * L2;
    AE_mat = nan(length(L_list), length(theta_list));  % 存储 AE 值

    for i = 1:length(L_list)
        L = L_list(i);
        for j = 1:length(theta_list)
            theta = theta_list(j);

            fprintf('正在计算：θ = %3.f°, L5 = %.2f, L = %.1f\n', rad2deg(theta), L5_ratios(r), L);

            % 检查几何可行性
            syms x y Lsym L1sym L2sym L5sym thetasym real
            exprLHS = (2*Lsym*L1sym*cos(thetasym) - 2*Lsym*L2sym - L1sym^2 + L2sym^2)^2;
            exprRHS = 4*L5sym^2 * (L1sym^2 - 2*L1sym*L2sym*cos(thetasym) + L2sym^2);
            isValid = double(subs(exprLHS <= exprRHS, ...
                [Lsym, L1sym, L2sym, L5sym, thetasym], ...
                [L, L1, L2, L5, theta]));

            if ~isValid
                continue;
            end

            % 解方程组
            syms x y
            eq1 = (x - L*sin(theta))^2 + (y + L1 - L*cos(theta))^2 == L5^2;
            eq2 = x^2 + y^2 == (x - L2*sin(theta))^2 + (y + L1 - L2*cos(theta))^2;
            sol = solve([eq1, eq2], [x, y], 'Real', true);

            x_all = double(sol.x);
            y_all = double(sol.y);

            % 点坐标
            A = [0, 0];
            B = [0, -L1];
            C = [L2*sin(theta), -L1 + L2*cos(theta)];

            % BC 线斜率
            xB = B(1); yB = B(2);
            xC = C(1); yC = C(2);
            slope_BC = (yC - yB) / (xC - xB);

            % AC 线斜率
            xA = A(1); yA = A(2);
            slope_AC = (yC - yA) / (xC - xA);

            % 逐个筛选合法解
            valid_idx = [];
            for k = 1:length(x_all)
                xE = x_all(k);
                yE = y_all(k);
                y_BC = yB + slope_BC * (xE - xB);
                y_AC = yA + slope_AC * (xE - xA);
                if xE > 0 && yE > y_BC && yE < y_AC
                    valid_idx(end+1) = k;
                end
            end

            if isempty(valid_idx)
                continue;
            end

            % 多解中选 x 最大的
            [~, best_k] = max(x_all(valid_idx));
            idx = valid_idx(best_k);
            E = [x_all(idx), y_all(idx)];

            AE = norm(E - A);
            AE_mat(i,j) = AE;
        end
    end

    % 去除 theta < 34° 的部分
    theta_mask = theta_deg_list >= 34;
    AE_mat(:, ~theta_mask) = NaN;

    % 去除 AE < 17.68 的部分
    AE_mat(AE_mat < 16.38) = NaN;

    % 计算有效面积（在L-θ平面上的投影面积）
    delta_L = 0.25;
    delta_theta = 1;  % 度数
    valid_area = sum(~isnan(AE_mat), 'all') * delta_L * delta_theta;
    area_result(r) = valid_area;

    % 绘图
    [L_grid, theta_grid] = meshgrid(L_list, 180-theta_deg_list);
    figure;
    surf(L_grid, theta_grid, AE_mat', 'EdgeColor','none'); colorbar
    xlabel('L'); ylabel('\theta (deg)'); zlabel('AE');
    title(['AE vs L and \theta, L5 = ', num2str(L5_ratios(r)), '·L2']);
    view(135, 30);
end

% 输出面积汇总表格
T = table(L5_ratios(:), area_result, 'VariableNames', {'L5_over_L2', 'RemainingArea'});
disp(T);
% 保存表格为 CSV 文件
writetable(T, 'AE_Area_vs_L5ratio.csv');
disp('表格已保存为 AE_Area_vs_L5ratio.csv');
