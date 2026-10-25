%% =============================================================
%  Построение графиков по результатам вычислительных экспериментов
%  Файлы: ConditionalDepend.csv, DimensionDepend.csv
% =============================================================
clear; clc; close all;

%% ---------- 1. Чтение данных -----------------------------------
%dir = "C:\Users\2pm4\source\repos\Lex0z\Solution-of-a-system-of-linear-algebraic-equations2\out\build\x64-Debug"
cond_data = readmatrix('../data/ConditionalDepend.csv');  % [cond, err, errd, nev, nev_norm]
dim_data  = readmatrix('../data/DimensionDepend.csv');    % [cond, dimension, iter]

cond     = cond_data(:,1);
err      = cond_data(:,2);
err_norm = cond_data(:,3);
nev      = cond_data(:,4);
nev_norm = cond_data(:,5);

dim  = dim_data(:,2);
iter = dim_data(:,3);

%% ---------- 2. График: err и nev от cond -----------------------
figure(1); clf;
loglog(cond, err_norm, 'o-', 'LineWidth', 1.6, 'MarkerSize', 9); hold on;
loglog(cond, nev_norm, 's-', 'LineWidth', 1.6, 'MarkerSize', 9);


% опорная линия cond * eps_M (теоретическая оценка сверху)
epsM = 1e-16;
loglog(cond, cond*epsM, 'k--', 'LineWidth', 1.2);
% loglog(cond, err, 'o--', 'LineWidth', 1.6, 'MarkerSize', 9); 


grid on; box on;
set(gca, 'XScale','log', 'YScale','log', 'FontSize', 12);
xlabel('Число обусловленности  cond(A)', 'FontSize', 14);
ylabel('Ошибка / невязка',               'FontSize', 14);
title('Влияние числа обусловленности на точность', 'FontSize', 14);
legend({'$\frac{\|x - x^*\|}{\|x^*\|}$', ...
        '$\frac{\|Ax - b\|}{\|b\|}$', ...
        '$\mathrm{cond}(A)\cdot\varepsilon_M$'}, ...
       'Location','best', 'FontSize', 12, 'Interpreter','latex');

%% ---------- 3. График: err и nev от cond -----------------------
figure(2); clf;
loglog(cond, err, 'o-', 'LineWidth', 1.6, 'MarkerSize', 9); hold on;
loglog(cond, nev, 's-', 'LineWidth', 1.6, 'MarkerSize', 9);


% опорная линия cond * eps_M (теоретическая оценка сверху)
epsM = 1e-16;
loglog(cond, cond*epsM, 'k--', 'LineWidth', 1.2);
% loglog(cond, err, 'o--', 'LineWidth', 1.6, 'MarkerSize', 9); 


grid on; box on;
set(gca, 'XScale','log', 'YScale','log', 'FontSize', 12);
xlabel('Число обусловленности  cond(A)', 'FontSize', 14);
ylabel('Ошибка / невязка',               'FontSize', 14);
title('Влияние числа обусловленности на точность', 'FontSize', 14);
legend({'$\|x - x^*\|$', ...
        '$\|Ax - b\|$', ...
        '$\mathrm{cond}(A)\cdot\varepsilon_M$'}, ...
       'Location','best', 'FontSize', 12, 'Interpreter','latex');

%% ---------- 4. График: iter от dim -----------------------------
figure(3); clf;
loglog(dim, iter, 'o-', 'LineWidth', 1.6, 'MarkerSize', 9); hold on;

% опорная линия ~ n^3, подобранная по последней точке
c = iter(end) / dim(end)^3;
n_ref = logspace(log10(min(dim)), log10(max(dim)), 100);
loglog(n_ref, c*n_ref.^3, 'k--', 'LineWidth', 1.2);

grid on; box on;
set(gca, 'XScale','log', 'YScale','log', 'FontSize', 12);
xlabel('Размерность системы  n', 'FontSize', 14);
ylabel('Число операций',         'FontSize', 14);
title('Влияние размерности на объём вычислений', 'FontSize', 14);
legend({'Эксперимент', '$\sim n^3$'}, ...
       'Location','northwest', 'FontSize', 12, 'Interpreter','latex');

%% ---------- 4. Числа для отчёта -------------------------------
% таблица iter/n^3
fprintf('--- iter / n^3 ---\n');
for i = 1:length(dim)
    fprintf('n = %4d | iter = %10d | iter/n^3 = %.4f\n', ...
            dim(i), iter(i), iter(i)/dim(i)^3);
end

% наклон iter ~ n^p
p = polyfit(log(dim), log(iter), 1);
fprintf('\nНаклон log(iter) vs log(n)   = %.4f  (ожидается ~3)\n', p(1));

% наклон err_norm ~ cond^q
mask = err_norm > 0 & cond > 0;
q1 = polyfit(log(cond(mask)), log(err_norm(mask)), 1);
fprintf('Наклон log(err_norm) vs log(cond) = %.4f  (ожидается ~1)\n', q1(1));

% наклон err_norm ~ cond^q
mask3 = err > 0 & cond > 0;
q3 = polyfit(log(cond(mask3)), log(err(mask3)), 1);
fprintf('Наклон log(err) vs log(cond) = %.4f  (ожидается ~1)\n', q3(1));

% наклон nev_norm ~ cond^q (должен быть ~0 — невязка не зависит от cond)
mask2 = nev_norm > 0 & cond > 0;
q2 = polyfit(log(cond(mask2)), log(nev_norm(mask2)), 1);
fprintf('Наклон log(nev_norm) vs log(cond) = %.4f  (ожидается ~0)\n', q2(1));

% наклон nev_norm ~ cond^q (должен быть ~0 — невязка не зависит от cond)
mask4 = nev > 0 & cond > 0;
q4 = polyfit(log(cond(mask4)), log(nev(mask4)), 1);
fprintf('Наклон log(nev) vs log(cond) = %.4f  (ожидается ~1)\n', q4(1));