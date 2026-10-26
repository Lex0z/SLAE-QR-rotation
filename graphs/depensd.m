%% Численный анализ QR-разложения методом вращений Гивенса
clear; clc; close all;

%% 1. Данные эксперимента по обусловленности
cond_data = readmatrix('../data/ConditionalDepend.csv');
cond     = cond_data(:, 1);
err_abs      = cond_data(:, 2); % фактическая ошибка
residual_abs = cond_data(:, 4); % невязка
err      = cond_data(:, 3); % относительная фактическая ошибка
residual = cond_data(:, 5); % относительная невязка

figure(1);
loglog(cond, err, 'o-', 'LineWidth', 1.6, 'MarkerSize', 7); hold on;
loglog(cond, residual_abs, 's-', 'LineWidth', 1.6, 'MarkerSize', 7); hold on;
loglog(cond, residual, 's-', 'LineWidth', 1.6, 'MarkerSize', 7);


epsM = 1e-16;
loglog(cond, cond*epsM, 'k--', 'LineWidth', 1.2);

grid on; box on;
xlabel('Число обусловленности cond(A)');
ylabel('Относительная величина');
title('Влияние числа обусловленности на точность решения');
legend({'Относительная ошибка решения', 'Абсолютная невязка', 'Относительная невязка', '$\mathrm{cond}(A)\cdot\varepsilon_M$'}, ...
       'Location', 'best', 'Interpreter','latex');

%% 2. Данные эксперимента по размерности
% CSV columns: cond, dimension, arithmetic operation count
dim_data = readmatrix('../data/DimensionDepend.csv');
dim  = dim_data(:, 2);
ops  = dim_data(:, 3);

figure(2);
loglog(dim, ops, 'o-', 'LineWidth', 1.6, 'MarkerSize', 7); hold on;
% Опорная кубическая зависимость, нормированная по последней точке.
reference = ops(end) / dim(end)^3 * dim.^3;
loglog(dim, reference, '--', 'LineWidth', 1.4);
grid on; box on;
xlabel('Размерность системы n');
ylabel('Количество учитываемых операций');
title('Зависимость количества операций от размерности');
legend({'Экспериментальные значения', 'Опорная зависимость n^3'}, ...
       'Location', 'northwest');

%% 3. Численные значения для раздела 8 отчёта
fprintf('Эксперимент 1: обусловленность\n');
for i = 1:numel(cond)
    fprintf('cond = %9.2e | rel.error = %.6e | rel.residual = %.6e\n', ...
        cond(i), err(i), residual(i));
end

fprintf('\nЭксперимент 2: размерность\n');
for i = 1:numel(dim)
    fprintf('n = %3d | operations = %10.0f | operations/n^3 = %.6e\n', ...
        dim(i), ops(i), ops(i) / dim(i)^3);
end
if all(ops > 0)
    fit = polyfit(log(dim), log(ops), 1);
    fprintf('\nНаклон log(operations) от log(n): %.4f (ожидается порядок около 3)\n', fit(1));
end
