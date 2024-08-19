%%

%%
%%
% Definir rutas de los archivos
fullpaths = {
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_8.1_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_8.2_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_8.3_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_8.4_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_8.5_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_8.6_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_8.7_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_8.8_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_8.9_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_9.0_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_9.1_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_9.2_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_9.3_1000_coupling/simulations.mat'
};

% Número de simulaciones por archivo
numSimulations = 1000;

% Inicializar variables
totalSimulations = numSimulations * length(fullpaths);
lengths = zeros(1, totalSimulations);
widths = zeros(1, totalSimulations);
meanslips = zeros(1, totalSimulations);
list = zeros(1, totalSimulations);

% Asignar valores a 'list' para cada grupo de simulaciones
for i = 1:length(fullpaths)
    list((i-1)*numSimulations + 1:i*numSimulations) = 8.0 + 0.1 * i;
end

% Procesar cada archivo de simulación
for i = 1:length(fullpaths)
    % Cargar datos de simulación
    data = load(fullpaths{i});
    simulations = data.simulations; % Ajustar según la estructura correcta
    
    % Indices para la ubicación en las matrices finales
    idx_start = (i-1) * numSimulations + 1;
    idx_end = i * numSimulations;
    
    % Extraer y vectorizar los datos necesarios
    lengths(idx_start:idx_end) = cellfun(@(sim) (sim.length(1)) * size(sim.lat, 1), simulations);
    widths(idx_start:idx_end) = cellfun(@(sim) (sim.width(1)) * size(sim.lat, 2), simulations);
    meanslips(idx_start:idx_end) = cellfun(@(sim) (mean(sim.slip(:))) , simulations);

end
figure()
subplot(2,1,1)
% Ajuste de regresión lineal
coefficients_lengths = polyfit(list, log10(lengths), 1);  % Ajuste lineal
fit_lengths = polyval(coefficients_lengths, list);  % Valores ajustados

% Dibujar la línea de regresión
scatter(list, log10(lengths), 'x');hold on;
plot(list, fit_lengths, 'b-', 'LineWidth', 1);

% Configuración del gráfico
title('Scatter Plot of Lengths vs Magnitude');
xlabel('Magnitude');
ylabel('Lengths');
grid on;
set(gca, 'FontSize', 12); % Cambiar tamaño de la fuente

% Crear la leyenda con la ecuación de la recta de regresión
legend({'Data', sprintf('Fit: %.2f + %.2f * Magnitude', coefficients_lengths(2), coefficients_lengths(1))}, 'Location', 'northwest');

% Scatter plot de anchos vs magnitud con regresión lineal
subplot(2,1,2)
scatter(list, log10(widths), 'x');
hold on;

% Ajuste de regresión lineal
coefficients_widths = polyfit(list, log10(widths), 1);  % Ajuste lineal
fit_widths = polyval(coefficients_widths, list);  % Valores ajustados

% Dibujar la línea de regresión
plot(list, fit_widths, 'b-', 'LineWidth', 1);

% Configuración del gráfico
title('Scatter Plot of Width vs Magnitude');
xlabel('Magnitude');
ylabel('Width');
grid on;
set(gca, 'FontSize', 12); % Cambiar tamaño de la fuente

% Crear la leyenda con la ecuación de la recta de regresión
legend({'Data', sprintf('Fit: %.2f + %.2f * Magnitude', coefficients_widths(2), coefficients_widths(1))}, 'Location', 'northwest');

% Configuración del gráfico
title('Scatter Plot of Width vs Magnitude');
xlabel('Magnitude');
ylabel('Mean Slip');
grid on;
set(gca, 'FontSize', 12); % Cambiar tamaño de la fuente

% Crear la leyenda con la ecuación de la recta de regresión
legend({'Data', sprintf('Fit: %.2f + %.2f * Magnitude', coefficients_slip(2), coefficients_slip(1))}, 'Location', 'northwest');
