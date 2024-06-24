fullpath1='/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_8.7_10000/simulations.mat';
fullpath2='/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_8.8_10000/simulations.mat';
fullpath3='/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_8.9_10000/simulations.mat';
fullpath4='/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_9.0_10000/simulations.mat';
fullpath5='/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_9.1_10000/simulations.mat';
fullpath6='/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_9.2_10000/simulations.mat';
fullpath7='/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_9.3_10000/simulations.mat';
list=[fullpath1 fullpath2 fullpath3 fullpath4 fullpath5 fullpath6 fullpath7];
lengths=[];
widths=[];
for i=1:7
    simulations=load(list(i));
    for k = 1:10000
    % Obtener los datos para la simulación i
    k=k+1;
    X_grid = simulations{1,k}.lon;   % Ajusta 'lon' según el campo correcto en tu estructura simulations
    Y_grid = simulations{1,k}.lat;   % Ajusta 'lat' según el campo correcto en tu estructura simulations
    Slip = simulations{1,k}.slip;    % Ajusta 'slip' según el campo correcto en tu estructura simulations
    length=sum(simulations{i,k}.length);
    width=sum(simulations{i,k}.width);
    lengths=[lengths width];
    widths=[widths length];
    end
end
list=ones(1,70000);
list(1:10000)=8.7;
list(10001:20000)=8.8;
list(20001:30000)=8.9;
list(30001:40000)=9.0;
list(40001:50000)=9.1;
list(50001:60000)=9.2;
list(60001:70000)=9.3;
%%
% Definir rutas de los archivos
fullpaths = {
    '/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_8.7_10000/simulations.mat',
    '/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_8.8_10000/simulations.mat',
    '/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_8.9_10000/simulations.mat',
    '/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_9.0_10000/simulations.mat',
    '/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_9.1_10000/simulations.mat',
    '/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_9.2_10000/simulations.mat',
    '/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_9.3_10000/simulations.mat'
};

% Número de simulaciones por archivo
numSimulations = 10000;

% Inicializar variables
totalSimulations = numSimulations * length(fullpaths);
lengths = zeros(1, totalSimulations);
widths = zeros(1, totalSimulations);
list = zeros(1, totalSimulations);

% Asignar valores a 'list' para cada grupo de simulaciones
for i = 1:length(fullpaths)
    list((i-1)*numSimulations + 1:i*numSimulations) = 8.6 + 0.1 * i;
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
end

% Crear scatter plot
%%

%%
%%
% Definir rutas de los archivos
fullpaths = {
    '/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_8.1_1000/simulations.mat',
    '/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_8.2_1000/simulations.mat',
    '/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_8.3_1000/simulations.mat',
    '/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_8.4_1000/simulations.mat',
    '/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_8.5_1000/simulations.mat',
    '/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_8.6_1000/simulations.mat',
    '/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_8.8_1000/simulations.mat',
    '/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_8.9_1000/simulations.mat',
    '/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_9.0_1000/simulations.mat',
    '/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_9.1_1000/simulations.mat',
    '/home/alex/Desktop/StochasticSlipGenerator/Output_data/Simulation_9.2_1000/simulations.mat'
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
    list((i-1)*numSimulations + 1:i*numSimulations) = 8.1 + 0.1 * i;
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
    subplot(2,2,1)
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
    subplot(2,2,2)
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
    % 
    subplot(2,2,3)
    scatter(list, log10(meanslips), 'x');
    hold on;
    
    % Ajuste de regresión lineal
    coefficients_slip = polyfit(list, log10(meanslips), 1);  % Ajuste lineal
    fit_slips = polyval(coefficients_slip, list);  % Valores ajustados
    
    % Dibujar la línea de regresión
    plot(list, fit_slips, 'b-', 'LineWidth', 1);

% Configuración del gráfico
title('Scatter Plot of Width vs Magnitude');
xlabel('Magnitude');
ylabel('Mean Slip');
grid on;
set(gca, 'FontSize', 12); % Cambiar tamaño de la fuente

% Crear la leyenda con la ecuación de la recta de regresión
legend({'Data', sprintf('Fit: %.2f + %.2f * Magnitude', coefficients_slip(2), coefficients_slip(1))}, 'Location', 'northwest');