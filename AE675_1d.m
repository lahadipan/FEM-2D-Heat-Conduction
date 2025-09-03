clc; clear;close all;
k = 50; % Thermal conductivity of steel (W/m-K)
Q_h = 1000; % Heat source (W/m^2)
T_top = 25; % Temperature at top edge

% Create 6x6 grid (5x5 divisions)
nx = 6; ny = 6;
x = linspace(0, 1, nx);
y = linspace(0, 1, ny);
[X, Y] = meshgrid(x, y);
nodes = [X(:), Y(:)];

% Generate elements (4 triangles per square)
elements = [];
for i = 1:ny-1
    for j = 1:nx-1
        % Node indices for current square
        A = (i-1)*nx + j;
        B = (i-1)*nx + j+1;
        C = i*nx + j;
        D = i*nx + j+1;
        O = (nodes(A,:)+nodes(B,:)+nodes(C,:)+nodes(D,:))/4;
        nodes = [nodes;O];
        O_pos = size(nodes,1);
        if isequal(O, [0.5, 0.5])
            cent_pos = O_pos;
        end
        
        % Create 4 triangles
        elements = [elements; 
                   A B O_pos;
                   B D O_pos;
                   D C O_pos;
                   C A O_pos];
    end
end


% Initialize global stiffness matrix and force vector
n_nodes = size(nodes, 1);
K = zeros(n_nodes, n_nodes);
F = zeros(n_nodes, 1);
F(cent_pos) = F(cent_pos)+Q_h;

% Assemble stiffness matrix
for e = 1:size(elements, 1)
    node_ids = elements(e, :);
    x = nodes(node_ids, 1);
    y = nodes(node_ids, 2);

    % Calculate element stiffness
    A = abs(polyarea(x, y));
    J = [x(1)-x(3) x(2)-x(3); y(1)-y(3) y(2)-y(3)];
    Q = [1 0 -1; 0 1 -1];
    B = transpose(inv(J))*Q;
    Ke = k * A * transpose(B) * B;

    % Assemble to global matrix
    K(node_ids, node_ids) = K(node_ids, node_ids) + Ke;
end

top_corner_nodes = find(nodes(:,2) == 1);
fixed_nodes = top_corner_nodes;
free_nodes = setdiff(1:n_nodes, fixed_nodes);

% Partition matrices
Kff = K(free_nodes, free_nodes);
Kfc = K(free_nodes, fixed_nodes);
Ff = F(free_nodes);

% Known temperatures
Tc = T_top * ones(length(fixed_nodes), 1);

% Solve for unknown temperatures
Tf = Kff \ (Ff - Kfc*Tc);

% Combine results
T = zeros(n_nodes, 1);
T(free_nodes) = Tf;
T(fixed_nodes) = Tc;

%% Visualization: Temperature Contour Plot
figure;
% Create grid for smoother contour plot
xi = linspace(0, 1, 100);
yi = linspace(0, 1, 100);
[XI,YI] = meshgrid(xi,yi);
TI = griddata(nodes(:,1), nodes(:,2), T, XI, YI, 'cubic');

contourf(XI, YI, TI, 20, 'LineColor', 'none');
title('Temperature Contour Plot');
xlabel('x (m)'); ylabel('y (m)');
colorbar;
colormap(jet);
axis equal;
hold on;
% Mark heat source location
plot(0.5, 0.5, 'rp', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
hold off;

%% Temperature Variation Along Centerline
centerline_y = 0.5;
tolerance = 1e-6;

% Find all nodes within tolerance of y=0.5
on_centerline = abs(nodes(:,2) - centerline_y) < tolerance;

% Extract and sort these nodes
centerline_nodes = find(on_centerline);
[centerline_x, idx] = sort(nodes(centerline_nodes,1));
centerline_T = T(centerline_nodes(idx));

% Create high-resolution interpolation
xi = linspace(0, 1, 500)'; % 500 points for smoothness
yi = centerline_y * ones(size(xi));

% Use scatteredInterpolant for better handling of unstructured data
F = scatteredInterpolant(nodes(:,1), nodes(:,2), T, 'natural');
TI = F(xi, yi);

figure;
plot(xi, TI, 'b-', 'LineWidth', 2);
hold on;
plot(centerline_x, centerline_T, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);

% Find and mark exact center point
[~, center_idx] = min(abs(xi - 0.5));
plot(xi(center_idx), TI(center_idx), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

title(sprintf('Temperature Variation Along y = %.1f m', centerline_y));
xlabel('x (m)'); ylabel('Temperature (°C)');
grid on;
xlim([0 1]);
legend('Interpolated', 'Node Values', 'Heat Source', 'Location', 'best');
hold off;