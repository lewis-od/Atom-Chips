%% Solve laplace's equation on a rectanguar region with 2 circles cut out
% Requires PDE Toolbox
clear all;

% Parameters
sigma = 5.96e7; % Conductivity of conductor [S*m^-1]
V0 = 1.0; % Voltage is V0 and -V0 at ends of conductor [V]
d = 0.1; % Thickness of conductor [m]
z = 0.5; % Distance from conductor [m]

params = [sigma, d, z];

% Specify geometry
R1 = [3, 4, -1, 1, 1, -1, 2, 2, -2, -2]; % 1x2 Rectangle
C1 = [1, 0, 0.75, 0.5]'; % Circle of radius 0.5 centreed at (0, 1)
C1 = [C1;zeros(length(R1) - length(C1),1)];
C2 = [1, 0, -0.75, 0.5]'; % Circle of radius 0.5 centreed at (0, -1)
C2= [C2;zeros(length(R1) - length(C2),1)];
R1 = R1';
gm = [R1, C1, C2];
sf = 'R1-C1-C2';
ns = char('R1', 'C1', 'C2')';

[phi, xq, yq, model, result] = calc_potential(V0, gm, ns, sf);

figure();
zero_index = size(xq);
zero_index = ceil(zero_index(2)/2);
for i = 1:4
    subplot(2, 2, i);
    params = [sigma, d, 0.1*i];
    [Bx, By] = calc_field(params, phi);
    B = sqrt(Bx.^2 + By.^2);
    yq_slice = yq(:, zero_index);
    B_slice = B(:, zero_index);
    plot(yq_slice, B_slice);
    title(['z=' num2str(params(3)) '; x=0']);
    xlabel('y');
    ylabel('B');
end

% Plot results
% plot_results(phi, Bx, By, xq, yq, z, d);
