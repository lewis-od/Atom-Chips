W = 1;
L = 100;
j = 10;
mu_0 = 4e-7 * pi;
B_bias = 1e-7;

x = linspace(-1.5, 1.5);
y = linspace(-1.5, 1.5);
z = 1.9192;
[x,y] = meshgrid(x,y);

Bx = eval_Bx(x,y,z, L/2, -W/2);
Bx = Bx + eval_Bx(x,y,z, -L/2, W/2);
Bx = Bx - eval_Bx(x,y,z, L/2, W/2);
Bx = Bx - eval_Bx(x,y,z, -L/2, -W/2);
Bx = ((mu_0*j)/(4*pi)) .* Bx;

Bz = eval_Bz(x,y,z, L/2, -W/2);
Bz = Bz + eval_Bz(x,y,z, -L/2, W/2);
Bz = Bz - eval_Bz(x,y,z, L/2, W/2);
Bz = Bz - eval_Bz(x,y,z, -L/2, -W/2);
Bz = ((mu_0*j)/(4*pi)) .* Bz;

B = sqrt(Bx.^2 + Bz.^2);
figure();
surf(x,y,B, 'EdgeColor', 'none');
xlabel('x');
ylabel('y')
view(2);
% hold on;
% Bx = Bx + B_bias;
% B = sqrt(Bx.^2 + Bz.^2);
% plot(z, B);
% legend({'No Bias', 'Bias'});