function [] = plot_results(phi, Bx, By, xq, yq, z, d)
%plot_results Plot the calculated potential and field

figure(1);
surf(xq, yq, phi, 'EdgeColor', 'none', 'FaceColor', 'interp');
title("Electric Potential", 'FontSize', 18);
xlabel('x', 'FontSize', 16);
ylabel('y', 'FontSize', 16);
c = colorbar();
c.Label.String = "\phi";
c.Label.FontSize = 16;
colormap jet;
axis equal;
view(2);

figure(2);
B = sqrt(Bx.^2 + By.^2);
surf(xq, yq, B, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap('jet');
c = colorbar();
c.Label.String = '|B|';
c.Label.FontSize = 16;
title(['B-Field for z=' num2str(z) ' and d=' num2str(d)], 'FontSize', 18);
xlabel('x', 'FontSize', 16);
ylabel('y', 'FontSize', 16);
view(2);

end

