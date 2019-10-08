function ff = drawnCountor(a_0, step_x, cylXY, N_cylinders)

%%%% отрисовка контуров цилиндров
x1 = [- a_0 :step_x/100: a_0+step_x/100];
y1 = sqrt(a_0^2 - x1.^2);
hold on
for ncyl = 1:N_cylinders
    plot(real(x1 - cylXY(ncyl,1)) / a_0, real(y1 - cylXY(ncyl,2)) / a_0, 'w', 'LineWidth',1);
    plot(real(x1 - cylXY(ncyl,1)) / a_0, real(-y1 - cylXY(ncyl,2)) / a_0, 'w', 'LineWidth',1)
end
hold off