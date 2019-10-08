function ff = videoReal(t, w0, x, nameX, y, nameY, z, nameZ)

figure
for j = 1:size(t,2)

    H_zz = z.* exp(i * w0 * t(j));

pcolor(x, y, real(H_zz))
% colormap hot
xlabel(nameX);
ylabel(nameY);
title(nameZ)
shading interp

    F = getframe;
end

movie(F)

% movie2avi(F,'film1.avi', 'compression', 'Indeo5');
