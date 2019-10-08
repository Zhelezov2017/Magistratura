function ff = pcolorReImModule(x, nameX, y, nameY, z, nameZ)

figure
axes('FontSize',16)
pcolor(x, y, real(z))
% xlabel(nameX,'FontSize',16);
% ylabel(nameY,'FontSize',16);
% title(['Re ' nameZ],'FontSize',16)
shading interp
colorbar('location','eastoutside')

% figure
% axes('FontSize',16)
% mesh(x, y, real(z))
% shading interp
% colorbar('location','eastoutside')
% 
% figure
% axes('FontSize',16)
% mesh(x, y, imag(z))
% shading interp
% colorbar('location','eastoutside')

%figure
%axes('FontSize',16)
%pcolor(x, y, imag(z))
% xlabel(nameX,'FontSize',16);
% ylabel(nameY,'FontSize',16);
% title(['Im ' nameZ],'FontSize',16)
%shading interp
%colorbar('location','eastoutside')



%figure
%axes('FontSize',16)
%pcolor(x, y, abs(z))
% colormap bone
%colormap hot
% colormap pink
% xlabel(nameX,'FontSize',16);
% ylabel(nameY,'FontSize',16);
% title(['|' nameZ '|'],'FontSize',16)
%shading interp
%colorbar('location','eastoutside')

% figure
% axes('FontSize',16)
% pcolor(x, y, abs(z))
% % xlabel(nameX,'FontSize',16);
% % ylabel(nameY,'FontSize',16);
% % title(['Re ' nameZ],'FontSize',16)
% shading interp
% colorbar('location','eastoutside')

% figure
% axes('FontSize',16)
% mesh(x, y, abs(z))
% % colormap bone
% % colormap hot
% % colormap pink
% % xlabel(nameX,'FontSize',16);
% % ylabel(nameY,'FontSize',16);
% % title(['|' nameZ '|'],'FontSize',16)
% shading interp
% colorbar('location','eastoutside')