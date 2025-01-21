figure;movegui('northeast');axis off;axis equal;set(gcf,'Renderer','OpenGL'); view3d rot;hold on;set(gcf,'color','white');
camorbit(0,0,'camera'); axis vis3d; view(-90,0);
h5 = scatter3(pcpts(pind&1,1),pcpts(pind&1,2), pcpts(pind&1,3),2,[0 0.2235 0.3705],'filled');
h6 = scatter3(pcpts(~pind,1),pcpts(~pind,2), pcpts(~pind,3),4,[0.8500 0.3250 0.0980],'filled');
% legend('Remained','Removed')
hold off


set(gcf,'Units','Inches','renderer','Painters');  
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
temp=['reduction',num2str(reducNum+1),'.pdf'];
saveas(gca,temp);