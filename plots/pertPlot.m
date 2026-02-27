function pertPlot(x,y,z,C,nPeriods)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Color','w')
hold on
axis equal                  

view(325, 40)

% colored line: 
surface([x x], [y y], [z z], [C C], ...
        'FaceColor','none', 'EdgeColor','interp', 'LineWidth', 1.6);

colormap("parula");               

cb = colorbar; 
cb.Label.String = 'Periods';
clim([0, max(C)])                   
title(['Orbit representation – num2str(nPeriods)', ' periods'], ...
      'Interpreter', 'latex');
warning('off')
end