function plotsun()
    Rsun_km = 696340;
    AU = 149597870.7;
    Rsun_AU = Rsun_km / AU;
    
    [Xs, Ys, Zs] = sphere(30);      % mesh per la sfera
    surf(Rsun_AU*Xs, Rsun_AU*Ys, Rsun_AU*Zs, ...
        'FaceColor',[1 1 0], ...    % giallo
        'EdgeColor','none', ...
        'FaceLighting','gouraud')
    
    lightangle(45,30)
end