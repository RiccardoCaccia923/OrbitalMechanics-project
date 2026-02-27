function plotEarth()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotEarth: this function creates a sphere covered with earth layout in
%            order to complete orbit plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate the sphere
[xT,yT,zT] = sphere(200);   
rT = astroConstants(23);    

% Load Earth image
earthImg = imread('EarthTexture.jpg');
earthImg = flipud(earthImg);      % fix nord axis as z positive

% Create surface with Earth texture
surface(xT*rT, yT*rT, zT*rT, ...
        'CData', earthImg, ...
        'FaceColor','texturemap', ...
        'EdgeColor','none')
axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')

end