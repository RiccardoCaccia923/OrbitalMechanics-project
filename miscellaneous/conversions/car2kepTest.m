clear
clc
mu = astroConstants(13);
% [~,~,s0] = orbitselect('POL');
s0 = [
   -18960.3685;
   -16881.0146;
     8353.2302;
       0.6993;
      -0.4907;
      -1.6214
];

kep = car2kep(mu,s0);
s = kep2car(mu,kep);

disp( s - s0)