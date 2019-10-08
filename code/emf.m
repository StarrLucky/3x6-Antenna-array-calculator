function Z_emf = emf(x,y,d,lambda_emf,H)
i = sqrt(-1);
 INT1 = quadl(@(ksi) i * 30*( exp(-i*2*pi/lambda_emf *sqrt(ksi.^2 -2*ksi*y+y^2+d^2))/sqrt(ksi.^2 -2*ksi*y +y^2+d^2) +  exp(-i*2*pi/lambda_emf*sqrt(ksi.^2 +2*ksi*y+y^2+d^2))/sqrt(ksi.^2 +2.*ksi*y +y^2+d^2) -  2*cos(2*pi*y/lambda_emf)*exp(-i*(2*pi/lambda_emf)*sqrt(ksi.^2 + d^2))/sqrt(ksi.^2 + d^2))*sin(2*pi/lambda_emf * (x-H+ksi)),H-x,H);
 INT2 = quadl(@(ksi) i*30*(exp( -i*2*pi/lambda_emf *sqrt(ksi.^2 - 2*ksi*y+y^2 +d^2))/sqrt(ksi.^2 - 2*ksi*y+y^2 +d^2)+exp(-i*(2*pi/lambda_emf)*sqrt(ksi.^2 +2*ksi*y+y^2 +d^2))/sqrt(ksi.^2 +2.*ksi*y+y^2 +d^2) - 2*cos(2*pi*y/lambda_emf)*exp(-i* (2*pi/lambda_emf)*sqrt(ksi.^2 +d^2))/sqrt(ksi.^2 +d^2) )*sin(2*pi/lambda_emf * (x+H-ksi)),H,H+x);

Z_emf = INT1+ INT2;

clear INT1 INT2;
end
