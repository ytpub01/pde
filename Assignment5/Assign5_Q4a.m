ezplot('(1+x+(1/2)*(x^2-y^2)+(1/6)*(x^3-3*x*y^2))^2+(y+x*y+(1/6)*(3*x^2*y-y^3))^2-1');
x=[-pi:.1:pi];
hold on
z=-1.25+1.25*exp(i*x);
plot(real(z), imag(z), 'red');
z=-0.62*((3/2)-2*exp(i*x)+(1/2)*exp(i*2*x));
plot(real(z), imag(z), 'black');