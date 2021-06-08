function [x_p y_p] = project_pt(xx,yy,aa,b)

x_p = (-b+(yy+xx/aa))/(aa+1/aa);
y_p = -x_p/aa+(yy+xx/aa);
