clear;clc;close;
x = -1;%randn([1,1]);
y = -1;%randn([1,1]);

s = [[x,y]];

a = -1;



eps = inf;
step = .01;
j = 0;

while eps > .00001 && j <2000
    j = j+1;
    
    
    b = [3*x^2 + a*y; 2*y + a*x];
    Hess = [6*x, a; a, 2];
    dxy = Hess\-b;

    x = x + step*dxy(1);
    y = x + step*dxy(2);

    eps = norm(b);
    s = [s;[x,y]];
end


m = 2;
[xx,yy]=meshgrid(linspace(-m,m,100),linspace(-m,m,100));



plot3(s(:,1),s(:,2),s(:,1).^3 + s(:,2).^2 + a.*s(:,1).*s(:,2),'r');hold on;
surf(xx,yy, xx.^3 + yy.^2 + a.*xx.*yy);

zlim([-1.5,1.5])