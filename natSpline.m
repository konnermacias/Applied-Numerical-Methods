% natSpline.m
function [mat]=natSpline(n, x, a)
    h(1) = x(2)-x(1); % use reference in book
    l(1) = 1;
    mu(1) = 0;
    z(1) = 0;
    for i=2:n % adjusted for MATLAB's indexing
        h(i) = x(i+1)-x(i);
        alpha = (3/h(i))*(a(i+1)-a(i)) - (3/h(i-1))*(a(i)-a(i-1));
        l(i) = 2*(x(i+1)-x(i-1))-h(i-1)*mu(i-1);
        mu(i)= h(i)/l(i);
        z(i) = (alpha-h(i-1)*z(i-1))/l(i);
    end
    l(n+1) = 1;
    z(n+1) = 0;
    c(n+1) = 0;
    for j=n:-1:1
        c(j) = z(j)-mu(j)*c(j+1);
        b(j) = (a(j+1)-a(j))/h(j)-(h(j)*(c(j+1)+2*c(j)))/3;
        d(j) = (c(j+1)-c(j))/(3*h(j));
    end
    
     mat = [a(1:n);b(1:n);c(1:n);d(1:n)]; % created matrix of all a,b,c,d vals
     hold on; % allows for multiple plots
     for i=1:n
        spl=mat(1:4,i);
        x1=x(i);
        x2=x(i+1);
        x_pts=linspace(x1,x2,20);
        f=spl(1) + spl(2)*(x_pts-x(i))+spl(3)*(x_pts-x(i)).^2+spl(4)*(x_pts-x(i)).^3;
        disp(spl);
        plot(x_pts,f);
        plot(x,a);
     end
end

% test.m
x_coord = [12,39,49,55,65,77,89,111,139,159,192,235,261,286,313,365,380,408,415,420,422,424,427,428,443,452,469,491,506,513,525,535,540,549,559,566,567,572,578,579];
y_coord = [317,255,211,177,153,141,134,129,128,131,133,138,141,148,148,132,128,115,98,85,81,78,77,56,34,24,12,9,13,20,32,45,70,75,91,120,148,171,192,205];
y_coord = 317 - y_coord; % flip image

% Break up into multiple splines!
hold on;
x_1to18 = x_coord(1:18);
y_1to18 = y_coord(1:18);
natSpline(17,x_1to18,y_1to18); % spline 1
x_18to24 = x_coord(18:24);
y_18to24 = y_coord(18:24);
natSpline(6,x_18to24,y_18to24); % spline 2
x_24to28 = x_coord(24:28);
y_24to28 = y_coord(24:28);
natSpline(4,x_24to28,y_24to28); % spline 3
x_28to36 = x_coord(28:36);
y_28to36 = y_coord(28:36);
natSpline(8,x_28to36,y_28to36); % spline 4
x_36to40 = x_coord(36:40);
y_36to40 = y_coord(36:40);
natSpline(4,x_36to40,y_36to40); % spline 5
            
        
        
