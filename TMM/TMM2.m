clear all
m = 1;
a = 0.14;
b = 0.18;
c = 2.3
alfa = pi*7/180

 X = [0.14 0.14];
 Y = [0.17 0]
% Y(2) = 0;
i = 3;
% R1 = sqrt(0.25^2 + 0.25^2);

for t = 0:pi/100:2*pi;
    X(i) = a * cos(t);
    Y(i) = b * sin(t);
    i = i + 1
end;
    for t = pi:pi/100:2*pi;
        X(i) = 0.1*cos(t) + 0.24;
        Y(i) = 0.1*sin(t);
        i = i + 1;
    end;
   for t = 0:pi/100:5*pi/3;
        X(i) = 0.34*cos(t);
        Y(i) = 0.34*sin(t);
        i = i + 1;
    end;    
% i = 1;
% for t = 0 : pi / 100 : 2 * pi ;
%     X = a * cos(t-alfa);
%     Y = b * sin(t);
%     Y(i) = (b * sin(t) + a * cos(t) / cos(alfa)) / (cos(alfa) + ((sin(alfa)) * sin(alfa) )/ cos(alfa));
%     X(i) = (a * cos(t) - Y(i) * sin(alfa)) / cos(alfa);
%     i = i + 1
% end;

% R = 0.5
% t = atan(-2/5);
% x0 = -R * cos(t);
% y0 = -R * sin(t);

% for k = 0:pi/100:2*pi;
%     X(i) = R * cos(k) + x0;
%     Y(i) = R * sin(k) + y0;
%     i = i + 1;
% end;


% k = 0.5
% j = 1
% for fi = 0:pi/100:100*pi/300
%        r(j) = k * fi;
%        X(i) = r(j) * cos(fi);
%        Y(i) = r(j) * sin(fi);
%        i = i + 1;
%        j = j + 1;
% end;
% 
% Rn = r(j-1) * cos(fi)
% R  = sqrt((X(i - 1) + 0.5)^2 + Y(i - 1)^2);
% for k = pi/6:pi/100:2*pi;
%     X(i) = R * cos(k) - 0.5;
%     Y(i) = R * sin(k);
%     i = i + 1;
% end;
% k = Rn:(Rn-0.01):-R
% for k = Rn:(Rn-0.01):-R
%     X(i) = k;
%     Y(i) = sqrt(R*R - k*k);
%     i = i + 1;
% end;

figure;
plot (X,Y);
grid on;
    

    
    






% for i = -17:17
%     q1 = i * pi / 18;
%     for j = -10:10
%         q2 = j * pi / 40;
%         for k = -18:8
%             q3 = k * pi / 72
%             X(m) = (-1) * sin(q1) * (c * cos(q2 + q3) - b * sin(q2));
%             Y(m) = cos(q1) * (c * cos(q2 + q3) - b * sin(q2));
%             Z(m) = c * sin(q2 + q3) + b * cos(q2) + a;
%             m = m + 1;
%         end
%     end
% end
% A = [X;Y;Z]
% mesh(A)
% %waterfall(A)
% figure;
% title 'X,Y and Z'
% xlabel ('x');
% ylabel ('y');
% zlabel ('z');
% plot3(X,Y,Z)
% grid on;
% figure;
% title 'X and Y'
% figure;
% plot (X,Y)
% title 'x and z'
% figure;
% plot(X,Z)
% title 'y and z'
% figure;
% plot(Y,Z)