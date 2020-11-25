function [ pFpx_nom, pFpu_nom, pFpwt_nom, phpx_nom ] = Problem_1(L,x_nom,u_nom)


pFpx_nom = [ 0, 0, -u_nom(1)*sin(x_nom(3)), 0, 0,           0;
             0, 0,  u_nom(1)*cos(x_nom(3)), 0, 0,           0;
             0, 0,             0,           0, 0,           0;
             0, 0,             0,           0, 0, -u_nom(3)*sin(x_nom(6));
             0, 0,             0,           0, 0,  u_nom(3)*cos(x_nom(6));
             0, 0,             0,           0, 0,           0              ];
                         
pFpu_nom = [   cos(x_nom(3)),                   0,                     0,        0;
               sin(x_nom(3)),                   0,                     0,        0;
             (1/L)*tan(u_nom(2)), (u_nom(1)/L)*(sec(u_nom(2))^2),      0,        0;
                      0,                        0,                cos(x_nom(6)), 0;
                      0,                        0,                sin(x_nom(6)), 0;
                      0,                        0,                     0,        1  ];
                                  
pFpwt_nom = eye(6);

dh1dx1 = (1/(1+((x_nom(5)-x_nom(2))/(x_nom(4)-x_nom(1)))^2))*((x_nom(5)-x_nom(2))/((x_nom(4)-x_nom(1))^2));
dh2dx1 = (0.5)*((((x_nom(1)-x_nom(4))^2)+((x_nom(2)-x_nom(5))^2))^(-.5))*(2*(x_nom(1)-x_nom(4))*(1));
dh3dx1 = (1/(1+((x_nom(2)-x_nom(5))/(x_nom(1)-x_nom(4)))^2))*((x_nom(5)-x_nom(2))/((x_nom(1)-x_nom(4))^2));
dh4dx1 = 0;
dh5dx1 = 0;

dh1dx2 = (1/(1+((x_nom(5)-x_nom(2))/(x_nom(4)-x_nom(1)))^2))*((-1)/((x_nom(4)-x_nom(1))));
dh2dx2 = (0.5)*((((x_nom(1)-x_nom(4))^2)+((x_nom(2)-x_nom(5))^2))^(-.5))*(2*(x_nom(2)-x_nom(5))*(1));
dh3dx2 = (1/(1+((x_nom(2)-x_nom(5))/(x_nom(1)-x_nom(4)))^2))*((1)/((x_nom(1)-x_nom(4))));
dh4dx2 = 0;
dh5dx2 = 0;

dh1dx3 = -1;
dh2dx3 = 0;
dh3dx3 = 0;
dh4dx3 = 0;
dh5dx3 = 0;

dh1dx4 = (1/(1+((x_nom(5)-x_nom(2))/(x_nom(4)-x_nom(1)))^2))*((x_nom(2)-x_nom(5))/((x_nom(4)-x_nom(1))^2));
dh2dx4 = (0.5)*((((x_nom(1)-x_nom(4))^2)+((x_nom(2)-x_nom(5))^2))^(-.5))*(2*(x_nom(1)-x_nom(4))*(-1));
dh3dx4 = (1/(1+((x_nom(2)-x_nom(5))/(x_nom(1)-x_nom(4)))^2))*((x_nom(5)-x_nom(2))/((x_nom(1)-x_nom(4))^2));
dh4dx4 = 1;
dh5dx4 = 0;

dh1dx5 = (1/(1+((x_nom(5)-x_nom(2))/(x_nom(4)-x_nom(1)))^2))*((1)/((x_nom(4)-x_nom(1))));
dh2dx5 = (0.5)*((((x_nom(1)-x_nom(4))^2)+((x_nom(2)-x_nom(5))^2))^(-.5))*(2*(x_nom(2)-x_nom(5))*(-1));
dh3dx5 = (1/(1+((x_nom(2)-x_nom(5))/(x_nom(1)-x_nom(4)))^2))*((-1)/((x_nom(4)-x_nom(1))));
dh4dx5 = 0;
dh5dx5 = 1;

dh1dx6 = 0;
dh2dx6 = 0;
dh3dx6 = -1;
dh4dx6 = 0;
dh5dx6 = 0;


phpx_nom = [ dh1dx1, dh1dx2, dh1dx3, dh1dx4, dh1dx5, dh1dx6;
             dh2dx1, dh2dx2, dh2dx3, dh2dx4, dh2dx5, dh2dx6;
             dh3dx1, dh3dx2, dh3dx3, dh3dx4, dh3dx5, dh3dx6;
             dh4dx1, dh4dx2, dh4dx3, dh4dx4, dh4dx5, dh4dx6;
             dh5dx1, dh5dx2, dh5dx3, dh5dx4, dh5dx5, dh5dx6  ];
         
         
end
    
