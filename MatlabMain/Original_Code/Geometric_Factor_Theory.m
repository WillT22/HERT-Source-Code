%Resets all variables and values in MATLAB
clear;
close all;
clc;
%% Geometric Factor- Theory
% Disc 1: First Detector
% Disc 2: Last Collimator Tooth
% Disc 3: First Collimator Tooth

L_12 = 6.0;%cm distance between first and last collimator teeth
L_23 = 0.3; %cm distance between last collimator tooth and first detector
r1 = 0.9; %cm radius of first detector
r2 = 0.9; %cm radius of last collimator tooth
r3 = 1.0; %cm radius of first collimator tooth

L_13 = L_12+L_23;

theta_m_12 = atan((r1+r2)/L_12);
theta_m_13 = atan((r1+r3)/L_13);
theta_m_23 = atan((r2+r3)/L_23);

theta_c_12 = atan(abs(r1-r2)/L_12);
theta_c_13 = atan(abs(r1-r3)/L_13);
theta_c_23 = atan(abs(r2-r3)/L_23);


G13 = 0.5*(pi^2)*((r1^2+r3^2+L_13^2)-(((r1^2+r3^2+L_13^2)^2-4*(r1^2)*(r3^2))^0.5));
G12 = 0.5*(pi^2)*((r1^2+r2^2+L_12^2)-((r1^2+r2^2+L_12^2)^2-4*(r1^2)*(r2^2))^0.5);
G23 = 0.5*(pi^2)*((r2^2+r3^2+L_23^2)-((r2^2+r3^2+L_23^2)^2-4*(r2^2)*(r3^2))^0.5);
 
if theta_c_12 >= theta_c_13
    G13 = 0.5*(pi^2)*((r1^2+r3^2+L_13^2)-(((r1^2+r3^2+L_13^2)^2-4*(r1^2)*(r3^2))^0.5));
    fprintf('Geometric_Factor_13 (Front Coll. to First Detector)= %7.5f  cm^2 sr\n \n',G13)
    FOV = 2*theta_m_13*180/pi
    
elseif theta_m_13 >= theta_m_12
    
    G12 = 0.5*(pi^2)*((r1^2+r2^2+L_12^2)-((r1^2+r2^2+L_12^2)^2-4*(r1^2)*(r2^2))^0.5);
    fprintf('Geometric_Factor_12 (First Tooth to Last Tooth) = %7.5f  cm^2 sr\n \n',G12)
    FOV = 2*theta_m_12*180/pi
    
elseif theta_c_12 >= theta_m_13
    G23 = 0.5*(pi^2)*((r2^2+r3^2+L_23^2)-((r2^2+r3^2+L_23^2)^2-4*(r2^2)*(r3^2))^0.5);
    fprintf('Geometric_Factor (Last Tooth to First Detector) = %7.5f  cm^2 sr\n \n',G23)
    
else
    theta_a = atan(((L_23*r1^2+L_12*r3^2-L_13*r2^2)^0.5)/(L_12*L_23*L_13));
    G13 = 0.5*(pi^2)*((r1^2+r3^2+L_13^2)-(((r1^2+r3^2+L_13^2)^2-4*(r1^2)*(r3^2))^0.5));
    
    Z12= Zij(theta_a,L_12,r1,r2,theta_c_12,theta_m_12,r2);
    Z13= Zij(theta_a,L_13,r1,r3,theta_c_13,theta_m_13,r2);
    Z23= Zij(theta_a,L_23,r2,r3,theta_c_23,theta_m_23,r2);
    
    G123=G13-pi^2*r2^2*sin(theta_a)^2+Z23+Z12-Z13;
    
    fprintf('Geometric_Factor (Three Disc Telescope) = %7.5f  cm^2 sr\n \n',G123)

end


% %Inner Ring Only
% for i = 1:8
%     rI(i)= 1; %cm
% end
% 
% 
% for i = 1:length(rI)
%     for j = 1:length(rI)
%         if i > j
%             L_inner(i,j) = 0/0;
%             theta_c_inner(i,j) = 0/0;
%             theta_m_inner(i,j) = 0/0;
%             
%         else
%             L_inner(i,j) = (j-i)*Length_inner;
%             theta_m_inner(i,j) = atand((rI(i)+rI(j))/L_inner(i,j));
%             theta_c_inner(i,j) = atand(abs((rI(i)-rI(j)))/L_inner(i,j));
%         end
%     end
% end
% 
% 
% for i = 1:length(rI)
%     theta_c_coll(i) = atand(abs((r_coll-rI(i)))/(l_coll+L_inner(1,i)));
%     theta_m_coll(i) = atand((r_coll+rI(i))/(l_coll+L_inner(1,i)));
% end
% 
% if max(min(theta_c_inner,[],1)) == theta_c_inner(1,8)
%     Gtest = 0.5*(pi^2)*((rI(8)^2+rI(1)^2+L_inner(1,8)^2)-(((rI(8)^2+rI(1)^2+L_inner(1,8)^2)^2-4*(rI(8)^2)*(rI(1)^2))^0.5));
%     %     x = rI(8)^2+rI(1)^2+L_inner(1,8)^2;
%     %     y = 4*(rI(8)^2)*(rI(1)^2);
%     %     Gtest1_inner = 0.5*(pi^2)*(x-(x^2-y)^0.5);
% end
% % Gtest1_inner = Gtest *30/180