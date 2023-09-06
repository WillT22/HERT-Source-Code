
rI= zeros(1,8);
theta_c_inner= zeros(length(rI));
L_inner = zeros(length(rI));
Length_inner = 0.25; %cm
theta_m_inner= zeros(length(rI));
r_coll=1;
l_coll=6.85;

%Inner Ring Only
for i = 1:8
    rI(i)= 1; %cm
end


for i = 1:length(rI)
    for j = 1:length(rI)
        if i > j
            L_inner(i,j) = 0/0;
            theta_c_inner(i,j) = 0/0;
            theta_m_inner(i,j) = 0/0;
            
        else
            L_inner(i,j) = (j-i)*Length_inner;
            theta_m_inner(i,j) = atand((rI(i)+rI(j))/L_inner(i,j));
            theta_c_inner(i,j) = atand(abs((rI(i)-rI(j)))/L_inner(i,j));
        end
    end
end


for i = 1:length(rI)
    for j = 1:length(rI)
        if i > j
            theta_m_inner(i,j) = 0/0;
            
        else
            theta_m_inner(i,j) = atand((rI(i)+rI(j))/L_inner(i,j));
        end
    end
end

for i = 1:length(rI)
    theta_c_coll(i) = atand(abs((r_coll-rI(i)))/(l_coll+L_inner(1,i)));
    theta_m_coll(i) = atand((r_coll+rI(i))/(l_coll+L_inner(1,i)));
end

if max(min(theta_c_inner,[],1)) == theta_c_inner(1,8)
    Gtest = 0.5*(pi^2)*((rI(8)^2+rI(1)^2+L_inner(1,8)^2)-(((rI(8)^2+rI(1)^2+L_inner(1,8)^2)^2-4*(rI(8)^2)*(rI(1)^2))^0.5));
    %     x = rI(8)^2+rI(1)^2+L_inner(1,8)^2;
    %     y = 4*(rI(8)^2)*(rI(1)^2);
    %     Gtest1_inner = 0.5*(pi^2)*(x-(x^2-y)^0.5);
end
% Gtest1_inner = Gtest *30/180
Gtest1_inner = 0.5*(pi^2)*((r_coll^2+rI(1)^2+l_coll^2)-(((r_coll^2+rI(1)^2+l_coll^2)^2-4*(r_coll^2)*(rI(1)^2))^0.5));
Gtest2_inner = 0.5*(pi^2)*((r_coll^2+rI(8)^2+(l_coll+L_inner(1,8))^2)-(((r_coll^2+rI(8)^2+(l_coll+L_inner(1,8))^2)^2-4*(r_coll^2)*(rI(8)^2))^0.5));
fprintf('Geometric_Factor (Front Coll. to Detector 1) = %7.5f  cm^2 sr\n \n',Gtest1_inner)
fprintf('Geometric_Factor_Inner (Front Coll. to Detector 15)= %7.5f  cm^2 sr\n \n',Gtest2_inner)

% fprintf('Geometric_Factor_Inner = %7.5f  cm^2 sr \n \n',Gtest)
% fprintf('Radius(cm)    Detector 1 Detector 3 Detector 5 Detector 7 Detector 9 Detector 11 Detector 13 Detector 15 \n')
% fprintf('             %7.2f     %7.2f    %7.2f    %7.2f    %7.2f    %7.2f     %7.2f     %7.2f \n',rI(1:8))
% fprintf('\n')
% fprintf('Length(cm)    Detector 1 Detector 3 Detector 5 Detector 7 Detector 9 Detector 11 Detector 13 Detector 15 \n')
% for i = (1:length(rI)-1)
%     fprintf('Detector %2.0f: %7.2f     %7.2f    %7.2f    %7.2f    %7.2f    %7.2f     %7.2f     %7.2f \n',1+(i-1)*2,L_inner(i,1:8))
% end
% fprintf('\n')
% fprintf('Theta_c(deg) Detector 1 Detector 3 Detector 5 Detector 7 Detector 9 Detector 11 Detector 13 Detector 15 \n')
% for i = (1:length(rI)-1)
%     fprintf('Detector %2.0f: %7.2f     %7.2f    %7.2f    %7.2f    %7.2f    %7.2f     %7.2f     %7.2f \n',1+(i-1)*2,theta_c_inner(i,1:8))
% end
% fprintf('Theta_m(deg) Detector 1 Detector 3 Detector 5 Detector 7 Detector 9 Detector 11 Detector 13 Detector 15 \n')
% for i = (1:length(rI)-1)
%     fprintf('Detector %2.0f: %7.2f     %7.2f    %7.2f    %7.2f    %7.2f    %7.2f     %7.2f     %7.2f \n',1+(i-1)*2,theta_m_inner(i,1:8))
% end


%Outer Ring Included

% 
% rO= zeros(1,9);
% theta_c_outer= zeros(length(rO));
% L_outer = zeros(length(rO));
% Length_outer = 0.25; %cm
% theta_m_outer= zeros(length(rO));
% 
% for i = 1:9
%     rO(i)= 2; %cm
% end
% 
% for i = 1:length(rO)
%     for j = 1:length(rO)
%         if i > j
%             L_outer(i,j) = 0/0;
%             theta_m_outer(i,j) = 0/0;
%             theta_c_outer(i,j) = 0/0;
%             
%         else
%             L_outer(i,j) = (j-i)*Length_outer;
%             theta_m_outer(i,j) = atand((rO(i)+rO(j))/L_outer(i,j));
%             theta_c_outer(i,j) = atand(abs((rO(i)-rO(j)))/L_outer(i,j));
%         end
%     end
% end


% if max(min(theta_c_outer,[],1)) == theta_c_outer(1,9)
    
%     Gtest1_inner = 0.5*(pi^2)*((r_coll^2+rO(1)^2+(l_coll+L_outer(1,8))^2)-(((r_coll^2+rO(1)^2+(l_coll+L_outer(1,8))^2)^2-4*(r_coll^2)*(rO(1)^2))^0.5));
    %     x = rO(9)^2+rO(1)^2+L_outer(1,9)^2;
    %     y = 4*(rO(9)^2)*(rO(1)^2);
    %     Gtest1_outer = 0.5*(pi^2)*(x-(x^2-y)^0.5)
% end
%Gtest1_outer = 0.5*(pi^2)*((rO(9)^2+rO(1)^2+9.375^2)-(((rO(9)^2+rO(1)^2+9.375^2)^2-4*(rO(9)^2)*(rO(1)^2))^0.5))
% Gtest1_outer = Gtest *30/180
% fprintf('Geometric_Factor_Outer = %7.5f  cm^2 sr\n \n',Gtest)
% fprintf('Radius(cm)    Detector 1 Detector 3 Detector 5 Detector 7 Detector 9 Detector 11 Detector 13 Detector 15  Detector 17 \n')
% fprintf('             %7.2f     %7.2f    %7.2f    %7.2f    %7.2f    %7.2f     %7.2f     %7.2f     %7.2f \n',rO(1:9))
% fprintf('\n')
% fprintf('Length(cm)    Detector 1 Detector 3 Detector 5 Detector 7 Detector 9 Detector 11 Detector 13 Detector 15  Detector 17 \n')
% for i = (1:length(rO))
%     fprintf('Detector %2.0f: %7.2f     %7.2f    %7.2f    %7.2f    %7.2f    %7.2f     %7.2f     %7.2f     %7.2f \n',1+(i-1)*2,L_outer(i,1:9))
% end
% fprintf('\n')
% fprintf('Theta_c(deg) Detector 1 Detector 3 Detector 5 Detector 7 Detector 9 Detector 11 Detector 13 Detector 15  Detector 17\n')
% for i = (1:length(rO))
%     fprintf('Detector %2.0f: %7.2f     %7.2f    %7.2f    %7.2f    %7.2f    %7.2f     %7.2f     %7.2f     %7.2f \n',1+(i-1)*2,theta_c_outer(i,1:9))
% end
% fprintf('Theta_m(deg) Detector 1 Detector 3 Detector 5 Detector 7 Detector 9 Detector 11 Detector 13 Detector 15  Detector 17\n')
% for i = (1:length(rO))
%     fprintf('Detector %2.0f: %7.2f     %7.2f    %7.2f    %7.2f    %7.2f    %7.2f     %7.2f     %7.2f     %7.2f\n',1+(i-1)*2,theta_m_outer(i,1:9))
% end


