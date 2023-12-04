%% Creating Function to find Max and Min Geometric Factors
function G3_whole = findG3whole(L_12,L_23,L_13,r1,r2,r3)
%Creating angles for simplifying criteria Eq. 14 in Sullivan
%Angles of incidence such that circles project onto each other depending on location of incidence of the particles
%maximum angle that the particle will always hit within both circles
theta_c_12 = atan(abs(r1-r2)/L_12); %first to last collimator tooth
theta_c_13 = atan(abs(r1-r3)/L_13); %first tooth to first detector
theta_c_23 = atan(abs(r2-r3)/L_23); %last tooth to detector

%maximum angle such that particles may hit both circles depending on incidence location
theta_m_12 = atan((r1+r2)/L_12); %first to last collimator tooth
theta_m_13 = atan((r1+r3)/L_13); %first tooth to first detector
theta_m_23 = atan((r2+r3)/L_23); %last tooth to detector

%geometric factor of first collimator tooth and first detector (high E particles)
G13 = 0.5*(pi^2)*((r1^2+r3^2+L_13^2)-((r1^2+r3^2+L_13^2)^2-4*(r1^2)*(r3^2))^0.5);
%geometric factor of last collimator tooth and first detector (low E particles)
G12 = 0.5*(pi^2)*((r1^2+r2^2+L_12^2)-((r1^2+r2^2+L_12^2)^2-4*(r1^2)*(r2^2))^0.5); 
%geometric factor of first to last collimator tooth
G23 = 0.5*(pi^2)*((r2^2+r3^2+L_23^2)-((r2^2+r3^2+L_23^2)^2-4*(r2^2)*(r3^2))^0.5); 

%Apply Simlifying Criteria
%if the 'always hits' critical angle is defined from the first tooth and the detector
if theta_c_12 >= theta_c_13
    G13 = 0.5*(pi^2)*((r1^2+r3^2+L_13^2)-(((r1^2+r3^2+L_13^2)^2-4*(r1^2)*(r3^2))^0.5));
    G3_whole = G13; %geometric factor is defined by the first tooth and first detector
    fprintf('Geometric_Factor_13 (First Tooth to First Detector)= %7.5f  cm^2 sr\n',G13)
    FOV = 2*theta_m_13*180/pi;

%otherwise, if the 'can hit' critical angle is defined by the collimator
elseif theta_m_13 >= theta_m_12
    G12 = 0.5*(pi^2)*((r1^2+r2^2+L_12^2)-((r1^2+r2^2+L_12^2)^2-4*(r1^2)*(r2^2))^0.5);
    G3_whole = G12; %geometric factor is defined by the collimator
    fprintf('Geometric_Factor_12 (First Tooth to Last Tooth) = %7.5f  cm^2 sr\n',G12)
    FOV = 2*theta_m_12*180/pi;
    
%otherwise, if the 'can hit' critical angle is defined by the last tooth and the detector
elseif theta_c_12 >= theta_m_13
    G23 = 0.5*(pi^2)*((r2^2+r3^2+L_23^2)-((r2^2+r3^2+L_23^2)^2-4*(r2^2)*(r3^2))^0.5);
    G3_whole = G23; %geometric factor is defined by the last tooth and first detector
    fprintf('Geometric_Factor_23 (Last Tooth to First Detector) = %7.5f  cm^2 sr\n',G23)

%otherwise, the geometric factor is defined by all three components
else
    theta_a = atan(((L_23*r1^2+L_12*r3^2-L_13*r2^2)^0.5)/(L_12*L_23*L_13));
    G13 = 0.5*(pi^2)*((r1^2+r3^2+L_13^2)-(((r1^2+r3^2+L_13^2)^2-4*(r1^2)*(r3^2))^0.5));
    
    %from Sullivan eq 16
    Z12= Zij(theta_a,L_12,r1,r2,theta_c_12,theta_m_12,r2);
    Z13= Zij(theta_a,L_13,r1,r3,theta_c_13,theta_m_13,r2);
    Z23= Zij(theta_a,L_23,r2,r3,theta_c_23,theta_m_23,r2);
    
    G123=G13-pi^2*r2^2*sin(theta_a)^2+Z23+Z12-Z13; %Sullivan eq 15
    G3_whole = G123;
    fprintf('Geometric_Factor (Three Disc Telescope) = %7.5f  cm^2 sr\n',G123)
end
end