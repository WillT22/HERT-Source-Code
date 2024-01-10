gapR = 0.108; % gap from back of detect to front of next detector
gapB = 0.10;  % gap from back of detect to front of next detector

depth = 0.15
d1_hd = 0.05 * depth; % half depth
d1_z = d1_hd;
d2_z = d1_z + depth + gapR;
d3_z = d2_z + depth + gapB;
d4_z = d3_z + depth + gapR;
d5_z = d4_z + depth + gapB;
d6_z = d5_z + depth + gapR;
d7_z = d6_z + depth + gapB;
d8_z = d7_z + depth + gapR;
d9_z = d8_z + depth + gapB;

disp(d9_z)