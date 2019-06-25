function [normal,lambda] = Estimate_Normal_CPS_MicroBRDF_Specular(h,~,I) 
    if size(h,1) > 3
        h = h.';
    end
    
    if size(I,1) > 1
        I = I.';
    end
    
    num = length(I); %number of light
    Isqrt = sqrt(I); Isqrt_ = mean(Isqrt); Isqrt_d = Isqrt/Isqrt_;
    
    hI = h.*repmat(sqrt(Isqrt),3,1); hIhIT_ = hI*hI.'/num;
    
    Mat0 = zeros(6,6); Vec0 = zeros(6,1);
    for i = 1:num
        Mat = Isqrt_d(i)*hIhIT_ - hI(:,i)*hI(:,i).';
        vec = [ Mat(1,1), Mat(1,2) + Mat(2,1), Mat(1,3) + Mat(3,1), Mat(2,2), Mat(2,3) + Mat(3,2), Mat(3,3)];
        Mat0 = Mat0 + vec.'*vec;
        Vec0 = Vec0 - 2*(1-Isqrt_d(i))*vec.';
    end
    
     [x,y,z] = Solver_CalibratedPS_Specular(Mat0,Vec0);
     error = inf; 
     for i = 1:length(x)
         vec_temp = [x(i)^2 x(i)*y(i) x(i)*z(i) y(i)^2 y(i)*z(i) z(i)^2];
         error_temp = vec_temp*Mat0*vec_temp.'+vec_temp*Vec0;
         if error_temp < error
             normal = [x(i); y(i); z(i)];
             w = (1 + normal.'*hIhIT_*normal)/Isqrt_;
             error = error_temp;
         end
     end
          
     %calculate lambda
     lambda = 1-norm(normal)^2/w;
     
     %surface normal
     normal = normal/norm(normal); 
 
    %double check the sign
    if sum(find(h.'*normal>0)) < sum(find(-h.'*normal>0)) 
        normal = -normal;
    end
    
end

function [x y z] = Solver_CalibratedPS_Specular(A,B)
    a1 = A(1,1); a2 = A(1,2); a3 = A(1,3); a4 = A(1,4); a5 = A(1,5); a6 = A(1,6);
                 a7 = A(2,2); a8 = A(2,3); a9 = A(2,4); a10 = A(2,5); a11 = A(2,6);
                              a12 = A(3,3); a13 = A(3,4); a14 = A(3,5); a15 = A(3,6);
                                            a16 = A(4,4); a17 = A(4,5); a18 = A(4,6);
                                                          a19 = A(5,5); a20 = A(5,6);
                                                                        a21 = A(6,6); 
    b1 = B(1); b2 = B(2); b3 = B(3); b4 = B(4); b5 = B(5); b6 = B(6);                                                                    
	% precalculate polynomial equations coefficients
	c(1) = 4*a1;
	c(2) = 6*a2;
	c(3) = 2*a7+4*a4;
	c(4) = 2*a9;
	c(5) = 6*a3;
	c(6) = 4*a5+4*a8;
	c(7) = 2*a10+2*a13;
	c(8) = 2*a12+4*a6;
	c(9) = 2*a11+2*a14;
	c(10) = 2*a15;
	c(11) = 2*b1;
	c(12) = b2;
	c(13) = b3;
	c(14) = 0;
	c(15) = 2*a2;
	c(16) = 2*a7+4*a4;
	c(17) = 6*a9;
	c(18) = 4*a16;
	c(19) = 2*a5+2*a8;
	c(20) = 4*a13+4*a10;
	c(21) = 6*a17;
	c(22) = 2*a11+2*a14;
	c(23) = 4*a18+2*a19;
	c(24) = 2*a20;
	c(25) = b2;
	c(26) = 2*b4;
	c(27) = b5;
	c(28) = 0;
	c(29) = 2*a3;
	c(30) = 2*a5+2*a8;
	c(31) = 2*a10+2*a13;
	c(32) = 2*a17;
	c(33) = 2*a12+4*a6;
	c(34) = 4*a14+4*a11;
	c(35) = 4*a18+2*a19;
	c(36) = 6*a15;
	c(37) = 6*a20;
	c(38) = 4*a21;
	c(39) = b3;
	c(40) = b5;
	c(41) = 2*b6;
	c(42) = 0;

	M = zeros(89, 116);
	ci = [68, 156, 244, 599, 687, 775, 1130, 1218, 1306, 1661, 1749, 2104, 2976, 3064, 3152, 3507, 3595, 3683, 4038, 4126, 4481, 5358, 5446, 5534, 5889, 5977, 6332, 7215, 7303, 7658, 8545];
	M(ci) = c(1);

	ci = [157, 245, 333, 688, 776, 864, 1219, 1307, 1395, 1750, 1838, 2193, 3065, 3153, 3241, 3596, 3684, 3772, 4127, 4215, 4570, 5447, 5535, 5623, 5978, 6066, 6421, 7304, 7392, 7747, 8634];
	M(ci) = c(2);

	ci = [246, 334, 422, 777, 865, 953, 1308, 1396, 1484, 1839, 1927, 2282, 3154, 3242, 3330, 3685, 3773, 3861, 4216, 4304, 4659, 5536, 5624, 5712, 6067, 6155, 6510, 7393, 7481, 7836, 8723];
	M(ci) = c(3);

	ci = [335, 423, 511, 866, 954, 1042, 1397, 1485, 1573, 1928, 2016, 2371, 3243, 3331, 3419, 3774, 3862, 3950, 4305, 4393, 4748, 5625, 5713, 5801, 6156, 6244, 6599, 7482, 7570, 7925, 8812];
	M(ci) = c(4);

	ci = [691, 779, 867, 1222, 1310, 1398, 1664, 1752, 1840, 2106, 2194, 2460, 3599, 3687, 3775, 4041, 4129, 4217, 4483, 4571, 4837, 5892, 5980, 6068, 6334, 6422, 6688, 7660, 7748, 8014, 8901];
	M(ci) = c(5);

	ci = [780, 868, 956, 1311, 1399, 1487, 1753, 1841, 1929, 2195, 2283, 2549, 3688, 3776, 3864, 4130, 4218, 4306, 4572, 4660, 4926, 5981, 6069, 6157, 6423, 6511, 6777, 7749, 7837, 8103, 8990];
	M(ci) = c(6);

	ci = [869, 957, 1045, 1400, 1488, 1576, 1842, 1930, 2018, 2284, 2372, 2638, 3777, 3865, 3953, 4219, 4307, 4395, 4661, 4749, 5015, 6070, 6158, 6246, 6512, 6600, 6866, 7838, 7926, 8192, 9079];
	M(ci) = c(7);

	ci = [1314, 1402, 1490, 1756, 1844, 1932, 2109, 2197, 2285, 2462, 2550, 2727, 4133, 4221, 4309, 4486, 4574, 4662, 4839, 4927, 5104, 6337, 6425, 6513, 6690, 6778, 6955, 8016, 8104, 8281, 9168];
	M(ci) = c(8);

	ci = [1403, 1491, 1579, 1845, 1933, 2021, 2198, 2286, 2374, 2551, 2639, 2816, 4222, 4310, 4398, 4575, 4663, 4751, 4928, 5016, 5193, 6426, 6514, 6602, 6779, 6867, 7044, 8105, 8193, 8370, 9257];
	M(ci) = c(9);

	ci = [1848, 1936, 2024, 2201, 2289, 2377, 2465, 2553, 2641, 2729, 2817, 2905, 4578, 4666, 4754, 4842, 4930, 5018, 5106, 5194, 5282, 6693, 6781, 6869, 6957, 7045, 7133, 8283, 8371, 8459, 9346];
	M(ci) = c(10);

	ci = [5586, 5674, 5762, 6028, 6116, 6204, 6381, 6469, 6557, 6734, 6822, 6999, 7337, 7425, 7513, 7690, 7778, 7866, 8043, 8131, 8308, 8562, 8650, 8738, 8915, 9003, 9180, 9440, 9528, 9705, 9969];
	M(ci) = c(11);

	ci = [5675, 5763, 5851, 6117, 6205, 6293, 6470, 6558, 6646, 6823, 6911, 7088, 7426, 7514, 7602, 7779, 7867, 7955, 8132, 8220, 8397, 8651, 8739, 8827, 9004, 9092, 9269, 9529, 9617, 9794, 10058];
	M(ci) = c(12);

	ci = [6120, 6208, 6296, 6473, 6561, 6649, 6737, 6825, 6913, 7001, 7089, 7177, 7782, 7870, 7958, 8046, 8134, 8222, 8310, 8398, 8486, 8918, 9006, 9094, 9182, 9270, 9358, 9707, 9795, 9883, 10147];
	M(ci) = c(13);

	ci = [7455, 7543, 7631, 7808, 7896, 7984, 8072, 8160, 8248, 8336, 8424, 8512, 8672, 8760, 8848, 8936, 9024, 9112, 9200, 9288, 9376, 9452, 9540, 9628, 9716, 9804, 9892, 9974, 10062, 10150, 10236];
	M(ci) = c(14);

	ci = [80, 168, 256, 611, 699, 787, 1142, 1230, 1318, 1673, 1761, 2116, 2985, 3073, 3161, 3516, 3604, 3692, 4047, 4135, 4490, 5364, 5452, 5540, 5895, 5983, 6338, 7218, 7306, 7661, 8546];
	M(ci) = c(15);

	ci = [169, 257, 345, 700, 788, 876, 1231, 1319, 1407, 1762, 1850, 2205, 3074, 3162, 3250, 3605, 3693, 3781, 4136, 4224, 4579, 5453, 5541, 5629, 5984, 6072, 6427, 7307, 7395, 7750, 8635];
	M(ci) = c(16);

	ci = [258, 346, 434, 789, 877, 965, 1320, 1408, 1496, 1851, 1939, 2294, 3163, 3251, 3339, 3694, 3782, 3870, 4225, 4313, 4668, 5542, 5630, 5718, 6073, 6161, 6516, 7396, 7484, 7839, 8724];
	M(ci) = c(17);

	ci = [347, 435, 523, 878, 966, 1054, 1409, 1497, 1585, 1940, 2028, 2383, 3252, 3340, 3428, 3783, 3871, 3959, 4314, 4402, 4757, 5631, 5719, 5807, 6162, 6250, 6605, 7485, 7573, 7928, 8813];
	M(ci) = c(18);

	ci = [703, 791, 879, 1234, 1322, 1410, 1676, 1764, 1852, 2118, 2206, 2472, 3608, 3696, 3784, 4050, 4138, 4226, 4492, 4580, 4846, 5898, 5986, 6074, 6340, 6428, 6694, 7663, 7751, 8017, 8902];
	M(ci) = c(19);

	ci = [792, 880, 968, 1323, 1411, 1499, 1765, 1853, 1941, 2207, 2295, 2561, 3697, 3785, 3873, 4139, 4227, 4315, 4581, 4669, 4935, 5987, 6075, 6163, 6429, 6517, 6783, 7752, 7840, 8106, 8991];
	M(ci) = c(20);

	ci = [881, 969, 1057, 1412, 1500, 1588, 1854, 1942, 2030, 2296, 2384, 2650, 3786, 3874, 3962, 4228, 4316, 4404, 4670, 4758, 5024, 6076, 6164, 6252, 6518, 6606, 6872, 7841, 7929, 8195, 9080];
	M(ci) = c(21);

	ci = [1326, 1414, 1502, 1768, 1856, 1944, 2121, 2209, 2297, 2474, 2562, 2739, 4142, 4230, 4318, 4495, 4583, 4671, 4848, 4936, 5113, 6343, 6431, 6519, 6696, 6784, 6961, 8019, 8107, 8284, 9169];
	M(ci) = c(22);

	ci = [1415, 1503, 1591, 1857, 1945, 2033, 2210, 2298, 2386, 2563, 2651, 2828, 4231, 4319, 4407, 4584, 4672, 4760, 4937, 5025, 5202, 6432, 6520, 6608, 6785, 6873, 7050, 8108, 8196, 8373, 9258];
	M(ci) = c(23);

	ci = [1860, 1948, 2036, 2213, 2301, 2389, 2477, 2565, 2653, 2741, 2829, 2917, 4587, 4675, 4763, 4851, 4939, 5027, 5115, 5203, 5291, 6699, 6787, 6875, 6963, 7051, 7139, 8286, 8374, 8462, 9347];
	M(ci) = c(24);

	ci = [5598, 5686, 5774, 6040, 6128, 6216, 6393, 6481, 6569, 6746, 6834, 7011, 7346, 7434, 7522, 7699, 7787, 7875, 8052, 8140, 8317, 8568, 8656, 8744, 8921, 9009, 9186, 9443, 9531, 9708, 9970];
	M(ci) = c(25);

	ci = [5687, 5775, 5863, 6129, 6217, 6305, 6482, 6570, 6658, 6835, 6923, 7100, 7435, 7523, 7611, 7788, 7876, 7964, 8141, 8229, 8406, 8657, 8745, 8833, 9010, 9098, 9275, 9532, 9620, 9797, 10059];
	M(ci) = c(26);

	ci = [6132, 6220, 6308, 6485, 6573, 6661, 6749, 6837, 6925, 7013, 7101, 7189, 7791, 7879, 7967, 8055, 8143, 8231, 8319, 8407, 8495, 8924, 9012, 9100, 9188, 9276, 9364, 9710, 9798, 9886, 10148];
	M(ci) = c(27);

	ci = [7467, 7555, 7643, 7820, 7908, 7996, 8084, 8172, 8260, 8348, 8436, 8524, 8681, 8769, 8857, 8945, 9033, 9121, 9209, 9297, 9385, 9458, 9546, 9634, 9722, 9810, 9898, 9977, 10065, 10153, 10237];
	M(ci) = c(28);

	ci = [267, 711, 799, 1154, 1242, 1330, 1685, 1773, 2128, 3082, 3170, 3525, 3613, 3701, 4056, 4144, 4499, 5370, 5458, 5546, 5901, 5989, 6344, 7221, 7309, 7664, 8547];
	M(ci) = c(29);

	ci = [356, 800, 888, 1243, 1331, 1419, 1774, 1862, 2217, 3171, 3259, 3614, 3702, 3790, 4145, 4233, 4588, 5459, 5547, 5635, 5990, 6078, 6433, 7310, 7398, 7753, 8636];
	M(ci) = c(30);

	ci = [445, 889, 977, 1332, 1420, 1508, 1863, 1951, 2306, 3260, 3348, 3703, 3791, 3879, 4234, 4322, 4677, 5548, 5636, 5724, 6079, 6167, 6522, 7399, 7487, 7842, 8725];
	M(ci) = c(31);

	ci = [534, 978, 1066, 1421, 1509, 1597, 1952, 2040, 2395, 3349, 3437, 3792, 3880, 3968, 4323, 4411, 4766, 5637, 5725, 5813, 6168, 6256, 6611, 7488, 7576, 7931, 8814];
	M(ci) = c(32);

	ci = [890, 1334, 1422, 1688, 1776, 1864, 2130, 2218, 2484, 3705, 3793, 4059, 4147, 4235, 4501, 4589, 4855, 5904, 5992, 6080, 6346, 6434, 6700, 7666, 7754, 8020, 8903];
	M(ci) = c(33);

	ci = [979, 1423, 1511, 1777, 1865, 1953, 2219, 2307, 2573, 3794, 3882, 4148, 4236, 4324, 4590, 4678, 4944, 5993, 6081, 6169, 6435, 6523, 6789, 7755, 7843, 8109, 8992];
	M(ci) = c(34);

	ci = [1068, 1512, 1600, 1866, 1954, 2042, 2308, 2396, 2662, 3883, 3971, 4237, 4325, 4413, 4679, 4767, 5033, 6082, 6170, 6258, 6524, 6612, 6878, 7844, 7932, 8198, 9081];
	M(ci) = c(35);

	ci = [1513, 1868, 1956, 2133, 2221, 2309, 2486, 2574, 2751, 4239, 4327, 4504, 4592, 4680, 4857, 4945, 5122, 6349, 6437, 6525, 6702, 6790, 6967, 8022, 8110, 8287, 9170];
	M(ci) = c(36);

	ci = [1602, 1957, 2045, 2222, 2310, 2398, 2575, 2663, 2840, 4328, 4416, 4593, 4681, 4769, 4946, 5034, 5211, 6438, 6526, 6614, 6791, 6879, 7056, 8111, 8199, 8376, 9259];
	M(ci) = c(37);

	ci = [2047, 2313, 2401, 2489, 2577, 2665, 2753, 2841, 2929, 4684, 4772, 4860, 4948, 5036, 5124, 5212, 5300, 6705, 6793, 6881, 6969, 7057, 7145, 8289, 8377, 8465, 9348];
	M(ci) = c(38);

	ci = [5785, 6140, 6228, 6405, 6493, 6581, 6758, 6846, 7023, 7443, 7531, 7708, 7796, 7884, 8061, 8149, 8326, 8574, 8662, 8750, 8927, 9015, 9192, 9446, 9534, 9711, 9971];
	M(ci) = c(39);

	ci = [5874, 6229, 6317, 6494, 6582, 6670, 6847, 6935, 7112, 7532, 7620, 7797, 7885, 7973, 8150, 8238, 8415, 8663, 8751, 8839, 9016, 9104, 9281, 9535, 9623, 9800, 10060];
	M(ci) = c(40);

	ci = [6319, 6585, 6673, 6761, 6849, 6937, 7025, 7113, 7201, 7888, 7976, 8064, 8152, 8240, 8328, 8416, 8504, 8930, 9018, 9106, 9194, 9282, 9370, 9713, 9801, 9889, 10149];
	M(ci) = c(41);

	ci = [7654, 7920, 8008, 8096, 8184, 8272, 8360, 8448, 8536, 8778, 8866, 8954, 9042, 9130, 9218, 9306, 9394, 9464, 9552, 9640, 9728, 9816, 9904, 9980, 10068, 10156, 10238];
	M(ci) = c(42);


	%Mr = rref(M);  % replace me with a MEX
    index1 = [1:59, 61:78, 82:90, 97:99];
    index2 = [116 115 114 113 112 111 110 109 108 107 106 105 104 103 102 101 100 96 95 94 93 92 91 81 80 79 60];
    Mr = M(:,index1)\M(:,index2);

	A = zeros(27);
	%amcols = [116 115 114 113 112 111 110 109 108 107 106 105 104 103 102 101 100 96 95 94 93 92 91 81 80 79 60];
	A(1, 4) = 1;
	A(2, 7) = 1;
	A(3, 9) = 1;
	A(4, 10) = 1;
	A(5, 13) = 1;
	A(6, 15) = 1;
	A(7, 16) = 1;
	A(8, :) = -Mr(89, :);
	A(9, :) = -Mr(88, :);
	A(10, :) = -Mr(87, :);
	A(11, 20) = 1;
	A(12, 22) = 1;
	A(13, 23) = 1;
	A(14, :) = -Mr(85, :);
	A(15, :) = -Mr(84, :);
	A(16, :) = -Mr(83, :);
	A(17, :) = -Mr(81, :);
	A(18, 26) = 1;
	A(19, :) = -Mr(76, :);
	A(20, :) = -Mr(75, :);
	A(21, :) = -Mr(73, :);
	A(22, :) = -Mr(72, :);
	A(23, :) = -Mr(71, :);
	A(24, :) = -Mr(58, :);
	A(25, :) = -Mr(56, :);
	A(26, :) = -Mr(55, :);
	A(27, :) = -Mr(31, :);

	[V D] = eig(A);
	sol =  V([4, 3, 2],:)./(ones(3, 1)*V(1,:));

	if (find(isnan(sol(:))) > 0)
		
		x = [];
		y = [];
		z = [];
	else
		
		I = find(not(imag( sol(1,:) )));
		x = sol(1,I);
		y = sol(2,I);
		z = sol(3,I);
	end
end
