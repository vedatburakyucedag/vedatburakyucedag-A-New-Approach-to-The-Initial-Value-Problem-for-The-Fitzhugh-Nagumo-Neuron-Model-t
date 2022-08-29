clc; clear;

v1 = load('original.mat');
v2 = load('adams_basforth.mat');
v3 = load('new_version.mat');
v4 = load('predictor_corrector.mat');

n = 10000;

figure(2); clf(2);

aa = 0; 
bb = 0; 
cc = 0; 
dd = 0; 

er1_v2 = 0; 
er1_v3 = 0; 
er1_v4 = 0; 

aa2 = 0; 
bb2 = 0; 
cc2 = 0; 
dd2 = 0;

er2_v2 = 0;
er2_v3 = 0;
er2_v4 = 0;

for i = 1:n
    
	aa = aa + v1.v(i);
    bb = bb + v2.v(i);
    cc = cc + v3.v(i);
    dd = dd + v4.v(i);
    
    er1_v2 = er1_v2 + (v2.v(i) - v1.v(i))^2;
    er1_v3 = er1_v3 + (v3.v(i) - v1.v(i))^2;
    er1_v4 = er1_v4 + (v4.v(i) - v1.v(i))^2;
    
    aa2 = aa2 + (v1.v(i))^2;
    bb2 = bb2 + (v2.v(i))^2;
    cc2 = cc2 + (v3.v(i))^2;
    dd2 = dd2 + (v4.v(i))^2;
    
    
   er2_v2 = er2_v2 + abs(v2.v(i) - v1.v(i));
   er2_v3 = er2_v3 + abs(v3.v(i) - v1.v(i));
   er2_v4 = er2_v4 + abs(v4.v(i) - v1.v(i));
        
end

root_mean_square_error_er1_v2 = sqrt(er1_v2/n)
root_mean_square_error_er1_v3 = sqrt(er1_v3/n)
root_mean_square_error_er1_v4 = sqrt(er1_v4/n)

norm_root_mean_square_error_er1_v2 = root_mean_square_error_er1_v2/(max(v1.v) - min(v1.v))
norm_root_mean_square_error_er1_v3 = root_mean_square_error_er1_v3/(max(v1.v) - min(v1.v))
norm_root_mean_square_error_er1_v4 = root_mean_square_error_er1_v4/(max(v1.v) - min(v1.v))

Relative_Spike_Energy_Error_er1_v2 = (aa2 - bb2)/(aa2) % perfect
Relative_Spike_Energy_Error_er1_v3 = (aa2 - cc2)/(aa2) % perfect
Relative_Spike_Energy_Error_er1_v4 = (aa2 - dd2)/(aa2) % perfect

mean_absolute_error_v2 = er2_v2/n
mean_absolute_error_v3 = er2_v3/n
mean_absolute_error_v4 = er2_v4/n




    
