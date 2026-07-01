%% code to find required reflux ratio using the bottom and top concentrations.
clear; 
clc; 
close all;

xD = 0.9;   
xB = 0.320;   
N = 4;        
alpha = 2.30; %relative volatility       

obj_func = @(R) mccabe_thiele_sim(R, xD, xB, N, alpha);

try
    R_req = fzero(obj_func, [0.1, 100]);
    
    fprintf('xD: %.4f\n', xD);
    fprintf('xB: %.4f\n', xB);
    fprintf('Required R: %.2f\n', R_req);
catch
    fprintf('Error: Cannot reach xD=%.4f from xB=%.4f with %d stages.\n', xD, xB, N);
end

function err = mccabe_thiele_sim(R, xD, target_xB, N, alpha)
    x = xD; 
    y = xD; 
    
    for i = 1:N
        x = y / (alpha - y * (alpha - 1));
        y = (R / (R + 1)) * x + (xD / (R + 1));
    end
    
    err = x - target_xB;
end