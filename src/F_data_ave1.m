function [A1, A2] = F_data_ave1(CNT, num_ave, A1, A2)
    
    A1(CNT,1) = mean(A1(CNT,2:num_ave+1));
    A2(CNT,1) = mean(A2(CNT,2:num_ave+1));
end
