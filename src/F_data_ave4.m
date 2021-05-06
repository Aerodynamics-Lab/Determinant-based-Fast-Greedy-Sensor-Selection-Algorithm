function [A1, A2, A3, A4] = F_data_ave4(CNT, num_ave, NumCrossVal, CV, A1, A2, A3, A4)
    
    A1(CNT, NumCrossVal+1) = mean(A1(CNT,CV+1:num_ave+CV));
    A2(CNT, NumCrossVal+1) = mean(A2(CNT,CV+1:num_ave+CV));
    A3(CNT, NumCrossVal+1) = mean(A3(CNT,CV+1:num_ave+CV));
    A4(CNT, NumCrossVal+1) = mean(A4(CNT,CV+1:num_ave+CV));    
end