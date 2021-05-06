function [A1, A2, A3, A4] = F_data_ave3(CNT, NumCrossVal, A1, A2, A3, A4)
    
    A1(CNT,1) = mean(A1(CNT,2:NumCrossVal+1)); %time
    A2(CNT,1) = mean(A2(CNT,2:NumCrossVal+1)); %det
    A3(CNT,1) = mean(A3(CNT,2:NumCrossVal+1)); %mean_error
    A4(CNT,1) = std(A3(CNT,2:NumCrossVal+1));  %std_error
        
end
