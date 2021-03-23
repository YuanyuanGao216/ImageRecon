function T = Make_T_matrix(Proc_data, SS_flag, shortSepChLst, condition)
nB = size(Proc_data.tbasis,2);
if condition > nB
    fprintf('condition is exceding\n')
end
T_HbO_brain = squeeze(Proc_data.dA(:,condition:nB:end,1)); 
T_HbR_brain = squeeze(Proc_data.dA(:,condition:nB:end,2)); 
t_HbO_brain = squeeze(Proc_data.tbasis(:,:,1));
t_HbR_brain = squeeze(Proc_data.tbasis(:,:,2)); 

if SS_flag == 0
    time = Proc_data.time;
    fs = 1/(time(2)-time(1));
    frq1=(1:2)/2/fs;
    T_scalp = zeros(length(time),2*length(frq1));

    for idx=1:length(frq1)
        T_scalp(:,1+(idx-1)*2)=sin(time*2*pi*frq1(idx))';
        T_scalp(:,2+(idx-1)*2)=cos(time*2*pi*frq1(idx))';
    end
    T.T_HbO_scalp = T_scalp;
    T.T_HbR_scalp = T_scalp;
elseif isempty(shortSepChLst)
    dc = Proc_data.dc;
    dc_HbO = dc.dataTimeSeries(:,1:3:end);
    dc_HbR = dc.dataTimeSeries(:,2:3:end);

    dc_HbO_ss = mean(dc_HbO,2);
    dc_HbR_ss = mean(dc_HbR,2);
    
    T.T_HbO_scalp = dc_HbO_ss;
    T.T_HbR_scalp = dc_HbR_ss;
else
    dc = Proc_data.dc;
    dc_HbO = dc.dataTimeSeries(:,1:3:end);
    dc_HbR = dc.dataTimeSeries(:,2:3:end);

    dc_HbO_ss = dc_HbO(:,shortSepChLst);
    dc_HbR_ss = dc_HbR(:,shortSepChLst);
    
    T.T_HbO_scalp = dc_HbO_ss;
    T.T_HbR_scalp = dc_HbR_ss;
end

T.T_HbO_brain = T_HbO_brain;
T.T_HbR_brain = T_HbR_brain;

T.t_HbO_brain = t_HbO_brain;
T.t_HbR_brain = t_HbR_brain;
