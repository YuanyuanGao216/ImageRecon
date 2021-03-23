function Proc_data = process_in_homer(file,rhoSD_ssThresh)

acquired = SnirfClass(file);
data = acquired.data;
probe = acquired.probe;

mlActMan =  cell(1,1);
mlActMan{1,1} = ones(length(data.measurementList),1);

tIncMan =  cell(1,1);
tIncMan{1,1} = ones(size(acquired.data.time));

stim = acquired.stim;
Aaux = [];
tIncAuto = [];
rcMap = [];

%% homer3 GLM processing

mlActAuto = hmrR_PruneChannels(data,probe,mlActMan,tIncMan,[0  10000000],4,[0  45]);
dod = hmrR_Intensity2OD(data);
dod = hmrR_BandpassFilt(dod,0,0.5);
dc = hmrR_OD2Conc(dod,probe,[6  6  1]);
[dcAvg,dcAvgStd,nTrials,dcNew,dcResid,dcSum2,beta,R,dA,tbasis,tHRF] = hmrR_GLM_time_basis(dc,stim,probe,mlActAuto,Aaux,tIncAuto,rcMap,[-2  20],1,2,[0.1 3.0 10.0 1.8 3.0 10.0],rhoSD_ssThresh,1,3);



Proc_data.dod = dod;
Proc_data.dA = dA;
Proc_data.time = dc.time;
Proc_data.tHRF = tHRF;
Proc_data.stim = stim;
Proc_data.tbasis = tbasis;
Proc_data.mlActAuto = mlActAuto;
Proc_data.dc = dc;