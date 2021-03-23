clear all
close all
clc
% delete(findall(0));
if isempty(getenv('COMPUTERNAME'))% on server
    mask_threshold = -4;
    threshold_brain = 5;
    threshold_scalp = 20;
    sigma_brain = 5;
    sigma_scalp = 20;
else % on my own computer
    mask_threshold = -2;
    threshold_brain = 15;
    threshold_scalp = 80;
    sigma_brain = 15;
    sigma_scalp = 80;
end
%% tool folder manage
add_tool_path(['toolbox',filesep, 'snirf_homer3'])
add_tool_path(['toolbox',filesep, 'Homer3-master',filesep, 'Homer3-master'])

%% load snirf data and process in Homer3
acquired_path = ['data',filesep];
file = '2020-12-25_006_probe_correct.snirf';
Proc_data = process_in_homer([acquired_path, file]);

%% Calculate L matrix
rhoSD_ssThresh = -1000;
[Adot ,Adot_scalp, E] = Get_A_dot(acquired_path,rhoSD_ssThresh);

%% Calculate T matrix
T = Make_T_matrix(Proc_data);
% here for scalp time series, now I am using sin and cos. Can i use short separation as time basis?

%% Calculate G matrix
At_file = 'atlasViewer.mat';

M = Make_mask([acquired_path,At_file],mask_threshold, Adot, Adot_scalp);
%%
alpha = 0.001;
beta = 0.001;
A = Make_A_matrix(Adot, Adot_scalp, E,  M, alpha, beta);
G = Make_G_matrix([acquired_path,At_file], M, threshold_brain, threshold_scalp, sigma_brain, sigma_scalp);

%% H matrix

OD = Proc_data.dod.dataTimeSeries;
Y = reshape(OD',[],1); 
device = 'gpu'; % device could be 'gpu' or 'cpu'
n_batch_brain = 1;
n_batch_scalp = 1;
fprintf('start gpu computation\n')
tic
[HTH,HTY] = make_H_matrix(A, G, T, Y, device, n_batch_brain, n_batch_scalp);% parallel computing toolbox is needed to run on gpu
toc
% fprintf('start cpu computation\n')
% tic
% [HTH,HTY] = make_H_matrix(A, G, T, Y, 'cpu',n_batch);
% toc
save('HTH.mat', 'HTH')
save('HTY.mat', 'HTY')
%% do inverse
tic
beta = (HTY)'/(HTH + alpha * eigs(HTH,1) * eye(size(HTH,1))); 
toc
%% approximate the conc Conc = GT*beta
tic
Conc = approximate_Conc(G,T,beta, device, n_batch_brain);
toc
save('Conc.mat','Conc')
%% visualization

visualize(Conc,acquired_path,M,Proc_data);