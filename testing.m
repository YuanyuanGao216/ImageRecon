% N = [2452356, ...
%     1566, ...
%     8441688, ...
%     20004, ...
%     8292, ...
%     8292, ...
%     4000800,...
%     82936584,...
%     82936584,...
%     165873168,...
%     414600,...
%     414600,...
%     20004,...
%     8292,...
%     8292,...
%     4000800,...
%     82936584,...
%     82936584,...
%     165873168,...
%     414600,...
%     414600];
% 
% GB = N*8/(1024)^3;
% GB = GB';
% sumGB = sum(GB)

% N = [2452356, ...
%     1566, ...
%     8441688, ...
%     20004, ...
%     8292, ...
%     8292, ...
%     4000800,...
%     165873168,...
%     414600,...
%     414600,...
%     20004,...
%     8292,...
%     8292,...
%     4000800,...
%     165873168,...
%     414600,...
%     414600];
% 
% GB = N*8/(1024)^3;
% GB = GB';
% sumGB = sum(GB)

% N = [349922400, ...
%     8441688, ...
%     82936584, ...
%     82936584, ...
%     8292, ...
%     8292,...
%     4000800,...
%     165873168,...
%     414600,...
%     ];
% 
% GB = N*8/(1024)^3;
% GB = GB';
% sumGB = sum(GB)
% 
% N = [57600, ...
%     240, ...
%     78800, ...
%     1576, ...
%     4146, ...
%     4146,...
%     315200,...
%     6534096,...
%     6534096,...
%     13068192,...
%     414600,...
%     414600,...
%     1576, ...
%     4146, ...
%     4146,...
%     315200,...
%     6534096,...
%     6534096,...
%     13068192,...
%     414600,...
%     414600,...
%     ];
% 
% GB = N*8/(1024)^3;
% GB = GB';
% sumGB = sum(GB)

% H_GPU = load('H_GPU.mat','H');
% H_CPU = load('H_CPU.mat','H');

% for i_H = 1:70
% %     fprintf('i_H is %d \n',i_H)
%     
%     H_GPU_1 = H_GPU.H(:,i_H);
%     H_CPU_1 = H_CPU.H(:,i_H);
% 
%     [row,col] = find(H_GPU_1 - H_CPU_1 > 1e-3,5);
%     if ~isempty(row)
%         for j = 1:length(row)
%             fprintf('col is %d, row is %d',i_H, row(j))
%             fprintf('H_GPU is %f, H_CPU is %f\n',H_GPU_1(row(j)), H_CPU_1(row(j)))
%         end
%     end
% end
% H1 = [H_GPU_1,H_CPU_1];
acquired = SnirfClass(['data_ss',filesep,'2021-01-13_001.snirf']);
acquired.probe.sourcePos2D = acquired.probe.sourcePos2D*10;
acquired.probe.detectorPos2D = acquired.probe.detectorPos2D*10;
acquired.probe.sourcePos3D = acquired.probe.sourcePos3D*10;
acquired.probe.detectorPos3D = acquired.probe.detectorPos3D*10;
acquired.Save(['data_ss',filesep,'2021-01-13_001_probe_correct.snirf'])