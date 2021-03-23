function [Adot ,Adot_scalp, E, channels, shortSepChLst] = Get_A_dot(acquired_path,rhoSD_ssThresh,mlActAuto)

% load At file
A_file = 'fw\Adot.mat';
A_scalp_file = 'fw\Adot_scalp.mat';

Adot        = load([acquired_path,A_file]);
Adot_scalp  = load([acquired_path,A_scalp_file]);

Adot = Adot.Adot;
Adot_scalp = Adot_scalp.Adot_scalp;

At_file = 'atlasViewer.mat';
atlasViewer = load([acquired_path,At_file]);
warning off

imgrecon = atlasViewer.imgrecon;

subjData    = imgrecon.subjData;
SD          = subjData.SD;

% use only active channels
% ml = SD.MeasList;
% if isfield(SD, 'MeasListAct') == 1
%     activeChLst = find(ml(:,4)==1 & SD.MeasListAct==1);
% end
% use only active channels
ml = SD.MeasList;
activeChLst = find(ml(:,4)==1 & mlActAuto{1,1}==1);

% get long separation channels only for reconstruction
lst = find(ml(:,4)==1 & SD.MeasListAct==1);
rhoSD = zeros(length(lst),1);
posM = zeros(length(lst),3);
for iML = 1:length(lst)
    rhoSD(iML) = sum((SD.SrcPos(ml(lst(iML),1),:) - SD.DetPos(ml(lst(iML),2),:)).^2).^0.5;
    posM(iML,:) = (SD.SrcPos(ml(lst(iML),1),:) + SD.DetPos(ml(lst(iML),2),:)) / 2;
end
% rhoSD_ssThresh = 0.8;
longSepChLst = lst(find(rhoSD>=rhoSD_ssThresh));
shortSepChLst = lst(find(rhoSD<rhoSD_ssThresh));
lstLS_all = [longSepChLst; longSepChLst+size(ml,1)/2]; % both wavelengths

if isempty(lstLS_all)
    menu(sprintf('All channels meet short separation threshold.\nYou need some long separation channels for image recon.\nPlease lower the threshold and retry.'), 'Okay');
    return;
end

[channels,~]=intersect(activeChLst,longSepChLst);

Adot = Adot(channels,:,:);
Adot_scalp = Adot_scalp(channels,:,:);

% Adot = Adot(activeChLst,:,:);
% Adot = Adot(longSepChLst,:,:);
% 
% Adot_scalp = Adot_scalp(activeChLst,:,:);
% Adot_scalp = Adot_scalp(longSepChLst,:,:);



% put A matrix together and combine with extinction coefficients
E = GetExtinctions([SD.Lambda(1) SD.Lambda(2)]);
E = E/10; %convert from /cm to /mm  E raws: wavelength, columns 1:HbO, 2:HbR
% Amatrix = [squeeze(Anew(:,:,1))*E(1,1) squeeze(Anew(:,:,1))*E(1,2);
%            squeeze(Anew(:,:,2))*E(2,1) squeeze(Anew(:,:,2))*E(2,2)]; % 100, 59132
