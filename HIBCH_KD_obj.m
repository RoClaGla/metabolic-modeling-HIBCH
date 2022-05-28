% Uncomment first run
%initCobraToolbox();
%changeCobraSolver('mosek','all');
%loadMod = readCbModel('iHepatocytes2322.xml');
model = loadMod;

% All units are mmol/gDW/h, millimoles per gram dry weight per hour

% Model to assess before and after KD/KU without medium
rawMod = model;

% Knock down HIBCH associated reactions: (both have lb = 0)
down = 0.5;
KDmod = changeRxnBounds(model, 'HMR_4741', 1000*(1-down), 'u');
KDmod = changeRxnBounds(KDmod, 'HMR_3755',1000*(1-down), 'u');

% Knock up HIBCH associated reactions: (both have lb = 0)
%up = 0.5;
%KUmod = changeRxnBounds(mod, 'HMR_4741', 1000*(1+up), 'u');
%KUmod = changeRxnBounds(KUmod, 'HMR_3755', 1000*(1+up), 'u');

% Knock down HIBCH associated reactions: (both have lb = 0)
rawModKD = changeRxnBounds(rawMod, 'HMR_4741', 1000*(1-down), 'u');
rawModKD = changeRxnBounds(rawModKD, 'HMR_3755', 1000*(1-down), 'u');

% 
%rawModKU = changeRxnBounds(rawMod, 'HMR_4741', 1000*(1+up), 'u');
%rawModKU = changeRxnBounds(rawModKU, 'HMR_3755', 1000*(1+up), 'u');

% Objectives (picked by HIBCH_KD_full)
objRxns = importdata('rxnsData\ObjRxns.txt');
objInds = objRxns.data;
objRxnNames = objRxns.textdata;

% active uptake reactions -- file contains names and a matrix of indices
% first column -- indices in terms of iHepatocytes model
% second       -- indices in terms of the uptake reactions index in
%                 the exchange reaction indicies
activeUptakeRxns = importdata('rxnsData\ActiveUptakeRxns.txt');
activeUptakeRxnInds = activeUptakeRxns.data(:,1);
activeUptakeRxnNames = activeUptakeRxns.textdata;

% All exchange reactions (among which are medium components)
exchRxns  = importdata('rxnsData\ExchRxns.txt');
exchRxnInds  = exchRxns.data;
exchRxnNames = exchRxns.textdata;

% Active exchange reactions in terms of indicies in exchRxnInds
actExchRxnInds = activeUptakeRxns.data(:,2);
               

% Set all exchange reaction uptakes to 0 (lower bound = 0) for medium
for i=1:length(exchRxnInds)
    model = changeRxnBounds(model, exchRxnNames(i),0,'l');
    KDmod = changeRxnBounds(KDmod, exchRxnNames(i),0,'l');
    %KUmod = changeRxnBounds(KUmod, exchRxnNames(i),0,'l');
end
% Set medium components to -1000 (lower bound = -1000)
for i=1:length(activeUptakeRxnInds)
    model = changeRxnBounds(model, model.rxns(activeUptakeRxnInds(i)),-1000,'l');
    KDmod = changeRxnBounds(KDmod, model.rxns(activeUptakeRxnInds(i)),-1000,'l');
    %KUmod = changeRxnBounds(KUmod,  mod.rxns(activeUptakeRxnInds(i)),-1000,'l');
end

% Make two more models in which we force flux
%   - one with lb = 500, ub = 1000 for SMFCA blood pool -> FAs (palmitate/oleate)  (ind 4988)
%   - one with lb = 500, ub = 1000 for NEFA blood pool --> FAs (EPA/DHA) (ind 4987)

SMFCAind = 4988;
SMFCAmod = changeRxnBounds(model, model.rxns(SMFCAind), 500,'l');
SMFCAKDmod = changeRxnBounds(KDmod, model.rxns(SMFCAind), 500,'l');

SMFCAmod_raw = changeRxnBounds(rawMod, model.rxns(SMFCAind), 500,'l');
SMFCAKDmod_raw = changeRxnBounds(rawModKD, model.rxns(SMFCAind), 500,'l');

NEFAind = 4987;
NEFAmod = changeRxnBounds(model, model.rxns(NEFAind), 500,'l');
NEFAKDmod = changeRxnBounds(KDmod, model.rxns(NEFAind), 500,'l');

NEFAmod_raw = changeRxnBounds(rawMod, model.rxns(NEFAind), 500,'l');
NEFAKDmod_raw = changeRxnBounds(rawModKD, model.rxns(NEFAind), 500,'l');


% Objectives across the two different setups
  % with medium before, after KD, and after KU
obj_before = length(objInds)*[];
obj_afterKD = length(objInds)*[];
obj_SMFCA_before = length(objInds)*[];
obj_SMFCA_afterKD = length(objInds)*[];
obj_NEFA_before = length(objInds)*[];
obj_NEFA_afterKD = length(objInds)*[];
%obj_afterKU = length(objInds)*[];

  % without medium before, after KD, and after KU
obj_raw_before = length(objInds)*[];
obj_raw_afterKD = length(objInds)*[];
obj_raw_SMFCA_before = length(objInds)*[];
obj_raw_SMFCA_afterKD = length(objInds)*[];
obj_raw_NEFA_before = length(objInds)*[];
obj_raw_NEFA_afterKD = length(objInds)*[];
%obj_raw_afterKU = length(objInds)*[];

% Exchange reactions and other reactions
  % with medium
exchRxn_before = length(objInds)*[];
exchRxn_afterKD = length(objInds)*[];
exchRxn_SMFCA_before = length(objInds)*[];
exchRxn_SMFCA_afterKD = length(objInds)*[];
exchRxn_NEFA_before = length(objInds)*[];
exchRxn_NEFA_afterKD = length(objInds)*[];
%exchRxn_afterKU = length(objInds)*[];

  % without medium
exchRxn_raw_before = length(objInds)*[];
exchRxn_raw_afterKD = length(objInds)*[];
exchRxn_raw_SMFCA_before = length(objInds)*[];
exchRxn_raw_SMFCA_afterKD = length(objInds)*[];
exchRxn_raw_NEFA_before = length(objInds)*[];
exchRxn_raw_NEFA_afterKD = length(objInds)*[];
%exchRxn_raw_afterKU = length(objInds)*[];

% Loop through reactions that change under HIBCH KD in the hepatocyte model
% print as you go.

% across the exchange reactions!
fID = fopen('results\HIBCHKD-obj.txt','w');
fIDr = fopen('results\HIBCHKD-obj-raw.txt','w');
% Exchange reactions
fID_eb = fopen('results\HIBCH-exchRxns-before.txt','w');
fKD_ea = fopen('results\OHIBCH-exchRxns-afterKD.txt','w');
%fKU_ea = fopen('results\medium-HIBCH-exchRxns-afterKU.txt','w');
fID_ebr = fopen('results\HIBCH-exchRxns-before-raw.txt','w');
fKD_ear = fopen('results\HIBCH-exchRxns-afterKD-raw.txt','w');
%fKU_ear = fopen('results\medium-HIBCH-exchRxns-afterKU-raw.txt','w');

fID_NEFA = fopen('results\HIBCHKD-obj-NEFA.txt','w');
fIDr_NEFA = fopen('results\HIBCHKD-obj-NEFA-raw.txt','w');
% Exchange reactions
fID_eb_NEFA = fopen('results\HIBCH-exchRxns-before-NEFA.txt','w');
fKD_ea_NEFA = fopen('results\HIBCH-exchRxns-afterKD-NEFA.txt','w');
%fKU_ea = fopen('results\medium-HIBCH-exchRxns-afterKU.txt','w');
fID_ebr_NEFA = fopen('results\HIBCH-exchRxns-before-NEFA-raw.txt','w');
fKD_ear_NEFA = fopen('results\HIBCH-exchRxns-afterKD-NEFA-raw.txt','w');

fID_SMFCA = fopen('results\HIBCHKD-obj-SMFCA.txt','w');
fIDr_SMFCA = fopen('results\HIBCHKD-obj-SMFCA-raw.txt','w');
% Exchange reactions
fID_eb_SMFCA = fopen('results\HIBCH-exchRxns-before-SMFCA.txt','w');
fKD_ea_SMFCA = fopen('results\OHIBCH-exchRxns-afterKD-SMFCA.txt','w');
%fKU_ea = fopen('results\medium-HIBCH-exchRxns-afterKU.txt','w');
fID_ebr_SMFCA = fopen('results\HIBCH-exchRxns-before-SMFCA-raw.txt','w');
fKD_ear_SMFCA = fopen('results\HIBCH-exchRxns-afterKD-SMFCA-raw.txt','w');
%fKU_ear = fopen('results\medium-HIBCH-exchRxns-afterKU-raw.txt','w');

fprintf(fID, 'ObjName, Before, After KD\n');
fprintf(fID_eb, 'ObjName, Before\n');
fprintf(fKD_ea, 'ObjName, After KD\n');
%fprintf(fKU_ea, 'ObjName, After KU\n');
fprintf(fIDr, 'ObjName, Before, After KD\n');
fprintf(fID_ebr, 'ObjName, Before\n');
fprintf(fKD_ear, 'ObjName, After KD\n');
%fprintf(fKU_ear, 'ObjName, After KU\n');

fprintf(fID_NEFA, 'ObjName, Before, After KD\n');
fprintf(fID_eb_NEFA, 'ObjName, Before\n');
fprintf(fKD_ea_NEFA, 'ObjName, After KD\n');
%fprintf(fKU_ea, 'ObjName, After KU\n');
fprintf(fIDr_NEFA, 'ObjName, Before, After KD\n');
fprintf(fID_ebr_NEFA, 'ObjName, Before\n');
fprintf(fKD_ear_NEFA, 'ObjName, After KD\n');
%fprintf(fKU_ear, 'ObjName, After KU\n');

fprintf(fID_SMFCA, 'ObjName, Before, After KD\n');
fprintf(fID_eb_SMFCA, 'ObjName, Before\n');
fprintf(fKD_ea_SMFCA, 'ObjName, After KD\n');
%fprintf(fKU_ea, 'ObjName, After KU\n');
fprintf(fIDr_SMFCA, 'ObjName, Before, After KD\n');
fprintf(fID_ebr_SMFCA, 'ObjName, Before\n');
fprintf(fKD_ear_SMFCA, 'ObjName, After KD\n');
%fprintf(fKU_ear, 'ObjName, After KU\n');

% Set all objectives to 0
model.c(1:length(model.c)) = 0;
KDmod.c(1:length(model.c)) = 0;
SMFCAmod.c(1:length(model.c)) = 0;
SMFCAKDmod.c(1:length(model.c)) = 0;
NEFAmod.c(1:length(model.c)) = 0;
NEFAKDmod.c(1:length(model.c)) = 0;
%KUmod.c(1:length(mod.c)) = 0;
rawMod.c(1:length(model.c)) = 0;
rawModKD.c(1:length(model.c)) = 0;
SMFCAmod_raw.c(1:length(model.c)) = 0;
SMFCAKDmod_raw.c(1:length(model.c)) = 0;
NEFAmod_raw.c(1:length(model.c)) = 0;
NEFAKDmod_raw.c(1:length(model.c)) = 0;
%rawModKU.c(1:length(mod.c)) = 0;

for j=1:length(objInds)
    % reset previous objective function (if applicable; i.e., j>1)
    if(j>1)
        model.c(objInds(j-1)) = 0;
        KDmod.c(objInds(j-1)) = 0;
        SMFCAmod.c(objInds(j-1)) = 0;
        SMFCAKDmod.c(objInds(j-1)) = 0;
        NEFAmod.c(objInds(j-1)) = 0;
        NEFAKDmod.c(objInds(j-1)) = 0;
        %KUmod.c(objInds(j-1)) = 0;
        rawMod.c(objInds(j-1)) = 0;
        rawModKD.c(objInds(j-1)) = 0;
        SMFCAmod_raw.c(objInds(j-1)) = 0;
        SMFCAKDmod_raw.c(objInds(j-1)) = 0;
        NEFAmod_raw.c(objInds(j-1)) = 0;
        NEFAKDmod_raw.c(objInds(j-1)) = 0;
        %rawModKU.c(objInds(j-1)) = 0;        
    end
    % Set current objective
    model.c(objInds(j)) = 1;
    KDmod.c(objInds(j)) = 1;
    SMFCAmod.c(objInds(j)) = 1;
    SMFCAKDmod.c(objInds(j)) = 1;
    NEFAmod.c(objInds(j)) = 1;
    NEFAKDmod.c(objInds(j)) = 1;
    %KUmod.c(objInds(j)) = 1;
    
    rawMod.c(objInds(j)) = 1;
    rawModKD.c(objInds(j)) = 1;
    %rawModKU.c(objInds(j)) = 1;
    
    % run simulations, one for original model, one for KD model
    FBAsol   = optimizeCbModel(model,'max');
    FBAsolKD = optimizeCbModel(KDmod,'max');
    FBAsolSMFCA   = optimizeCbModel(SMFCAmod,'max');
    FBAsolSMFCAKD = optimizeCbModel(SMFCAKDmod,'max');
    FBAsolNEFA   = optimizeCbModel(NEFAmod,'max');
    FBAsolNEFAKD = optimizeCbModel(NEFAKDmod,'max');
    %FBAsolKU = optimizeCbModel(KUmod,'max');
    
    FBA_raw_sol   = optimizeCbModel(rawMod,'max');
    FBA_raw_solKD = optimizeCbModel(rawModKD,'max');
    FBA_raw_solSMFCA   = optimizeCbModel(SMFCAmod_raw,'max');
    FBA_raw_solSMFCAKD = optimizeCbModel(SMFCAKDmod_raw,'max');
    FBA_raw_solNEFA   = optimizeCbModel(NEFAmod_raw,'max');
    FBA_raw_solNEFAKD = optimizeCbModel(NEFAKDmod_raw,'max');
    %FBA_raw_solKU = optimizeCbModel(rawModKU,'max');
    
    % fetch data for objectives
    obj_before(j) = FBAsol.v(objInds(j));
    obj_afterKD(j)  = FBAsolKD.v(objInds(j));
    obj_SMFCA_before(j) = FBAsolSMFCA.v(objInds(j));
    obj_SMFCA_afterKD(j)  = FBAsolSMFCAKD.v(objInds(j));
    obj_NEFA_before(j) = FBAsolNEFA.v(objInds(j));
    obj_NEFA_afterKD(j)  = FBAsolNEFAKD.v(objInds(j));
    %obj_afterKU(j)  = FBAsolKU.v(objInds(j));
    
    obj_raw_before(j) = FBA_raw_sol.v(objInds(j));
    obj_raw_afterKD(j)  = FBA_raw_solKD.v(objInds(j));
    obj_raw_SMFCA_before(j) = FBA_raw_solSMFCA.v(objInds(j));
    obj_raw_SMFCA_afterKD(j)  = FBA_raw_solSMFCAKD.v(objInds(j));
    obj_raw_NEFA_before(j) = FBA_raw_solNEFA.v(objInds(j));
    obj_raw_NEFA_afterKD(j)  = FBA_raw_solNEFAKD.v(objInds(j));
    %obj_raw_afterKU(j)  = FBA_raw_solKU.v(objInds(j));
    
    % store flux results from exchange reactions as the jth row of exchRxn_    
    exchRxn_before(j,:) = FBAsol.v(exchRxnInds);
    exchRxn_afterKD(j,:)  = FBAsolKD.v(exchRxnInds);
    exchRxn_SMFCA_before(j,:) = FBAsolSMFCA.v(exchRxnInds);
    exchRxn_SMFCA_afterKD(j,:)  = FBAsolSMFCAKD.v(exchRxnInds);
    exchRxn_NEFA_before(j,:) = FBAsolNEFA.v(exchRxnInds);
    exchRxn_NEFA_afterKD(j,:)  = FBAsolNEFAKD.v(exchRxnInds);
    %exchRxn_afterKU(j,:)  = FBAsolKU.v(exchRxnInds);
    
    % store flux results from exchange reactions as the jth row of exchRxn_    
    exchRxn_raw_before(j,:) = FBA_raw_sol.v(exchRxnInds);
    exchRxn_raw_afterKD(j,:)  = FBA_raw_solKD.v(exchRxnInds);
    exchRxn_raw_SMFCA_before(j,:) = FBA_raw_solSMFCA.v(exchRxnInds);
    exchRxn_raw_SMFCA_afterKD(j,:)  = FBA_raw_solSMFCAKD.v(exchRxnInds);
    exchRxn_raw_NEFA_before(j,:) = FBA_raw_solNEFA.v(exchRxnInds);
    exchRxn_raw_NEFA_afterKD(j,:)  = FBA_raw_solNEFAKD.v(exchRxnInds);
    %exchRxn_raw_afterKU(j,:)  = FBA_raw_solKU.v(exchRxnInds);
    
    % print for model with medium
    A = [string(objRxnNames(j)), obj_before(j), obj_afterKD(j)];%, obj_afterKU(j)];
    formatSpec = '%s, %.f, %.f\n';%, %.f\n';
    fprintf(fID,formatSpec, A);
    
    formatSpec = ['%s,', repmat('%.f,',1,length(exchRxnInds)-1),'%.f','\n'];
    fprintf(fID_eb,formatSpec, [string(model.rxns(j)), exchRxn_before(j,:)]);
    fprintf(fKD_ea,formatSpec, [string(model.rxns(j)), exchRxn_afterKD(j,:)]);
    %fprintf(fKU_ea,formatSpec, [j, exchRxn_afterKU(j,:)]);
    
    % print for model without medium (_raw_)
    A = [string(objRxnNames(j)), obj_raw_before(j), obj_raw_afterKD(j)];%, obj_raw_afterKU(j)];
    formatSpec = '%s, %.f, %.f\n';%, %.f\n';
    fprintf(fIDr,formatSpec, A);
    
    formatSpec = ['%s,', repmat('%.f,',1,length(exchRxnInds)-1),'%.f','\n'];
    fprintf(fID_ebr,formatSpec, [string(objRxnNames(j)), exchRxn_raw_before(j,:)]);
    fprintf(fKD_ear,formatSpec, [string(objRxnNames(j)), exchRxn_raw_afterKD(j,:)]);
    %fprintf(fKU_ear,formatSpec, [string(objRxnNames(j)), exchRxn_raw_afterKU(j,:)]);
    
    % print for model with medium and NEFA 
    A = [string(objRxnNames(j)), obj_NEFA_before(j), obj_NEFA_afterKD(j)];%, obj_raw_afterKU(j)];
    formatSpec = '%s, %.f, %.f\n';%, %.f\n';
    fprintf(fID_NEFA,formatSpec, A);
    
    formatSpec = ['%s,', repmat('%.f,',1,length(exchRxnInds)-1),'%.f','\n'];
    fprintf(fID_ebr_NEFA,formatSpec, [string(objRxnNames(j)), exchRxn_NEFA_before(j,:)]);
    fprintf(fKD_ear_NEFA,formatSpec, [string(objRxnNames(j)), exchRxn_NEFA_afterKD(j,:)]);
    %fprintf(fKU_ear,formatSpec, [string(objRxnNames(j)), exchRxn_raw_afterKU(j,:)]);
    
    % print for model without medium (_raw_) and NEFA 
    A = [string(objRxnNames(j)), obj_raw_NEFA_before(j), obj_raw_NEFA_afterKD(j)];%, obj_raw_afterKU(j)];
    formatSpec = '%s, %.f, %.f\n';%, %.f\n';
    fprintf(fIDr_NEFA,formatSpec, A);
    
    formatSpec = ['%s,', repmat('%.f,',1,length(exchRxnInds)-1),'%.f','\n'];
    fprintf(fID_ebr_NEFA,formatSpec, [string(objRxnNames(j)), exchRxn_raw_NEFA_before(j,:)]);
    fprintf(fKD_ear_NEFA,formatSpec, [string(objRxnNames(j)), exchRxn_raw_NEFA_afterKD(j,:)]);
    %fprintf(fKU_ear,formatSpec, [string(objRxnNames(j)), exchRxn_raw_afterKU(j,:)]);
    
    % print for model with medium and SMFCA
    A = [string(objRxnNames(j)), obj_SMFCA_before(j), obj_SMFCA_afterKD(j)];%, obj_raw_afterKU(j)];
    formatSpec = '%s, %.f, %.f\n';%, %.f\n';
    fprintf(fID_SMFCA,formatSpec, A);
    
    formatSpec = ['%s,', repmat('%.f,',1,length(exchRxnInds)-1),'%.f','\n'];
    fprintf(fID_ebr_SMFCA,formatSpec, [string(objRxnNames(j)), exchRxn_SMFCA_before(j,:)]);
    fprintf(fKD_ear_SMFCA,formatSpec, [string(objRxnNames(j)), exchRxn_SMFCA_afterKD(j,:)]);
    %fprintf(fKU_ear,formatSpec, [string(objRxnNames(j)), exchRxn_raw_afterKU(j,:)]);
    
    % print for model without medium (_raw_) and NEFA 
    A = [string(objRxnNames(j)), obj_raw_SMFCA_before(j), obj_raw_SMFCA_afterKD(j)];%, obj_raw_afterKU(j)];
    formatSpec = '%s, %.f, %.f\n';%, %.f\n';
    fprintf(fIDr_SMFCA,formatSpec, A);
    
    formatSpec = ['%s,', repmat('%.f,',1,length(exchRxnInds)-1),'%.f','\n'];
    fprintf(fID_ebr_SMFCA,formatSpec, [string(objRxnNames(j)), exchRxn_raw_SMFCA_before(j,:)]);
    fprintf(fKD_ear_SMFCA,formatSpec, [string(objRxnNames(j)), exchRxn_raw_SMFCA_afterKD(j,:)]);
    %fprintf(fKU_ear,formatSpec, [string(objRxnNames(j)), exchRxn_raw_afterKU(j,:)]);
    
    if(mod(j,25)==0)
        j
    end
end

fclose(fID);
fclose(fID_eb);
fclose(fKD_ea);
%fclose(fKU_ea);
fclose(fIDr);
fclose(fID_ebr);
fclose(fKD_ear);
%fclose(fKU_ear);

fclose(fID_NEFA);
fclose(fID_eb_NEFA);
fclose(fKD_ea_NEFA);
%fclose(fKU_ea);
fclose(fIDr_NEFA);
fclose(fID_ebr_NEFA);
fclose(fKD_ear_NEFA);
%fclose(fKU_ear);

fclose(fID_SMFCA);
fclose(fID_eb_SMFCA);
fclose(fKD_ea_SMFCA);
%fclose(fKU_ea);
fclose(fIDr_SMFCA);
fclose(fID_ebr_SMFCA);
fclose(fKD_ear_SMFCA);
%fclose(fKU_ear);