% Uncomment first run
initCobraToolbox();
changeCobraSolver('mosek','all');
loadMod = readCbModel('iHepatocytes2322.xml');
model = loadMod;

% All units are mmol/gDW/h, millimoles per gram dry weight per hour

% Model to assess before and after KD/KU without medium
rawMod = model;

% Knock down HIBCH associated reactions: (both have lb = 0)
down = 0.5;
KDmod = changeRxnBounds(model, 'HMR_4741', 1000*(1-down), 'u');
KDmod = changeRxnBounds(KDmod, 'HMR_3755',1000*(1-down), 'u');

% Knock down HIBCH associated reactions: (both have lb = 0)
rawModKD = changeRxnBounds(rawMod, 'HMR_4741', 1000*(1-down), 'u');
rawModKD = changeRxnBounds(rawModKD, 'HMR_3755', 1000*(1-down), 'u');


% Objectives (picked by HIBCH_KD_full)
objRxns = importdata('ObjRxns.txt');
objInds = objRxns.data;
objRxnNames = objRxns.textdata;

% active uptake reactions -- file contains names and a matrix of indices
% first column -- indices in terms of iHepatocytes model
% second       -- indices in terms of the uptake reactions index in
%                 the exchange reaction indicies
activeUptakeRxns = importdata('ActiveUptakeRxns.txt');
activeUptakeRxnInds = activeUptakeRxns.data(:,1);
activeUptakeRxnNames = activeUptakeRxns.textdata;

% All exchange reactions (among which are medium components)
exchRxns  = importdata('ExchRxns.txt');
exchRxnInds  = exchRxns.data;
exchRxnNames = exchRxns.textdata;

% Active exchange reactions in terms of indicies in exchRxnInds
actExchRxnInds = activeUptakeRxns.data(:,2);
               

% Set all exchange reaction uptakes to 0 (lower bound = 0) for medium
for i=1:length(exchRxnInds)
    model = changeRxnBounds(model, exchRxnNames(i),0,'l');
    KDmod = changeRxnBounds(KDmod, exchRxnNames(i),0,'l');
end
% Set medium components to -1000 (lower bound = -1000)
for i=1:length(activeUptakeRxnInds)
    model = changeRxnBounds(model, model.rxns(activeUptakeRxnInds(i)),-1000,'l');
    KDmod = changeRxnBounds(KDmod, model.rxns(activeUptakeRxnInds(i)),-1000,'l');
end


% Objectives across the two different setups
  % with medium before, after KD, and after KU
obj_before = length(objInds)*[];
obj_afterKD = length(objInds)*[];

  % without medium before, after KD, and after KU
obj_raw_before = length(objInds)*[];
obj_raw_afterKD = length(objInds)*[];

% Exchange reactions and other reactions
  % with medium
exchRxn_before = length(objInds)*[];
exchRxn_afterKD = length(objInds)*[];

  % without medium
exchRxn_raw_before = length(objInds)*[];
exchRxn_raw_afterKD = length(objInds)*[];

% Loop through reactions that change under HIBCH KD in the hepatocyte model
% print as you go.

% across the exchange reactions!
fID = fopen('results\HIBCHKD-obj.txt','w');
fIDr = fopen('results\HIBCHKD-obj-raw.txt','w');
% Exchange reactions
fID_eb = fopen('results\HIBCH-exchRxns-before.txt','w');
fKD_ea = fopen('results\OHIBCH-exchRxns-afterKD.txt','w');
fID_ebr = fopen('results\HIBCH-exchRxns-before-raw.txt','w');
fKD_ear = fopen('results\HIBCH-exchRxns-afterKD-raw.txt','w');

fprintf(fID, 'ObjName, Before, After KD\n');
fprintf(fID_eb, 'ObjName, Before\n');
fprintf(fKD_ea, 'ObjName, After KD\n');
fprintf(fIDr, 'ObjName, Before, After KD\n');
fprintf(fID_ebr, 'ObjName, Before\n');
fprintf(fKD_ear, 'ObjName, After KD\n');

% Set all objectives to 0
model.c(1:length(model.c)) = 0;
KDmod.c(1:length(model.c)) = 0;
rawMod.c(1:length(model.c)) = 0;
rawModKD.c(1:length(model.c)) = 0;


for j=1:length(objInds)
    % reset previous objective function (if applicable; i.e., j>1)
    if(j>1)
        model.c(objInds(j-1)) = 0;
        KDmod.c(objInds(j-1)) = 0;
        rawMod.c(objInds(j-1)) = 0;
        rawModKD.c(objInds(j-1)) = 0;     
    end
    % Set current objective
    model.c(objInds(j)) = 1;
    KDmod.c(objInds(j)) = 1;
    
    rawMod.c(objInds(j)) = 1;
    rawModKD.c(objInds(j)) = 1;
    
    % run simulations, one for original model, one for KD model
    FBAsol   = optimizeCbModel(model,'max');
    FBAsolKD = optimizeCbModel(KDmod,'max');
    
    FBA_raw_sol   = optimizeCbModel(rawMod,'max');
    FBA_raw_solKD = optimizeCbModel(rawModKD,'max');
    
    % fetch data for objectives
    obj_before(j) = FBAsol.v(objInds(j));
    obj_afterKD(j)  = FBAsolKD.v(objInds(j));
    
    obj_raw_before(j) = FBA_raw_sol.v(objInds(j));
    obj_raw_afterKD(j)  = FBA_raw_solKD.v(objInds(j));
    
    % store flux results from exchange reactions as the jth row of exchRxn_    
    exchRxn_before(j,:) = FBAsol.v(exchRxnInds);
    exchRxn_afterKD(j,:)  = FBAsolKD.v(exchRxnInds);
    
    % store flux results from exchange reactions as the jth row of exchRxn_    
    exchRxn_raw_before(j,:) = FBA_raw_sol.v(exchRxnInds);
    exchRxn_raw_afterKD(j,:)  = FBA_raw_solKD.v(exchRxnInds);
    
    % print for model with medium
    A = [string(objRxnNames(j)), obj_before(j), obj_afterKD(j)];
    formatSpec = '%s, %.f, %.f\n';%, %.f\n';
    fprintf(fID,formatSpec, A);
    
    formatSpec = ['%s,', repmat('%.f,',1,length(exchRxnInds)-1),'%.f','\n'];
    fprintf(fID_eb,formatSpec, [string(model.rxns(j)), exchRxn_before(j,:)]);
    fprintf(fKD_ea,formatSpec, [string(model.rxns(j)), exchRxn_afterKD(j,:)]);
    
    % print for model without medium (_raw_)
    A = [string(objRxnNames(j)), obj_raw_before(j), obj_raw_afterKD(j)];%, obj_raw_afterKU(j)];
    formatSpec = '%s, %.f, %.f\n';%, %.f\n';
    fprintf(fIDr,formatSpec, A);
    
    formatSpec = ['%s,', repmat('%.f,',1,length(exchRxnInds)-1),'%.f','\n'];
    fprintf(fID_ebr,formatSpec, [string(objRxnNames(j)), exchRxn_raw_before(j,:)]);
    fprintf(fKD_ear,formatSpec, [string(objRxnNames(j)), exchRxn_raw_afterKD(j,:)]);
end

fclose(fID);
fclose(fID_eb);
fclose(fKD_ea);
fclose(fIDr);
fclose(fID_ebr);
fclose(fKD_ear);
