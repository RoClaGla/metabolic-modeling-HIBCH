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
% Implement medium for mod; leave rawMod alone
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

% Objective value before and after KD
% Model with medium implementation
obj_before = length(model.c)*[];
obj_afterKD = length(model.c)*[];
%obj_afterKU = length(objInds)*[];

% Model without medium implementation
obj_raw_before = length(model.c)*[];
obj_raw_afterKD = length(model.c)*[];
%obj_raw_afterKU = length(objInds)*[];

% Exchange reactions with medium
exchRxn_before = length(model.c)*[];
exchRxn_afterKD = length(model.c)*[];
%exchRxn_afterKU = length(objInds)*[];

% without medium
exchRxn_raw_before = length(model.c)*[];
exchRxn_raw_afterKD = length(model.c)*[];
%exchRxn_raw_afterKU = length(objInds)*[];

% Loop through reactions that change under HIBCH KD in the hepatocyte model
% print as you go; you need one file per reaction? No, print means and stds
% across the exchange reactions!
fID = fopen('results\Full-HIBCHKD-obj.txt','w');
fIDr = fopen('results\Full-HIBCHKD-obj-raw.txt','w');

% Exchange reactions
fID_eb = fopen('results\Full-HIBCHKD-exchRxns-before.txt','w');
fKD_ea = fopen('results\Full-HIBCHKD-exchRxns-afterKD.txt','w');
%fKU_ea = fopen('results\HIBCHKD-exchRxns-afterKU.txt','w');

fID_ebr = fopen('results\Full-HIBCHKD-exchRxns-before-raw.txt','w');
fKD_ear = fopen('results\Full-HIBCHKD-exchRxns-afterKD-raw.txt','w');
%fKU_ear = fopen('results\HIBCHKD-exchRxns-afterKU-raw.txt','w');

% Interesting reactions
fID_ib = fopen('results\Full-HIBCHKD-intRxns-before.txt','w');
fKD_ia = fopen('results\Full-HIBCHKD-intRxns-afterKD.txt','w');
%fKU_ia = fopen('results\HIBCHKD-intRxns-afterKU.txt','w');

fID_ibr = fopen('results\Full-HIBCHKD-intRxns-before-raw.txt','w');
fKD_iar = fopen('results\Full-HIBCHKD-intRxns-afterKD-raw.txt','w');
%fKU_iar = fopen('results\HIBCHKD-intRxns-afterKU-raw.txt','w');

fprintf(fID, 'ObjName, Before, After KD\n');

fprintf(fID_eb, 'ObjName, Before\n');
fprintf(fKD_ea, 'ObjName, After KD\n');
%fprintf(fKU_ea, 'ObjName, After KU\n');

fprintf(fID_ib, 'ObjName, Before\n');
fprintf(fKD_ia, 'ObjName, After KD\n');
%fprintf(fKU_ia, 'ObjName, After KU\n');

fprintf(fIDr, 'ObjName, Before, After KD\n');

fprintf(fID_ebr, 'ObjName, Before\n');
fprintf(fKD_ear, 'ObjName, After KD\n');
%fprintf(fKU_ear, 'ObjName, After KU\n');


fprintf(fID_ibr, 'ObjName, Before\n');
fprintf(fKD_iar, 'ObjName, After KD\n');
%fprintf(fKU_iar, 'ObjName, After KU\n');

% Set all objectives to 0
model.c(1:length(model.c)) = 0;
KDmod.c(1:length(model.c)) = 0;
%KUmod.c(1:length(mod.c)) = 0;
rawMod.c(1:length(model.c)) = 0;
rawModKD.c(1:length(model.c)) = 0;
%rawModKU.c(1:length(mod.c)) = 0;

for j=1:length(model.c)
    % reset previous objective function (if applicable; i.e., j>1)
    if(j>1)
        model.c(j-1) = 0;
        KDmod.c(j-1) = 0;
        %KUmod.c(j-1) = 0;
        rawMod.c(j-1) = 0;
        rawModKD.c(j-1) = 0;
        %rawModKU.c(j-1) = 0;        
    end
    % Set current objective
    model.c(j) = 1;
    KDmod.c(j) = 1;
    %KUmod.c(j) = 1;
    
    rawMod.c(j) = 1;
    rawModKD.c(j) = 1;
    %rawModKU.c(j) = 1;
    
    % run simulations, one for original model, one for KD model
    FBAsol   = optimizeCbModel(model,'max');
    FBAsolKD = optimizeCbModel(KDmod,'max');
    %FBAsolKU = optimizeCbModel(KUmod,'max');
    
    FBA_raw_sol   = optimizeCbModel(rawMod,'max');
    FBA_raw_solKD = optimizeCbModel(rawModKD,'max');
    %FBA_raw_solKU = optimizeCbModel(rawModKU,'max');
    
    % fetch data for objectives
    obj_before(j) = FBAsol.v(j);
    obj_afterKD(j)  = FBAsolKD.v(j);
    %obj_afterKU(j)  = FBAsolKU.v(j);
    
    obj_raw_before(j) = FBA_raw_sol.v(j);
    obj_raw_afterKD(j)  = FBA_raw_solKD.v(j);
    %obj_raw_afterKU(j)  = FBA_raw_solKU.v(j);
    
    % store flux results from exchange reactions as the jth row of exchRxn_    
    exchRxn_before(j,:) = FBAsol.v(exchRxnInds);
    exchRxn_afterKD(j,:)  = FBAsolKD.v(exchRxnInds);
    %exchRxn_afterKU(j,:)  = FBAsolKU.v(exchRxnInds);
    
    % store flux results from exchange reactions as the jth row of exchRxn_    
    exchRxn_raw_before(j,:) = FBA_raw_sol.v(exchRxnInds);
    exchRxn_raw_afterKD(j,:)  = FBA_raw_solKD.v(exchRxnInds);
    %exchRxn_raw_afterKU(j,:)  = FBA_raw_solKU.v(exchRxnInds);
    
    % print for model with medium
    A = [string(model.rxns(j)), obj_before(j), obj_afterKD(j)];%, obj_afterKU(j)];
    formatSpec = '%s, %.f, %.f\n';%, %.f\n';
    fprintf(fID,formatSpec, A);
    
    formatSpec = ['%s,', repmat('%.f,',1,length(exchRxnInds)-1),'%.f','\n'];
    fprintf(fID_eb,formatSpec, [string(model.rxns(j)), exchRxn_before(j,:)]);
    fprintf(fKD_ea,formatSpec, [string(model.rxns(j)), exchRxn_afterKD(j,:)]);
    %fprintf(fKU_ea,formatSpec, [j, exchRxn_afterKU(j,:)]);
    
    % print for model without medium
    A = [string(model.rxns(j)), obj_raw_before(j), obj_raw_afterKD(j)];%, obj_raw_afterKU(j)];
    formatSpec = '%s, %.f, %.f\n';%, %.f\n';
    fprintf(fIDr,formatSpec, A);
    
    formatSpec = ['%s,', repmat('%.f,',1,length(exchRxnInds)-1),'%.f','\n'];
    fprintf(fID_ebr,formatSpec, [string(model.rxns(j)), exchRxn_raw_before(j,:)]);
    fprintf(fKD_ear,formatSpec, [string(model.rxns(j)), exchRxn_raw_afterKD(j,:)]);
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

% Medium model
% Sort objectives by the changes induced by KD
[objValueKD,objIndexKD] = sort(obj_afterKD-obj_before,'descend');
% Fetch the nonzero elements
orderedObjValueKD = length(objValueKD(objValueKD>1|objValueKD<-1))*[];
orderedObjIndexKD = length(objValueKD(objValueKD>1|objValueKD<-1))*[];
k = 1; % running index
for i=1:length(objValueKD)
   if(objValueKD(i)>1||objValueKD(i)<-1)
       orderedObjValueKD(k) = objValueKD(i);
       orderedObjIndexKD(k) = objIndexKD(i);
       k = k + 1;
   end
end
orderedObjValueKD;
orderedObjIndexKD;

% RAW MODEL
% Sort objectives by the changes induced by KD
[objValueKD_raw,objIndexKD_raw] = sort(obj_raw_afterKD-obj_raw_before,'descend');
% Fetch the nonzero elements
orderedObjValueKD_raw = length(objValueKD(objValueKD_raw>1|objValueKD_raw<-1))*[];
orderedObjIndexKDo_raw = length(objValueKD(objValueKD_raw>1|objValueKD_raw<-1))*[];
k = 1; % running index
for i=1:length(objValueKD_raw)
   if(objValueKD_raw(i)>1||objValueKD_raw(i)<-1)
       orderedObjValueKD_raw(k) = objValueKD_raw(i);
       orderedObjIndexKDo_raw(k) = objIndexKD_raw(i);
       k = k + 1;
   end
end
orderedObjValueKD_raw;
orderedObjIndexKDo_raw;

fID = fopen('rxnsData\ObjRxns.txt','w');
for i=1:length(orderedObjIndexKD)
    formatspec = '%s, %s\n';
    fprintf(fID, formatspec, [string(model.rxns(orderedObjIndexKD(i))), orderedObjIndexKD(i) ]);
end
fclose(fID);