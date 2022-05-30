# metabolic-modeling-HIBCH
Perform FBA on human hepatocytes (liver cells) after HIBCH gene knock down.

Build: MATLAB R2018a with Cobra Toolbox, and MOSEK 9.3

Requires: MATLAB, Cobra Toolbox (https://opencobra.github.io/cobratoolbox/stable/index.html), MOSEK 9.3 (https://docs.mosek.com/latest/releasenotes/index.html)

Model: iHepatocytes2322 https://www.ebi.ac.uk/biomodels/model/download/MODEL1402200003.3?filename=MODEL1402200003_url.xml

Files:
- HIBCH_KD_full.m knocks down HIBCH related reactions and performs FBA on the full model
- HIBCH_KD_obj.m knocks down HIBCH related reactions and performs FBA on the reactions output by HIBCH_KD_full.m
- ObjRxns.txt: output from HIBCH_KD_full.m
- ExchRxns.txt: uptake reactions in the model (manually curated)
- ActiveUptakeRxns.txt: active uptake reactions (simulating a medium)
- run-full.sh runs FBA on the full model and identifies the reactions affected by knock down
- run-obj.sh runs only the reactions affected by knock down
