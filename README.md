# metabolic-modeling-HIBCH
Perform FBA on human hepatocytes (liver cells) after HIBCH gene knock down.

Requirements
====
Bash for wrapper scripts, `MATLAB` with `Mosek 9.3` and `Cobra toolbox`
* Cobra Toolbox (https://opencobra.github.io/cobratoolbox/stable/index.html)
*  MOSEK 9.3 (https://docs.mosek.com/latest/releasenotes/index.html)

Wrapper scripts
====
* `run-full.sh` -- wrapper script to run flux balance analysis to metabolically prioritize each reaction in the genome-scale metabolic reconstruction of human liver cells
  * outputs `ObjRxns.txt` for use with `run-obj.sh` -- this is committed to the repository for convenience.
* `run-obj.sh`  -- wrapper script to run flux balance analysis and output only the reactions -- termed objectives -- that were affected by HIBCH knock down

Code
====
* HIBCH_KD_full.m knocks down HIBCH related reactions and performs FBA on the full model
* HIBCH_KD_obj.m knocks down HIBCH related reactions and performs FBA on the reactions output by HIBCH_KD_full.m

Model iHepatocytes2322 https://www.ebi.ac.uk/biomodels/model/download/MODEL1402200003.3?filename=MODEL1402200003_url.xml
====

Additional files
====
* `ObjRxns.txt` is output from `HIBCH_KD_full.m` listing the reactions that were affected by HIBCH knock down
* `ExchRxns.txt` are the uptake reactions in the model (manually curated)
* `ActiveUptakeRxns.txt` are the active uptake reactions (simulating a medium)
