## RooUnfold
Install RooUnfold on your local computer https://gitlab.cern.ch/RooUnfold/RooUnfold

Hint: Use 

```bash 
git clone https://gitlab.cern.ch/RooUnfold/RooUnfold.git
``` 
instead of :7999 option






After install , add to your  ~/.bashrc
```bash
source /path_to/RooUnfold/build/setup.sh
``` 


And to ~/.rootlogon.C(example here) the following:
```cpp
{
...
  gSystem->Load("libRooUnfold");
...
}
```

## Mc Jet Tree
* Embedding jet tree file is on SDCC `/gpfs01/star/pwg/prozorov/HFjets/myJetFramework/output_jets.root`, you should copy it to your directory.

* If interested, the macro used to create this file is in `StHIOverlayAngularities.cxx(.h)` in StRoot directory


## Unfolding macros
* Macro to create response Matrix is in `/unfold/createResponseMatrixAngularity3CentralityBins.C`

* The equiprobale binning for the Mc sample (and Unfolded real data) is calculated in `/unfold/fillTestHistsImproved.C`






