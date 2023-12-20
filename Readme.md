Install RooUnfold on your local computer https://gitlab.cern.ch/RooUnfold/RooUnfold

Hint: Use 

```bash 
git clone https://gitlab.cern.ch/RooUnfold/RooUnfold.git
``` 
instead of :7999 option

After installing everything add to your ~/.rootlogon.C(example here) the following:
```cpp
{
...
  gSystem->Load("libRooUnfold");
...
}
```


Embedding jet tree fils is on SDCC `/gpfs01/star/pwg/prozorov/HFjets/myJetFramework/output_jets.root`
If interested, the macro used to create this file is in `StHIOverlayAngularities.cxx(.h)` in StRoot directory


Macro to create response Matrix is in `/unfold/createResponseMatrixAngularity3CentralityBins.C`

The equiprobale binning for the Mc sample (and Unfolded real data) is calculated in `/unfold/fillTestHistsImproved.C`




