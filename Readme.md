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
