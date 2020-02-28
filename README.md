# EnrichR
Functional Enrichment Analysis and Network Construction 
#### PS: EnrichR will be changed to richR shortly
## Description  
__EnrichR__ is a package can be used for functional enrichment analysis and network construction based on enrichment analysis results. It supported almost all species pubished by ENSEMBL and included with Bioconductor. Now the _EnrichR_ provide function to direct download annotation dataset from the MsigDB. 
## Dependencies  
R>2.15
## Installation
```   
library(devtools)    
install_github("guokai8/EnrichR")
```
## Getting started
```
library(EnrichR)
```  
More detail please see [vignette](https://github.com/guokai8/EnrichR/wiki)
```    
vignette("EnrichR")
```   
## Some useful commands
If you want tranform ID from one type another type("SYMBOL"->"ENSEMBL")
``` 
idconvert(keys=vector_of_symbols,species="human",fkeytype= "SYMBOL",tkeytype="ENSEMBL")
```  
If you want have more details include input data and enrichment results
```  
getdetail(enrichres,input_data)
```  
## Contact information
I still working on this package and will add more functions here.   
For any questions please contact guokai8@gmail.com  
