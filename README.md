# DoubleRobGam
*Ilaria Prosdocimi*   (prosdocimi.ilaria@gmail.com) 

An R package to obtain (robust) smooth estimates for both mean and dispersion functions in extended Generalised Additive Models. The methods are presented in the PhD thesis:

Prosdocimi, I. (2010). Smooth and robust estimation of mean and dispersion functions in regression models. PhD thesis, KULeuven, available [here](https://lirias.kuleuven.be/handle/123456789/280610).


The thesis itself was mostly based on the following papers: 

* Gijbels, I. and Prosdocimi, I. (2012). Flexible Mean and Dispersion Function Estimation in Extended Generalized Additive Models, *Communications in Statistics - Theory and Methods*, **41**, DOI:10.1080/03610926.2012.654881

* Croux, C., Gijbels, I. and Prosdocimi, I. (2012).  Robust Estimation of Mean and Dispersion Functions in Extended Generalized Additive Models. *Biometrics*, **268**, 31-44. doi:10.1111/j.1541-0420.2011.01630.x.

* Gijbels I. and  Prosdocimi I. (2011). Smooth estimation of mean and dispersion function in extended Generalized Additive Models with application to Italian Induced Abortion data. *Journal of Applied Statistics*, **38**. DOI:10.1080/02664763.2010.550039

* Gijbels, I., Prosdocimi, I. and  Claeskens, G. (2010). Nonparametric estimation of mean and dispersion functions in extended Generalized Linear Models. *Test*, **19**. DOI:10.1007/s11749-010-0187-1

The coding effort has mostly been done by me (Ilaria Prosdocimi) - the co-authors of the above mentioned papers should not be held responsible for any bugs present in the code. 


### Disclaimer 
The code has been developed some years ago and I never managed to do all the upgrades I had planned. I have nevertheless tried to put some effort in documenting the functions and decided to make the code available as is. I am aware of some issues and I expect users to find even more bugs and bottlenecks. I am happy to receive pull requests and get informed about bugs, although I can not guarantee I will have much time to spend improving these functions. 

### Usage 
To use the library in R, you would need to install it from GitHub via the `devtools` library: 

```
# the remotes package is needed to be able to load the package from GitHub
# install.packages("remotes")
remotes::install_github("ilapros/DoubleRobGam")
library(DoubleRobGam)
``` 

The package vignette gives some examples of how to use the `DoubleGam` and `DoubleRobGam` functions. See the vignette by typing:

```
browseVignettes("DoubleRobGam")
``` 

