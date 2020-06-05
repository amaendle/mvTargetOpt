# Install the package

To install the package dirctly from GitLab enter the following in the R-console:

```
install.packages("devtools")
mypersonalaccesstoken <- "abc123???" # put yout personal access token from GitLab here
devtools::install_gitlab(repo="maendle/mvtargetopt", auth_token=mypersonalaccesstoken, host = "https://gitlab.informatik.uni-bremen.de/")
```

# References

MDPI and ACS Style
Bader, A.; Toenjes, A.; Wielki, N.; Mändle, A.; Onken, A.-K.; Hehl, A.; Meyer, D.; Brannath, W.; Tracht, K. Parameter Optimization in High-Throughput Testing for Structural Materials. Materials 2019, 12, 3439.

# Acknowledgments

Financial support of subproject P02 ‘Heuristic, Statistical and Analytical Experimental Design’ of the Collaborative Research Center SFB 1232 "Farbige Zustände" by the German Research Foundation (DFG) is gratefully acknowledged.