To install the package dirctly from GitLab enter the following in the R-console:

```
install.packages("devtools")
mypersonalaccesstoken <- "abc123???" # put yout personal access token from GitLab here
devtools::install_gitlab(repo="maendle/mvtargetopt", auth_token=mypersonalaccesstoken, host = "https://gitlab.informatik.uni-bremen.de/")
```

To use it using some JSON files you have to define the corresponding paths as parameter:

```
mvTargetOpt::usejson(processdat="processparameter.json", 
                     formdat="formdata_20200211.json", 
                     maindat="http://test.sfb1232.de:85/all_metric_dl", 
                     outjson="out.json")
```

If the chosen input is not appropriate there will be an error (as console output and in out.json).