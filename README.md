To install the package dirctly from GitLab enter the following in the R-console:

install.packages("devtools")
mypersonalaccesstoken <- "abc123???" # put yout personal access token from GitLab here
devtools::install_gitlab(repo="maendle/mvtargetopt", auth_token=mypersonalaccesstoken, host = "https://gitlab.informatik.uni-bremen.de/")