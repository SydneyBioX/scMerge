# profvis::profvis({
#   BiocSingular::runExactSVD(Y0)
#   BiocSingular::runExactSVD(Y0, fold = 5)
#   
#   BiocSingular::runRandomSVD(Y0)
#   BiocSingular::runRandomSVD(Y0, fold = 5)
# })
# 
# Y0_da = as(Y0, "HDF5Array")
# 
# profvis::profvis({
#   BiocSingular::runExactSVD(Y0_da)
#   BiocSingular::runExactSVD(Y0_da, fold = 5)
#   
#   BiocSingular::runRandomSVD(Y0_da)
#   BiocSingular::runRandomSVD(Y0_da, fold = 5)
# })
# 
# Y0_sa = as(Y0, "dgCMatrix")
# class(Y0_sa)
# 
# profvis::profvis({
#   BiocSingular::runExactSVD(Y0_sa)
#   BiocSingular::runExactSVD(Y0_sa, fold = 5)
#   
#   BiocSingular::runRandomSVD(Y0_sa)
#   BiocSingular::runRandomSVD(Y0_sa, fold = 5)
# })
# 
# 
# 
# pryr::object_size(Y0)
# pryr::object_size(Y0_da)
# pryr::object_size(Y0_sa)
# 
# class(Y0 %*% t(Y0))
# class(Y0_da %*% t(Y0_da))
# class(Y0_sa %*% t(Y0_sa))
# class(eigenMatMult(Y0_da,t(Y0_da)))
# class(eigenMatMult(Y0_sa,t(Y0_sa)))