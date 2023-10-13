
BayesPrism <- function(bk.dat,sc.dat,cell.state.labels,cell.type.labels,ncore=50){
  library(ggplot2)
  library(reshape2)
  library(BayesPrism)
  library(Seurat)
  set.seed(12345)
  
  bk.stat <- plot.bulk.outlier(
    bulk.input=bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID 
    sc.input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
    cell.type.labels=cell.type.labels,
    species="hs", #currently only human(hs) and mouse(mm) annotations are supported
    return.raw=TRUE
    #pdf.prefix="gbm.bk.stat" specify pdf.prefix if need to output to pdf
  )
  
 
  plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                   bulk.input = bk.dat
                   #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
  )
  

 
  
  #Construct a prism object
  myPrism <- new.prism(
    reference=sc.dat.filtered.pc, 
    mixture=bk.dat,
    input.type="count.matrix", 
    cell.type.labels = cell.type.labels, 
    cell.state.labels = cell.state.labels,
    key=NULL,
    outlier.cut=0.01,
    outlier.fraction=0.1,
  )
  
  #Run BayesPrism
  bp.res <- run.prism(prism = myPrism, n.cores=ncore)
  
  return(bp.res)
}

