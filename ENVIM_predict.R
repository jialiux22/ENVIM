
ENVIM_predict = function(microbio.train,
                 microbio.test,
                 metab.train,
                 seed = 1234,
                 outputdirectory = getwd(),
                 fold_rf = 10,
                 fold_ENVIM = 10,
                 alpha_range = seq(0,1,0.1)){
  
  
  library(glmnet, quietly = TRUE)
  library(MASS, quietly = TRUE)
  library(caret, quietly = TRUE)
  library(GenABEL, quietly = TRUE)
  library(caTools, quietly = TRUE)
  
  print("Welcome to use ENVIM algorithm")
  
  ## set up output directory  
  setwd(outputdirectory)
  
  ## set up variable of result
  time = c()
  mse.train = c()
  train.cor = c()
  test.pred = data.frame(numeric(nrow(microbio.test)))
  rownames(test.pred) = rownames(microbio.test)
  name.gene.train = matrix(NA,ncol = ncol(metab.train),nrow = 10000)
  name.gene.test = matrix(NA,ncol = ncol(metab.train),nrow = 10000)
  colnames(name.gene.train) = colnames(metab.train)
  colnames(name.gene.test) = colnames(metab.train)
  
  print("Start to running ENVIM algorithm")
  
  ## running algorithm for each metabolite
  
  for (col.n in 1:ncol(metab.train)) {
    
    t = Sys.time()
    y.train = metab.train[,col.n]
    
    names(y.train) = rownames(metab.train)
    
    ## use Box-Cox transformation for training metabolite
    b.train = boxcox(y.train~1,plotit = F)
    lambda.train = b.train$x[which.max(b.train$y)]
    if (lambda.train!=0){
      y.train = (y.train^(lambda.train) - 1)/lambda.train
    }
    
    if (lambda.train==0){
      y.train = log(y.train)}
    
    ## microbial data
    otu_t = data.frame(microbio.train)
    
    
    ## merge training metabolite with microbial data
    merge.train.sig = merge(y.train,otu_t,by=0,sort = F)
    
    rownames(merge.train.sig) = merge.train.sig$Row.names
    
    merge.train.sig = merge.train.sig[,-1]
    
    ## microbial data transformation (ranked transformation)
    merge.train.sig[,2:ncol(merge.train.sig)] =  apply(merge.train.sig[,2:ncol(merge.train.sig)],2,rntransform)
    
    
    ## run 10-fold cross-validated random forest model to find variable importance of microbial features
    set.seed(seed)
    
    control = trainControl(method="cv", number=fold_rf)
    
    model = train(x ~., data=merge.train.sig, method="rf", trControl=control)
    
    importance <- varImp(model)
    
    imp_mat = as.matrix(importance$importance)
    
    ## find variable importance level (times of 10)
    tune.level = sort(as.numeric(unique(floor(imp_mat/10))),decreasing = T) * 10
    
    upper = 1
    while (length(colnames(merge.train.sig)[c(1,which(imp_mat>=tune.level[upper])+1)])<3) {
      upper = upper + 1
    }
    
    tune.level = tune.level[tune.level<=tune.level[upper]]
    
    
    ##run 10-fold cross validated ENVIM
    accu.mat = data.frame(rep(tune.level,each = length(alpha_range)))
    
    accu.mat[,2] = rep(alpha_range,length(tune.level))
    
    
    accu.lambda = c()
    
    accu.error = c()
    
    ## Assign an fold id (keep reproducible)
    set.seed(seed)
    fold.id = sample(rep(seq(fold_ENVIM), length = nrow(merge.train.sig)))
    
    
    for (level in tune.level) { 
      
      
      train = merge.train.sig[,colnames(merge.train.sig)[c(1,which(imp_mat>=level)+1)]]
      
      
      for (i in alpha_range) {
        
        fit = cv.glmnet(data.matrix(train[,-1]),train[,1],alpha=i,foldid = fold.id)
        
        accu.error = c(accu.error,fit$cvm[which(fit$lambda == fit$lambda.1se)])
        
        accu.lambda =  c(accu.lambda, fit$lambda.1se)
        
      }
      
    }
    
    
    accu.mat[,3] = accu.lambda
    
    accu.mat[,4] = accu.error
    
    colnames(accu.mat) = c("level_imp","alpha","lambda.1se",'err.1se')
    
    ## Find parameter of ENVIM (alpha, lambda, and level) (minimum mean square error) 
    min.alpha = accu.mat$alpha[which.min(accu.mat$err.1se)]
    
    min.lambda = accu.mat$lambda.1se[which.min(accu.mat$err.1se)]
    
    min.level = accu.mat$level_imp[which.min(accu.mat$err.1se)]
    
    
    ## Find training model with the parameter such that the model has the minimum mean square error
    train = merge.train.sig[,colnames(merge.train.sig)[c(1,which(imp_mat>=min.level)+1)]]
    
    mod = glmnet(data.matrix(train[,-1]),train[,1],alpha=min.alpha, lambda = min.lambda)
    
    coefficient = coef(mod)
    
    name.coef = rownames(coefficient)
    
    weight = matrix(coefficient)
    
    rownames(weight) =name.coef
    
    ## training stage
    
    ## training weight matrix
    if (length(name.coef[weight!=0][-1]) != 0) {
      
      val.weight.train = paste0(name.coef[weight!=0][-1],'_(',signif(as.numeric(weight[weight!=0])[-1],4),")")
      
      val.weight.train = val.weight.train[order(abs(as.numeric(weight[weight!=0])[-1]),decreasing = T)]
      
      name.gene.train[1:length(name.coef[weight!=0][-1]),col.n] = val.weight.train
    }
    
    predict.y = predict(mod, data.matrix(train[,-1]))
    
    if (lambda.train!=0){
      predict.y = exp(log(1 + lambda.train * predict.y)/lambda.train)
      
      measure.y = exp(log(1 + lambda.train * train[,1])/lambda.train)
    }
    
    if (lambda.train==0){
      predict.y = exp(predict.y)
      measure.y = exp(train[,1])
    }
    
    ## training correlation
    train.cor[col.n] = cor(predict.y,measure.y, method = "spearman")
    
    ## training mean square error
    mse.train[col.n] = mean((predict.y - measure.y)^2)
    
    ##testing stage
    test.x = microbio.test
    
    test.x =  apply(test.x,2,rntransform)
    
    inter.gene = intersect(name.coef,colnames(test.x))
    
    inter.weight = weight[c("(Intercept)",inter.gene),]
    
    ## testing weight matrix
    if (length(names(inter.weight)[inter.weight!=0][-1]) != 0) {
      
      val.weight.test = paste0(names(inter.weight)[inter.weight!=0][-1],'_(',signif(as.numeric(inter.weight[inter.weight!=0])[-1],4),")")
      
      val.weight.test = val.weight.test[order(abs(as.numeric(inter.weight[inter.weight!=0])[-1]),decreasing = T)]
      
      name.gene.test[1:length(names(inter.weight)[inter.weight!=0][-1]),col.n] = val.weight.test
    }
    
    ## find overlapped gene between training microbial features and testing microbial features
    test.x1 = test.x[, inter.gene] 
    
    test.x1 = as.matrix(cbind(rep(1,nrow(microbio.test)),test.x1))
    
    pred.test = test.x1 %*% inter.weight
    
    ## transformed testing metabolite back to the original scale
    if (lambda.train!=0){
      pred.test = exp(log(1 + lambda.train * pred.test)/lambda.train)
    }
    
    if (lambda.train==0){
      pred.test = exp(pred.test)}
    
    ## testing prediction
    test.pred[,(col.n+1)] = pred.test
    
    time[col.n] = round(difftime(Sys.time(), t, units = "mins"),1)
    
    print(paste0("Finish running metabolite ",colnames(metab.train)[col.n], " (with running time ",round(difftime(Sys.time(), t, units = "mins"),1)," mins)"))
    
  }
  
  test.pred = test.pred[,-1]
  
  rownames(test.pred) = rownames(microbio.test)
  
  colnames(test.pred) = colnames(metab.train)
  
  name.gene.train = name.gene.train[1:max(colSums(!is.na(name.gene.train))),]
  
  name.gene.test = name.gene.test[1:max(colSums(!is.na(name.gene.test))),]
  
  write.csv(name.gene.train,"weight.train.csv",quote = F)
  
  write.csv(name.gene.test,"weight.test.csv",quote = F)
  
  
  ## summary table
  sum.table <- data.frame(colnames(metab.train))
  colnames(sum.table)[1] = "Compound_name"
  sum.table$Train_spearman_cor = train.cor
  sum.table$Train_mse = mse.train
  sum.table$Time_mins = time
  
  ## save output 
  write.csv(name.gene.train,"weight.train.csv",quote = F)
  write.csv(name.gene.test,"weight.test.csv",quote = F)
  write.csv(sum.table,"summary.csv",quote = F)
  write.csv(test.pred,"Pred_test.csv",quote = F)
  
  print(paste0("Weight.train.csv, Weight.test.csv, Summary.csv, and Pred_test.csv are stored in ",outputdirectory))
  
}

setwd("~/Desktop/ENVIM")
mtr = read.csv("micro.train.csv",row.names = 1)
mte = read.csv("micro.test.csv",row.names = 1)

metab = read.csv("metab.train.csv",row.names = 1)
setwd("~/Desktop/ENVIM/aa")

ENVIM(microbio.train = mtr,microbio.test = mte,metab.train =metab )
