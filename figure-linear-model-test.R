library(data.table)
library(ggplot2)
library(dplyr)

testFold.vec <- Sys.glob("neuroblastoma-data/data/*/cv/*/testFolds/*")
testFold.vec <- Sys.glob(one.testFold.path)
one.testFold.path <- file.path("C:/Users/jonat/Documents/ROC-Curve-Analysis/neuroblastoma-data/data/ATAC_JV_adipose/cv/equal_labels/testFolds/4")
testFold.path <- one.testFold.path

#OneFold <- function(testFold.path)
{
  #Filename and directory initializations
  cv.path <- dirname(dirname(testFold.path))
  folds.csv <- file.path(cv.path, "folds.csv")
  cv.type <- basename(cv.path)
  test.fold <- basename(testFold.path)
  data.dir <- dirname(dirname(cv.path))
  data.name <- basename(data.dir)
  data.list <- list()
  
  
  #Read the "inputs.csv", "outputs.csv", and "evaluations.csv", for each data 
  for(f in c("inputs", "outputs", "evaluation")){
    f.csv.xz <- file.path(data.dir, paste0(f, ".csv.xz"))
    if(file.exists(f.csv.xz)){
      system(paste("unxz", f.csv.xz))
    }
    f.csv <- file.path(data.dir, paste0(f, ".csv"))
    f.dt <- data.table::fread(f.csv)
    data.list[[f]] <- f.dt
  }
  
  
  ## replace positive fp/fn at end with 0 to avoid AUM=Inf.
  data.list[["evaluation"]][min.log.lambda==-Inf & 0<fn, fn := 0]
  data.list[["evaluation"]][max.log.lambda==Inf & 0<fp, fp := 0]
  
  
  ## read folds.csv for the specific data type
  folds.dt <- data.table::fread(folds.csv)
  #initialize my own validation folds
  n.val.folds <- 4
  cv.folds <- 1:n.val.folds
  #set test folds 
  folds.dt[fold == test.fold, set := "test"]
  folds.dt[fold != test.fold, set:= "subtrain"]
  val.fold.assignments <- rep(cv.folds, l=nrow(folds.dt[fold != test.fold]))
  
  iteration.dt.list <- list()
  for(val.fold in cv.folds)
  {
  folds.dt[set != "test",] <- folds.dt[set != "test",] %>% 
    mutate(set = ifelse(val.fold.assignments == val.fold, "validation", "subtrain"))
  seqs.train <- folds.dt[["sequenceID"]]
  X.all <- scale(data.list$inputs[, -1])
  rownames(X.all) <- data.list$inputs$sequenceID
  X.finite <- X.all[, apply(is.finite(X.all), 2, all)]
  set.list <- list()
  for(s in unique(folds.dt$set)){
    set.list[[s]] <- rownames(X.finite) %in% folds.dt[s==set, sequenceID]
  }
  X.list <- lapply(set.list, function(i)X.finite[i, ])
  neg.t.X.subtrain <- -t(X.list[["subtrain"]])
  y.train <- data.list[["outputs"]][
    seqs.train,
    cbind(min.log.lambda, max.log.lambda),
    on="sequenceID"]
  keep <- apply(is.finite(y.train), 1, any)
  X.train <- X.finite[seqs.train, ]
  init.fun.list <- list(
    IntervalRegressionCV=function(){
      fit <- penaltyLearning::IntervalRegressionCV(
        X.train[keep, ],
        y.train[keep, ])  
      fit[["param.mat"]]
    }
  )
  
  for(seed in 1:4)for(init.name in names(init.fun.list)){
    init.fun <- init.fun.list[[init.name]]
    set.seed(seed)
    int.weights <- init.fun()
    weight.vec <- int.weights[-1]
    intercept <- int.weights[1]
    computeAUM <- function(w, i, is.set){
      pred.pen.vec <- (X.finite %*% w) + i
      pred.dt <- data.table(
        sequenceID=rownames(pred.pen.vec),
        pred.log.lambda=as.numeric(pred.pen.vec))
      set.dt <- pred.dt[is.set]
      penaltyLearning::ROChange(
        data.list$evaluation, set.dt, "sequenceID")
    }
    for(iteration in 1:50){
      summary.dt.list <- list()
      set.roc.list <- list()
      for(set in names(set.list)){
        set.roc.list[[set]] <- computeAUM(weight.vec, intercept, set.list[[set]])
        summary.dt.list[[set]] <- with(set.roc.list[[set]], data.table(
          set,
          thresholds[threshold=="predicted"],
          auc,
          aum))
      }
      summary.dt <- do.call(rbind, summary.dt.list)
      iteration.dt.list[[paste(seed, init.name, iteration, val.fold)]] <- data.table(
        seed, init.name, iteration, paste(val.fold), summary.dt)
      cat(sprintf(
        "it=%d seed=%d init=%s val=%d\n",
        iteration, seed, init.name, val.fold))
      g.dt <- set.roc.list[["subtrain"]][["aum.grad"]]
      ## If aum.grad has some problems with no changes in error then
      ## they may be missing.
      g.vec <- rep(0, ncol(neg.t.X.subtrain))
      names(g.vec) <- colnames(neg.t.X.subtrain)
      g.vec[
        g.dt[["sequenceID"]]
      ] <- g.dt[["lo"]]
      direction.vec <- neg.t.X.subtrain %*% g.vec
      take.step <- function(s){
        weight.vec + s*direction.vec
      }
      
      #line search
      set.aum.list <- list()
      for(step.size in 10^seq(-10, 0, by=0.5)){
        new.weight.vec <- take.step(step.size)
        for(set in "subtrain"){
          set.roc <- computeAUM(new.weight.vec, 0, set.list[[set]])
          set.aum.list[[paste(step.size, set)]] <- data.table(
            step.size, set, aum=set.roc$aum,
            intercept=set.roc$thresholds[
              threshold=="min.error", (max.thresh+min.thresh)/2])
        }#line search
      }
      set.aum <- do.call(rbind, set.aum.list)
      best.dt <- set.aum[, .SD[min(aum)==aum], by=set]
      weight.vec <- take.step(best.dt[["step.size"]])
      intercept <- best.dt[["intercept"]]
    }#iteration
  }#seed/init.name
  
  }#validation fold
  cv.dt <- data.table(
    do.call(rbind, iteration.dt.list),
    data.name, cv.type, test.fold)
  
  mean.dt <- cv.dt[, .(mean.aum=mean(aum)), by=.(seed, iteration, set)]
  val.cv.dt <- cv.dt[set == "validation"]
  val.mean.dt <- mean.dt[set == "validation"]
  val.min.dt <- val.mean.dt[, .(iteration = which.min(mean.aum), min.aum = min(mean.aum)), by=.(seed) ]
  
  ggplot(data = val.cv.dt) +
    geom_line(mapping = aes(x = iteration, y = aum, color = val.fold)) +
    geom_line(data = val.mean.dt, mapping = aes(x = iteration, y = mean.aum)) +
    geom_point(data = val.min.dt, mapping = aes(x = iteration, y = min.aum)) +
    facet_wrap(.~seed)
  
  selected.iter <- round(mean(val.min.dt$iteration))
  
  
  folds.dt[fold == test.fold, set := "test"]
  folds.dt[fold != test.fold, set:= "subtrain"]
  iteration.dt.list <- list ()
  iteration.dt.list <- list()

    seqs.train <- folds.dt[["sequenceID"]]
    X.all <- scale(data.list$inputs[, -1])
    rownames(X.all) <- data.list$inputs$sequenceID
    X.finite <- X.all[, apply(is.finite(X.all), 2, all)]
    set.list <- list()
    for(s in unique(folds.dt$set)){
      set.list[[s]] <- rownames(X.finite) %in% folds.dt[s==set, sequenceID]
    }
    X.list <- lapply(set.list, function(i)X.finite[i, ])
    neg.t.X.subtrain <- -t(X.list[["subtrain"]])
    y.train <- data.list[["outputs"]][
      seqs.train,
      cbind(min.log.lambda, max.log.lambda),
      on="sequenceID"]
    keep <- apply(is.finite(y.train), 1, any)
    X.train <- X.finite[seqs.train, ]
    init.fun.list <- list(
      IntervalRegressionCV=function(){
        fit <- penaltyLearning::IntervalRegressionCV(
          X.train[keep, ],
          y.train[keep, ])  
        fit[["param.mat"]]
      }
    )
    
    for(seed in 1:4)for(init.name in names(init.fun.list)){
      init.fun <- init.fun.list[[init.name]]
      set.seed(seed)
      int.weights <- init.fun()
      weight.vec <- int.weights[-1]
      intercept <- int.weights[1]
      computeAUM <- function(w, i, is.set){
        pred.pen.vec <- (X.finite %*% w) + i
        pred.dt <- data.table(
          sequenceID=rownames(pred.pen.vec),
          pred.log.lambda=as.numeric(pred.pen.vec))
        set.dt <- pred.dt[is.set]
        penaltyLearning::ROChange(
          data.list$evaluation, set.dt, "sequenceID")
      }
      for(iteration in 1:selected.iter){
        summary.dt.list <- list()
        set.roc.list <- list()
        for(set in names(set.list)){
          set.roc.list[[set]] <- computeAUM(weight.vec, intercept, set.list[[set]])
          summary.dt.list[[set]] <- with(set.roc.list[[set]], data.table(
            set,
            thresholds[threshold=="predicted"],
            auc,
            aum))
        }
        summary.dt <- do.call(rbind, summary.dt.list)
        iteration.dt.list[[paste(seed, init.name, iteration)]] <- data.table(
          seed, init.name, iteration, summary.dt)
        cat(sprintf(
          "it=%d seed=%d init=%s\n",
          iteration, seed, init.name))
        g.dt <- set.roc.list[["subtrain"]][["aum.grad"]]
        ## If aum.grad has some problems with no changes in error then
        ## they may be missing.
        g.vec <- rep(0, ncol(neg.t.X.subtrain))
        names(g.vec) <- colnames(neg.t.X.subtrain)
        g.vec[
          g.dt[["sequenceID"]]
        ] <- g.dt[["lo"]]
        direction.vec <- neg.t.X.subtrain %*% g.vec
        take.step <- function(s){
          weight.vec + s*direction.vec
        }
        
        #line search
        set.aum.list <- list()
        for(step.size in 10^seq(-10, 0, by=0.5)){
          new.weight.vec <- take.step(step.size)
          for(set in "subtrain"){
            set.roc <- computeAUM(new.weight.vec, 0, set.list[[set]])
            set.aum.list[[paste(step.size, set)]] <- data.table(
              step.size, set, aum=set.roc$aum,
              intercept=set.roc$thresholds[
                threshold=="min.error", (max.thresh+min.thresh)/2])
          }#line search
        }
        set.aum <- do.call(rbind, set.aum.list)
        best.dt <- set.aum[, .SD[min(aum)==aum], by=set]
        weight.vec <- take.step(best.dt[["step.size"]])
        intercept <- best.dt[["intercept"]]
      }#iteration
    }#seed/init.name
  
  post.cv.dt <- data.table(
    do.call(rbind, iteration.dt.list),
    data.name, cv.type, test.fold)
  post.cv.dt$seed <- paste(post.cv.dt$seed)
  
  test.post.cv.dt <- post.cv.dt[set == "test"]
  
  best.linear.dt <- test.post.cv.dt[, .(aum = min(aum), type = "best.linear"), by = seed ]
  selected.dt <- test.post.cv.dt[, .(aum = aum[selected.iter], type = "selected"), by = seed]
  initial.dt <- test.post.cv.dt[, .(aum = aum[1], type = "initial"), by = seed]
  
  test.aum.dt <- rbind(best.linear.dt, selected.dt, initial.dt)
  
  ggplot(data = test.aum.dt) +
    geom_point(aes(x = aum, y = type, color = seed))
  
}


all.it.list <- list()
for(testFold.i in seq_along(testFold.vec)){
  fdir <- testFold.vec[testFold.i]
  out.csv <- file.path(fdir, "linear-model-aum.csv")
  all.it.list[[testFold.i]] <- if(file.exists(out.csv)){
    data.table::fread(out.csv)
  }else{
    cat(sprintf("%4d / %4d %s\n", testFold.i, length(testFold.vec), fdir))
    iteration.dt <- OneFold(fdir)
    data.table::fwrite(iteration.dt, out.csv)
    iteration.dt
  }
}
all.it <- do.call(rbind, all.it.list)

subtrain.it <- all.it[set=="subtrain"]
subtrain.it[, diff := c(NA, diff(aum)), by=.(init.name, data.name, test.fold, seed)]
subtrain.it[, .(init.name, data.name, test.fold, iteration, aum, diff)]
subtrain.it[diff>1e-6]




validation.it <- all.it[set=="validation"]

valid.best.ids <- all.it[
  set=="validation",
  .SD[which.min(aum), .(iteration)],
  by=.(data.name, test.fold, init.name, seed)]
test.best.ids <- all.it[
  set=="test",
  .SD[which.min(aum), .(iteration)],
  by=.(data.name, test.fold, init.name, seed)]

## model selection.
test.it1 <- all.it[set=="test" & iteration==1]
test.selected <- all.it[set=="test"][valid.best.ids, on=names(valid.best.ids)]
test.best <- all.it[set=="test"][test.best.ids, on=names(test.best.ids)]

## compare with best predictions (no linear model).
best.compare <- best.aum[
  test.best,
  .(data.name, test.fold, init.name, seed, aum, best.aum),
  on=.(data.name, test.fold)]
best.compare[, aum.diff := aum-best.aum]
ggplot()+
  geom_point(aes(
    aum.diff, init.name),
    shape=1,
    data=best.compare)+
  facet_grid(. ~ data.name + test.fold, scales="free", labeller=label_both)+
  theme_bw()+
  scale_x_log10()+
  theme(panel.spacing=grid::unit(0, "lines"))
best.compare[, .(
  min.diff=min(aum.diff),
  max.diff=max(aum.diff)
), by=.(data.name, test.fold, init.name)]

best.pred <- best.aum[
  unique(test.best[, .(data.name, test.fold)]),
  on=.(data.name, test.fold)]
test.show <- rbind(
  data.table(iterations="initial", test.it1),
  data.table(iterations="best.linear", test.best),
  data.table(iterations="selected", test.selected))
ifac <- function(x)factor(
  x, c("initial", "selected", "best.linear", "best.pred"))
test.show[, Iterations := ifac(iterations)]
best.pred[, Iterations := ifac("best.pred")]


test.show[, neg.auc := -auc]
test.show.tall <- melt(
  test.show[init.name=="IntervalRegressionCV"],
  measure.vars=c("neg.auc", "error.percent", "aum"))
test.iCV <- dcast(
  test.show.tall,
  data.name + test.fold + variable + seed ~ iterations)

