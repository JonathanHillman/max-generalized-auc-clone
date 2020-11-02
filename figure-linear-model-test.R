

testFold.vec <- Sys.glob("neuroblastoma-data/data/*/cv/*/testFolds/*")

testFold.vec <- c(#file.path("neuroblastoma-data/data/ATAC_JV_adipose/cv/equal_labels/testFolds/4"),
                  file.path("neuroblastoma-data/data/detailed/cv/R-3.6.0-profileSize/testFolds/3"),
                  file.path("neuroblastoma-data/data/H3K27ac-H3K4me3_TDHAM_BP/cv/equal_labels/testFolds/2"),
                  file.path("neuroblastoma-data/data/H3K4me3_TDH_other/cv/equal_labels/testFolds/1"),
                  file.path("neuroblastoma-data/data/H3K4me3_XJ_immune/cv/equal_labels/testFolds/1"),
                  file.path("neuroblastoma-data/data/H3K4me3_XJ_immune/cv/equal_labels/testFolds/2"),
                  file.path("neuroblastoma-data/data/H3K4me3_XJ_immune/cv/equal_labels/testFolds/4"),
                  file.path("neuroblastoma-data/data/systematic/cv/R-3.6.0-profileSize/testFolds/1"),
                  file.path("neuroblastoma-data/data/detailed/cv/chrom/testFolds/15"),
                  file.path("neuroblastoma-data/data/detailed/cv/R-3.6.0-chrom/testFolds/15"))

one.testFold.path <- file.path("neuroblastoma-data/data/ATAC_JV_adipose/cv/equal_labels/testFolds/2")
#testFold.path <- one.testFold.path

OneFold <- function(testFold.path)
{
  # install.packages("data.table")
  # install.packages("ggplot2")
  # install.packages("dplyr")
  library(data.table)
  library(ggplot2)
  library(dplyr)
  
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
      #system(file.path(data.dir, file.path(paste("unxz", xz.file)))
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
  #cv.folds <- 1:n.val.folds
  cv.folds <- unique(folds.dt[fold != test.fold]$fold)
  
  #set test folds 
  folds.dt[fold == test.fold, set := "test"]
  folds.dt[fold != test.fold, set:= "subtrain"]
  
  cv.dt.list <- list()
  for(val.fold in cv.folds)
  {
  folds.dt[set != "test",] <- folds.dt[set != "test",] %>% 
    mutate(set = ifelse(fold == val.fold, "validation", "subtrain"))
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
  
  
  
  fit.model <- function(X.train, y.train, set.list, num.iterations, seed)
  {
    iteration.dt.list <- list()
    for (init.name in names(init.fun.list))
    {
      init.fun <- init.fun.list[[init.name]]
      set.seed(seed)
      int.weights <- init.fun()
      
      weight.vec <- int.weights[-1]
      intercept <- int.weights[1]
      computeAUM <- function(w, i, is.set) {
        pred.pen.vec <- (X.finite %*% w) + i
        pred.dt <- data.table(
          sequenceID = rownames(pred.pen.vec),
          pred.log.lambda = as.numeric(pred.pen.vec)
        )
        set.dt <- pred.dt[is.set]
        penaltyLearning::ROChange(data.list$evaluation, set.dt, "sequenceID")
      }
      for (iteration in 1:num.iterations) {
        summary.dt.list <- list()
        set.roc.list <- list()
        for (set in names(set.list)) {
          set.roc.list[[set]] <-
            computeAUM(weight.vec, intercept, set.list[[set]])
          summary.dt.list[[set]] <-
            with(set.roc.list[[set]],
                 data.table(set,
                            thresholds[threshold == "predicted"],
                            auc,
                            aum))
        }
        summary.dt <- do.call(rbind, summary.dt.list)
        iteration.dt.list[[paste(seed, init.name, iteration)]] <-
          data.table(seed, init.name, iteration, summary.dt)
        cat(
          sprintf(
            "it=%d seed=%d init=%s\n",
            iteration,
            seed,
            init.name
          )
        )
        g.dt <- set.roc.list[["subtrain"]][["aum.grad"]]
        ## If aum.grad has some problems with no changes in error then
        ## they may be missing.
        g.vec <- rep(0, ncol(neg.t.X.subtrain))
        names(g.vec) <- colnames(neg.t.X.subtrain)
        g.vec[g.dt[["sequenceID"]]] <- g.dt[["lo"]]
        direction.vec <- neg.t.X.subtrain %*% g.vec
        take.step <- function(s) {
          weight.vec + s * direction.vec
        }
        
        #line search
        set.aum.list <- list()
        for (step.size in 10 ^ seq(-10, 0, by = 0.5)) {
          new.weight.vec <- take.step(step.size)
          for (set in "subtrain") {
            set.roc <- computeAUM(new.weight.vec, 0, set.list[[set]])
            set.aum.list[[paste(step.size, set)]] <-
              data.table(step.size,
                         set,
                         aum = set.roc$aum,
                         intercept = set.roc$thresholds[threshold == "min.error", (max.thresh +
                                                                                     min.thresh) / 2])
          }#line search
        }
        set.aum <- do.call(rbind, set.aum.list)
        best.dt <- set.aum[, .SD[min(aum) == aum], by = set]
        weight.vec <- take.step(best.dt[["step.size"]])
        intercept <- best.dt[["intercept"]]
      }#iteration
    }
    output.dt <- data.table(do.call(rbind, iteration.dt.list),
                            data.name,
                            cv.type,
                            test.fold)
    
    output.dt
    
  }
  
  for(curr.seed in 1:4)
  {
    cv.dt.list[[paste(curr.seed, val.fold)]] <- data.table( fit.model(X.train, 
                                                                 y.train, 
                                                                 set.list,
                                                                 num.iterations = 50, 
                                                                 seed = curr.seed),
                                                       val.fold = paste(val.fold))
    
  }#seed
  cat(sprintf("Validation fold %d complete\n", val.fold))
  }#validation fold
  
  cv.dt <- do.call(rbind, cv.dt.list)
  
  cv.csv <- file.path(testFold.path, "linear-model-cv-aum.csv")
  data.table::fwrite(cv.dt, cv.csv)
  
  for(curr.set in unique(cv.dt$set))
  {
  mean.dt <- cv.dt[, .(mean.aum=mean(aum)), by=.(seed, iteration, set)]
  set.cv.dt <- cv.dt[set == curr.set]
  min.set.cv.dt <- set.cv.dt[, .(iteration = which.min(aum), min.aum = min(aum)), by=.(seed, val.fold)]
  set.mean.dt <- mean.dt[set == curr.set]
  set.min.mean.dt <- set.mean.dt[, .(iteration = which.min(mean.aum), min.aum = min(mean.aum)), by=.(seed) ]
  
  out.png.path <- file.path(testFold.path, paste0("linear-model-", curr.set, "-aum-line-graph.png"))
  png(out.png.path)
  
  print(ggplot(data = set.cv.dt) +
    geom_line(mapping = aes(x = iteration, y = aum, color = val.fold)) +
    geom_point(data = min.set.cv.dt, mapping = aes(x = iteration, y = min.aum, color = val.fold)) +
    geom_line(data = set.mean.dt, mapping = aes(x = iteration, y = mean.aum), size = 1) +
    geom_point(data = set.min.mean.dt, mapping = aes(x = iteration, y = min.aum), size = 3) +
    facet_wrap(.~seed) +
    ggtitle(paste(curr.set, "AUM with for each validation fold (mean in bold)", collapse = "")))
  
  dev.off()
  
  if(curr.set == "validation")
  { 
    selected.iter.dt <- data.table( iteration = round(set.min.mean.dt$iteration), seed = set.min.mean.dt$seed)
  }
  }
  
  
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
    
    
    post.cv.dt.list <- list()
    for(curr.seed in 1:4)
    {
      post.cv.dt.list[[paste(curr.seed)]] <- fit.model(X.train, 
                                                  y.train,
                                                  set.list,
                                                  num.iterations = selected.iter.dt[seed == curr.seed]$iteration, 
                                                  seed = curr.seed)
    }#seed
      
    post.cv.dt <- do.call(rbind, post.cv.dt.list)
  
    post.cv.csv <- file.path(testFold.path, "linear-model-post-cv-aum.csv")
    data.table::fwrite(post.cv.dt, post.cv.csv)
      
  test.post.cv.dt <- post.cv.dt[set == "test"]
  
  best.linear.dt <- test.post.cv.dt[, .(aum = min(aum),
                                        type = "best.linear"),
                                    by = .(seed)]
  selected.dt <-
    test.post.cv.dt[, .(aum = last(aum),
                        type = "selected"), 
                    by = seed]
  initial.dt <-
    test.post.cv.dt[, .(aum = first(aum),
                        type = "initial"),  
                    by = seed]
  
  test.aum.dt <- rbind(best.linear.dt, selected.dt, initial.dt)
  
  out.png.path <- file.path(testFold.path, "linear-model-test-aum-comparison.png")
  png(out.png.path)
  
  print(ggplot(data = test.aum.dt) +
    geom_point(aes(x = aum, y = type, color = seed)))
}

OneFold(one.testFold.path)


for(curr.testFold in testFold.vec)
{
  OneFold(curr.testFold)
}


