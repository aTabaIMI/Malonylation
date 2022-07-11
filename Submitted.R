rm(list=ls()) #cleaning the memory froma variables
#loading the required libraries
library(xlsx) 
library(readxl)
library(ggplot2)
library(dplyr)
library(plotrix)
library(ROCit)
library(pROC)
library(PRROC)
library(nnet)
library(neuralnet)
library(e1071)
library(caTools)
library(caret)
library(class)
library(reticulate)
library(keras)
virtualenv_create("myenv")
use_virtualenv("myenv")
library("xgboost")
library(Ckmeans.1d.dp)
library("randomForest")

#Replace the following path with the ones corresponding to your local storage and according to the fasta data
pathC1 <- "C:/Users/Administrator/Documents/UseCases/ProteinsSequences/PSSM/Class_1/Class_1/Class_1/pssm/Class_1"
pathC2 <- "C:/Users/Administrator/Documents/UseCases/ProteinsSequences/PSSM/Class_2/Class_2/Class_2/pssm/Class_2"
pathC3 <- "C:/Users/Administrator/Documents/UseCases/ProteinsSequences/PSSM/Class_3/Class_3/Class_3/pssm/Class_3"
pathC4 <- "C:/Users/Administrator/Documents/UseCases/ProteinsSequences/PSSM/Class_4/Class_4/Class_4/pssm/Class_4"

P1 <- "C:/Users/Administrator/Documents/UseCases/ProteinsSequences/PSSM/Class_1/Class_1/Class_1/pssm/Class_1"
P2 <- "C:/Users/Administrator/Documents/UseCases/ProteinsSequences/PSSM/Class_2/Class_2/Class_2/pssm/Class_2"
P3 <- "C:/Users/Administrator/Documents/UseCases/ProteinsSequences/PSSM/Class_3/Class_3/Class_3/pssm/Class_3"
P4 <- "C:/Users/Administrator/Documents/UseCases/ProteinsSequences/PSSM/Class_4/Class_4/Class_4/pssm/Class_4"
PATH2POSITIONS <- c(P1,P2,P3,P4)
PathModel <- "C:/Users/Administrator/Documents/UseCases/ProteinsSequences/PSSM/Class_1/Class_1/Class_1/pssm/"

Q1 <- "C:/Users/Administrator/Documents/UseCases/ProteinsSequences/PSSM/Class_1/Class_1/Class_1.xlsx"
Q2 <- "C:/Users/Administrator/Documents/UseCases/ProteinsSequences/PSSM/Class_2/Class_2/Class_2.xlsx"
Q3 <- "C:/Users/Administrator/Documents/UseCases/ProteinsSequences/PSSM/Class_3/Class_3/Class_3.xlsx"
Q4 <- "C:/Users/Administrator/Documents/UseCases/ProteinsSequences/PSSM/Class_4/Class_4/Class_4.xlsx"
PATH2NAMES <- c(Q1,Q2,Q3,Q4)
PATH2PHYSCHEIMFTR <- "C:/Users/Administrator/Documents/UseCases/ProteinsSequences/PCP.xlsx"

##Reading the source sequences files
#Replace the following path with the ones corresponding to your local storage and according to the fasta data
path <- "C:/Users/Administrator/Documents/UseCases/ProteinsSequences/NewSource.txt"
MyFullDataOriginal <- readLines( path)
Len <- length(MyFullDataOriginal)
MyFullData <- MyFullDataOriginal[1]
l <- 1

for( i in 1:Len ){ #removing empty lines from sequences
  if(MyFullDataOriginal[i]!=""){
    MyFullData[l]<- MyFullDataOriginal[i]
    l <- l+1}
}

#Reading the labels from Excel file
#Replace the following path with the ones corresponding to your local storage and according to the fasta data
Path2Label <- "C:/Users/Administrator/Documents/UseCases/ProteinsSequences/PositiveSites.xlsx"
MyOriginalLabels <- read_excel(Path2Label)
MyVerifiedLabels <- MyOriginalLabels
MyNames <- vector(mode = "character", length = length(MyOriginalLabels$Codes))
MySites <- vector(mode = "numeric", length = length(MyOriginalLabels$Site))
##Extracing the "K" value locations from sequences
l <- 1
pre <- 0
Index <- 0 #initializing a matrix showing "K" positions
Index[1]<-0
MyRepeatName <- "null"
MyRepeatName[1] <- "null"
SeqID <- 0
j<-0
Lenn <- length(MyFullData)
MyTempSequence <- vector(mode = "character", length = 1602)
TemplateTable <- NULL #the templates found and stored 
MySeqName <- NULL
MySeqName[1] <- NULL
##merging the sequence lines into one line for each sequence
for( i in 1:Lenn ){
  if((strsplit(MyFullData[i],split='')[[1]][1]==">") && (strsplit(MyFullData[i],split='')[[1]][11]=="|") ){#Name line
    SeqID <- SeqID+1
    MySeqName[SeqID] <- paste0(strsplit(MyFullData[i],split='')[[1]][5],strsplit(MyFullData[i],split='')[[1]][6],
                               strsplit(MyFullData[i],split='')[[1]][7],strsplit(MyFullData[i],split='')[[1]][8],
                               strsplit(MyFullData[i],split='')[[1]][9],strsplit(MyFullData[i],split='')[[1]][10])
    
    pre <-0}
  else if((strsplit(MyFullData[i],split='')[[1]][1]==">") && (strsplit(MyFullData[i],split='')[[1]][11]=="-")){
    SeqID <- SeqID+1
    MySeqName[SeqID] <- paste0(strsplit(MyFullData[i],split='')[[1]][5],strsplit(MyFullData[i],split='')[[1]][6],
                               strsplit(MyFullData[i],split='')[[1]][7],strsplit(MyFullData[i],split='')[[1]][8],
                               strsplit(MyFullData[i],split='')[[1]][9],strsplit(MyFullData[i],split='')[[1]][10],
                               strsplit(MyFullData[i],split='')[[1]][11],strsplit(MyFullData[i],split='')[[1]][12])
    
    pre <-0}else{#merging the lines of each sequence into a single sequence
      MyTempSequence[SeqID] <- paste0(MyTempSequence[SeqID], MyFullData[i])
      for(j in 1:length(strsplit(MyFullData[i],split = '')[[1]])){
        if(strsplit(MyFullData[i],split='')[[1]][j]=="K"){#caution: capital K
          Index[l] <- pre+j
          MyRepeatName[l] <- MySeqName[SeqID]
          l <- l+1}}
      pre <- pre+j
    }
}
#Extracting the labels (verified)
MyUniqCodeExtrac <- unique(MyRepeatName)
MyUniqCodeOriginal <- unique(MyOriginalLabels$Codes)



for(i in 1:1602){
  MyOriginalIndex <-which(MyOriginalLabels$Codes==MyUniqCodeExtrac[i])
  MyRepeatIndex <- which(MyRepeatName==MyUniqCodeExtrac[i])
  for(j in 1:length(MyOriginalIndex)){
    for(k in 1:length(MyRepeatIndex))
      if(MyOriginalLabels$Site[MyOriginalIndex[j]]==Index[MyRepeatIndex[k]]){
        MyVerifiedLabels$Site[MyOriginalIndex[j]] <- Index[MyRepeatIndex[k]]
      }
  }
}

#Extracting the lengths of each sequence
Length <- vector(mode = 'numeric' , length = 1602)
for(i in 1:1602)
  Length[i] = length(strsplit(MyTempSequence[i],split = '')[[1]])


#Extracting the feature values corresponding to PSSM features and physiochemistry properties:
#building empty feature matrix:
SubSeqLen <- 5 #the window sizes of 5, 7, 9,11,13,15,17,.... can be used
NrSelPhyChem <- 16 #an integer number between 0 and 16, 0 when no features than PSSM is involved.
Temp <- vector(mode = "numeric", length = (20*SubSeqLen+NrSelPhyChem*SubSeqLen+20+400+1) ) #PSSM + Phisiochemistry properties and another 20+400 is for sequencing features


FeatureMatrix <- data.frame(t(Temp))
FeatureNames <- paste0("F",c(1:(20*SubSeqLen+NrSelPhyChem*SubSeqLen+20+400))) #20*SubSeqLen is for PSSM and the rest for phisiochemis features and another 20+400 is for sequencing features
FeatureNames <- c(FeatureNames , "Label")
names(FeatureMatrix) <- FeatureNames
FeatureMatrix[,(20*SubSeqLen+NrSelPhyChem*SubSeqLen+1):(20*SubSeqLen+NrSelPhyChem*SubSeqLen+20)] <- 0
FeatureMatrix[,(20*SubSeqLen+NrSelPhyChem*SubSeqLen+21):(20*SubSeqLen+NrSelPhyChem*SubSeqLen+21+400)] <- 0
Counter <- 0
#Physiochemistry feature values:
MyPCMFeatures <- as.matrix(read_excel(PATH2PHYSCHEIMFTR, col_names = FALSE))
AminoAcidDic <- c("A", "R", "N", "D",  "C",  "Q", "E",  "G", "H", "I", "L", "K", "M", "F","P", "S", "T", "W", "Y", "V")

#Building PSSM feature values, PCP, and amino acid frequency per each sequence of associated size 
for(PathCount in 1:4) {
  FileName <- PATH2NAMES[PathCount]
  MyClassNames <- read_excel(FileName, col_names = FALSE)
  MyListC <- list.files(PATH2POSITIONS[PathCount])
  for (i in 1:length(MyListC)) {
    #extracting positive positions
    PositiveIndex <-
      MyVerifiedLabels$Site[MyVerifiedLabels$Codes == MyClassNames[[1]][i]]
    MyFullPositions <-
      readLines(paste0(PATH2POSITIONS[PathCount], "/", "Class_", PathCount, "_", i, ".csv"))
    #extracting negative positions:
    TotalIndex <- c(1:Length[which(MySeqName == MyClassNames[[1]][i])])
    TotalIndex <- TotalIndex[-PositiveIndex]
    NegativeIndex <-
      sample(TotalIndex,
             length(PositiveIndex),
             replace = FALSE,
             prob = NULL)
    #extracting feature values from both negative and positive positions
    for (j in 1:length(PositiveIndex)) {
      #writing to files
      if (PositiveIndex[j] > (trunc(SubSeqLen/2)+4) &&
          (PositiveIndex[j] +4+trunc(SubSeqLen/2)) <= (length(MyFullPositions)-6)) {
        Counter <- Counter + 1
        FeatureMatrix[Counter,(20*SubSeqLen+NrSelPhyChem*SubSeqLen+1):(20*SubSeqLen+NrSelPhyChem*SubSeqLen+20)] <- 0
        FeatureMatrix[Counter,(20*SubSeqLen+NrSelPhyChem*SubSeqLen+21):(20*SubSeqLen+NrSelPhyChem*SubSeqLen+21+400)] <- 0
        for (k in 0:(SubSeqLen-1)) {
          if (k != (SubSeqLen-1)) {
          #Extracting the features corresponding to the specified length
          Temp <-
            strsplit (t(MyFullPositions[(PositiveIndex[j] + k- trunc(SubSeqLen/2)+4)]), split = ' ') #Feeding positive index
          Temp_ <-
            strsplit (t(MyFullPositions[(PositiveIndex[j] + k- trunc(SubSeqLen/2)+5)]), split = ' ') #Feeding positive index
          FeatureMatrix [Counter, (1 + 20 * k):(20 * (k + 1))] <-
            Temp[[1]][-which(Temp[[1]] == "")][3:22] #PSSM values
          for(NPSSM in 0:(SubSeqLen-1)){ #Non PSSM features (physiochemistry features)
            FeatureMatrix [Counter , (20*SubSeqLen+NPSSM*NrSelPhyChem+1):(20*SubSeqLen+(NPSSM+1)*NrSelPhyChem)]  <-
              as.numeric(MyPCMFeatures[which(MyPCMFeatures[,1]==Temp[[1]][-which(Temp[[1]] == "")][2]), 2:(NrSelPhyChem+1)])
          }#sequencing features
          SeqIndex <- which(AminoAcidDic==Temp[[1]][-which(Temp[[1]] == "")][2])
          FeatureMatrix [Counter , (20*SubSeqLen+NrSelPhyChem*SubSeqLen+SeqIndex)] <-
            FeatureMatrix [Counter , (20*SubSeqLen+NrSelPhyChem*SubSeqLen+SeqIndex)]+(1/SubSeqLen)
          SeqIndex_ <- which(AminoAcidDic==Temp[[1]][-which(Temp_[[1]] == "")][2])
          FeatureMatrix [Counter , (20*SubSeqLen+NrSelPhyChem*SubSeqLen+20+((SeqIndex-1)*20+SeqIndex_))] <-
            FeatureMatrix [Counter , (20*SubSeqLen+NrSelPhyChem*SubSeqLen+20+((SeqIndex-1)*20+SeqIndex_))]+(1/(SubSeqLen*SubSeqLen))
          }
          else {
            Temp <-
              strsplit (t(MyFullPositions[(PositiveIndex[j] + k- trunc(SubSeqLen/2)+4)]), split = ' ') #Feeding positive index
            FeatureMatrix [Counter, (1 + 20 * k):(20 * (k + 1))] <-
              Temp[[1]][-which(Temp[[1]] == "")][3:22] #PSSM values
            for(NPSSM in 0:(SubSeqLen-1)){ #Non PSSM features (physiochemistry features)
              FeatureMatrix [Counter , (20*SubSeqLen+NPSSM*NrSelPhyChem+1):(20*SubSeqLen+(NPSSM+1)*NrSelPhyChem)]  <-
                as.numeric(MyPCMFeatures[which(MyPCMFeatures[,1]==Temp[[1]][-which(Temp[[1]] == "")][2]), 2:(NrSelPhyChem+1)])
            }#sequencing features
            SeqIndex <- which(AminoAcidDic==Temp[[1]][-which(Temp[[1]] == "")][2])
            FeatureMatrix [Counter , (20*SubSeqLen+NrSelPhyChem*SubSeqLen+SeqIndex)] <-
              FeatureMatrix [Counter , (20*SubSeqLen+NrSelPhyChem*SubSeqLen+SeqIndex)]+(1/SubSeqLen)
            
          }
        }
        FeatureMatrix[Counter, (20*SubSeqLen+NrSelPhyChem*SubSeqLen+20+400+1)] <- 1
        paste0(MyClassNames[[1]][i],"__Positive") #"1" is a label for positive positions
      }
    }
    #Feeding negative index
    for (j in 1:length(NegativeIndex)) {
      #writing to files
      if (NegativeIndex[j] > (trunc(SubSeqLen/2)+4) &&
          (NegativeIndex[j] +4+ trunc(SubSeqLen/2)) < (length(MyFullPositions)-6 )) {
        Counter <- Counter + 1
        FeatureMatrix[Counter,(20*SubSeqLen+NrSelPhyChem*SubSeqLen+1):(20*SubSeqLen+NrSelPhyChem*SubSeqLen+20)] <- 0
        FeatureMatrix[Counter,(20*SubSeqLen+NrSelPhyChem*SubSeqLen+21):(20*SubSeqLen+NrSelPhyChem*SubSeqLen+21+400)] <- 0
        for (k in 0:(SubSeqLen-1)) {
          if (k != (SubSeqLen-1)) {
          #Extracting the features corresponding to the length
          Temp <-
            strsplit (t(MyFullPositions[(NegativeIndex[j] + k- trunc(SubSeqLen/2)+4)]), split = ' ') #Feeding positive index
          Temp_ <-
            strsplit (t(MyFullPositions[(NegativeIndex[j] + k- trunc(SubSeqLen/2)+5)]), split = ' ') #Feeding positive index
          FeatureMatrix [Counter, (1 + 20 * k):(20 * (k + 1))] <-
            Temp[[1]][-which(Temp[[1]] == "")][3:22]
          for(NPSSM in 0:(SubSeqLen-1)){ #Non PSSM features (physiochemistry features)
            FeatureMatrix [Counter , (20*SubSeqLen+NPSSM*NrSelPhyChem+1):(20*SubSeqLen+(NPSSM+1)*NrSelPhyChem)]  <-
              as.numeric(MyPCMFeatures[which(MyPCMFeatures[,1]==Temp[[1]][-which(Temp[[1]] == "")][2]), 2:(NrSelPhyChem+1)])
          }#sequencing features
          SeqIndex <- which(AminoAcidDic==Temp[[1]][-which(Temp[[1]] == "")][2])
          FeatureMatrix [Counter , (20*SubSeqLen+NrSelPhyChem*SubSeqLen+SeqIndex)] <-
            FeatureMatrix [Counter , (20*SubSeqLen+NrSelPhyChem*SubSeqLen+SeqIndex)]+(1/SubSeqLen)
          SeqIndex_ <- which(AminoAcidDic==Temp[[1]][-which(Temp[[1]] == "")][2])
          FeatureMatrix [Counter , (20*SubSeqLen+NrSelPhyChem*SubSeqLen+20+((SeqIndex-1)*20+SeqIndex_))] <-
            FeatureMatrix [Counter , (20*SubSeqLen+NrSelPhyChem*SubSeqLen+20+((SeqIndex-1)*20+SeqIndex_))]+(1/(SubSeqLen*SubSeqLen))
          }
          else {
            #Extracting the features corresponding to the length
            Temp <-
              strsplit (t(MyFullPositions[(NegativeIndex[j] + k- trunc(SubSeqLen/2)+4)]), split = ' ') #Feeding positive index
            FeatureMatrix [Counter, (1 + 20 * k):(20 * (k + 1))] <-
              Temp[[1]][-which(Temp[[1]] == "")][3:22]
            for(NPSSM in 0:(SubSeqLen-1)){ #Non PSSM features (physiochemistry features)
              FeatureMatrix [Counter , (20*SubSeqLen+NPSSM*NrSelPhyChem+1):(20*SubSeqLen+(NPSSM+1)*NrSelPhyChem)]  <-
                as.numeric(MyPCMFeatures[which(MyPCMFeatures[,1]==Temp[[1]][-which(Temp[[1]] == "")][2]), 2:(NrSelPhyChem+1)])
            }
            SeqIndex <- which(AminoAcidDic==Temp[[1]][-which(Temp[[1]] == "")][2])
            FeatureMatrix [Counter , (20*SubSeqLen+NrSelPhyChem*SubSeqLen+SeqIndex)] <-
              FeatureMatrix [Counter , (20*SubSeqLen+NrSelPhyChem*SubSeqLen+SeqIndex)]+(1/SubSeqLen)
            }
        }
        FeatureMatrix[Counter, (20*SubSeqLen+NrSelPhyChem*SubSeqLen+20+400+1)] <- 0
        paste0(MyClassNames[[1]][i],"__Negative") #"0" is a label for positive positions
      }
    } #uncomment the following to save the results in a .csv file
    #FileName <- paste0(P1,MyClassNames[[1]][i], "_PositiveFeature.csv")
    #write.table(t(MyFullPositions[(PositiveIndex[j]+3)]), file=FileName, sep="     ", col.names = FALSE, append=TRUE)
  }
}
#writing into a file:
##write.table(FeatureMatrix,file ="FeatureDataLen_35.csv" , sep = "  ", col.names = FALSE )

PositiveFeatureValuesIndexes <-which(FeatureMatrix$Label==1)
PositiveFeatureValues <- FeatureMatrix[PositiveFeatureValuesIndexes,]
#spliting the set into train and test sets


#converting to numerical values
FeatureMatrix <- as.data.frame(sapply(FeatureMatrix,as.numeric))
MySample <- sample.split(FeatureMatrix$Label , SplitRatio=0.8)
trainingSet <- subset(FeatureMatrix,MySample==TRUE)#because of cross validation #subset(FeatureMatrix,MySample==TRUE)
trainingSet <- FeatureMatrix #in case of cross validation
testSet <- subset(FeatureMatrix,MySample==FALSE)
trainingSet$Label <- as.factor(trainingSet$Label)
testSet$Label <- as.factor(testSet$Label)
#normalizing the extracted feature set
normal <-function(x) { ifelse(x==0,0, (x -mean(x))/(sd(x)))}
trainingSet[,-(20*SubSeqLen+NrSelPhyChem*SubSeqLen+20+400+1)] <- lapply(trainingSet[,-(20*SubSeqLen+NrSelPhyChem*SubSeqLen+20+400+1)],normal)
testSet[,-(20*SubSeqLen+NrSelPhyChem*SubSeqLen+20+400+1)] <- lapply(testSet[,-(20*SubSeqLen+NrSelPhyChem*SubSeqLen+20+400+1)],normal)


##cross fold validation model in DNN:
folds = createFolds(trainingSet$Label, k = 10)

cv = lapply(folds, function(x) { # start of function
  # in the next two lines we will separate the Training set into it's 10 pieces
  trainingFold = trainingSet[-x, ] # training fold =  training set minus (-) it's sub test fold
  testFold = trainingSet[x, ] # here we describe the test fold individually
  # now apply (train) the classier on the training_fold
  MyDeepModel <- keras_model_sequential() 
  MyDeepModel %>% 
    layer_dense(units = 128, activation = "relu", input_shape = c(600)) %>%
    layer_batch_normalization() %>%
    layer_dropout(0.5) %>%
    layer_dense(units = 64, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(0.2) %>%
    layer_dense(units = 32, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(0.2) %>%
    layer_dense(units = 16, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(0.2) %>%
    layer_dense(units = 32, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(0.2) %>%
    layer_dense(units = 64, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(0.5) %>%
    layer_dense(units = 128, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dense(units = 2, activation = "softmax") %>%
    compile(
      loss = 'mse', 
      optimizer = optimizer_sgd(lr = 0.01, decay = 1e-6, momentum = 0.94, nesterov = TRUE),
      metrics = c('accuracy')
    )
  
  train_x <- matrix(unlist(trainingFold[,-dim(trainingSet)[2]]), nrow = dim(trainingFold)[1], ncol=dim(trainingFold)[2]-1)
  train_y <- matrix(as.numeric(unlist(trainingFold[,dim(trainingFold)[2]])), nrow = dim(trainingFold)[1], ncol=1)
  train_y <- train_y-1
  train_y <- to_categorical(train_y)
  
  test_x <- matrix(unlist(testFold[,-dim(testFold)[2]]), nrow = dim(testFold)[1], ncol=dim(testFold)[2]-1)
  test_y <- matrix(as.numeric(unlist(testFold[,dim(testFold)[2]])), nrow = dim(testFold)[1], ncol=1)
  test_y <- test_y-1
  test_y <- to_categorical(test_y)
  
  PathModel <- "C:/Users/Administrator/Documents/UseCases/ProteinsSequences/PSSM/Class_1/Class_1/Class_1/pssm/my_best_model.hd5"
  CallModel <- callback_model_checkpoint(
    PathModel,
    monitor = "val_accuracy",
    verbose = 0,
    save_best_only = TRUE,
    save_weights_only = TRUE,
    mode = "max",
    save_freq = "epoch"
  )
  
  MyKerasModel <- MyDeepModel %>% fit( train_x,train_y ,batch_size=16,
                                       epochs=70, callbacks = list(CallModel),verbose=0)
  print("###")
  print(max(MyKerasModel$metrics$val_accuracy))
  print("###")
  Predict <-predict(MyDeepModel, test_x, probability = TRUE, type="prob")
  CM<-confusionMatrix(as.factor(apply(Predict, 1, which.max)-1), as.factor(apply(test_y,1,which.max)-1))
  print(CM)
  ROCit_obj <- rocit(score=Predict[,2],class=as.numeric(apply(test_y,1,which.max)-1))
  plot(ROCit_obj, type="l", lwd=6,lty = c(4,4),col=c(2,4), YIndex = FALSE, values = TRUE, grid=FALSE)
  PR<-pr.curve(Predict[,1],Predict[,2], curve = TRUE, rand.compute = TRUE)
  plot(PR$curve[,1],PR$curve[,2], col=2,lwd=2, xlim=c(0,1), ylim=c(0.5,1), auc.main=FALSE, rand.plot = TRUE, type = 'l', ylab = 'Precision', xlab = 'Recall')
  lines(x = seq(0,1,by=0.1), y=rep(0.5,11),type = 'l',col=4,lwd=2)
  plot(PR$curve[,1],PR$curve[,2], col=2,lwd=2, xlim=c(0,1), ylim=c(0.0,1), auc.main=FALSE, rand.plot = TRUE, type = 'l', ylab = 'Precision', xlab = 'Recall')
  lines(x = seq(0,1,by=0.1), y=rep(0.5,11),type = 'l',col=4,lwd=2)
  AUC <- ciAUC(ROCit_obj, level = 0.9)
  print(AUC)
  return(CM)
  
})

##################################################################################
#SVM model
svmfit = svm(Label ~ ., data = trainingSet, scale = TRUE, probability=TRUE )

Predict <-predict(svmfit, testSet, probability = TRUE)
confusionMatrix(Predict, testSet$Label)

###Neural Network
NN_Model<-nnet(as.numeric(trainingSet$Label)-1 ~.,data=trainingSet ,size=5, decay=9e-5, maxit=400,MaxNWts=19000)

Predict <-predict(NN_Model, testSet)
confusionMatrix(as.factor(as.numeric(Predict > 0.5)), testSet$Label)

Predict <- as.numeric(Predict > 0.5)
Pluspredict <- Predict[which(Predict>0.5)]
PlusTest <- testSet[-which(testSet$Label==1),]


###############

#KNN and SVM cross validation model 


folds = createFolds(trainingSet$Label, k = 10)

cv = lapply(folds, function(x) { # start of function
  # in the next two lines we will separate the Training set into it's 10 pieces
  trainingFold = trainingSet[-x, ] # training fold =  training set minus (-) it's sub test fold
  testFold = trainingSet[x, ] # here we describe the test fold individually
  # now apply (train) the classier on the training_fold
  classifier = svm(formula = Label ~ .,
                   data = trainingFold,
                   type = 'C-classification',
                   kernel = 'radial')
  #classifier = knn(trainingFold,testFold,cl=trainingFold[,141],k=13)
  # next step in the loop, we calculate the predictions and cm and we equate the accuracy
  # note we are training on training_fold and testing its accuracy on the test_fold
  y_pred = predict(classifier, newdata = testFold[-(20*SubSeqLen+NrSelPhyChem*SubSeqLen+20+400+1)])
  cm = table(testFold[, (20*SubSeqLen+NrSelPhyChem*SubSeqLen+20+400+1)], y_pred)
  accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})


Predict <-predict(svmfit, testSet)
confusionMatrix(Predict, testSet$Label)
ROCit_obj <- rocit(score=as.numeric(Predict)-1,class=testSet$Label)
plot(ROCit_obj, type="l", lwd=6,lty = c(4,4),col=c(2,4), YIndex = FALSE, values = TRUE, grid=FALSE)
PR<-pr.curve(testSet$Label ,as.numeric(Predict)-1, curve = TRUE, rand.compute = TRUE)
plot(PR,col=2, xlim=c(0,1), ylim=c(0.5,1), auc.main=TRUE, rand.plot = TRUE)
AUC <- ciAUC(ROCit_obj, level = 0.9)
print(AUC)
###Boosting methods
#pure XGB models
bstDense <- xgboost(data = as.matrix(trainingSet[,-(20*SubSeqLen+NrSelPhyChem*SubSeqLen+20+400+1)]), label = as.numeric(as.character(trainingSet$Label)), max.depth =16, eta = 0.007, nthread = 4,nrounds = 1200, objective = "binary:logistic", gamma=0.4,subsample=0.5,tree_method="hist", booster = "gbtree", eval_metric="error")
Predict <- predict(bstDense, as.matrix(testSet[,-(20*SubSeqLen+NrSelPhyChem*SubSeqLen+20+400+1)]))
Prediction <- as.numeric(Predict > 0.5)
confusionMatrix(as.factor(Prediction) , testSet$Label)
xgb_imp <- xgb.importance( model = bstDense)
(FeaturePlot <- xgb.ggplot.importance(xgb_imp, measure = "Frequency", rel_to_first = FALSE, top_n=90))
FeaturePlot+ggplot2::ylab("Frequency")
sum(FeaturePlot$data$Frequency)




dtrain <- xgb.DMatrix(data = as.matrix(trainingSet[,-(20*SubSeqLen+NrSelPhyChem*SubSeqLen+20+400+1)]), label=as.numeric(as.character(trainingSet$Label)))
dtest <- xgb.DMatrix(data = as.matrix(testSet[,-(20*SubSeqLen+NrSelPhyChem*SubSeqLen+20+400+1)]), label=as.numeric(as.character(testSet$Label)))
watchlist <- list(train=dtrain, test=dtest)


bst <- xgb.cv(data=dtrain,nfold = 10,  max.depth=16, eta=0.4, nthread = 4,gamma=0.1, watchlist=watchlist, booster = "gbtree",objective = "binary:logistic", eval_metric=c(  "error") ,nrounds = 220, early_stopping_rounds = 160)
xgb.save(bst, 'Myxgb.model')
bst <- xgb.load('Myxgb.model')

Predict <- predict(bst, as.matrix(testSet[,-(20*SubSeqLen+NrSelPhyChem*SubSeqLen+20+1)]))
Prediction <- as.numeric(Predict > 0.5)
ROCit_obj <- rocit(score=Prediction,class=testSet$Label)
plot(ROCit_obj)
Prediction <- as.factor(as.numeric(Predict > 0.5))
ConfusionMatrix <- confusionMatrix(Prediction, testSet$Label)
ConfusionMatrix
###############manual cross validation:
#XGB cross validation model in a cross validation settings
folds = createFolds(trainingSet$Label, k = 10)

cv = lapply(folds, function(x) { # start of function
  # in the next two lines we will separate the Training set into it's 10 pieces
  trainingFold = trainingSet[-x, ] # training fold =  training set minus (-) it's sub test fold
  testFold = trainingSet[x, ] # here we describe the test fold individually
  # now apply (train) the classier on the training_fold
  dtrain <- xgb.DMatrix(data = as.matrix(trainingFold[,-(20*SubSeqLen+NrSelPhyChem*SubSeqLen+20+400+1)]), label=as.numeric(as.character(trainingFold$Label)))
  dtest <- xgb.DMatrix(data = as.matrix(testFold[,-(20*SubSeqLen+NrSelPhyChem*SubSeqLen+20+400+1)]), label=as.numeric(as.character(testFold$Label)))
  watchlist <- list(train=dtrain, test=dtest)
  classifier = xgb.train(data=dtrain, 
                         max.depth=17, eta=0.007, 
                         nthread = 8, 
                         gamma=0.3,
                         lambda=1.3,
                         subsample=0.5,
                         tree_method="hist",
                         #watchlist=watchlist, 
                         objective = "binary:logistic",
                         eval_metric=c(  "error") ,
                         nrounds = 1200, 
                         verbose = 0)
  # next step in the loop, we calculate the predictions and cm and we equate the accuracy
  # note we are training on training_fold and testing its accuracy on the test_fold
  y_pred = predict(classifier, newdata = dtest)
  cm = table(testFold$Label, as.numeric(y_pred>0.5))
  CM <-confusionMatrix(testFold$Label, as.factor(as.numeric(y_pred>0.5)))
  print(CM)
  #accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  ROCit_obj <- rocit(score=y_pred,class=testFold$Label)
  plot(ROCit_obj, type="l", lwd=6,lty = c(4,4),col=c(2,4), YIndex = FALSE, values = TRUE, grid=FALSE)
  PR<-pr.curve(testFold$Label ,y_pred, curve = TRUE, rand.compute = TRUE)
  plot(rev(PR$curve[,3]),PR$curve[,2],lwd=2,col=2, xlim=c(0,1), ylim=c(0.5,1), auc.main=TRUE, rand.plot = TRUE, xlab='Recall',ylab='Precision',type='l')
  lines(x = seq(0,1,by=0.1), y=rep(0.5,11),type = 'l',col=4)
  AUC <- ciAUC(ROCit_obj, level = 0.9)
  print(AUC)
  plot(ROCit_obj)
  return(CM)
})
#calculating the statistics from cross fold model
MeanAcc <- 0
MeanSen <- 0
MeanSp <- 0
MeanKp <- 0
AccMax <- 0
SenMax <- 0
SPMax <- 0
KPMax <- 0
for (i in 1:length(folds)) {
  Te <- (cv[[i]]$table[1,1]+cv[[i]]$table[2,2])/(cv[[i]]$table[1,1]+cv[[i]]$table[2,2]+cv[[i]]$table[2,1]+cv[[i]]$table[1,2])
  ifelse(AccMax<Te, AccMax <- Te, AccMax <- AccMax ) 
  MeanAcc <- MeanAcc+ Te
  Te <- cv[[i]]$byClass[1]
  ifelse(SenMax<Te, SenMax <- Te, SenMax <- SenMax ) 
  MeanSen <- MeanSen+ Te
  Te <- cv[[i]]$byClass[2]
  ifelse(SPMax<Te, SPMax <- Te, SPMax <- SPMax ) 
  MeanSp <- MeanSp+ Te
  Te <- cv[[i]]$overall[2]
  ifelse(KPMax<Te, KPMax <- Te, KPMax <- KPMax ) 
  MeanKp <- MeanKp+ Te
}
MeanAcc <- MeanAcc/length(folds)
MeanSen <- MeanSen/length(folds)
MeanSp <- MeanSp/length(folds)
MeanKp <- MeanKp/length(folds)
print(MeanAcc)
print(MeanSen)
print(MeanSp)
print(MeanKp)
print(AccMax)
print(SenMax)
print(SPMax)
print(KPMax)
#########################################
#pure Random forest model which can be used directly in a cross validation format, set the lower and upper bound of the tree numbers
for(NrTree in 80:100){#To check with different sizes of tree
  fitModelRF<- randomForest(Label ~., data = trainingSet , ntree=NrTree, proximity=TRUE)
  #fitModelRF<-train(Label ~., data = trainingSet,method="rf", proximity=TRUE)
  Predict <-predict(fitModelRF, testSet,type="prob")
  print(confusionMatrix(as.factor(as.numeric(Predict[,2]>0.5)), testSet$Label))
  ROCit_obj <- rocit(score=as.numeric(Predict[,2]),class=testSet$Label)
  plot(ROCit_obj, type="l", lwd=6,lty = c(4,4),col=c(2,4), YIndex = FALSE, values = TRUE, grid=FALSE)
  PR<-pr.curve(testSet$Label ,Predict[,2], curve = TRUE, rand.compute = TRUE)
  plot(PR$curve[,3],rev(PR$curve[,2]),col=2,lwd=2, xlim=c(0,1), ylim=c(0.5,1), auc.main=TRUE, rand.plot = TRUE, type = 'l', ylab = 'Precision', xlab = 'Recall')
  lines(x = seq(0,1,by=0.1), y=rep(0.5,11),type = 'l',col=4)
  AUC <- ciAUC(ROCit_obj, level = 0.9)
  print(AUC)
  plot(ROCit_obj)
  print(NrTree)}


