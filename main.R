library(shiny)
library(bslib)
library(shinyjs)
library(shinybusy)
library(knitr)
library(rmarkdown)
library(here)
library(tools)
library(spsComps)
library(nestedcv)
library(caret)
library(shapviz)
library(probably)
library(ggplot2)
library(tidyverse)
library(kableExtra)
library(Boruta)
library(nnet)
library(randomForest)

ui <- fluidPage(
  titlePanel("Food authentication using foodomics - Metabolomics & Machine Learning"),
  sidebarPanel(
    radioButtons("method", "Choose Method:",c("Random Forest"="rf","Linear Discriminant Analysis"="lda","Neural Network"="nnet")),
    radioButtons("preproc", "Preprocessing Method:",c("Center and Scale"="cs","None"="none")),
    radioButtons("filter", "Perform Feature Reduction:",c("Yes"="yes","No"="no")),
    radioButtons("mpred", "Make Predictions on Test Data:",c("Yes"="yes","No"="no")),
    actionButton("predict","Predict Quality"),
  use_busy_spinner(spin="fading-circle")),
  mainPanel(
    h4("Predict Olive Oil Quality"),
    h5(em("Summary")),
    p("This demo app (written in R) uses metabolomics data from mass spectrometry to predict if any olive oil sample is adulterated with vegetable oil or not. The training dataset comprises of untargeted MS data from pure (un-adulterated) extra-virgin olive oil samples and olive oil samples adulterated with vegetable oils. The machine learning algorithm (for classification) and whether to use feature reduction and preprocessing methods are available choices. Additionally, this app can be used to make predictions on example test dataset (one adulterated and one pure olive oil sample)."),
    h5(em("Details")),
    p("Model evaluation is performed using nested cross validation with the inner loop for hyperparameter tuning (number of inner cv folds = 5) and outer loop for estimating performance of the model (number of outer cv folds = 10). Feature selection uses Boruta algorithm and SHAP (SHapley Additive ExPlanation) method is used for model interpretation and explain predictions on test data. Prediction summary includes predicted class and corresponding class probabilities. SHAP force plot is used to explain how each feature affected the prediction outcome for test data."),
    h5(em("R Packages used")),
    p("The app uses these R packages: shiny, MALDIquant, nestedcv, caret, shapviz, probably, ggplot2, tidyverse, kableExtra, rsconnect"),
    navset_card_tab(
    #   # Tab 1
    nav_panel("Model Summary",br(),tableOutput("modelsum")),
    #   # Tab 2
    nav_panel("Model Plots",br(),plotOutput("plot"),uiOutput("drocpr")),
    #   # Tab 3
    nav_panel("Prediction",br(),textOutput("filler"),tableOutput("pred")),
    #   # Tab 4
    nav_panel("Explainer",br(),textOutput("sfiller"),plotOutput("splot"),br(),uiOutput("dshap"))
  )
  )
)
server <- function(input,output) {
  observeEvent(input$predict,{
  show_spinner()
  # Set up CV parameters
  icv.fold = 5 # inner cv folds
  ocv.fold = 10 # outer cv folds
  
  # Read training data
  train.data = read.csv("Train_Adult_normalized.csv",stringsAsFactors = T)
  
  # Preprocess training data (centering and scaling) 
  train.data.proc = predict(preProcess(train.data[,-c(1,2)]),train.data[,-c(1,2)],method=c("center","scale"))
  
  # Set training parameters for inner loop (hyper-parameter tuning)
  # Set summaryFunction=defaultSummary for accuracy mteric and summaryFunction=mnLogLoss for minimum log loss
  train_control = trainControl(method="cv",number=icv.fold,classProbs=TRUE,summaryFunction = mnLogLoss,savePredictions="all")
  
  # Custom function to preprocess data
  preproc <- function(x,y){caret::preProcess(x,method=c("center","scale"))}
  
  # Nested cross validation (with boruta filter)
  if(input$filter == "yes")
  {
    if (input$preproc == "cs"){
    nested.cv = nestedcv::nestcv.train(y=train.data$Label,x=train.data[,-c(1,2)],method = input$method,filterFUN = boruta_filter,modifyX=preproc,modifyX_useY=TRUE,outer_method = "cv",n_outer_folds = ocv.fold,trControl = train_control,finalCV = T)
    }
    else{
      nested.cv = nestedcv::nestcv.train(y=train.data$Label,x=train.data[,-c(1,2)],method = input$method,filterFUN = boruta_filter,outer_method = "cv",n_outer_folds = ocv.fold,trControl = train_control,finalCV = T)
    }
    }
  else {
    if (input$preproc == "cs"){
    nested.cv = nestedcv::nestcv.train(y=train.data$Label,x=train.data[,-c(1,2)],method = input$method,modifyX=preproc,modifyX_useY=TRUE,outer_method = "cv",n_outer_folds = ocv.fold,trControl = train_control,finalCV = T)
    } else{
      nested.cv = nestedcv::nestcv.train(y=train.data$Label,x=train.data[,-c(1,2)],method = input$method,outer_method = "cv",n_outer_folds = ocv.fold,trControl = train_control,finalCV = T)  
    }
    }

  #Outer CV predictions (with class probabilities)
  outer.cv.pred = nested.cv$output
  # For pr curve
  pred.cal = as_tibble(outer.cv.pred) %>% mutate(.pred_class=predy,class=testy,.pred_Adulterated=1-predyp,.pred_Pure=predyp) %>% dplyr::select(.pred_class,class,.pred_Adulterated,.pred_Pure)
  # Outer CV results summary (Bal.Acc,Acc,AUC)
  outer.cv.summary = nested.cv$summary
  # Outer CV results confusion matrix and statistics
  outer.cv.cmat = confusionMatrix(data=outer.cv.pred$predy,reference=outer.cv.pred$testy,mode="everything")
  # Outer CV performance metrics
  outer.cv.met = outer.cv.cmat$byClass # sensitivity = outer.cv.met["Sensitivity"]
  metrics.all = round(append(outer.cv.met,outer.cv.summary$metrics["AUC"]),3)
  # Put performance metric in a table
  metrics.out = cbind.data.frame(names(metrics.all),value=metrics.all)
  rownames(metrics.out) = NULL
  colnames(metrics.out) = c("Performance Metrics","Value")
  roc.data = cbind.data.frame(x=1-nested.cv$roc$specificities,y=nested.cv$roc$sensitivities,t=nested.cv$roc$thresholds)
  
  #Put model summary in table
  mtable = kableExtra::kable(metrics.out,align='c') %>% kableExtra::kable_styling(bootstrap_options = c("striped","bordered"),full_width = F) %>% column_spec(1,"15em") %>% column_spec(2,"10em")
  
  hide_spinner()
  
   # Read new data (unseen data - usually customer data)
  new.data = read.csv("New_Adult_normalized.csv",stringsAsFactors = T)
  
  if(input$mpred == "yes"){
  # Prepare new data for predictions
  # Preprocessing transformation from training data is applied to new data
  new.data.proc = predict(preProcess(train.data[nested.cv$final_vars]),newdata=new.data[nested.cv$final_vars],method=c("center","scale"))
  
  # Make predictions on new data
  new.pred = predict(nested.cv,new.data.proc)
  new.pred.prob = predict(nested.cv,new.data.proc,type="prob")
  new.pred.final <- cbind.data.frame(Samples = new.data$Samples,.pred_class=as.factor(new.pred),new.pred.prob)
  colnames(new.pred.final) = c("Samples","Predicted Class","Class Probability (Adulterated)","Class Probability (Pure)")
  
  #Put prediction summary in table
  ptable = kableExtra::kable(new.pred.final,align='c') %>% column_spec(2,bold=T,color=ifelse(new.pred.final$Samples==new.pred.final$`Predicted Class`,"green","red")) %>% kableExtra::kable_styling(bootstrap_options = c("striped","bordered"),full_width = F) %>% column_spec(1,"15em") %>% column_spec(2,"10em")
  
  # SHAP values for new data (Explanation for individual predictions)
  shapvalues.new.data = fastshap::explain(object=nested.cv,X=train.data.proc,newdata = new.data.proc,pred_wrapper = pred_train,adjust=T,nsim=3)
  shap.plot = shapviz::shapviz(shapvalues.new.data,X=new.data.proc)
  
  # Predictions summary and plots
  output$filler <- NULL
  output$pred <- function(){ptable}
  splot1 <- shapviz::sv_force(shap.plot,row_id = 1) + ggtitle(paste("Sample:",new.pred.final$Samples[1],"|","Prediction:",new.pred[1]))
  splot2 <- shapviz::sv_force(shap.plot,row_id = 2) + ggtitle(paste("Sample:",new.pred.final$Samples[2],"|","Prediction:",new.pred[2]))
  pred.shap.plot = ggpubr::ggarrange(splot1,splot2,nrow=2)
  output$splot <- renderPlot({pred.shap.plot})
  
  output$dshap <- renderUI({
    downloadButton("downloadshap", "Download SHAP curves")})
  
  output$downloadshap <- downloadHandler(
    filename = "Prediction-SHAP.png",
    content = function(file) {
      ggsave(file,plot=pred.shap.plot)
    })
  
  } else {
    output$filler <- renderText({"Please set 'Make Predictions on test Data' to 'Yes' to see results here"})
    output$sfiller <- renderText({"Please set 'Make Predictions on test Data' to 'Yes' to see results here"})
    output$pred <- NULL
    output$splot <- NULL
    output$dshap <- NULL
  }
  
  # Prepare results for ui page
  # Model summary and plots
  output$modelsum <- function(){mtable}
  roc.plot = ggplot(roc.data,aes(x = x, y = y)) + geom_path() + geom_abline(lty = 3) + coord_equal() + theme_bw() + xlab("1-Specificity") + ylab("Sensitivity") + annotate(geom="text",x=0.75,y=0.1,label=paste("AUC:",round(nested.cv$roc$auc,3))) + ggtitle("ROC Curve")
  pr.plot = autoplot(yardstick::pr_curve(pred.cal,truth=class,.pred_Adulterated)) + ggtitle("Precision-Recall Curve")
  roc.pr.plot = ggpubr::ggarrange(roc.plot,pr.plot,nrow=1)
  output$plot <- renderPlot({roc.pr.plot})
  
  # Download model and prediction plots
   output$drocpr <- renderUI({
    downloadButton("downloadrocpr", "Download curves")})
   
   # Rshiny download tab
   output$dreport <- renderUI({
     downloadButton("downloadrep", "Download report")
   })
   
  # Download handler
  output$downloadrocpr <- downloadHandler(
    filename = "ROC-PR_curves.png",
    content = function(file) {
      ggsave(file,plot=roc.pr.plot)
    })
  
  })
  }
  
shinyApp(ui, server)