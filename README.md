# Predict Olive Oil Quality
## Summary
Metabolomics data from mass spectrometry (MS) is combined with machine learning (ML) to predict if any olive oil sample is adulterated with vegetable oil or not. The training dataset comprises of untargeted MS data from pure (un-adulterated) extra-virgin olive oil samples and olive oil samples adulterated with vegetable oils. This Rshiny app evaluates ML models built on metabolomics training data and makes predictions on example test data. 
## Details
Model evaluation is performed using nested cross validation with the inner loop for hyperparameter tuning (number of inner cv folds = 5) and outer loop for estimating performance of the model (number of outer cv folds = 10). Feature selection uses Boruta algorithm and SHAP (SHapley Additive ExPlanation) method is used for model interpretation and explain predictions on test data. Prediction summary includes predicted class and corresponding class probabilities. SHAP force plot is used to explain how each feature affected the prediction outcome for test data.
## How to use Rshiny app?
The machine learning algorithm (for classification) and whether to use feature reduction and preprocessing methods are available choices. Additionally, this app can be used to make predictions on example test dataset (one adulterated and one pure olive oil sample). Please make appropriate selections and click 'Predict Quality'.
## App output
ML model performance, predictions on test data and explanations using SHAP plots are outputs from the app.
## App Access
Check the link for Rshiny app: [Predict Quality](https://nandhinidev.shinyapps.io/PredQual) 
