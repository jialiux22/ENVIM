# ENVIM

ENVIM (elastic net variable importance model) is a software to predict metabolite abundance by using microbial features (metagenomics data or metatranscriptomics data) via ENVIM. There are two function in the folder: ENVIM and ENVIM_predict.

For using both function, you need to install following packages in R in order to run the code successfully:
1.glmnet
2.MASS
3.caret
4.GenABEL
5.caTools
6.optparse

ENVIM function is to use microbial training data and metabolite training data to predict metabolite testing data and compare the result between observed metabolite test data and predicted metabolite testing data.

Input 1, 2, 3, and 4 should be normalized into relative abundance and in CSV file. The column names for all four files are microbial gene names (for files 1 and 2) metabolites name (for file 3 and 4). The row names are sample names. 

Sample files of inputs and outputs are provided into the same directory.

Inputs:
1. Microbial training data
2. Microbial testing data
3. Metabolite training data
4. Metabolite testing data
5. Seed number (default: 1234)
6. Number of fold for random forest model (default: 10)
7. Number of fold for ENVIM (default: 10)
8. Output directory (default: same directory as the code)

Output:

1.Summary table (including training spearman correlation, training mean square error, testing spearman correlation, testing mean square error, and running time in minutes)

2.Training weight matrix

3.Testing weight matrix

We can run the ENVIM in R or bash.

Run software on R:
We will run ENVIM.R file in R. 
Sample code: 
      source(paste0(getwd(),"/ENVIM.R"))

      ENVIM(microbio.train = microbio.train,

      microbio.test = microbio.test,
      
      metab.train = metab.train,
      
      metab.test = metab.test,
      
      seed = 1234,
      
      outputdirectory = out,
      
      fold_rf = 10,
      
      fold_ENVIM = 10)

Run software on bash:
We will run ENVIM_predict.R on bash.
Sample code: 

      Rscript ENVIM_predict.R -a micro.train.csv -b micro.test.csv -c metab.train.csv -d metab.test.csv -e 1234 -f 10 -g 10 -o your_output_directory


ENVIM_predict function is only to predict metabolite testing data by using the model of microbial training data and metabolite training data while we do not have metabolite testing data.

Input 1, 2, and 3 should be normalized into relative abundance and in CSV file. The column names for all four files are microbial gene names (for files 1 and 2) metabolites name (for file 3 and 4). The row names are sample names. 

Sample files of inputs and outputs are provided into the same directory.

Inputs:
1. Microbial training data
2. Microbial testing data
3. Metabolite training data
4. Seed number (default: 1234)
5. Number of fold for random forest model (default: 10)
6. Number of fold for ENVIM (default: 10)
7. Output directory (default: same directory as the code)

Output:
1. Summary table (including training spearman correlation, training mean square error, and running time in minutes)
2. Training weight matrix
3. Testing weight matrix
4. Predicting metabolite data

Sample code in R:
source(paste0(getwd(),"/ENVIM.R"))


ENVIM_predict(microbio.train = microbio.train,
      microbio.test = microbio.test,
      metab.train = metab.train,
      seed = 1234,
      outputdirectory = out,
      fold_rf = 10,
      fold_ENVIM = 10)
