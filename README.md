# TEDProject

This is our implementation for the paper:

TED+: Towards Discovering Top-𝑘 Edge-Diversified Patterns in a Graph Database

Ted is implemented with Java (JDK1.8).

# Environments

JDK 1.8

Eclipse 2019

# Dataset

1. AIDS antiviral dataset      https://wiki.nci.nih.gov/display/NCIDTPdata/AIDS+Antiviral+Screen+Data

2. PubChem dataset             ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/

3. eMolecule dataset           https://www.emolecules.com/info/plus/download-database

# Example to run the codes

We first illustrate how to run TED with an exapmle.  Then, examples of its variants (BASE, DSS, PRM, and TEDLite,  see Section 7.1) are presented.

Let the dataset be AIDS40K, inputfile "AIDS40k", outputfile "AIDS40k_result.txt".

1. TED 

Step 1: Parameter Setting in main/Arguments.java. 

        1） hyperparameter setting. 
        
             inFilePath  =  "AIDS40k";  outFilePath =  "AIDS40k_result.txt"; minSup  =  1;  swapcondition  = "swap1"; maxNodeNum = 11;
             
             numberofpatterns = 5;      numberofgraphs   = 40000;  isPESIndex  = true;
             
        2)  parameter setting for optimization strategies
        
             Since all optimization strategies are used by TED,  we set hasPRM = true;  hasDSS = true;  hasInitialPatternGenerator  = true
             
             In addition, we set isLightVersion = false.
             
Step 2: Run the main class, main/TEDSMain.java. The results can be found in the outputfile.

2. TED's variants

Step 1:   Unless specified otherwise, the parameter settings are the same as above.

         1) PRM:  Set hasInitialPatternGenerator  = false; 2) DSS:  Set hasInitialPatternGenerator  = false and hasPRM = false; 
         
         3) BASE: Set hasInitialPatternGenerator  = false, hasPRM = false, and hasDSS = false;
         
         4) TEDLite:  Set isLightVersion = true.

Step 2:  Run the main class, main/TEDSMain.java. The results can be found in the outputfile.  
