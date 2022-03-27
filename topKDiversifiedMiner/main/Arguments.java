package main;

import java.util.Scanner;

public class Arguments {
	  private static Arguments arguments;
      private String[] args;
      String inFilePath  =  "AIDS40k"; 
      String outFilePath =  "AIDS40k_result.txt";
     // String inFilePath  =  "pubchem23238"; 
     //  String inFilePath  =  "pubchem1000000clean.txt"; 
     // String outFilePath =  "pubchem23238_result.txt";
     // String inFilePath  =  "emolecul10000"; 
      //  String outFilePath =  "emolecul10000_result.txt";
      public long minSup        =  1;
      //long minNodeNum = 0;
      //long maxNodeNum = Long.MAX_VALUE;
      
      String  swapcondition  = "swap1";  //"swap1", "swap2", "swapalpha"
      Double  swapAlpha    = 0.99;
      long minNodeNum = 2;
      long maxNodeNum = 11;
      String strategy  = "topk";
       // String strategy  = "greedy";
      Integer numberofpatterns = 5;
      Integer numberofgraphs   = 40000;
      
      Boolean hasPRM = true;
      Boolean isPESIndex  = true;
      Boolean isSimpleIndex = !isPESIndex;
      
      
      Boolean hasDSS = true;
      Boolean hasInitialPatternGenerator  = true;
      
     // double maintainTime = 0;
      
      Boolean isLightVersion = false;
      Integer ReadNumInEachBatch = 100000;
      Double  AvgE = 43.783167;         ///For Pubchem
      Double  SampleEdgeNum = 500000.0; /// For Pubchem
      

      Boolean isSimplified = false;
      
      
      private Arguments(String[] args) {
          this.args = args;
      }

      static Arguments getInstance(String[] args) {
          arguments = new Arguments(args);
          if (args.length > 0) {
              arguments.initFromCmd();
          } else {
            //arguments.initFromRun();
          }
          
          
          return arguments;
      }
      private void initFromCmd() {

      }
      private void initFromRun() {
          try (Scanner sc = new Scanner(System.in)) {
          	
              System.out.println("Please input the file path of data set: ");
              inFilePath = sc.nextLine();
              
              System.out.println("Please set the minimum support: ");
              minSup = sc.nextLong();
              
              outFilePath = inFilePath + "_result.txt";
          }
      }
}
