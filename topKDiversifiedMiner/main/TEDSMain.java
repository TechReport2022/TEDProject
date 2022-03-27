package main;

//import org.apache.commons.cli.*;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

public class TEDSMain {
    public static void main(String[] args) throws IOException {
    	Long Time1 = System.currentTimeMillis();	
        Arguments arguments = Arguments.getInstance(args);
        File inFile = new File(arguments.inFilePath);
        File outFile = new File(arguments.outFilePath);
        try (FileReader reader = new FileReader(inFile)) {
            try (FileWriter writer = new FileWriter(outFile)) {
            	if(arguments.isSimplified) {
            		TEDSSimplifiedProcessor processor = new TEDSSimplifiedProcessor();
                 	System.out.println("Start TEDSSimplifiedProcessor...");
                 	processor.run(reader, writer, arguments);
                 	System.out.println("Finished TEDSSimplifiedProcessor...");
            		
            	}else {
                	if(arguments.isLightVersion) {
                		TEDSLightProcessor processor = new TEDSLightProcessor();
                     	System.out.println("Start TEDSLightProcessor...");
                     	processor.run(reader, writer, arguments);
                     	System.out.println("Finished TEDSLightProcessor...");
                	}else {
                		TEDSProcessor processor = new TEDSProcessor();
                     	System.out.println("Start TEDSProcessor...");
                     	processor.run(reader, writer, arguments);
                     	System.out.println("Finished TEDSProcessor...");
                	}
            	}
            	
            }
        }
        Long Time2 = System.currentTimeMillis();
		System.out.println("Running Time(s) : " + (Time2 - Time1)*1.0 /1000);
    }
}
