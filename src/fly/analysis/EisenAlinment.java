/**
 * 
 */
package fly.analysis;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Iterator;

import org.biojava.bio.BioException;
import org.biojava.bio.program.gff.GFFEntrySet;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.GFFTools;
import org.biojava.utils.ParserException;

/**
 * @author hk3
 *
 */
public class EisenAlinment {

	/**
	 * @param args 
	 * 07-09-2011 Hashem Koohy
	 * 
	 * this is to check what fraction of the fly genome  is covered by the Eisen Alignment blocks
	 */
	
	
	private String fileName = null;
	private String outputFileName = null;
	private String aGffFile       = null;
	
	
	
	private void setFileName(String fn){
		this.fileName = fn;
	}/*setFileName*/
	
	private void setOutputFileName(String ofn){
		this.outputFileName = ofn;
	}/*setOutputFileName*/
	
	private void setAGffFile(String agff){
		this.aGffFile= agff;
	}/*setAGffFile*/
	
	private String getFileName(){
		return fileName;
	}/*getFileName*/
	
	private String getOutputFileName(){
		return outputFileName;	
	}/*getOutputFileName*/
	
	private String getAGffFile(){
		return aGffFile;
	}/*getAGffFile*/
	
	
	public static void main(String[] args) throws FileNotFoundException, ParserException, BioException, IOException {
		// TODO Auto-generated method stub
	
		EisenAlinment app = new EisenAlinment();
		
		 app.setFileName("/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ALIGNMENTS_EISENLAB/X_List_of_Alignment_Files.txt");
		 app.setOutputFileName("/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ALIGNMENTS_EISENLAB/ALIGNMENT_WISENLAB_STATISTICS/X_lengths_of_alignment_blocks.txt");
		 app.printoutThisBlockInformation(app.getFileName(), app.getOutputFileName());
		 
		 
		 app.setAGffFile("/Users/hk3/Desktop/Main/Composure_Droshophila_Model/Unionof_DnaseI_CodingRegions_and_RedFlyEnhancers/UnionOfbdtnpDnaseAccAllStages.gff");
		// app.printoutDnasIInformation(app.getAGffFile(), "4");
		
	}/*main*/
	
	
	
	public void printoutThisBlockInformation( String inputFile, String outputFile){
		try{
			//input file initiation
			FileInputStream fstream = new FileInputStream(inputFile);
			  DataInputStream in = new DataInputStream(fstream);
			  BufferedReader br = new BufferedReader(new InputStreamReader(in));
			  
			  //output file initiation
			  FileWriter fw = new FileWriter(outputFile);
			  
			  
			  String strLine;
			  int numberOfAlignments = 0;
			  int overLengthofAlignments = 0;
			  while ((strLine = br.readLine()) != null)   {
				  numberOfAlignments++;
				  
				  String alignmetSpecifications [] = strLine.split("_");
				  String Chr = alignmetSpecifications[0];
				  int start = Integer.valueOf(alignmetSpecifications[1]);
				  int end   = Integer.valueOf(alignmetSpecifications[2]);
				  int oneLength = end-start;
				  fw.write(String.valueOf(oneLength));
				  fw.write(",");
				  if(oneLength <0 ){
					  System.out.println("line " + numberOfAlignments + " looks abnormal!");
				  }
				  overLengthofAlignments = overLengthofAlignments+oneLength;
				  // Print the content on the console
				 // System.out.println (strLine);
				  
				  }
			  System.out.println("Number of Alignment Blocks in over this chromosome is: " + numberOfAlignments);
			  System.out.println("Over length of alingment blocks over this chro is " + overLengthofAlignments);
				  //Close the input stream
			  	fw.flush();
			  	fw.close();
				  in.close();
		}catch(Exception e){
			System.err.println("Error: " + e.getMessage());
		}
		
		
	}/*printoutThisBlockInformation*/
	
	
	public void printoutDnasIInformation(String DnasIinputFileName, String chr) throws FileNotFoundException, ParserException, BioException, IOException{//note: it is assumed that the input file is a gff file which includes union of DNasI records
		File DnaseInputFile = new File(DnasIinputFileName);
		GFFEntrySet  gff    = GFFTools.readGFF(DnaseInputFile);
		
		String outputFilePath = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ALIGNMENTS_EISENLAB/ALIGNMENT_WISENLAB_STATISTICS/";
		String outFileName    = outputFilePath+chr+"_UnionOfDNasI.txt";
		
		FileWriter fw  = new FileWriter(outFileName);
		
		int recordCounter =0;
		int sumOfRecordLengths = 0;
		for(Iterator<?>  gffi = gff.lineIterator(); gffi.hasNext();){
			Object o = gffi.next();
			if(o instanceof GFFRecord){
				GFFRecord gffrec = (GFFRecord) o;
				String thisChro = gffrec.getSeqName();
				if (thisChro.equals(chr)){
					recordCounter++;
					int oneRecordLength = gffrec.getEnd()- gffrec.getStart();
					if(oneRecordLength<0 ){
						System.out.println("this record sounds crazy! because has negative length! " + gffrec.getSeqName() + "   " + gffrec.getStart() + "  " + gffrec.getEnd());
					}
					fw.write(String.valueOf(oneRecordLength));				
					fw.write(",");				
					sumOfRecordLengths= sumOfRecordLengths+ oneRecordLength;
				}
			}
		}
		fw.flush();
		System.out.println("number of records is " + recordCounter + " and total length of records is " + sumOfRecordLengths);
	}/*printoutDnasIInformation*/
	
	

}/*EisenAlinment*/
