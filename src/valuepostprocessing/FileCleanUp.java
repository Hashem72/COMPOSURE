/**
 * 
 */
package valuepostprocessing;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.Iterator;

import org.biojava.bio.program.gff.GFFEntrySet;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.GFFTools;
import org.biojava.bio.program.gff.GFFWriter;

/**
 * @author hk3
 * 15-04-2011
 */
public class FileCleanUp {

	/**
	 * @param args
	 */
	
	private String inputGffFileName = null;
	private String outputGffFileName = null;
	
	//SETTERS
	public void setInputGffFileName(String ign){
		this.inputGffFileName = ign;
	}/*setInputGffFileName*/
	
	public void setOutputGffFileName(String ogf){
		this.outputGffFileName = ogf;
	}/*setOutputGffFileName*/
	
	//GETTERS
	public String getInputGffFileName(){
		return inputGffFileName;
	}/*getInputGffFileName*/
	
	public String getOutputGffFileName(){
		return outputGffFileName;
	}/*getOutputGffFileName*/
	
	
	//note: this is to print out only exons and cds from the input file but not starting and ending codons:
	public void removeStartingCodonRecordsFromAGffFile(String  inputFileName, String outputFileName) throws IOException, Exception{
		
			File gffinputFile  = new File(inputFileName);
			File gffoutputFile = new File (outputFileName);
			
	        PrintWriter pw = new PrintWriter(new FileWriter(gffoutputFile));
	        GFFWriter gffw = new GFFWriter(pw);
	        

			
			GFFEntrySet gff = GFFTools.readGFF(new BufferedReader(new FileReader(gffinputFile)));
			for(Iterator<?> gffi = gff.lineIterator(); gffi.hasNext();){
				Object o = gffi.next();
				if(o instanceof GFFRecord){
					GFFRecord gffRec = (GFFRecord) o;
					String oneRecordFeature = gffRec.getFeature();
					
					boolean Chr   = false;
					String chr    = gffRec.getSeqName();
					
					//Only "2L", "2R", "3L", "3R", "4", "U" and "X" is accepted as chromosomes
					if( chr.equals("2L") || chr.equals("2R") || chr.equals("3L") || chr.equals("3R") || chr.equals("4") || chr.equals("U") || chr.equals("X")){
						Chr = true;
					}
					
					
					//accept only "exon"s and "CDS"s but not strting and ending codons which are there in input file.
					if( Chr== true &&   ( oneRecordFeature.equals("exon") || oneRecordFeature.equals("CDS")) ){
						gffw.recordLine(gffRec);
					}
				}
			}
			pw.flush();
	}/*removeStartingCodonRecordsFromAGffFile*/
	
	
	
	public static void main(String[] args)  throws IOException, Exception{
		// TODO Auto-generated method stub
		
		FileCleanUp aFile = new FileCleanUp();
		
		 String inputFileName = args[0];
		 String outputFileName = args[1];
		 
		 aFile.setInputGffFileName(outputFileName);
		 aFile.setOutputGffFileName(outputFileName);
		 
		 aFile.removeStartingCodonRecordsFromAGffFile(inputFileName, outputFileName);
		 

	}/*main*/

}/*FileCleanUp*/
