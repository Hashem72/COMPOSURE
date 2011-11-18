/**
 * 
 */
package valuepostprocessing;

import java.io.File;
import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.regex.Pattern;
import java.io.FileReader;
import java.io.Reader;
import java.util.regex.Matcher;
import java.util.regex.Pattern;



import java.io.FileNotFoundException;

import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.biojava.bio.BioException;
import org.biojava.bio.program.gff.GFFEntrySet;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.GFFTools;
import org.biojava.bio.program.gff.GFFWriter;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.biojava.bio.seq.StrandedFeature.Strand;
import org.biojava.utils.ParserException;
import org.bjv2.util.IntSet;
import org.bjv2.util.IntSets;


/**
 * @author hk3 30-03-2011
 * this is meant to get union of overlapped blocks between entries of a  gff file
 *
 */
public class UnionOfFeaturesInAGffFile {

	/**
	 * @param args
	 */
	private String gffGenesInputFile     = null;
	private String gffGenesOutputFile    = null;
	private String gffDnaseIInputFile    = null;
	private String gffDnaseIOutputFile   = null;
	
	private String gffEnhancersInputFile  = null;
	private String gffEnhancersOutputFile = null;
	
	private String [] dMelChromosomes = {"22"};//{"2L", "3L","2R","3R","X", "4"};
	
	

	
	//setters
	
	public void setGffInputFile (String gif){
		this.gffGenesInputFile = gif;
	}/*setGffInputFile*/
	
	public void setGffOutputFile(String gof){
		this.gffGenesOutputFile = gof;
	}/*setGffOutputFile*/
	
	public void setGffDnaseIInputFile(String dif){
		this.gffDnaseIInputFile= dif;
	}/*setGffDnaseIInputFile*/
	
	public void setGffDnaseIOutputFile(String dof){
		this.gffDnaseIOutputFile= dof;
	}/*setGffDnaseIOutputFile*/
	
	public void setGffEnhancersInputFile(String eif){
		this.gffEnhancersInputFile = eif;
	}/*setGffEnhancersInputFile*/
	
	public void setGffEnhancersOutputFile(String eof){
		this.gffEnhancersOutputFile = eof;
	}/*setGffEnhancersOutputFile*/
	
	public void setDmelChromosomes (String [] Chro){
		this.dMelChromosomes = Chro;
	}/*setDmelChromosomes*/
	

	
	
	
	//getters
	public String getGffInputFile(){
		return gffGenesInputFile;
	}/*getGffInputFile*/
	
	public String getGffOutpuFile(){
		return gffGenesOutputFile;
	}/*getGffOutpuFile*/
	
	public String getGffDnaseIInputFile(){
		return gffDnaseIInputFile;
	}/*getGffDnaseIInputFile*/
	
	public String getGffDnaseIOutputFile(){
		return gffDnaseIOutputFile;
	}/*getGffDnaseIOutputFile*/
	
	public String getGffEnhancersInputFile(){
		return gffEnhancersInputFile;
	}/*getGffEnhancersInputFile*/
	
	public String getGffEnhancersOutputFile(){
		return gffEnhancersOutputFile;
	}/**/
	
	public String [] getDmelChromosomes (){
		return dMelChromosomes;
	}/*getDmelChromosomes*/
	
	public int getNumberOfLinesOfInAFile(String fileName) throws IOException {
		int lineCounter = 0;
		try{
		Reader reader = new FileReader(fileName);
		BufferedReader br = new BufferedReader(reader);
		String line;
		while((line = br.readLine()) != null){
			lineCounter++;
		}
		}
		catch(FileNotFoundException fe){
			fe.fillInStackTrace();
		}
		return lineCounter;
	}/*getNumberOfLinesOfInAFile*/
	
	public ArrayList<String> getLinesOfAFile(String fileName) throws IOException{
		ArrayList<String> lines = new ArrayList<String>();
		
		
		Reader      reader  = new FileReader(fileName);
		BufferedReader br        = new BufferedReader(reader);
		
		int NumberOfLines = getNumberOfLinesOfInAFile(fileName);
		for(int i= 0; i<NumberOfLines; i++){
			String line = br.readLine();
			lines.add(line);
		}
		
			return lines;
		
	}/*getLinesOfAFile*/
	
	
	public IntSet getUnionOfRecordsFromAGFFFile(String inputFile) throws IOException{
		IntSet union;
		
		FileInputStream fis = new FileInputStream(new File(inputFile));
		int empty = fis.available();
		if(empty ==0 ){
			union = null;
			System.out.println("Warning: File " + inputFile + " was empty!");
		}
		else
		{
		
		
		
			ArrayList<String> lines = getLinesOfAFile(inputFile);
			String firstLine = lines.get(0);
			String [] firstRecord = firstLine.split("\t");
			int firstStart = Integer.parseInt(firstRecord[3]);
			int firstEnd  = Integer.parseInt(firstRecord[4]);
			union = IntSets.range(firstStart, firstEnd);
		
			for(int i = 0; i<lines.size(); i++){
				String oneLine = lines.get(i);
				String [] oneRecord = oneLine.split("\t");
				
				int start          = Integer.parseInt(oneRecord[3]);
				int end            = Integer.parseInt(oneRecord[4]);
				IntSet thisPositions = IntSets.range(start, end);
				union = getUnionOfTwoIntSets(union,thisPositions);
			}
		}


		return union;
	}/*getUnionOfRecordsFromAGFFFile*/
	
	public void getUnionOfRecordsFromAGFFFileAndWriteThemIntoAGFFFile(String inputFile, String outputFile, String [] Chromosomes) throws IOException, ParserException, BioException{
		
    	File gffoutputFile      = new File(outputFile);	    			
		PrintWriter pw          = new PrintWriter(new FileWriter(gffoutputFile, true));
		GFFWriter gffw          = new GFFWriter(pw);
		
        File    aGffFile = new File(inputFile);
		GFFEntrySet gff = GFFTools.readGFF(aGffFile);

		

		for(String Chr:Chromosomes){
			System.out.println("Looking for union of features over chormosome  " + Chr);
			String tempFileName  = "temp.gff";
			File   tempFile      = new File(tempFileName);
			PrintWriter tpw      = new PrintWriter(new FileWriter(tempFile));
			GFFWriter tgffw      = new GFFWriter(tpw);
			 for (Iterator<?> gffi = gff.lineIterator(); gffi.hasNext(); ){
				 Object o = gffi.next();
				 if(o instanceof GFFRecord){
					 GFFRecord gro = (GFFRecord) o;
					 if(gro.getSeqName().equals(Chr)){
						 tgffw.recordLine(gro);
					 }
				 }/*if*/
			 }/*for_Iter*/
			tpw.flush();
		//}/*for Chr*/

		
		
		IntSet union = getUnionOfRecordsFromAGFFFile(tempFileName);
		
		
		if(union == null){
			SimpleGFFRecord rDummy    = new SimpleGFFRecord();
			rDummy.setSeqName("dummy");
			rDummy.setSource("dummy");
			rDummy.setFeature("dummy");
			rDummy.setStrand(null);
			rDummy.setStart(0);
			rDummy.setEnd(1);
			//gffw.recordLine(rDummy);	
		}
		else
		{
		
		
		int start = -1;
		int end   = -1;
		int dummyEnd = union.getMax()+5;//this dummy five is offset ( to make sure whole sequence is covered)
		int dummyStart = union.getMin()-5;
		
		for(int i= dummyStart; i<dummyEnd; i++){
			SimpleGFFRecord r    = new SimpleGFFRecord();
			boolean isMin = isMin(i, union);
			if (isMin){
				 start = i;				
			}
			boolean isMax = isMax(i,union);
				if(isMax){
					end = i;
					if(start == -1){
						start = union.getMin();
					}
				 //System.out.println("start = " + start + " end = " + end);
					String dataType = getDataType(inputFile); 
					r.setSeqName(Chr);
					r.setSource("extended");
					r.setFeature(dataType);
					r.setStrand(null);
				 //r.setScore(0);
					r.setStart(start);
					r.setEnd(end);
					gffw.recordLine(r);
					}
				}
			}
		pw.flush();
		deleteFile(tempFileName);
		}/*for Chr*/
	}/*getUnionOfRecordsFromAGFFFileAndWriteThemIntoAGFFFile*/

	
	public boolean isMin(int i, IntSet A){
		boolean isMin = false;
		if(A.contains(i) && !(A.contains(i-1))){
			isMin = true;
		}
		return isMin;
	}/*isMin*/
	
	public boolean isMax(int i , IntSet A){
		boolean isMax = false;
		if(A.contains(i) && !(A.contains(i+1))){
			isMax = true;
		}
		return isMax;
	}/*isMax*/
	

	

	
	public IntSet getUnionOfTwoIntSets(IntSet A, IntSet B){//it seems to me that the IntSets.union doesnt work properly. when one set contains another, it doesnt return mathematical union of two sets. This is t get rid of that.
		
		IntSet  C = null;
		if(IntSets.contain(A, B)){
			C = A;
		}
		else if(IntSets.contain(B, A)){
			C= B;
		}
		else{
			C = IntSets.union(A,B);
		}
	return C;		
	}/*getUnionOfTwoIntSets*/
	
	
	public void testUnion(){
		IntSet A  = IntSets.range(10, 20);
		IntSet B   = IntSets.range(15, 25);
		IntSet C = null;
		if(IntSets.contain(A, B)){
			C = A;
		}
		else if(IntSets.contain(B, A)){
			C= B;
		}
		else{
			C = IntSets.union(A,B);
		}
		System.out.println("Union: ");
		for(int ii:C){
			System.out.print(ii + ", ");
		}
		
		IntSet inter = IntSets.intersection(A,B);
		System.out.println("Intersection: ");
		for(int jj :inter){
			System.out.print( jj + ", ");
		}
	}/*testUnion*/
	
	public void getUnion(){
		IntSet A  = IntSets.range(3, 12);
		IntSet B   = IntSets.range(3, 12);
		
		IntSet Union = IntSets.union(A,B);
		IntSet Intersection = IntSets.intersection(A,B);
		for(int i:Intersection){
			System.out.print(i+ ", ");
		}
		
	}/*getUnion*/
	
	public String getDataType (String inputFileName){
		String dataType = null;
		String [] patterns = {"Genes", "Enhancers","DNaseI", "dmel_transcripts_and_exons", "bdtnpDnaseAccAllStages", "redfly_download","Chaines", "ChainesInRealCoordinages"};
		boolean found = false;
		for(String x:patterns){
			Pattern p = Pattern.compile(x);
			Matcher matcher = p.matcher(inputFileName);

			if(matcher.find()){
				dataType = x;
				found = true;
			}
		}
		if(!found){
			System.out.println("Bad data type!");
		}
		
		
		if(  dataType.equals("Genes")  || dataType.equals("dmel_transcripts_and_exons") ){
			dataType= "exons";
		}
		if(dataType.equals("Enhancers")  || dataType.equals("redfly_download")){
			dataType = "EnhancersRedFly";
		}
		if(dataType.equals("DNaseI") || dataType.equals("bdtnpDnaseAccAllStages")   ){
			dataType= "DNaseIUCSC";
		}
		if(dataType.equals("Chaines") || dataType.equals("ChainesInRealCoordinages")){
			dataType = "ChainFromPD";
		}
		//System.out.println("data type = "+ dataType);
		return dataType;
		
	}/*getDataType*/
	public void deleteFile(String file){
		File f = new File(file);
		boolean exists = f.exists();
		if(exists){
			boolean success = f.delete();
			if(!success){
				System.out.println("Deletion failed.");
				System.exit(0);
			}
		}
	}/*deleteFile*/

	

	
	
	public static void main(String[] args) throws ParserException,BioException,IOException {
		// TODO Auto-generated method stub
		
		//String gffGenesInputFile     = args[0];
		//String gffGenesOutputFile    = args[1];
		String gffDnaseIInputFile    = args[2];
		String gffDnaseIOutputFile   = args[3];
		//String gffEnhancersInputFile = args[4];
		//String gffEnhancesOutputFile = args[5];
		
		UnionOfFeaturesInAGffFile aGFF = new UnionOfFeaturesInAGffFile();
		
		String dMelChr [] = aGFF.getDmelChromosomes();
		//aGFF.setGffInputFile(gffGenesInputFile);
		//aGFF.setGffOutputFile(gffGenesOutputFile);
		aGFF.setGffDnaseIInputFile(gffDnaseIInputFile);
		aGFF.setGffDnaseIOutputFile(gffDnaseIOutputFile);
		//aGFF.setGffEnhancersInputFile(gffEnhancersInputFile);
		//aGFF.setGffEnhancersOutputFile(gffEnhancesOutputFile);
		
		//delete old files if exists
		//aGFF.deleteFile(gffGenesOutputFile);
		aGFF.deleteFile(gffDnaseIOutputFile);
		//aGFF.deleteFile(gffEnhancesOutputFile);
		
		//aGFF.getOverlaps(gffGenesInputFile, gffGenesOutputFile);
		//aGFF.getUnionOfRecordsFromAGFFFileAndWriteThemIntoAGFFFile(gffGenesInputFile, gffGenesOutputFile,dMelChr);
		aGFF.getUnionOfRecordsFromAGFFFileAndWriteThemIntoAGFFFile(gffDnaseIInputFile, gffDnaseIOutputFile,dMelChr);
		//aGFF.getUnionOfRecordsFromAGFFFileAndWriteThemIntoAGFFFile(gffEnhancersInputFile, gffEnhancesOutputFile, dMelChr );
	
		
		
	}/*main*/

}/*UnionOfFeaturesInAGffFile*/




