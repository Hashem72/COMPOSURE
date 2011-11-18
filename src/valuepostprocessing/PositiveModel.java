/**
 * 
 */
package valuepostprocessing;



import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.BioException;
import org.biojava.bio.program.gff.GFFEntrySet;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.GFFTools;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ParserException;


import org.bjv2.util.IntSet;
import org.bjv2.util.IntSets;

import series.SeriesIO;


import net.derkholm.nmica.matrix.Matrix2D;
import net.derkholm.nmica.matrix.SimpleMatrix2D;
import net.derkholm.nmica.utils.CliTools;



/**
 * @author hk3
 * this class is to make DnaseI model by reading all DnaseI records in chromosome and learn the parameters.
 * 25-05-2011
 */
public class PositiveModel {

	/**
	 * @param args
	 */
	
	private String blockId = null;
	private String pathToGffFiles = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/CodingRegions_Mapped_to_Alignments"; // this dir contaings DnaseI mapped to Alignments
	private String compDir        = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/Eisenlab_Composure_Output_Model_20051212_With_Window_Length_350_SER_temp"; //dir that contains comp files
	private String aliDir         = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ALIGNMENTS_EISENLAB/2L";
	private String outputFileName ="/Users/hk3/Desktop/Main/dummy.txt"; //"/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/MarkovModels/CodingRegionsPositive.txt";// positive model as a matrix will be written into this file
	private int minLength         = 200; // it is assuemed enhancers are longer than 100 and therefore sequences shorter than this threshold is not informative for training the model hence they will  be ignored
	
	
	//SETTERS
	public void setBlockId(String Id){
		this.blockId = Id;
	}/*setBlockId*/
	
	public void setMinLength(int minLen){
		this.minLength = minLen;
	}/*minLen*/
	
	public void setPathToGffFiles (String path){
		this.pathToGffFiles = path;
	}/*setPathToGffFiles*/
	
	public void setCompDir(String cd){
		this.compDir  = cd;
	}/*setCompDir*/
	
	public void setAliDir(String ad){
		this.aliDir = ad;
	}/*setAliDir*/
	
	//GETTERS
	public String  getBlockId(){
		return blockId;
	}/*getBlockId*/
	
	public int getMinLength(){
		return minLength;
	}/*getMinLength*/
	
	public String getPathToGffFiles(){
		return pathToGffFiles;
	}/*getPathToGffFiles*/
	
	public String getCompDir(){
		return compDir;
	}/*getCompDir*/
	
	public String getAliDir(){
		return aliDir;
	}/*getAliDir*/
	
	public void setOutputFileName(String of){
		this.outputFileName = of;
	}/*setOutputFileName*/
	
	
	public String getOutputFileName(){
		return outputFileName;
	}/*getOutputFileName*/
	
	
	public static int getIndexOfMaximumElementInAnArray(double [] arr){
		int index = -1;
		double max = Double.NEGATIVE_INFINITY;
		for(int i =0; i<arr.length; i++){
			if(arr[i] >= max){
				max = arr[i];
			}
		}
		System.out.println("max is " + max );
		for(int i = 0; i<arr.length; i++){
			if (max == arr[i]){
				index = i;
			}
		}
		return index;
		
	}/*getIndexOfMaximumElementInAnArray*/
	
	public void writeAMatrixIntoAFile(String fileName, Matrix2D m) throws IOException{
		File outputFile = new File(fileName);
		int numCol   = m.columns();
		int numRows  = m.rows();
		try{
			FileWriter fstream  = new FileWriter(outputFile);
			BufferedWriter bw = new BufferedWriter(fstream); 
			
			for(int i = 0; i < numRows; i++){
				String oneLine ="";
				for (int j = 0; j < numCol; j++){
				oneLine = oneLine+ m.get(i,j) + "\t";
				}
				bw.write(oneLine);
				bw.newLine();
			
			}
			bw.close();			
		}catch(Exception e){
			System.err.println("Error: " + e.getMessage());
		}	
	}/*writeAMatrixIntoAFile*/
	
	
	
	public String getStringOfSymbolsForOneRecord(GFFRecord gffr, Matrix2D m, Sequence seq){
		String str= "";
		int recordStart = gffr.getStart();
		int rcordEnd    = gffr.getEnd();
		//System.out.println("st = " + recordStart + " end = " + rcordEnd);
		IntSet range   = IntSets.range(recordStart, rcordEnd);
		Symbol gap = seq.getAlphabet().getGapSymbol();
		//System.out.println("seq length = " + seq.length());
		
		for(int i: range){
			if( (i!= seq.length()) && (seq.symbolAt(i+1) != gap)){// first condition is because m gets i ranging from 0 up to  m.rows (not m.row itself) which is equal to seq.length 
				int  maxIndex = -1;
				ArrayList<Double> oneRow = new ArrayList<Double>();
				
				for(int j =0;j<3;j++){
					
					oneRow.add(m.get(i, j));
				}
				
				Object obj_max = Collections.max(oneRow);
				maxIndex = oneRow.indexOf(obj_max);
				if(maxIndex ==0){
					str = str+"R";
				}
				else if(maxIndex == 1){
					str = str+"G";
				}
				else if (maxIndex == 2){
					str = str+"B";
				}
				else{
					System.err.println("Out of bound!");
				}
			}			
		}
		
		
		//System.out.println("st = " +  str);
		return str;
	}/*getStringOfSymbolsForOneRecord*/
	
	public ArrayList<String> getSringOfSymbolsForAllRecordsInAGffFile(File gffFile, Matrix2D serMatrix, Sequence seq, int lengthThreshold) 
	throws FileNotFoundException, ParserException, BioException, IOException{// serMatrix is the data from the ser file that has been used to map this gff file into alignment
		ArrayList<String> listOfStringOfSymbols = new ArrayList<String>();
		
		GFFEntrySet records  = GFFTools.readGFF(gffFile);
		 for(Iterator<?> gffi= records.lineIterator(); gffi.hasNext();){
			 Object o = gffi.next();
			 if(o instanceof GFFRecord){
				 GFFRecord gffr = (GFFRecord) o;
				 String str = getStringOfSymbolsForOneRecord(gffr,serMatrix, seq);
				 if(str.length()>=lengthThreshold){
					 listOfStringOfSymbols.add(str);					 
				 }
				 else{
					 System.out.println("A sequence ignored because after removing gaps its length became shorter than length threshold!");
				 }
				// System.out.println("Str = " + str);
			 }
		 }		
		return listOfStringOfSymbols;
	}/*getSringOfSymbolsForAllRecordsInAGffFile*/
	
	public ArrayList<String> getAllStringOfSymbolsAssociatedToAllSerFiles(String serDirName, String pathToGffFiles, int lengthThreshold) throws Exception{
		ArrayList<String> stringOfSymbols = new ArrayList<String>();
		
		File serDir = new File( serDirName);
		Pattern fnPattern = Pattern.compile("(([^_]+)_([0-9]+)_([0-9]+)).ser");
		
		
		
		int missingFileIndicator =0;
		for(File aSerFile: serDir.listFiles()){	
			if(aSerFile.getName().startsWith(".")){
				continue;
			}
			//System.out.println("file name is " + aSerFile.getName());
			//System.out.println("aSerFile" + aSerFile.getName());
			Matrix2D m          = SeriesIO.readSeries(aSerFile);			
			Matcher fnMatcher = fnPattern.matcher(aSerFile.getName());
			String blockId;
			
			if(fnMatcher.matches()){
				blockId = fnMatcher.group(1);
				}
			else {
				System.out.println("Bad file name detected: " + aSerFile.getName());
				continue;
				}
			String aliFileName = getAliDir()+"/"+blockId+".fa";
			File    aliFile    = new File(aliFileName);
			if(!aliFile.exists()){
				System.out.println("Warning! Could not find alignment file for block " + blockId );
				continue;
			}
			Sequence firstSeq = SeqIOTools.readFastaDNA(new BufferedReader(new FileReader(aliFile))).nextSequence();
			String gffFileName = pathToGffFiles+"/"+blockId+".gff";
			File  gffFile      = new File(gffFileName);
			//System.out.println("block = "+ blockId + "gff file " + gffFileName);
			boolean gffExists = gffFile.exists();
			
			if(!gffExists){
				//System.out.println("There is no gff file for  block " + blockId);
				missingFileIndicator++;
				continue;
			}
			ArrayList<String> listOfStringsForOneBlock = new ArrayList<String>();
			
			
			listOfStringsForOneBlock = getSringOfSymbolsForAllRecordsInAGffFile(gffFile,m,  firstSeq, lengthThreshold);
			stringOfSymbols.addAll(listOfStringsForOneBlock);
			
			
			}
		if(missingFileIndicator >=1){
			System.out.println("Note:  there were " +  missingFileIndicator + " blocks withiout any hit from this feature." );
		}
		return stringOfSymbols;
		
	}/*getAllStringOfSymbolsAssociatedToAllSerFiles*/
	
	public Matrix2D getTransitionMatrixFromAListOfStrings(ArrayList<String> listOfStrings){// note it is assumed that this ArrayList contains strings like "RGBBGRRRGGGBBB" ie strings alphabet R, G and B only
		Matrix2D  transitionMatirx = new SimpleMatrix2D(3,3);
		for(String str:listOfStrings){
			Matrix2D oneMatrix  = new SimpleMatrix2D(3,3);
			oneMatrix           = getFrequencyOfDiSymbolsFromAStringOfSymbols(str);
			//update transition matrix:
			transitionMatirx    = addTwoMatrices(transitionMatirx,oneMatrix);
		}
		return transitionMatirx;
	}/*getTransitionMatrixFromAListOfStrings*/
	
	public Matrix2D getFrequencyOfDiSymbolsFromAStringOfSymbols(String str){//note:here  is like str = "RRRRGGGBBBRGBR" where R stands for red, G stands for green and B stands for blue, given such a string
		// we want to see how many times for example R is followd by G. becasue the alphabet is {R,G, B} this function will fill a 3 by 3 matrix.
		Matrix2D freqMatrix = new SimpleMatrix2D(3,3); 
		int strLength = str.length();
		//R follows
		int RfR = 0; 
		int RfG = 0; // means  number of times letter  G followed  letter R ie occurrences of  "GR" in the sequence
		int RfB = 0;
		
		//G follows
		int GfR = 0;
		int GfG = 0;
		int GfB = 0;
		
		//B follows
		int BfR = 0;
		int BfG = 0;
		int BfB = 0;
		
		
		for(int i = 0; i< strLength-1; i++){
			String substr = str.substring(i, i+2);
			//System.out.println("i = " + i + " substr =" + substr);
			
			if(substr.equals("RR")){
				RfR++;
			}
			else if(substr.equals("GR")){
				RfG++;
			}
			else if(substr.equals("BR")){
				RfB++;
			}
			else if(substr.equals("RG")){
				GfR++;
			}
			else if(substr.equals("GG")){
				GfG++;
			}
			else if(substr.equals("BG")){
				GfB++;
			}
			else if(substr.equals("RB")){
				BfR++;
			}
			else if(substr.equals("GB")){
				BfG++;
			}
			else if(substr.equals("BB")){
				BfB++;
			}
			else{
				System.err.print("Unknown symbol!");
			}
			
		}
		
		//fill in the matrxi accordingly. note that 
		freqMatrix.set(0, 0, RfR);
		freqMatrix.set(0, 1, GfR);
		freqMatrix.set(0, 2, BfR);
		freqMatrix.set(1, 0, RfG);
		freqMatrix.set(1, 1, GfG);
		freqMatrix.set(1, 2, BfG);
		freqMatrix.set(2, 0, RfB);
		freqMatrix.set(2, 1, GfB);
		freqMatrix.set(2, 2, BfB);	
		return freqMatrix;
		
	}/*getFrequencyOfDiSymbolsFromAStringOfSymbols*/
	
	public Matrix2D convertFrequencyToProbability(Matrix2D A){
		int numOfCols = A.columns();
		int numOfRows = A.rows();
		for(int i =0; i< numOfRows; i++){
			double sumOverOnRow = 0;
			
			//get sum over one row
			int c = 0;
				while (c < numOfCols ){
				sumOverOnRow = sumOverOnRow + A.get(i, c);
				c++;			
			}
			
			//divide each element of that row by the sum of the row
			for (int j= 0; j<numOfCols; j++ ){
				A.set(i, j, A.get(i, j)/sumOverOnRow);
			}
		}
		return A;
	}/*convertFrequencyToProbability*/
	
	public Matrix2D addTwoMatrices(Matrix2D A, Matrix2D B){
		
		int numberOfRwosInA  = A.rows();
		int numberOfRowsInB  = B.rows();
		int numberOfColsInA  = A.columns();
		int numberOfColsInB  = B.columns();
		
		if ( numberOfRwosInA !=  numberOfRowsInB   && numberOfColsInA != numberOfColsInB){
			System.err.println("Matrices have got different dimensions and is not possible to add them up!");
			System.exit(0);
		}
		SimpleMatrix2D  sum = new SimpleMatrix2D(numberOfRwosInA, numberOfRowsInB);
		for(int i= 0; i< numberOfRwosInA; i++){
			for (int j = 0; j< numberOfColsInA; j++){
				double v = A.get(i, j)+ B.get(i, j);
				sum.set(i, j, v);
			}
		}
		
		return sum;
		
	}/*addTwoMatrices*/
	
	public String getBlockIdFromCompFileName (String compFileName){
		String [] specifications =  compFileName.split("\\.");
		String blockBaseName = specifications[0];
		return blockBaseName;
	}/*getBlockIdFromCompFileName*/
	
	
	
	
	
	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		
		PositiveModel app = new PositiveModel();
		args = CliTools.configureBean(app, args);
		app.run(args);

	}/*main*/
	
	private void run(String [] args) throws Exception{
		
		String serDir     = getCompDir();
		String pathToGffs = getPathToGffFiles();
		int enhancerLengthThreshold = getMinLength();
		ArrayList<String> listOfStrings = getAllStringOfSymbolsAssociatedToAllSerFiles(serDir,pathToGffs, enhancerLengthThreshold);
		
		
		
		System.out.println("size of listOfStrings is  " + listOfStrings.size());
		Matrix2D         frequencyMatrixFromAllBlocks = getTransitionMatrixFromAListOfStrings(listOfStrings);
		Matrix2D         transtionMatrixFromAllBlocks = convertFrequencyToProbability(frequencyMatrixFromAllBlocks);
		
		for(int i = 0; i <3; i++){
			for(int j = 0 ; j<3; j++){
				System.out.print(transtionMatrixFromAllBlocks.get(i, j) + "  ");
			}
			System.out.println();
		}
		
		writeAMatrixIntoAFile(getOutputFileName(),transtionMatrixFromAllBlocks);
		
	}/*run*/

}/*PositiveModel*/
