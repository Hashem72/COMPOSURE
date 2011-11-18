/**
 * 
 */
package valuepostprocessing;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.BioException;
import org.biojava.bio.program.gff.GFFEntrySet;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.GFFTools;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.utils.ParserException;

import series.SeriesIO;

import cern.jet.random.engine.RandomGenerator;


import net.derkholm.nmica.matrix.Matrix2D;
import net.derkholm.nmica.matrix.SimpleMatrix2D;
import net.derkholm.nmica.utils.CliTools;
import net.derkholm.nmica.utils.ConfigurationException;

/**
 * @author hk3 26-05-2011
 * this is to generate the 3 by 3 transition matrix for the negative set. Our positive set was DnaseI or even RedFly enhancers. Our negative set includes the same number of sequences
 * each of which has a length equal to one of the sequences in the positive set, and it has been ranmdomly picked from  a block which is in turn randomly chosen.
 *
 */
public class NegativeModel {

	/**
	 * @param args
	 */
	private String serDir         = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/Eisenlab_Composure_Output_Model_20051212_With_Window_Length_350_SER";
	private String aliDir         = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ALIGNMENTS_EISENLAB/2L";
	private String pathToGffFiles = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/Enhancers_Mapped_to_Alignments";
	private String chr            = "2L";
	private String outputFileName = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/MarkovModels/EnhancersNegative.txt";// this is to output negative model
	private int minLength         = 200; // sequence shorter (after removing gaps) than this will not be accepted for training  the model.
	
	
	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		
		NegativeModel app = new NegativeModel();
		args = CliTools.configureBean(app, args);
		app.run(args);


	}/*main*/
	
	private void run(String [] args) throws Exception{
		
			
		String serDir                    = getSerDir();
		String pathToGffs                = getPathToGffFiles();
		int minThresholdForLength        = getMinLength();
		PositiveModel aPositiveSetObject = new PositiveModel();
		String pathToAlidir              = getAliDir();
		
		ArrayList<Integer> listOfLengthOfAllGSRecords       = getlistOfLengthsOfAllGoldStandardRecords(serDir,pathToGffs);
		ArrayList<String> listOfStringInNegativeSet         = makeNegativeSetOfStrings(listOfLengthOfAllGSRecords,pathToGffs, serDir, pathToAlidir,minThresholdForLength );		
		Matrix2D   frequencyMatrixFromAllBlocksNegativeSet  = aPositiveSetObject.getTransitionMatrixFromAListOfStrings(listOfStringInNegativeSet);
		Matrix2D   transitionMatrixFromAllBlocksNegativeSet = aPositiveSetObject.convertFrequencyToProbability(frequencyMatrixFromAllBlocksNegativeSet);
		aPositiveSetObject.writeAMatrixIntoAFile(getOutputFileName(), transitionMatrixFromAllBlocksNegativeSet);
		
	}/*run*/
	
	////////// METHODS 
	
	public void setOutputFileName(String of){
		this.outputFileName = of;
	}/*setOutputFileName*/
	
	public void setSerDir(String sd){
		this.serDir = sd;
	}/*setSerDir*/
	
	public void  setMinLength(int minLen){
		this.minLength = minLen;
	}/*setMinLength*/
	
	public void setAliDir(String ad){
		this.aliDir = ad;
	}/*setAliDir*/
	
	public void setPathToGffFiles(String pgff){
		this.pathToGffFiles = pgff;
	}/*setPathToGffFiles*/
	
	public void setChr(String chr){
		this.chr = chr;
	}/*setChr*/
	
	public String getSerDir(){
		return serDir;
	}/*getSerDir*/
	
	public int getMinLength(){
		return minLength;
	}/*getMinLength*/
	
	public String getAliDir(){
		return aliDir;
	}/*getAliDir*/
	
	public String getOutputFileName(){
		return outputFileName;
	}/*getOutputFileName*/
	
	public String getPathToGffFiles(){
		return pathToGffFiles;
	}/*getPathToGffFiles*/
	
	public String  getChr(){
		return chr;
	}/*getChr*/
	
	public String getAliFileNameFromSerFileName(String serFileName, String pathToAliDir){// it is assumed that ser file is of chr_start_end.ser format
		if(serFileName.startsWith(".")){
			System.err.print("Bad ser file name!");
		}
		String []  serFileSpecifications = serFileName.split("\\.");
		String blockId = serFileSpecifications[0];	
		String aliFileName = pathToAliDir+"/"+blockId+".fa";
		return aliFileName;
	}/*getAliFileNameFromSerFileName*/
	
	public ArrayList<String> makeNegativeSetOfStrings(ArrayList<Integer> listOfLengthsInPositiveSet, String pathToGffFiles, String serDirName, String pathToAlidir, int lengthThreshold)
	throws Exception{
		ArrayList<String> negativeSetOfStrings = new ArrayList<String>();
		
		File serDir = new File( serDirName);
		File[] listOfSerFiles = serDir.listFiles(); 
		String aliFileName = null;
		Sequence firstSeq = null;
		
		int numberOfSerFiles = listOfSerFiles.length;
		Random generator = new Random();
		
		for(int lengthOfOneStr: listOfLengthsInPositiveSet){
			//int lengthOfOneStr = str.length();
			Matrix2D m = null;
			int lengthOfThisBlock = -1;
			while(lengthOfThisBlock <= 2*lengthOfOneStr ){
				int aRandomNumber = generator.nextInt( numberOfSerFiles);
				//System.out.println("0 < " +aRandomNumber + " < " +numberOfSerFiles );
				if( !listOfSerFiles[aRandomNumber].getName().startsWith(".")){
					File serFile = listOfSerFiles[aRandomNumber];
					aliFileName= getAliFileNameFromSerFileName(serFile.getName(),pathToAlidir );
					
					 m = SeriesIO.readSeries(serFile);
					lengthOfThisBlock = m.rows();
					}
			}
			int diffBetweenLengths = lengthOfThisBlock- lengthOfOneStr;
			int strStart           = generator.nextInt(diffBetweenLengths);
			int strEnd             = strStart+lengthOfOneStr;
			String chr             = getChr();
			SimpleGFFRecord r = new SimpleGFFRecord();
			r.setSeqName(chr);
			r.setStart(strStart);
			r.setEnd(strEnd);
			PositiveModel aPosObject = new PositiveModel();
			//System.out.println("str = " + strStart + " end = " + strStart + " blockLength = "+ lengthOfThisBlock);
			File aliFile = new File(aliFileName);
			if(! aliFile.exists()){
				System.out.println("ali file  " + aliFile.getName() + " not exists!");
				continue;
			}
			firstSeq = SeqIOTools.readFastaDNA(new BufferedReader(new FileReader(aliFile))).nextSequence();
			String oneNegativeStr = aPosObject.getStringOfSymbolsForOneRecord(r, m, firstSeq);
			int lengthOfStr = oneNegativeStr.length();
			if( lengthOfStr >= lengthThreshold){
				negativeSetOfStrings.add(oneNegativeStr); 
			}
			else{
				System.out.println("one sequence ignored because after removing  gaps its length became shorter than enhancer length threshold!");
			}
						
		}
		
		return negativeSetOfStrings;
	}/*makeNegativeSetOfStrings*/
	
	public ArrayList<Integer> getlistOfLengthsOfAllGoldStandardRecords(String pathToSerFiles, String  pathToGffFiles) 
	throws FileNotFoundException, ParserException, BioException, IOException{
		ArrayList<Integer> listOfAllLengths = new ArrayList<Integer>();
		
		File serDir = new File( pathToSerFiles);
		Pattern fnPattern = Pattern.compile("(([^_]+)_([0-9]+)_([0-9]+)).ser");
		int numberOfmissingFiles = 0;
		
		for(File aSerFile: serDir.listFiles()){
			if(aSerFile.getName().startsWith(".")){
				continue;
			}
			Matcher fnMatcher = fnPattern.matcher(aSerFile.getName());
			String blockId;
			
			if(fnMatcher.matches()){
				blockId = fnMatcher.group(1);
				}
			else {
				System.out.println("Bad file name detected: " + aSerFile.getName());
				continue;
				}
			String gffFileName = pathToGffFiles+"/"+blockId+".gff";
			File  gffFile      = new File(gffFileName);
			boolean gffExists = gffFile.exists();
			
			if(!gffExists){
				//System.out.println("There is no gff file for  block " + blockId);
				numberOfmissingFiles++;
				continue;
			}
			GFFEntrySet  gff    = GFFTools.readGFF(gffFile);
			for(Iterator<?> gffi= gff.lineIterator(); gffi.hasNext();){
				Object o = gffi.next();
				if(o instanceof GFFRecord){
					GFFRecord grec = (GFFRecord) o;
					int oneStart = grec.getStart();
					int oneEnd   = grec.getEnd();
					int oneLength = oneEnd-oneStart;
					listOfAllLengths.add(oneLength);
				}
			}
		}
		

		
		return listOfAllLengths;
	}/*getlistOfLengthsOfAllGoldStandardRecords*/


}/*NegativeModel*/
