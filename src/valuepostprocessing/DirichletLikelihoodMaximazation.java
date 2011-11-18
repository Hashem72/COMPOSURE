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
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import net.derkholm.nmica.matrix.Matrix2D;
import net.derkholm.nmica.utils.CliTools;
import net.derkholm.nmica.utils.ConfigurationException;

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
import umontreal.iro.lecuyer.probdistmulti.DirichletDist;

/**
 * @author hk3
 *
 */
public class DirichletLikelihoodMaximazation {

	/**
	 * @param args
	 * 20-06-2011
	 * Given a Training data set X in a matrix form (N by 3 : where N is the number of data and 3 is because I want to look at 3 profiles: Red, Green and Blue at each base)
	 * this is to find alpha_1, alpha_2 and alpha_3 that maximize P(X|alpha) where alpha = (alpha_1,alpha_2,alpha_3) 
	 * 
	 * 
	 */
	
	private String blockId = null;
	private String pathToGffFiles = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/TrainingDataSet/JunkState"; // this dir contaings training data set as gff files
	private String serDir        = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/Eisenlab_Composure_Output_Model_20051212_With_Window_Length_350_SER_temp"; //dir that contains ser files
	private String aliDir         = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ALIGNMENTS_EISENLAB/2L";
	private String chr            = "2L";
	private String posOutputFile  = "/Users/hk3/Desktop/pos.txt";
	private String negOutputFile  = "/Users/hk3/Desktop/neg.txt";
	private int minLengthThreshold = 2;//after removing gaps sequences shorter than this are not used for training the model.
	private String pathToDnaseTem  = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/TESTING_SET";

	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		DirichletLikelihoodMaximazation positiveApp = new DirichletLikelihoodMaximazation();
		args = CliTools.configureBean(positiveApp, args);
		positiveApp.run(args);
		
		
		
	}/*main*/
	
	private void run(String []args) throws Exception{
		
		String pathToSerDir = getPathToSerDir();
		String pathToGffs   = getPathToGffFiles();
		String pathToAli    = getAliDir();
		
		
		ArrayList<ArrayList<Double>> trainingDataProfileAsArryList = getProfilesFromTrainingDataSet(pathToSerDir,pathToGffs);
		System.out.println("Number of traning data is: " + trainingDataProfileAsArryList.size());
		double [][] trainingDataAsMatrix   = converTo2DMatrix(trainingDataProfileAsArryList);
		double [] posAlpha = DirichletDist.getMLE(trainingDataAsMatrix, trainingDataProfileAsArryList.size(), 3);
		for(int i = 0; i<posAlpha.length; i++){
			System.out.println(posAlpha[i]);
		}
		
		System.out.println("******************************");
		
		ArrayList<Integer> listOfLengthsInPositiveSet = getListOfLengthsOfAllRecordsInTrainingDataSet(pathToSerDir,pathToGffs);
		ArrayList<ArrayList<Double>> negativeDataSetProfilesAsArryList = getNegativeDataSet(listOfLengthsInPositiveSet,pathToSerDir,pathToGffs);
		double [][]negativeDataAsMatrix = converTo2DMatrix(negativeDataSetProfilesAsArryList);
		double [] negAlpha = DirichletDist.getMLE(negativeDataAsMatrix, negativeDataSetProfilesAsArryList.size(), 3);
		
		//these lines are just for a temp test and must be deleted later on
		/*double [] apos= {0.22525924730402658,0.17246189124006153,0.15114089753516846};
		double [] aneg = {0.15115037756222185,0.1920604160826608,0.20991774344689473};
		calculateScoresAndPrintThemIntoOutputFiles(pathToSerDir, pathToGffs,pathToAli, apos,aneg,posOutputFile,negOutputFile);		
		//end of temp, delet above mentioned line and then un-comment the following line
		 * 
		 */
		
		//calculateScoresAndPrintThemIntoOutputFiles(pathToSerDir, pathToGffs,pathToAli, posAlpha,negAlpha,posOutputFile,negOutputFile);
		calculateScoresAndPrintThemIntoOutputFiles(pathToSerDir, pathToDnaseTem,pathToAli, posAlpha,negAlpha,posOutputFile,negOutputFile);
		
		for(int i = 0;i<negAlpha.length;i++){
			System.out.println(negAlpha[i]);
		}
		

		
	}/*run*/
	
	public void setBlockId(String Id){
		this.blockId = Id;
	}/*setBlockId*/
	
	public void setChr(String chro){
		this.chr = chro;
	}/*setChr*/
	
	public void setPathToGffFiles(String gffDir){
		this.pathToGffFiles = gffDir;
	}/*setPathToGffFiles*/
	
	public void setPathToSerDir(String compDir){
		this.serDir = compDir;
	}/*setPathToCompDir*/
	
	public void setAliDir(String aliDir){
		this.aliDir = aliDir;
	}/*setAliDir*/
	
	public void setMinLengthThreshold(int mlt){
		this.minLengthThreshold = mlt;
	}/*minLengthThreshold*/
	
	
	public int getMinLengthThreshold(){
		return minLengthThreshold;
	}/*getMinLengthThreshold*/
	
	public String getBlockId(){
		return blockId;
	}/*getBlockId*/
	
	public String getPathToGffFiles(){
		return pathToGffFiles;
	}/*getPathToGffFiles*/
	
	public String getPathToSerDir(){
		return serDir;
	}/*getPathToCompDir*/
	
	public String getAliDir(){
		return aliDir;
	}/*getAliDir*/
	
	public String getChr(){
		return chr;
	}/*getChr*/
	
	
	
	public ArrayList<ArrayList<Double>> getProfilesFromOneRecord(GFFRecord gffr, Matrix2D m, Sequence seq){
		ArrayList<ArrayList<Double>> profiles = new ArrayList<ArrayList<Double>>();
		int recordStart = gffr.getStart();
		int recordEnd   = gffr.getEnd();
		IntSet range = IntSets.range(recordStart, recordEnd);
		Symbol gap = seq.getAlphabet().getGapSymbol();
		
		
		for(int i:range){
			if(  (i!= seq.length()) && (seq.symbolAt(i+1) != gap) ){//first condition is because m gets i ranging from 0 up to  m.rows (not m.row itself) which is equal to seq.length 
				ArrayList<Double> oneRow = new ArrayList<Double>();
				//System.out.println("i= " + i);
				for(int j = 0; j<3; j++){
					oneRow.add(m.get(i, j));
				}
				
				profiles.add(oneRow);
			}
		}
		return profiles;
	}/*getProfilesFromOneRecord*/
	
	public ArrayList<ArrayList<Double>> getProfilesFromAllRecordsInAGffFile(File gffFile, Matrix2D serMatrix, Sequence seq)
	throws FileNotFoundException, ParserException, BioException, IOException{
		ArrayList<ArrayList<Double>> profilesForAllRecords = new ArrayList<ArrayList<Double>>();
		int lengthThreshold       = getMinLengthThreshold();
		
		GFFEntrySet records = GFFTools.readGFF(gffFile);
		for(Iterator<?> gffi = records.lineIterator(); gffi.hasNext();){
			Object o = gffi.next();
			if(o instanceof GFFRecord){
				GFFRecord gffr = (GFFRecord) o;
				ArrayList<ArrayList<Double>> profilesForOneRecord = getProfilesFromOneRecord(gffr,serMatrix,seq );
				if(profilesForOneRecord.size() < lengthThreshold){
					System.out.println("a sequence with too short length ignored!");
					continue;
				}
				profilesForAllRecords.addAll(profilesForOneRecord);
			}
		}
		return profilesForAllRecords;
	}/*getProfilesFromAllRecordsInAGffFile*/
	
	public ArrayList<ArrayList<Double>> getProfilesFromTrainingDataSet(String serDirName, String pathToTrainingGffFiles) throws Exception{
		ArrayList<ArrayList<Double>> profilesFromTrainingSet = new ArrayList<ArrayList<Double>>();
		
		File serDir = new File(serDirName);
		Pattern fnPattern = Pattern.compile("(([^_]+)_([0-9]+)_([0-9]+)).ser");
		
		int missingFileIndicator =0;
		for(File aSerFile: serDir.listFiles()){
			if(aSerFile.getName().startsWith(".")){
				continue;
			}
			Matrix2D m = SeriesIO.readSeries(aSerFile);
			Matcher fnMatcher = fnPattern.matcher(aSerFile.getName());
			String blockId;
			if(fnMatcher.matches()){
				blockId = fnMatcher.group(1);
			}
			else{
				System.out.println("Bad file name detected: " + aSerFile.getName());
				continue;
			}
			String aliFileName = getAliDir()+"/"+blockId+".fa";
			File aliFile = new File(aliFileName);
			if(!aliFile.exists()){
				System.out.println("Warning! Could not find alignment file for block " + blockId );
				continue;
			}
			Sequence firstSeq = SeqIOTools.readFastaDNA(new BufferedReader(new FileReader(aliFile))).nextSequence();
			String gffFileName = getPathToGffFiles()+"/"+blockId+".gff";
			File gffFile = new File(gffFileName);
			boolean gffExists = gffFile.exists();
			if(!gffExists){
				missingFileIndicator++;
				continue;
			}
			
			//checkData(m, gffFile);
			ArrayList<ArrayList<Double>> profilesFromOneFile = getProfilesFromAllRecordsInAGffFile(gffFile,m,firstSeq);
			profilesFromTrainingSet.addAll(profilesFromOneFile);
		}
		if(missingFileIndicator >=1){
			System.out.println("Note:  there were " +  missingFileIndicator + " blocks withiout any hit from this feature." );
		}
		
		return profilesFromTrainingSet;		
	}/*getProfilesFromTrainingDataSet*/
	
	public double [][] converTo2DMatrix(ArrayList<ArrayList<Double>> dataArraysOfArrays){//this it to convert Arrays of Arrays of profiles into a 2d matrix because Dirichlet MLE accepts data only as a 2d matrix.
		int numberOfDataset = dataArraysOfArrays.size();
		double  matrixOfData [][] = new double[numberOfDataset][3];
		for(int i = 0; i<numberOfDataset; i++){
			for(int j = 0; j<3; j++){
				double oneElement = dataArraysOfArrays.get(i).get(j);
				matrixOfData[i][j]= oneElement;
			}
		}
		return matrixOfData;
		
		
	}/*converTo2DMatrix*/
	
	public double getLogDensityForOneSequence(ArrayList<ArrayList<Double>> profilesFromOnSequence, double [] alpha, boolean lengthNormalized){
		double density = 0;
		int numberOfCol = profilesFromOnSequence.size();
		for(int i = 0; i<numberOfCol; i++){
			double sum = sum(profilesFromOnSequence.get(i));
			if(sum ==0){
				System.out.println("a position ignored because red green and blue profiles all are equal to zero");
				continue;
			}
			double [] oneCol = new double[3];
			double oneColSum = 0;
			double oneDensity = 0;
			for(int j = 0; j<3; j++){
				oneCol[j] = profilesFromOnSequence.get(i).get(j);
				oneColSum=oneColSum+oneCol[j];
			}
			
			oneDensity = log2(DirichletDist.density(alpha, oneCol));
			
			
			density = density+ oneDensity;
		}
		if(lengthNormalized){
			density = density/numberOfCol;
		}
		return density;
	}/*getLogDensityForOneSequence*/
	
	public double getLogRatioScoreForOneSequence(ArrayList<ArrayList<Double>> profilesFromOnSequence, double [] positiveAlpha, double [] negativeAlpha, boolean lengthNormalized){
		if(profilesFromOnSequence.isEmpty()){
			System.err.println("******** empty profile was passed to function ********");
		}
		double scoreFromPosAlpha = getLogDensityForOneSequence(profilesFromOnSequence, positiveAlpha,lengthNormalized);
		double scoreFromNegAlpha = getLogDensityForOneSequence(profilesFromOnSequence,negativeAlpha,lengthNormalized);
		double logRatioScore  = scoreFromPosAlpha-scoreFromNegAlpha; //note because this scores are in log form, the ration changes to subtraction
		return logRatioScore;
		
	}/*getLogRatioScoreForOneSequence*/
	
	public double log2(double x){
		return(Math.log(x)/Math.log(2));
	}/*log2*/
	
	public ArrayList<Integer> getListOfLengthsOfAllRecordsInTrainingDataSet(String pathToSerFiles, String  pathToTrainingGffFiles)
	throws FileNotFoundException, ParserException, BioException, IOException{
		ArrayList<Integer> listOfAllLengths = new ArrayList<Integer>();
		
		File serDir = new File( pathToSerFiles);
		Pattern fnPattern = Pattern.compile("(([^_]+)_([0-9]+)_([0-9]+)).ser");
		
		for(File aSerFile: serDir.listFiles()){
			if(aSerFile.getName().startsWith(".")){
				continue;
			}
			Matcher fnMatcher = fnPattern.matcher(aSerFile.getName());
			String blockId;
			
			if(fnMatcher.matches()){
				blockId = fnMatcher.group(1);
			}
			else{
				System.out.println("Bad file name detected: " + aSerFile.getName());
				continue;
			}
			String gffFileName = pathToTrainingGffFiles+"/"+blockId+".gff";
			File gffFile = new File(gffFileName);
			boolean gffExists = gffFile.exists();
			if(!gffExists){
				System.out.println("There is no gff file for  block " + blockId);
				continue;
			}
			GFFEntrySet gff = GFFTools.readGFF(gffFile);
			for(Iterator<?> gffi = gff.lineIterator(); gffi.hasNext();){
				Object o = gffi.next();
				if(o instanceof GFFRecord){
					GFFRecord grec = (GFFRecord) o;
					int oneStart = grec.getStart();
					int oneEnd = grec.getEnd();
					int oneLength = oneEnd-oneStart;
					listOfAllLengths.add(oneLength);
				}
			}
				
			
		}
		
		return listOfAllLengths;		
	}/*getListOfLengthsOfAllRecordsInTrainingDataSet*/
	
	public ArrayList<ArrayList<Double>> getNegativeDataSet(ArrayList<Integer> listOfLength, String pathToSer, String pathToPosTrainingGffFiles) throws Exception{
		ArrayList<ArrayList<Double>> negDataSet = new ArrayList<ArrayList<Double>>();
		int lengthThreshold = getMinLengthThreshold();
		
		File serDir = new File(pathToSer);
		File [] listOfSerFiles = serDir.listFiles();
		String aliFileName = null;
		Sequence firstSeq = null;
		
		int numberOfSerFiles  = listOfSerFiles.length;
		Random generator      = new Random();
		
		for(int lenthOfOneRecord: listOfLength){
			Matrix2D m = null;
			int lengthOfThisBlock = -1;
			while(lengthOfThisBlock <= 2*lenthOfOneRecord){
				int aRandomNumber = generator.nextInt( numberOfSerFiles);
				if(!listOfSerFiles[aRandomNumber].getName().startsWith(".")){
					File serFile = listOfSerFiles[aRandomNumber];
					aliFileName = getAliFileNameFromSerFileName(serFile.getName(),getAliDir() );
					m = SeriesIO.readSeries(serFile);
					lengthOfThisBlock = m.rows();
				}
			}
			int diffBetweenLength = lengthOfThisBlock-lenthOfOneRecord;
			int dummyStart = generator.nextInt(diffBetweenLength);
			int dummyEnd   = dummyStart+lenthOfOneRecord;
			String chr     = getChr();
			SimpleGFFRecord r = new SimpleGFFRecord();
			r.setSeqName(chr);
			r.setStart(dummyStart);
			r.setEnd(dummyEnd);
			File aliFile = new File(aliFileName);
			if(!aliFile.exists()){
				System.out.println("ali file  " + aliFile.getName() + " not exists!");
				continue;
			}
			firstSeq = SeqIOTools.readFastaDNA(new BufferedReader(new FileReader(aliFile))).nextSequence();
			ArrayList<ArrayList<Double>> dataForOneRecord = getProfilesFromOneRecord(r,m,firstSeq);
			if(dataForOneRecord.size() < lengthThreshold){
				System.out.println("A too short sequence ignored! ");
				continue;
			}
			negDataSet.addAll(dataForOneRecord);
		}
		
		return negDataSet;
		
	}/*getNegativeDataSet*/
	
	
	public String getAliFileNameFromSerFileName(String serFileName, String pathToAliDir){// it is assumed that ser file is of chr_start_end.ser format
		if(serFileName.startsWith(".")){
			System.err.print("Bad ser file name!");
		}
		String []  serFileSpecifications = serFileName.split("\\.");
		String blockId = serFileSpecifications[0];	
		String aliFileName = pathToAliDir+"/"+blockId+".fa";
		return aliFileName;
	}/*getAliFileNameFromSerFileName*/
	
	public void calculateScoresAndPrintThemIntoOutputFiles(String pathToSer, String pathToTrainingGffFiles,String pathToAliDir, double [] posAlpha, double [] negAlpha, String posOutputFile, String negOutputFile)
	throws Exception{
		ArrayList<Double> posScores = new ArrayList<Double>();
		ArrayList<Double> negScores = new ArrayList<Double>();
		int lengthThreshold   = getMinLengthThreshold();
		
		File gffDir = new File(pathToTrainingGffFiles);
		Pattern fnPattern = Pattern.compile("(([^_]+)_([0-9]+)_([0-9]+)).gff");
		
		for(File aGffFile: gffDir.listFiles()){
			if(aGffFile.getName().startsWith(".")){
				continue;
			}
			Matcher fnMatcher = fnPattern.matcher(aGffFile.getName());
			String blockId;	
			if(fnMatcher.matches()){
				blockId = fnMatcher.group(1);
			}
			else{
				System.out.println("Bad file name detected: " + aGffFile.getName());
				continue;
			}
			String aliFileName = getAliDir()+"/"+blockId+".fa";
			File aliFile       = new File(aliFileName);
			String serFileName = getPathToSerDir()+"/"+blockId+".ser";
			File   serFile     = new File(serFileName);
			if( !aliFile.exists() && !serFile.exists()   ){
				System.out.println("One of the ser file or ali file for block " + blockId + " not exists!");
				continue;
			}
			Matrix2D m = SeriesIO.readSeries(serFile);
			int blockLength = m.rows();
			Sequence firstSeq = SeqIOTools.readFastaDNA(new BufferedReader(new FileReader(aliFile))).nextSequence();
			
			GFFEntrySet records = GFFTools.readGFF(aGffFile);
			for(Iterator<?> gffi = records.lineIterator(); gffi.hasNext();){
				Object o = gffi.next();
				if(o instanceof GFFRecord){
					GFFRecord gffr = (GFFRecord) o;
					int oneStart = gffr.getStart();
					int oneEnd   = gffr.getEnd();
					int oneLength = oneEnd-oneStart;
					ArrayList<ArrayList<Double>> profilesForOneRecord = getProfilesFromOneRecord(gffr,m, firstSeq );
					if(profilesForOneRecord.size()< lengthThreshold){
						System.out.println("A too short sequence ignored!");
						continue;
					}
					double oneLogRatioScore = getLogRatioScoreForOneSequence(profilesForOneRecord, posAlpha,negAlpha,true);
					posScores.add(oneLogRatioScore);
					//negative set:
					Random generator      = new Random();
					int upperBoundForRand = blockLength-oneLength-1;
					int aRandomNumber = generator.nextInt( upperBoundForRand);
					SimpleGFFRecord rDummy = new SimpleGFFRecord();
					rDummy.setStart(aRandomNumber);
					rDummy.setEnd(aRandomNumber+oneLength);
					rDummy.setSeqName("2L");
					ArrayList<ArrayList<Double>> profilesForOneNegRecord = getProfilesFromOneRecord(rDummy,m, firstSeq );
					if(profilesForOneNegRecord.size()< lengthThreshold){
						System.out.println("A too short sequence ignored!");
						continue;
					}
					
					double oneLogRatioScoreFromNeg = getLogRatioScoreForOneSequence(profilesForOneNegRecord,posAlpha,negAlpha,true);
					negScores.add(oneLogRatioScoreFromNeg);
				}
			}
					
		}
		writeArraysIntoAFile(posScores,posOutputFile);
		writeArraysIntoAFile(negScores,negOutputFile);
	}/*calculateScoresAndPrintThemIntoOutputFiles*/
	
	public void writeArraysIntoAFile(ArrayList<Double> someDoubles, String OutputFileName){
		
		try{
			FileWriter fw = new FileWriter(OutputFileName);
			BufferedWriter bw = new BufferedWriter(fw);
			for(double d:someDoubles){
				bw.write(Double.toString(d)+ ",");
			}
			bw.close();
			
		}
		catch(IOException e){
			System.out.println("IOException : " + e);
		}
	}/*writeArraysIntoAFile*/
	
	public double sum(ArrayList<Double> list){
		if (list == null || list.size()<1){
			return 0;
		}
		double sum = 0;
		for(Double d:list){
			sum = sum+d;
		}
		return sum;
	}/*sum*/
	
	public void checkData(Matrix2D m, File aGffFile ){//it seems to me that some of rows in m are consiting only from zeros, this is to check that.
		int numRows = m.rows();
		for(int i = 0; i< numRows; i++){
			double red = m.get(i, 0);
			double green = m.get(i, 1);
			double blue = m.get(i, 2);
			if(red ==0 || green == 0 || blue ==0){
				System.out.println("i = " + i +  "red = " + red +  " green = " +green + " blue = " + blue +"File = " +  aGffFile.getName());
			}
		}
	}

	
}/*DirichletLikelihoodMaximazation*/
