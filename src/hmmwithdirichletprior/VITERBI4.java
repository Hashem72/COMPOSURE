/**
 * 
 */
package hmmwithdirichletprior;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import net.derkholm.nmica.matrix.Matrix2D;

import org.biojava.bio.Annotation;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dp.DP;
import org.biojava.bio.dp.EmissionState;
import org.biojava.bio.dp.IllegalTransitionException;
import org.biojava.bio.dp.MarkovModel;
import org.biojava.bio.dp.ScoreType;
import org.biojava.bio.dp.SimpleEmissionState;
import org.biojava.bio.dp.SimpleMarkovModel;
import org.biojava.bio.dp.State;
import org.biojava.bio.dp.StatePath;
import org.biojava.bio.dp.onehead.SingleDP;
import org.biojava.bio.dp.onehead.SingleDPMatrix;
import org.biojava.bio.program.gff.GFFWriter;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleAlphabet;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.bjv2.util.cli.CliTools;

import series.SeriesIO;
import umontreal.iro.lecuyer.probdistmulti.DirichletDist;


/**
 * @author hk3  25-11-201
 * This is in slight modification of VITERBI3. The key difference is to change it such that to get a ser file as an input and run the code only for this file at a time. Therefore I would be able to 
 * use a job array and make use of all cpus in farm to run it over any number of ser files. Previous VITERBI3  was taking a ser directory and then running over all ser files in that directory and therefore I had
 * memory problem.  
 */
public class VITERBI4 {

	/**
	 * @param args
	 */
	private String serFileName;
	//private String serDir; // = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/Eisenlab_Composure_Output_Model_20051212_With_Window_Length_350_SER_temp";
	private String aliDir; // = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ALIGNMENTS_EISENLAB/2L";
	private String chainDir;  // = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/Chains_From_MarkovModel_ChainFinder_In_Value_Coordinates"; 
	private int blockLengthThreshold; // = 2000; // some blocks are really short, those shorter than this threshold will be ignored
	private double scoreThreshold; 
	
	
	//private String serDir = "/nfs/th_group/hk3/ANALYSIS_EISENLAB/2L/Eisenlab_Composure_Output_Model_20051212_With_Window_Length_350_SER_temp";
	//private String aliDir = "/nfs/th_group/hk3/ALIGNMENTS_EISENLAB/2L";
	//private String chainDir = "/nfs/th_group/hk3/ANALYSIS_EISENLAB/2L/Chains_From_MarkovModel_ChainFinder_In_Value_Coordinates"; 

	
	
	private int    minLengthThreshold = 20; //  adjacent peaks with a distance shorter than this will be connected to each other
	private int    minPeakLengths       = 100; // peaks shorter than this wont be accepted
	
	private  String [] states =               {"M",  "r", "E", "G", "g", "e","R", "J"};//M:Mixed, r:RedFollowedbyMix, E:RedToGreenEdge, G:Green, g:GreenFollowedbyMix, e:GreenToRedEdge, R:red
	private  double [] startProbabilities  =  {0.7,  0.15, 0.02,   0.02,  0.05, 0.02,  0.02, 0.02};//note: order is important, the first element corresponds to "P" and the second element corresponds to "N"
	private  double [][] observationProfiles = null; // note this is a 3 by T matrix
	private  double [][] transitionProbability =
	{ 
			//M,        r,           E,       G,      g,           e,       R,     J  
			{0.9999,     1e-05,        0,       0,      9e-05,     0,       0,       0}, //M: Mix {0.9999,     1e-05,        0,       0,      9e-05,     0,       0,       0}
			{0,        0.997,         0.003,    0,      0,           0,       0,     0}, //r: Red after M {0,        0.997,         0.003,    0,      0,           0,       0,     0} 
			{0,        0,            0,       1,      0,           0,       0,     0}, //E: Red Green Edge {0,        0,            0,       1,      0,           0,       0,     0} 
			{0.00075,   0,            0,       0.998,   0,           0.00075,  0,     0.0005}, //G: Green {0.00075,   0,            0,       0.998,   0,           0.00075,  0,     0.0005}, 
			{0,        0,            0,       0,      0.998,        0.002,    0,     0}, //g: Green after M {0,        0,            0,       0,      0.998,        0.002,    0,     0} 
			{0,        0,            0,       0,      0,           0,       1,     0}, //e: Green to Red Edge {0,        0,            0,       0,      0.998,        0.002,    0,     0} 
			{0.00125,  0,            0.00125, 0,      0,           0,       0.997,  0.0005},  //R : Red {0.00125,  0,            0.00125, 0,      0,           0,       0.997,  0.0005} 
			{0,        0.499995,     0,       0,      0.499995,    0,       0,    0.00001 } //J {0,        0.499995,     0,       0,      0.499995,    0,       0,    0.00001 }
	};

	
	private boolean filterMeaningLessPicks = false;
	
	
	public void setSerFileName(String fn){
		this.serFileName = fn;
	}/*setSerFileName*/
	
	public void setBlockLengthThreshold(int l){
		this.blockLengthThreshold = l;
	}/*setBlockLengthThreshold*/
	
	public void setScoreThreshold(double s){
		this.scoreThreshold = s;
	}/*setScoreThreshold*/
	
	
	public void setAliDir(String ad){
		this.aliDir = ad;
	}/*setAliDir*/
	
	public void setChainDir(String cd){
		this.chainDir = cd;
	}/*setChainDir*/
	
	public void setStates (String[] states){
		this.states = states;
	}/*setStates*/
	
	public void setStartProbabilities(double [] sProb){
		this.startProbabilities = sProb;
	}/*setStartProbabilities*/
	
	public void setTransionProbabilities(double [][] traProb){
		this.transitionProbability = traProb;
	}/*setTransionProbabilities*/
	
	
	public String getSerFileName(){
		return serFileName;
	}/*getSerFileName*/
	
	
	
	public int getBlockLengthThreshold(){
		return blockLengthThreshold;
	}/*getBlockLengthThreshold*/
	
	//public String getSerDir(){
	//	return serDir;
	//}/*getSerDir*/
	
	public String getAliDir(){
		return aliDir;
	}/*getAliDir*/
	
	public String [] getStates(){
		return states;
	}/*getStates*/
	
	public double[] getStartProbabilities(){
		return startProbabilities;
	}/*getStartProbabilities*/
	
	public double [][] getTransitionProbability(){
		return transitionProbability;
	}/*getTransitionProbability*/
	
	public boolean getFilterMeaningLessPicks(){
		return filterMeaningLessPicks;
	}/*getFilterMeaningLessPicks*/
	
	public double getScoreThreshold(){
		return scoreThreshold;
	}/*getScoreThreshold*/
	
	public int getMinLengthThreshold(){
		return minLengthThreshold;
	}/*minLengthThreshold*/
	
	public int getMinPeakLengths(){
		return minPeakLengths;
	}/*getMinPeakLengths*/
	
	public String getChainDir(){
		return chainDir;
	}/*getChainDir*/
	

	
public static void main(String [] args) throws Exception{
	VITERBI4 app = new VITERBI4();
	args = CliTools.configureBean(app, args);
	app.run(args);
	
}/*main*/
	
	
	
	public  void run(String[] args) throws Exception {
		// TODO Auto-generated method stub
		VITERBI4 app = new VITERBI4();
		
		
		//String serDirName     = app.getSerDir();
		app.setSerFileName(this.serFileName);
		app.setAliDir(this.aliDir);
		app.setChainDir(this.chainDir);
		app.setBlockLengthThreshold(this.blockLengthThreshold);
		app.setScoreThreshold(this.scoreThreshold);
		
		String oneSerFileName    = app.getSerFileName();
		String aliDirName     = app.getAliDir();
		String [] stateNames  = app.getStates();
		double [] strProb     = app.getStartProbabilities();
		double [][]tranProb   = app.getTransitionProbability();
		app.checkTransitionMatrixNormality(tranProb);
		
		
		
		double scoreThr                 = app.getScoreThreshold();
		int    minLenForTrhoughs        = app.getMinLengthThreshold();
		int    minLenForPeakDis         = app.getMinPeakLengths();
		int lengthThresholdForBlocks  = app.getBlockLengthThreshold();
		
		double backGroundAlpha           []  = {0.7, 0.8, 1};     //{0.15106438458196936,0.20973711740901407,0.22900628496706044}  {0.8, 0.9, 1}
		double greenStateAlpha           []  = {0.9,2.3,1};       //{0.1713002602078479,11.729438447100865,1.3763120488468623}; {0.9,2.3,1}
		double greenFollowingMixedAlpha  []  = {0.8,2.3,0.9};     //{0.16327302239493097,9.219505194760753,1.8283941303615243} {0.8,2.3,0.9}
		double greenToRedEdgeAlpah       []  = {1.3,1.4,1.2};     //{0.4799631840073914,0.4944929532552769,0.2857067163030646}; {1.3,1.4,1.2}
		double redStateAlpha             []  = {2.5,1,0.8};       //{71.25302335371902,0.201959003605188,0.20299272598532336} {5,1.1,0.9}
		double redToGreenEdgeAlpha       []  = {5,1.1,0.9};     //{0.39728078185044136,0.3710506533397517,0.2490053255191496}{1.4,1.3,1.2}
		//double redFollowingMixedAlpha  []  = {10,0.1,0.1};      //{4.5200496933818615,0.1772646945273629,0.17979061871698085};
		double junctionStateAlpha        []  = {0.7, 0.8, 1};
			
		LinkedHashMap<String, double[]> statesAndDirPar = new LinkedHashMap<String, double[]>();
		//statesAndDirPar.put("P", x);
		
		statesAndDirPar.put("M", backGroundAlpha);
		statesAndDirPar.put("r", redStateAlpha);
		statesAndDirPar.put("E", redToGreenEdgeAlpha);
		statesAndDirPar.put("G", greenStateAlpha);
		statesAndDirPar.put("g", greenFollowingMixedAlpha);
		statesAndDirPar.put("e", greenToRedEdgeAlpah);
		statesAndDirPar.put("R", redStateAlpha);
		statesAndDirPar.put("J", junctionStateAlpha);
		
		
		//for(File aSerFile: serDir.listFiles()){
		//	String oneFileName = aSerFile.getName();
			if(oneSerFileName.startsWith(".")){
				System.err.println("Is this a configuration file?");
			}
			
			File aSerFile  = new File(oneSerFileName); 
			BLOCK block = new BLOCK(aSerFile, aliDirName);
			String blockId = block.blockId;
			int blockLength  = app.getBlockLength(blockId);
			if(blockLength < lengthThresholdForBlocks){
			//System.out.println(" found sequence :" + blockId + " with too short length :" + blockLength);
			System.err.println("got a block shorter than length threshold!");
			}
			Matrix2D m     = block.blockMatrix;
			Sequence seq   = block.blockSeq;
			ArrayList<Integer> rowsWithSumZero = block.missingData;
			SimpleAlphabet observedSeqAlphabet = block.blockObservedSeqAlphabet;
			SimpleSymbolList symbolList        = block.blockSimpleSymbolList;
			
			//System.out.println("seq " + blockId + " is processed");
			
			MarkovModel mm = BLOCK.makeMarkovModel(observedSeqAlphabet, tranProb, strProb, statesAndDirPar, "dirichletMM");
			 DP dp    = new SingleDP(mm);	
			 
			 
			 
			SymbolList [] symList = {symbolList};
			
			StatePath viterbiPath = dp.viterbi(symList, ScoreType.PROBABILITY);
			SymbolList symbolsInViterbi             = viterbiPath.symbolListForLabel(StatePath.STATES);
			ArrayList<String> viterbiPathAsAnArrayList           = app.symbolListToArrayList(symbolsInViterbi);
			System.out.println("veterbi path length for " + blockId + " is " +viterbiPath.length() );
			
			
			ArrayList<String> ViterbiPath = new ArrayList<String>();
			 
			
			
			
			for(int i=1; i<= symbolsInViterbi.length(); i++){
				Symbol oneSym = symbolsInViterbi.symbolAt(i);
				ViterbiPath.add(oneSym.getName());
				//System.out.print(oneSym.getName());
			}
			
			System.out.println();
			boolean isFiltering  = app.getFilterMeaningLessPicks();
			app.filterMeaningLessPicks(ViterbiPath, isFiltering);
			
			ArrayList<Double> viterbiTranslatedToNumbersAndMappedToAlingment = app.mapViterbiPathToAlignment(ViterbiPath, seq, rowsWithSumZero);
			String outputFileName  = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/ViterbiPaths/VITERBI2_"+blockId+".txt";
			app.printoutViterbiPath(viterbiTranslatedToNumbersAndMappedToAlingment, outputFileName);
			
		
			

	        SingleDPMatrix forwardMatrix = (SingleDPMatrix) dp.forwardMatrix(
	                new SymbolList[] {symbolList},
	                ScoreType.PROBABILITY
	            );
	        
	        double score = forwardMatrix.getScore();
	        //System.err.printf("Forward: %g%n", score);
	        System.out.println();
	        SingleDPMatrix backwardMatrix = (SingleDPMatrix) dp.backwardMatrix(
	        		new SymbolList[] {symbolList}, 
	        		ScoreType.PROBABILITY
	        		);
	        
	        ArrayList<Double> posteriorScores = app.getPosteriorDecodingScores(forwardMatrix,backwardMatrix);
	        //posteriorScores = app.fillInShortTroughs(posteriorScores, scoreThr, minLenForTrhoughs);
	        //posteriorScores = app.filterOutShortPeaks(posteriorScores, scoreThr, minLenForPeakDis);
	       // posteriorScores = app.getEnhancerRegions(posteriorScores, scoreThr);
	        ArrayList<Double> PDScoresMapedToalignment = app.mapPDScoresToAlignment(posteriorScores, seq, rowsWithSumZero);
	        String posteriorOutputFileName  = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/ViterbiPaths/Posterior_"+blockId+".txt";
	        app.printoutViterbiPath(PDScoresMapedToalignment, posteriorOutputFileName);
	        
	        //writte chains into gff files
	        
	        String chainDir = app.getChainDir();
	        app.wirteChainsIntoAGffFile(PDScoresMapedToalignment, blockId, chainDir, scoreThr);
	        
	       
		




	}/*run*/

	
	public static class BLOCK{
		public String      blockId;
		public Matrix2D    blockMatrix; //this either comp or ser matrix
		public Sequence    blockSeq;	//sequence from the alignment file
		public SimpleAlphabet blockObservedSeqAlphabet;
		public SimpleSymbolList   blockSimpleSymbolList;
		public ArrayList<Integer> missingData;
		
		public BLOCK(File serFile, String pathToAli) throws Exception{
			Pattern fnPattern  = Pattern.compile("(([^_]+)_([0-9]+)_([0-9]+)).ser");
			String serFileName = serFile.getName();
			Matcher fnMatcher  = fnPattern.matcher(serFileName);

			String blockId;
			if(fnMatcher.matches()){
				blockId = fnMatcher.group(1);
			}else{
				blockId = null;
				System.out.println("Bad file name!");
				System.exit(0);
			}
			this.blockId = blockId;
			
			Matrix2D m = SeriesIO.readSeries(serFile);
			this.missingData = getRowsThatSumupToZero(m);
			
			this.blockMatrix = m;
			String aliFileName = pathToAli+"/"+blockId+".fa";
			File aliFile       = new File(aliFileName);
			if(!aliFile.exists()){
				System.err.println(" Could not find alignment file for block " + blockId );
			}
			Sequence firstSeq  = SeqIOTools.readFastaDNA(new BufferedReader(new FileReader(aliFile))).nextSequence();
			this.blockSeq = firstSeq;
			
			//System.out.println("seq length is " + firstSeq.length());
			
			LinkedHashMap<SimpleAlphabet, SimpleSymbolList> alphabetAndSybmolList = getAlphabetAndSimpleSymbolList(m, firstSeq);
			int numberOfAlphabetAndSymbols = alphabetAndSybmolList.size();
			if(numberOfAlphabetAndSymbols>1){
				System.err.println("Expected one set of alphabet but obtainded more!");
			}
			SimpleAlphabet alphabet = new SimpleAlphabet();
			SimpleSymbolList  symbolList = null;
			for(Map.Entry<SimpleAlphabet, SimpleSymbolList> me:alphabetAndSybmolList.entrySet()){
				alphabet= me.getKey();
				symbolList = me.getValue();
			}
			this.blockObservedSeqAlphabet = alphabet;
			this.blockSimpleSymbolList = symbolList;
		}/*BLOCK-Constructor*/
		
		public ArrayList<Integer> getRowsThatSumupToZero(Matrix2D aMatrix){
			ArrayList<Integer> result = new ArrayList<Integer>();
			int numberOfRows = aMatrix.rows();
			for(int i=0; i<numberOfRows;i++){
				double redValue = aMatrix.get(i, 0);
				double greenValue = aMatrix.get(i, 1);
				double blueValue  = aMatrix.get(i, 2);
				double rowSum   = redValue+greenValue+blueValue;
				if(rowSum ==0.0){
					result.add(i);
				}
			}
			return result;
			
		}/*getRowsThatSumupToZero*/
		
		

		
		public LinkedHashMap<SimpleAlphabet, SimpleSymbolList> getAlphabetAndSimpleSymbolList(Matrix2D m, Sequence sequence) throws IllegalSymbolException{
			LinkedHashMap<SimpleAlphabet, SimpleSymbolList> alphabetAndSymbolList = new LinkedHashMap<SimpleAlphabet, SimpleSymbolList>();
			
			SimpleAlphabet alphabet = new SimpleAlphabet();
			List<AtomicSymbol> listOfSymbols = new ArrayList<AtomicSymbol>(); 
			alphabet.setName("ObservedSequenceAlphabet");
			int numberofRows = m.rows();
			int seqLength    = sequence.length();
			if(numberofRows != seqLength){
				System.err.print("It was assumed your sequence has a length equal to  the number of rows of the matrix, but found a case that is not true!");
			}
			
			Symbol gap = sequence.getAlphabet().getGapSymbol();
			
			for(int i = 0; i< numberofRows; i++){
				List<Symbol> oneListOfSymbol = new ArrayList<Symbol>(3);
				 
				
				//red is match, green is flanking and blue is background
				double redValue = m.get(i, 0);
				double greenValue = m.get(i, 1);
				double blueValue  = m.get(i, 2);
				double onesum = redValue+greenValue+blueValue;
				boolean isGap = false;
				isGap = ( sequence.symbolAt(i+1) == gap);
				if(onesum == 0){//sum of this value is supposed to be one, but for some rows it sums up to zero, this is to ignore those up to time we found out why these naouthy rows sums up to zero!
					//continue;
					//note these three lines is only a dummy solution for positions where red, green and blue are summed up to zero!
					redValue = 0.3333;
					greenValue = 0.3333;
					blueValue  = 1 - (redValue +greenValue);
				}
				
				
				//make one triplet symbol from three symbols
				Symbol redSymbol  = AlphabetManager.createSymbol(Double.toString(redValue));
				Symbol greenSymol = AlphabetManager.createSymbol(Double.toString(greenValue));
				Symbol blueSymbol = AlphabetManager.createSymbol(Double.toString(blueValue));
				oneListOfSymbol.add(redSymbol);
				oneListOfSymbol.add(greenSymol);
				oneListOfSymbol.add(blueSymbol);
				
				//now create symbol and add it to alphabet
				AtomicSymbol oneSym = (AtomicSymbol) AlphabetManager.createSymbol(Annotation.EMPTY_ANNOTATION, oneListOfSymbol, alphabet);
				alphabet.addSymbol(oneSym);
				listOfSymbols.add(oneSym);
				
			}
			
			
			SimpleSymbolList ssl = new SimpleSymbolList(alphabet,listOfSymbols);
			
			alphabetAndSymbolList.put(alphabet, ssl);
			return alphabetAndSymbolList;
		}/*getAlphabetAndSimpleSymbolList*/
		
		
		
		
		private static MarkovModel makeMarkovModel(
				SimpleAlphabet alphabet,  double [][] transitionMatrix, double [] startProbabilities, LinkedHashMap<String , 
				double []> statesAndCorresponingDirichletParameters, String modelName ) throws Exception{
			
			
			int [] advance = { 1 };
			SimpleAlphabet anAlphabet = alphabet;
			
			int numberOfStates = statesAndCorresponingDirichletParameters.size();
			ArrayList<String > stateNames = new ArrayList<String>();
			ArrayList<double []> arraysOfDirichletParameters = new ArrayList<double []>();
			
			for(Map.Entry<String, double[]> me:statesAndCorresponingDirichletParameters.entrySet() ){
				double oneDirichletPar [] = me.getValue();
				arraysOfDirichletParameters.add(oneDirichletPar);
				String oneState           = me.getKey();
				stateNames.add(oneState);
			}
						
			//Distribution initiation 
			Distribution 	[] dists          = new Distribution[numberOfStates];
			EmissionState   [] emissionStates = new SimpleEmissionState[numberOfStates];
			
			try{
				for(int i = 0; i< numberOfStates;i++){
					dists[i] = DistributionFactory.DEFAULT.createDistribution(anAlphabet);
					String oneState = stateNames.get(i);
					emissionStates[i] =  new SimpleEmissionState(oneState, Annotation.EMPTY_ANNOTATION,advance,dists[i] );
				}
			}
			catch(Exception e){
				throw new Exception("Can't create distributions!");
			}
			
			//create HMM
			SimpleMarkovModel mm  = new SimpleMarkovModel(1, anAlphabet, modelName );
			
			//add states to the model
			for(State s:emissionStates ){
				try{
					mm.addState(s);
				}
				catch(Exception e){
					throw new Exception("Can't add states to model!");
				}
			}
			
			
			//create transitions
			try{
				State magic = mm.magicalState();
				for(State i:emissionStates ){
					mm.createTransition(magic, i);
					mm.createTransition(i, magic);
					for(State j: emissionStates){
						mm.createTransition(i, j);
					}
				}
			}
			
			
			catch(Exception e){
				throw new Exception("Can't create transitions!");
			}
			

			//set up emission scores			
			for(Iterator<?> i = anAlphabet.iterator(); i.hasNext();){
				AtomicSymbol oneSym = (AtomicSymbol) i.next();
				//System.out.println(oneSym.getName());
				double [] symbolsInThisSymbolAsArrayOfDoubles = makeArrayOfDoublesFromASymbol(oneSym);
				for(int d =0 ; d< dists.length; d++){
					double dirichletPar [] = arraysOfDirichletParameters.get(d);
					double oneDensity         = DirichletDist.density(dirichletPar, symbolsInThisSymbolAsArrayOfDoubles);
					
					//System.out.print( oneDensity + "  ");
					
					oneDensity = oneDensity +100; //because some of these scores are zero, this one unit is to be able to take log! 100
					
					dists[d].setWeight(oneSym,oneDensity );
					
				}
				//System.out.println();
			}
			
				//set  transition scores
				Distribution transDist;
				
				//magical to others
				transDist = mm.getWeights(mm.magicalState());
				for(int i=0; i<emissionStates.length; i++){
					transDist.setWeight(emissionStates[i], startProbabilities[i]);					
				}
				
				//from each state to others
				for(int i =0; i<emissionStates.length;i++){
					transDist = mm.getWeights(emissionStates[i]);
					transDist.setWeight(mm.magicalState(), 0.0001);
						transDist.setWeight(emissionStates[i], transitionMatrix[i][i]);
					for(int j =0; j<emissionStates.length; j++){
						if(i !=j){
							transDist.setWeight(emissionStates[j], transitionMatrix[i][j]);
						}
					}
				}			
			//System.out.println("One HMM set up!");
			return mm;			
		}/*makeMarkovModel*/
		
		private static double []  makeArrayOfDoublesFromASymbol(AtomicSymbol as){
			List symbols = as.getSymbols();
			int numberOfSym = symbols.size();
			double x[] = new double[numberOfSym];
			if(numberOfSym != 3){
				System.err.println(" Number of symbols is expected to be 3 but it is not here!");
			}
			//double sum = 0;
			for(int i = 0; i<symbols.size(); i++){
		        Symbol sym = (Symbol)symbols.get(i);
		        double doubeValueOfSym =  Double.parseDouble(sym.getName());
		     //   sum = sum + Double.parseDouble(sym.getName());
		        x[i] = doubeValueOfSym;
			}
			return x;
		} /*makeArrayOfDoublesFromASymbol*/
		
		public double getSumOfElementsInADoubleArray(double []array){
			double sum = 0;
			for(int i = 0; i< array.length; i++){
				sum = sum + array[i];
			}
			return sum;
		}/*getSumOfElementsInADoubleArray*/		
	}/*BLOCK_Class*/
	
	
	public void checkTransitionMatrixNormality( double [][]traMatrix){
		int numberOfRows = traMatrix.length;
		int numberOfCol   = traMatrix[0].length;
		for(int r = 0; r<numberOfRows; r++){
			double oneRowSum =0;
			BigDecimal bdOneRowSum  = new BigDecimal(oneRowSum);
			for(int c=0; c<numberOfCol; c++){
				BigDecimal bdOneElement = new BigDecimal(traMatrix[r][c]);
				bdOneRowSum = bdOneRowSum.add(bdOneElement);
			}
			oneRowSum = bdOneRowSum.doubleValue();
			if(oneRowSum != 1){
				System.out.println("At row "+ r + " Trainsition probabilities are not summing up to " + oneRowSum  + " whereas assumed to be one!" );
				//System.exit(0);
			}
		}
	}/*checkTransitionMatrixNormality*/
	
	
	public int getBlockLength(String blockID){
		String blockSpecifications[] = blockID.split("_");
		int blockLength = Integer.parseInt(blockSpecifications[2])- Integer.parseInt(blockSpecifications[1]);
		return blockLength;
		
	}/*getBlockLength*/

	
	public ArrayList<Double> mapViterbiPathToAlignment(ArrayList<String> viterbiPath, Sequence seq, ArrayList<Integer> indicesOfRowsWithSumZero){
		ArrayList<Double> result = new ArrayList<Double>();
		
		int seqLength = seq.length();
		int nonGapPosition = 0;
		Symbol gap = seq.getAlphabet().getGapSymbol();
		for(int i =0; i< seqLength; i++){
			double oneScore = Double.MIN_VALUE;
			boolean isSumZero = false;
			isSumZero = (indicesOfRowsWithSumZero.contains(i));
			boolean isGap = false;
			isGap   = (seq.symbolAt(i+1) == gap); //note: these +1 in indices are because seq and viterbi path start from one and not zero!
			if(isSumZero){
				oneScore =0.0;
			}
			else {
				if(viterbiPath.get(nonGapPosition).equals("M")){
					oneScore  = 0.0;
				}
				else if(viterbiPath.get(nonGapPosition).equals("R")){
					oneScore = 2.0;
				}
				else if(viterbiPath.get(nonGapPosition).equals("r")){
					oneScore = 2.0;
				}
				else if(viterbiPath.get(nonGapPosition).equals("G")){
					oneScore = 2.0;
				}
				else if(viterbiPath.get(nonGapPosition).equals("g")){
					oneScore = 2.0;
				}
				else if(viterbiPath.get(nonGapPosition).equals("E")){
					oneScore = 1.0;
				}
				else if(viterbiPath.get(nonGapPosition).equals("e")){
					oneScore = 1.0;
				}
				else if(viterbiPath.get(nonGapPosition).equals("J")){
					oneScore =1.5;
				}
				else{
					System.err.println("Unknown charecter detected as a state name!");
				}
				
				nonGapPosition++;
			}
			result.add(oneScore);
		}
		return result;
		
	}/*mapViterbiPathToAlignment*/
	
	
	public ArrayList<Double> mapPosteriorPathToAlignment(ArrayList<String> posteriorPath, Sequence seq, ArrayList<Integer> indicesOfRowsWithSumZero){
		ArrayList<Double> result = new ArrayList<Double>();
		
		int seqLength = seq.length();
		int nonGapPosition = 0;
		Symbol gap = seq.getAlphabet().getGapSymbol();
		for(int i=0; i< seqLength; i++){
			double oneScore = Double.MIN_VALUE;
			boolean isSumZero = false;
			isSumZero = (indicesOfRowsWithSumZero.contains(i));
			boolean isGap = false;
			isGap   = (seq.symbolAt(i+1) == gap); //note: these +1 in indices are because seq  starts from one and not zero!
			if(isGap || isSumZero){
				oneScore =-1.0;
			}
			else
			{
				if(posteriorPath.get(nonGapPosition).equals("M") ){
					oneScore =0.0;
				}
				else if (posteriorPath.get(nonGapPosition).equals("R") || 
						posteriorPath.get(nonGapPosition).equals("G")  || 
						posteriorPath.get(nonGapPosition).equals("r")  || 
						posteriorPath.get(nonGapPosition).equals("g")  || 
						posteriorPath.get(nonGapPosition).equals("E")  || 
						posteriorPath.get(nonGapPosition).equals("e")  || 
						posteriorPath.get(nonGapPosition).equals("J") )
				{
			
					oneScore =1.0;
				}
				else
				{
					System.err.println("Unknown charecter detected as a state name!");
				}
				nonGapPosition++;		
			}
			result.add(oneScore);
		}
		return result;		
	}/*mapPosteriorPathToAlignment*/
	
	public static double getWeight(String s){
		double result = Double.MIN_VALUE; 
		if(s.equals("M")){
			result =0.0;
		}
/*		else if(s.equals("R") || s.equals("G") || s.equals("J") || s.equals("!-1")){
			result =1;
		}
		else if(s.equals("r") || s.equals("g")){
			result = 0.5;
		}
		else if (s.equals("e") || s.equals("E")){
			result = 0.9;
		}
*/		
		else{
			result =1.0;
		}
		return result;
	}/*getWeight*/
	
	public ArrayList<Double> mapPDScoresToAlignment(ArrayList<Double> pdScores, Sequence seq, ArrayList<Integer> indicesOfRowsWithSumZero){
		ArrayList<Double> result = new ArrayList<Double>();
		
		int seqLength  = seq.length();
		int nonGapPosition  =0;
		Symbol gap = seq.getAlphabet().getGapSymbol();
		for(int i =0; i<seqLength; i++){
			double oneScore = Double.MIN_VALUE;
			boolean isSumZero = false;
			isSumZero = (indicesOfRowsWithSumZero.contains(i));
			boolean isGap = false;
			isGap   = (seq.symbolAt(i+1) == gap); //note: these +1 in indices are because seq  starts from one and not zero!
			if(  isSumZero){
				oneScore =0.0;
			}
			else
			{
				oneScore = pdScores.get(nonGapPosition);
				nonGapPosition++;
			}
			result.add(oneScore);
		}		
		return result;
	}/*mapPDScoresToAlignment*/
	
	
	
	public ArrayList<String> symbolListToArrayList(SymbolList sl){
		ArrayList<String> result = new ArrayList<String>();
		int numberOfSymbols = sl.length();
		for(int i =1; i<= numberOfSymbols; i++ ){
			String onestate = sl.symbolAt(i).getName();
			result.add(onestate);
		}
		return result;
		
	}/*symbolListToArrayList*/
	
	public void printoutViterbiPath( ArrayList<Double> viterbiPath, String outputFileName){
		try{
			FileOutputStream fo = new FileOutputStream(outputFileName);
			PrintStream ps      = new PrintStream(fo);
			for(int i = 0; i< viterbiPath.size(); i++){
				String oneLine = Integer.toString(i+1) + "\t" + Double.toString(viterbiPath.get(i));
				ps.println(oneLine);
			}
			ps.close();
			fo.close();
		}catch(Exception e){
			
		}
	}/*printoutViterbiPath*/
	
	
	public ArrayList<String> getPosteriorDecodingPath(SingleDPMatrix fMatrix, SingleDPMatrix bMatrix){
		ArrayList<String> pdPath = new ArrayList<String>();
		State states[]   = fMatrix.states();
		double fScore    = fMatrix.getScore();
		double fMatrixScores [][] = fMatrix.scores;
		double bMatrixScores [][] = bMatrix.scores;
		
		int numberOfObservations = fMatrixScores.length;
		int numberOfStates       = states.length;
		for(int i=1; i<numberOfObservations-1; i++){
			ArrayList<Double> probOfStates = new ArrayList<Double>();
			for(int s =0; s<numberOfStates; s++){
				double probOfOneState = Math.exp(fMatrixScores[i][s]+ bMatrixScores[i][s]-fScore);
				probOfStates.add(probOfOneState);
			}
			Object max_obj    =  Collections.max(probOfStates);
			int maxValueIndex     =  probOfStates.indexOf(max_obj);
			State topScoredState = states[maxValueIndex];
			pdPath.add(topScoredState.getName());
			//System.out.print(topScoredState.getName());
		}
		//System.out.println();
		//System.out.println("\n number of rows is " + numberOfObservations);
		//System.out.println("\n length of posterior decoding path is " + pdPath.size());
		return pdPath;
		
	}/*getPosteriorDecodingPath*/
	
	public ArrayList<Double> getPosteriorDecodingScores(SingleDPMatrix fMatrix, SingleDPMatrix bMatrix){
		ArrayList<Double> pdScores = new ArrayList<Double>();
		State states[]             = fMatrix.states();
		double fScore              =fMatrix.getScore();
		double fMatrixScores[][]   = fMatrix.scores;
		double bMatrixScores[][]   = bMatrix.scores;
		
		int numberOfObservations  = fMatrixScores.length;
		int numberOfStates        = states.length;
		
		for(int i  =1; i<numberOfObservations-1; i++){
			double sum =0;
			double g = Double.MIN_VALUE;			
			for(int s=0; s<numberOfStates; s++){
				  g= getWeight(states[s].getName());
/*				
				if( states[s].getName().equals("M")){
					g=0;
				}
				else
				{
					g=1;
				}
				
*/				
				double probOfOneState = Math.exp(fMatrixScores[i][s]+ bMatrixScores[i][s]-fScore);
				sum                   = sum + (probOfOneState*g);	
			}	
			pdScores.add(sum);
		}	
		return pdScores;
	}/*getPosteriorDecodingScores*/
	
	public void filterMeaningLessPicks(ArrayList<String> path, boolean needsFiltering){
		
		if(needsFiltering){
			for(int i =4; i< path.size()-3; i++){
				
				if(path.get(i).equals("E")  && path.get(i+1).equals("G") && path.get(i+2).equals("M")){
					int endPosition = i;
					int startPosition =i;
					
					while(path.get(startPosition-1).equals("r") &&  (startPosition >1) ){
						startPosition--;
					}/*while*/
					//System.out.println("startPosi = "+ startPosition + " endPosi = "+ endPosition);
					if(path.get(startPosition-1).equals("M")){
						for(int p = startPosition-1; p <= endPosition+1; p++){
							path.set(p, "M");
						}
					}
					
				}/*if*/
				
				else if(path.get(i).equals("e")  && path.get(i+1).equals("R") && path.get(i+2).equals("M")){
					int endPosition = i;
					int startPosition =i;
					while(path.get(startPosition-1).equals("g") &&  (startPosition >1)){
						startPosition--;
					}/*while*/
					if(path.get(startPosition-1).equals("M")){
						for(int p = startPosition-1; p <= endPosition+1; p++){
							path.set(p, "M");
						}
					}
				}/*if*/
				
				else
				{
					continue;
				}
				
				
				
			}/*for i*/
		}/*if*/
	}/*filterMeaningLessPicks*/

	public ArrayList<Double> fillInShortTroughs(ArrayList<Double> scores, double threshold, int lengthThreshold){
		ArrayList<Double> result = new ArrayList<Double>();
		result = scores;
		int troughStart = 0;
		int troughEnd   = 0;
		double artificialScore = threshold + ( 1-threshold)/2;
		for(int i = 0; i<result.size()-1; i++){
			double thisScore = result.get(i);
			double nextScore = result.get(i+1);
			if(thisScore >threshold && nextScore <= threshold){
				troughStart = i+1;
				troughEnd  = i+1;
				while(result.get( troughEnd) <= threshold &&  troughEnd< result.size()-1 ){
					troughEnd++;
				}/*while*/
				troughEnd = troughEnd-1;
				int thisTroughLength = troughEnd-troughStart;
				if(thisTroughLength < lengthThreshold){
					//System.out.println("start =" + troughStart+ " end = " + troughEnd + " length is " + thisTroughLength);
					for(int x=troughStart; x<= troughEnd; x++ ){				
						result.set(x, artificialScore);
					}/*for*/
				}/*if*/
				i=troughEnd;
			}/*if*/
		}/*for*/
		return result;
	}/*filShortTroughs*/
	
	public ArrayList<Double> filterOutShortPeaks(ArrayList<Double> scores, double threshold, int minPeakLength){
		ArrayList<Double> result = new ArrayList<Double>();
		result                   = scores;
		
		int peakStart = 0;
		int peakEnd   = 0;
		double artificialScore = threshold - ( 1-threshold)/2;
		for(int i =0; i<result.size()-1; i++){
			double thisScore = result.get(i);
			double nextScore = result.get(i+1);
			if(thisScore < threshold && nextScore >= threshold){
				peakStart = i+1;
				peakEnd   = i+1;
				while( result.get(peakEnd) >= threshold && peakEnd < result.size()-1 ){
					peakEnd++;
				}
				peakEnd = peakEnd-1;
				int peakLength = peakEnd- peakStart;
				if(peakLength < minPeakLength ){
					//System.out.println("start =" + peakStart+ " end = " + peakEnd + " length is " + peakLength);
					for(int x= peakStart; x<= peakEnd; x++){
						result.set(x, artificialScore);
					}
				}
				i= peakEnd;
			}
		}
		
		return result;
	}/*filterOutShortPeaks*/
	
	
	public void wirteChainsIntoAGffFile(ArrayList<Double> PDScores, String blockId, String pathToChainsDir, double scoreThreshold) throws IOException{
		String blockSpecfications [] = blockId.split("_");
		String Chr                   = blockSpecfications[0];
		String gffoutputFileName = pathToChainsDir+"/"+blockId+".gff";
		PrintWriter pw = new PrintWriter(new FileWriter(gffoutputFileName));
		final GFFWriter gffw = new GFFWriter(pw);
		int numberOfPDScores = PDScores.size();
		int oneStart = 0;
		int oneEnd   = 0;
		int dummyPos = 0;
		
		
		for(int i =0; i<numberOfPDScores-1; i++){ 
			double thisScore = PDScores.get(i);
			double nextScore = PDScores.get(i+1);
			if(PDScores.get(0)>= scoreThreshold ){
				thisScore = thisScore-0.001;
			}
			if(thisScore  < scoreThreshold   &&  nextScore >= scoreThreshold){
				oneStart = i+1;
				dummyPos = i+1;
				double sum =0;
				while( PDScores.get(dummyPos) >= scoreThreshold && dummyPos< numberOfPDScores-1){
					sum = sum + PDScores.get(dummyPos+1);
					dummyPos++;
				}/*while*/
				//System.out.println("                      start = " + oneStart + " end = " + dummyPos);
				oneEnd = dummyPos;
				int oneChainLength = oneEnd- oneStart;
				SimpleGFFRecord oneRecord = new SimpleGFFRecord();
				oneRecord.setStart(oneStart+1);
				oneRecord.setEnd(oneEnd);
				oneRecord.setSeqName(Chr);
				oneRecord.setScore(sum/oneChainLength);//oneRecord.setScore(sum/oneChainLength+1);
				oneRecord.setSource("composure");
				oneRecord.setFeature("chainFromPD");
				
				
				
				if(oneChainLength> 5){
					//System.out.println(Chr + "\t" + oneStart + "\t" + oneEnd + "\t" +  oneChainLength);
					gffw.recordLine(oneRecord);
				}
				i = oneEnd-1;
				
			}/*if*/
		}/*for*/
		
		pw.flush();
	}/*wirteChainsIntoAGffFile*/
	
	
	public ArrayList<Double> getEnhancerRegions(ArrayList<Double> scores, double thershold){
		ArrayList<Double> result = new ArrayList<Double>();
		for(int i =0; i< scores.size(); i++){
			double thisScore = scores.get(i);
			if(thisScore > thershold){
				result.add(1.0);
			}
			else
			{
				result.add(0.0);
			}
			
		}
		return result;
	}/*getEnhancerRegions*/


	
	
}/*VITERBI4*/
