/**
 * 
 */
package hmmwithdirichletprior;

import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import net.derkholm.nmica.matrix.Matrix2D;

import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.Symbol;

import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.GFFWriter;
import org.biojava.bio.program.gff.SimpleGFFRecord;


import series.SeriesIO;
import umontreal.iro.lecuyer.probdistmulti.DirichletDist;

/**
 * @author hk3 04-07-2011
 * this is to model a HMM model with a dirichlet prior. Here I am going to implement the viterbi algorithm to fit my data.
 * 	PLEASE NOTE: Here in this version of VITERBI you see an implementation of the viterbi algorithm that has been done entirely by MYSELF. BUT you may see another version of VITERBI in this package which is implementation
 * viterbi algorithm based on BIOJAVA dp.
 */
public class VITERBI {

	/**
	 * @param args
	 */
	private String serDir = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/Eisenlab_Composure_Output_Model_20051212_With_Window_Length_350_SER_temp";
	private String aliDir = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ALIGNMENTS_EISENLAB/2L";
	private String chainDir = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/Chains_From_MarkovModel_ChainFinder_In_Value_Coordinates";
	
	
	private  String [] states =               {"M",  "r", "E", "G", "g", "e","R", "J"};//M:Mixed, r:RedFollowedbyMix, E:RedToGreenEdge, G:Green, g:GreenFollowedbyMix, e:GreenToRedEdge, R:red
	private  double [] startProbabilities  =  {0.8,  0.15, 0.0,   0.0,  0.05, 0.0,  0.0, 0.0};//note: order is important, the first element corresponds to "P" and the second element corresponds to "N"
	private  double [][] observationProfiles = null; // note this is a 3 by T matrix
	private  double [][] transitionProbability =/*{ 
		//M,    r,     E,     G,     g,     e,    R, 
		{0.7,   0.2,   0,      0,    0.1,   0,    0}, //M: Mix {0.7,  0.2,  0,     0,    0.1,  0,   0  } 
		{0,     0.7,   0.3,    0,    0,     0,    0}, //r: Red after M {0   , 0.8 , 0.2 , 0   , 0   , 0   , 0  } 
		{0,     0,     0,      1,    0,     0,    0}, //E: Red Green Edge {0   , 0   , 0   , 1   , 0   , 0   , 0  } 
		{0.1,   0,     0,      0.6,  0,    0.3, 0}, //G: Green {0.1,  0,    0,      0.6,  0,     0.3, 0  } 
		{0,     0,     0,      0,    0.6,   0.4,  0}, //g: Green after M {0   , 0   , 0   , 0   , 0.8 , 0.2 , 0  } 
		{0,     0,     0,      0,    0,     0,    1}, //e: Green to Red Edge {0   , 0   , 0   , 0   , 0   , 0   , 1  } 
		{0.1,   0,     0.2,   0,    0,     0,    0.7}  //R : Red {0.05 , 0   ,0.05  , 0   , 0   , 0   , 0.9} 
};*/


	
		//this is the transition matrix for eight state model to check if works better
	
	{ 
			//M,    r,     E,     G,     g,     e,    R,  J  
			{0.7,   0.2,   0,      0,    0.1,   0,    0, 0}, //M: Mix {0.7,  0.2,  0,     0,    0.1,  0,   0  } 
			{0,     0.7,   0.3,    0,    0,     0,    0, 0 }, //r: Red after M {0   , 0.8 , 0.2 , 0   , 0   , 0   , 0  } 
			{0,     0,     0,      1,    0,     0,    0, 0}, //E: Red Green Edge {0   , 0   , 0   , 1   , 0   , 0   , 0  } 
			{0.1,   0,     0,      0.64,  0,     0.25, 0, 0.9}, //G: Green {0.1,  0,    0,      0.6,  0,     0.3, 0  } 
			{0,     0,     0,      0,    0.6,   0.4,  0,  0}, //g: Green after M {0   , 0   , 0   , 0   , 0.8 , 0.2 , 0  } 
			{0,     0,     0,      0,    0,     0,    1,  0}, //e: Green to Red Edge {0   , 0   , 0   , 0   , 0   , 0   , 1  } 
			{0.1,   0,     0.15,   0,    0,     0,    0.74, 0.9},  //R : Red {0.05 , 0   ,0.05  , 0   , 0   , 0   , 0.9} 
			{0,       0.499995, 0,     0,    0.499995,  0,    0,    0.00001 } //J
	};
	
	
	private  double [][] emissionProbability =  null;
	private  double [] alpha = null; // this vector is dirichlet parameter
	
	//some post processing parameters
	private boolean fileterMeaninglessPeaks = false;
	private boolean conncetAdjacentPeaks    =   false;
	private int     minDistanceBetweenTwoPeaks = 10; // two peaks closer than this will be connected to each other
	private int     minChainLength             = 10; // discard shorter chains 
	
	//some debugging parameters
	private static boolean printoutViterviScores = true;
	

	
	public void setMinChainLength(int l){
		this.minChainLength = l;
	}/*setMinChainLength*/
	
	public void setChainDir(String cd){
		this.chainDir = cd;
	}/*setChainDir*/
	
	public void setPrintoutViterviScores(boolean b){
		VITERBI.printoutViterviScores =b;
	}/*setPrintoutViterviScores*/
	
	public void setFileterMeaninglessPeaks(boolean b){
		this.fileterMeaninglessPeaks = b;
	}/*setFilterUnsurroundedPeaks*/
	
	public void setConncetAdjacentPeaks(boolean b){
		this.conncetAdjacentPeaks = b;
	}/*setConncetAdjacentPeaks*/
	
	public void setMinDistanceBetweenTwoPeaks(int d){
		this.minDistanceBetweenTwoPeaks = d;
	}/*setMinDistanceBetweenTwoPeaks*/
	
	
	public  void setSerDir(String sd){
		this.serDir = sd;
	}/*setSerDir*/
	
	public void setAliDir(String ad){
		this.aliDir =ad;
	}/*setAliDir*/
	
	public void setAlpha(double [] alpha){
		this.alpha = alpha;
	}/*setAlpha*/
	
	public void setStates(String [] s){
		this.states = s;
	}/*setStates*/
	
	public void setStartProbabilities(double [] sp){
		this.startProbabilities = sp;
	}/*setStartProbabilities*/
	
	public void setObservationProfiles(double [][] op){
		this.observationProfiles = op;
	}/*setObservationProfiles*/
	
	public void setTransitionProbability(double [][] tp){
		this.transitionProbability = tp;
	}/*setTransitionProbability*/
	
	public void setEmissionProbability(double [][] ep){
		this.emissionProbability = ep;
	}/*setEmissionProbability*/
	
	
	public int getMinChainLength(){
		return minChainLength;
	}/*getMinChainLength*/
	
	public String getChainDir(){
		return chainDir;
	}/*getChainDir*/
	
	public double [][] getEmissionProbability(){
		return emissionProbability;
	}/*getEmissionProbability*/
	
	public static boolean getPrintoutViterviScores(){
		return printoutViterviScores;
	}/*getPrintoutViterviScores*/
	
	public boolean getFileterMeaninglessPeaks(){
		return fileterMeaninglessPeaks;
	}/*getFilterUnsurroundedPeaks*/
	
	public boolean getConncetAdjacentPeaks(){
		return conncetAdjacentPeaks;
	}/*getConncetAdjacentPeaks*/
	
	public int getMinDistanceBetweenTwoPeaks(){
		return minDistanceBetweenTwoPeaks;
	}/*getMinDistanceBetweenTwoPeaks*/
	
	
	
	public double[][] getTransitionProbability(){
		return transitionProbability;
	}/*getTransitionProbability*/
	
	public double [][] getObservationProfiles(){
		return observationProfiles;
	}/*getObservationProfiles*/
	
	public double [] getStartProbabilities(){
		return startProbabilities;
	}/*getStartProbabilities*/
	
	public String [] getStates(){
		return states;
	}/*getStates*/
	
	
	public String getSerDir(){
		return serDir;
	}/*getSerDir*/
	
	public String getAliDir(){
		return aliDir;
	}/*getAliDir*/
	
	public double [] getAlpha(){
		return alpha;
	}/*getAlpha*/
	
	
	public static class BLOCK{
		
		public String      blockId;
		public Matrix2D    blockMatrix; //this either comp or ser matrix
		public Sequence    blockSeq;	//sequence from the alignment file
		public double [][] observedProfiles;
		public double [][] emissionMatrx;
		public ArrayList<Integer> positionOfMissingData = null; // in some positions, the sum of three profiles is zero! I dont know why yet, but this is to keep track of those positions
		
		public BLOCK(File serFile, String pathToAli, LinkedHashMap<String, double[]> statesAndParameters) throws Exception{
			
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
			this.blockMatrix = m;
			String aliFileName = pathToAli+"/"+blockId+".fa";
			File aliFile       = new File(aliFileName);
			if(!aliFile.exists()){
				System.err.println(" Could not find alignment file for block " + blockId );
			}
			Sequence firstSeq  = SeqIOTools.readFastaDNA(new BufferedReader(new FileReader(aliFile))).nextSequence();
			this.blockSeq = firstSeq;
			
			HashMap<ArrayList<ArrayList<Double>>, ArrayList<Integer>> observedProfilesAndMissingPositions = new HashMap<ArrayList<ArrayList<Double>>, ArrayList<Integer>>();
			observedProfilesAndMissingPositions = getProfiles(m,firstSeq);
			int hmSize = observedProfilesAndMissingPositions.size();
			if(hmSize>1){
				System.err.print("found more than one observed profile matrix!");
			}
			double [][] observedData = null;
			for( Map.Entry<ArrayList<ArrayList<Double>>, ArrayList<Integer>> me: observedProfilesAndMissingPositions.entrySet() ){
				ArrayList<ArrayList<Double>> observedProfileAsList = me.getKey();
				ArrayList<Integer>          missingDataPositions   = me.getValue();
				observedData         = convertArrayListOfArrayListToMatrix(observedProfileAsList);
				this.observedProfiles = observedData;
				this.positionOfMissingData = missingDataPositions;
			}
			double [][] emissionData =  getEmissionProbabilities(observedData,statesAndParameters);
			this.emissionMatrx = emissionData;
			
		}/*BLOCK CONSTRUCTOR*/
		
		
		public double [][]getEmissionProbabilities(double [][] observedData, LinkedHashMap<String , double []> statesAndCorresponingDirichletParameters){
			int numberOfCols = observedData[1].length;
			int numberOfStates = statesAndCorresponingDirichletParameters.size();
			double [][] emissionMatrix = new double[numberOfStates][numberOfCols];
			
			for(int col =0; col < numberOfCols ; col++){
				double [] oneColAsAnArray = new double[3];
				for(int row =0; row < 3 ; row++){
					oneColAsAnArray[row] = observedData[row][col];
				}
				int stateNumber = 0;
				for(Map.Entry<String, double[]> me: statesAndCorresponingDirichletParameters.entrySet()){
					double oneDirichletPar [] = me.getValue();
					double oneDensity         = DirichletDist.density(oneDirichletPar, oneColAsAnArray);
					emissionMatrix[stateNumber][col] = oneDensity;
					stateNumber++;
					
				}
			}		
			return emissionMatrix; 
		}/*getEmissionProbabilities*/
		
		
		
		public HashMap<ArrayList<ArrayList<Double>>, ArrayList<Integer>> getProfiles( Matrix2D m, Sequence seq){
			HashMap<ArrayList<ArrayList<Double>>, ArrayList<Integer> > map = new HashMap<ArrayList<ArrayList<Double>>, ArrayList<Integer>>();
			ArrayList<ArrayList<Double>> obsProf = new ArrayList<ArrayList<Double>>();
			ArrayList<Integer> positionsWithZeroSum = new ArrayList<Integer>();
			
			//ArrayList<ArrayList<Double>> profiles = new ArrayList<ArrayList<Double>>();
			int numberOfRows = m.rows();
			Symbol gap = seq.getAlphabet().getGapSymbol();
			int seqLength = seq.length();
			if(numberOfRows != seqLength){
				System.err.println("Sequence length is not equal to number rows in the matrix " + " numberOfRows = " + " and seqLength = " + seqLength);
			}
			for(int i =0; i<numberOfRows; i++){
				if( (i!= seqLength) && (seq.symbolAt(i+1)) != gap ){
					ArrayList<Double> oneRow = new ArrayList<Double>();
					for(int j=0; j<3; j++){
						double oneElement = m.get(i, j);
						oneRow.add(oneElement);
					}
					if(sum(oneRow) == 0){
						positionsWithZeroSum.add(i);
					}else{
						obsProf.add(oneRow);
					}
				}
				
			}
			
			map.put(obsProf, positionsWithZeroSum);		
			return map;			
		}/*getProfiles*/
		
		
		public static double sum(ArrayList<Double> list){
			if (list == null || list.size()<1){
				return 0;
			}
			double sum = 0;
			for(Double d:list){
				sum = sum+d;
			}
			return sum;
		}/*sum*/
		
		public double[][] convertArrayListOfArrayListToMatrix(ArrayList<ArrayList<Double>> x){
			int sizeOfArrayList = x.size();
			 double [][] matrix =  new  double [3][sizeOfArrayList] ;
			 for(int i=0; i<sizeOfArrayList; i++){
				 for(int j = 0; j<3; j++){
					 double oneElement = x.get(i).get(j);
					 matrix[j][i] = oneElement;
				 }
			 }
			return matrix;
		}/*convertArrayListOfArrayListToMatrix*/
		
		public static class PAIR{
			private double first;
			private int    second;
			PAIR(double first, int second){
				this.first = first;
				this.second = second;
			}/*PAIR-Constructor*/
		}/*PAIR- class*/
		
		public static class PAIROFLISTS{
			private ArrayList<String>  first;
			private ArrayList<Double>  second ;
			PAIROFLISTS(ArrayList<String> x, ArrayList<Double> y){
				this.first = x;
				this.second = y;
			}
		}/*PAIROFLISTS*/
		
		PAIROFLISTS viterbi(String [] states, double [] startProb, double [][] observedData, double [][] emissionData, double [][] transitionMatrix){//return is a pair of arraylist of strings which includes the viterbi path and a double which is its score
			
			ArrayList<String> viterbiPath   = new ArrayList<String>();
			ArrayList<Double> viterbiScores = new ArrayList<Double>();
			int numberOfStates              = states.length;
			int numberOfObservations        = observedData[0].length; // note observedData is 3 by T data where T is the length of sequence (filtering out gaps and missing data) and 3 is because of red, green and blue states
			double [][] viterbiMatrix       = new double[numberOfStates][numberOfObservations];
			int    [][] argmaxIndices       = new int   [numberOfStates][numberOfObservations];
			
			
			ArrayList<Double> lastColInViterbiMatrix = new ArrayList<Double>();
			
			//initialisation
			for(int state = 0; state<numberOfStates; state ++){
				viterbiMatrix[state][0] = Math.log(startProb[state]) + Math.log(emissionData[state][0]);
				argmaxIndices[state][0] = -1;// note: this col ie the zero'th col,  is never accessed in calculating the path (pi's in Richard's book), therefore it is assigned as -1)
				
			}
			
			for(int obs= 1; obs <numberOfObservations; obs ++ ){
				for (int state =0; state < numberOfStates; state ++){
					PAIR maxValueAndItsArg = getArgMaxAndItsMaxValue(obs, state, viterbiMatrix,transitionMatrix);
					double oneScore = maxValueAndItsArg.first;
					int argmax       = maxValueAndItsArg.second;
					argmaxIndices[state][obs] = argmax;
					viterbiMatrix[state][obs]   = Math.log(emissionData[state][obs]) + oneScore;
					if(obs == numberOfObservations-1){
						lastColInViterbiMatrix.add(viterbiMatrix[state][obs]);
					}
				}

			}
			
			
			ArrayList<Integer> viterbiPathIndices = new ArrayList<Integer>();
			PAIR argAndScoreFromLastCol           = getMaxValueAndItsIndexInAGivenArrayList(lastColInViterbiMatrix);
			int  lastColArgmax                    = argAndScoreFromLastCol.second;
			viterbiPathIndices.add(0, lastColArgmax);
			
			
			
			for(int i = numberOfObservations-2; i>= 0; i--){
				int latestArgmax = viterbiPathIndices.get(0);
				int nextArgmax    =  argmaxIndices[latestArgmax][i+1];
				//update
				viterbiPathIndices.add(0, nextArgmax);
				//System.out.println(nextArgmax);
				
			}
			
			int numOfColsInVM = viterbiMatrix[0].length;
			int numOFRowsInVM = viterbiMatrix.length;
			
			for(int i =0; i<numOfColsInVM;i++){
				int oneIndex = viterbiPathIndices.get(i);
				viterbiPath.add(states[oneIndex]);
			}
			
			/*for(int i = 0; i<1000;i++){
				System.out.println("i= " + i + " " + states[viterbiPathIndices.get(i)]);
			}*/
			
			
			PAIROFLISTS viterbiPathAndScores = new PAIROFLISTS(viterbiPath,viterbiScores);
			
			
			
			return viterbiPathAndScores;
			
		}/*viterbi*/
		
		PAIR getArgMaxAndItsMaxValue(int observation, int state,  double [][] vitMatrix, double [][] transitionMatrix){
			int numberOfStates = transitionMatrix[1].length;
			ArrayList<Double> x =  new ArrayList<Double>(); // this will be used to store v_k(i-1)*a_k,l as described in Richard Durbin's book
			for(int row = 0; row<numberOfStates ; row ++){
			
				double oneScore = vitMatrix[row][observation-1] + Math.log(transitionMatrix[row][state]);
				x.add(oneScore);
			}
			PAIR aPair     = getMaxValueAndItsIndexInAGivenArrayList(x);
			return aPair;
		}/*getArgMaxAndItsMaxValue*/
		
		PAIR getMaxValueAndItsIndexInAGivenArrayList(ArrayList<Double> x){
			double maxValue = Double.MAX_VALUE;
			int maxValueIndex =  -1;
			Object max_obj    =  Collections.max(x);
			maxValueIndex     =  x.indexOf(max_obj);
			maxValue          =  (Double) max_obj;
			PAIR aPair = new PAIR(maxValue,maxValueIndex);
			return aPair;
		}/*getMaxValueAndItsIndexInAGivenArrayList*/
		
		
	}/*BLOCK*/
	
	public ArrayList<Double> mapViterbiPathToAlignment(Sequence seq, Matrix2D m, ArrayList<String> vitebiPath){//gaps are mapped to zero, Positive state is mapped to 1 and Negative State is mapped to 0.5
		ArrayList<Double> viterbiPathMappedToAlignment = new ArrayList<Double>();
		
		Symbol gap = seq.getAlphabet().getGapSymbol();
		int sequenceLength = seq.length();
		
		int nonGapPosition = 0;
		for(int i =0; i<sequenceLength; i++ ){
			double oneScore = Double.MIN_VALUE;
			ArrayList<Double> oneRowInMatrix = new ArrayList<Double>();
			for(int col =0; col<3; col++){
				oneRowInMatrix.add(m.get(i, col));
			}
			
			if( (seq.symbolAt(i+1) != gap   &&  (BLOCK.sum(oneRowInMatrix) != 0)  ) ){
				if(vitebiPath.get(nonGapPosition).equals("M")){
					oneScore = 0;
				}
				else if(vitebiPath.get(nonGapPosition).equals("R")){
					oneScore = 2;
				}
				else if(vitebiPath.get(nonGapPosition).equals("r")){
					oneScore = 2;
				}
				else if(vitebiPath.get(nonGapPosition).equals("G")){
					oneScore = 2;
				}
				else if(vitebiPath.get(nonGapPosition).equals("g")){
					oneScore = 2;
				}
				else if(vitebiPath.get(nonGapPosition).equals("E")){
					oneScore = 1;
				}
				else if(vitebiPath.get(nonGapPosition).equals("e")){
					oneScore = 1;
				}
				else if (vitebiPath.get(nonGapPosition).equals("J")){
					oneScore = 1.5;
				}
				
				else{
					System.err.println("Unkown state!");
				}
				nonGapPosition++;	
			}
			else{
				oneScore = -1;
			}
			
			viterbiPathMappedToAlignment.add(oneScore);
		}
		 boolean filterMeaninglessPeaks = getFileterMeaninglessPeaks(); 
		 
		 
		 if(filterMeaninglessPeaks){
			 viterbiPathMappedToAlignment = filterMeaninglessPeaks(viterbiPathMappedToAlignment);
		 }
		 
		 boolean connectAdjacenPeaks = getConncetAdjacentPeaks();
		 
		 int lengthThreshold = getMinDistanceBetweenTwoPeaks();
		 if(connectAdjacenPeaks){
			 viterbiPathMappedToAlignment = connectAdjacentPeaks(viterbiPathMappedToAlignment, lengthThreshold);
		 }
		
		return viterbiPathMappedToAlignment;
	}/*mapViterbiPathToAlignment*/
	
	public ArrayList<Double> filterMeaninglessPeaks(ArrayList<Double> scores){// peaks like ..mgggeREGGGmm.. or ..mrrrEGeRRRmm.. are not acceptable
		int numberOfScores = scores.size();
		for(int i = 1; i<numberOfScores-4; i++){
			
			if(  scores.get(i) ==1  && scores.get(i+2) ==1){
				int dummyStart = i-1;
				int dummyEnd   = i+3;
				scores.set(i, (double) 0);
				scores.set(i+2, (double) 0);
				scores.set(i+1, (double) 0  );
				while(scores.get(dummyStart) != 0){
					scores.set(dummyStart, (double) 0);
					dummyStart = dummyStart-1;
				}
				while(scores.get(dummyEnd) != 0){
					scores.set(dummyEnd, (double) 0);
					dummyEnd = dummyEnd+1;
				}
				i = dummyEnd;
			}
			
		}
		return scores;
	}/*filterMeaninglessPeaks*/
	
	
		
	
	public ArrayList<Double> connectAdjacentPeaks(ArrayList<Double> scores, int lengthThreshold){
		int numberOfScores = scores.size();
		
		int firstPeakPos  = 0;
		int secondPeakPos = 0;
		int pos           = 0; 
		
		System.out.println("numberOfScores = " + numberOfScores);
			
		for(int i = 0; i< numberOfScores-1; i++){
			if(scores.get(i) > 0 && scores.get(i+1) ==0 ){
				//System.out.println("i = " + i );
				firstPeakPos = i;
				pos = i;
				while(scores.get(pos+1) < 1 && pos< numberOfScores-2 ){
					pos++;
					//System.out.println("          " + pos);
				}
				secondPeakPos = pos+1;
				int distance = secondPeakPos- firstPeakPos;
				if(distance<lengthThreshold){
					for(int l=firstPeakPos ; l< secondPeakPos; l++){
						scores.set(l, (double) 2);
					}
					i = i + distance;
				}
			}
		}
		return scores;
	}/*connectAdjacentPeaks*/
	
	
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
	
	public void wirteChainsIntoAGffFile(ArrayList<Double> vitebiPath, String blockId, String pathToChainsDir) throws IOException{
		String blockSpecfications [] = blockId.split("_");
		String Chr                   = blockSpecfications[0];
		String gffoutputFileName = pathToChainsDir+"/"+blockId+".gff";
		PrintWriter pw = new PrintWriter(new FileWriter(gffoutputFileName));
		final GFFWriter gffw = new GFFWriter(pw);
		int lengthOfViterPath = vitebiPath.size();
		int oneStart = 0;
		int oneEnd   = 0;
		int dummyPos = 0;
		
		for(int i =0; i<lengthOfViterPath-1; i++){
			if(vitebiPath.get(i) > 0 && vitebiPath.get(i+1) > 0){
				oneStart = i;
				dummyPos = i;
				while( vitebiPath.get(dummyPos+1) != 0 && dummyPos< lengthOfViterPath-2){
					dummyPos++;
				}/*while*/
				oneEnd = dummyPos+1;
				int oneChainLength = oneEnd- oneStart;
				SimpleGFFRecord oneRecord = new SimpleGFFRecord();
				oneRecord.setStart(oneStart);
				oneRecord.setEnd(oneEnd);
				oneRecord.setSeqName(Chr);
				oneRecord.setScore(100);
				oneRecord.setSource("composure");
				oneRecord.setFeature("chainFromMarkovModel");
				
				
				int minChainLength = getMinChainLength();
				if(oneChainLength> minChainLength){
					gffw.recordLine(oneRecord);
				}
				i = i + oneChainLength;
				
			}/*if*/
		}/*for*/
		
		pw.flush();
	}/*wirteChainsIntoAGffFile*/

	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		
		VITERBI app = new VITERBI();
		
		String serDirName = app.getSerDir();
		String aliDirName = app.getAliDir();
		double [] startProbability         = app.getStartProbabilities();
		double [][] transitionProbability  = app.getTransitionProbability();
		String [] states                   = app.getStates();
		
		double backGroundAlpha     []  = {0.15106438458196936,0.20973711740901407,0.22900628496706044};
		double greenStateAlpha     []  = {0.1713002602078479,11.729438447100865,1.3763120488468623};//{0.1713002602078479,11.729438447100865,1.3763120488468623};
		double greenFollowingMixedAlpha [] = {0.16327302239493097,9.219505194760753,1.8283941303615243};//{0.16327302239493097,9.219505194760753,1.8283941303615243}
		double greenToRedEdgeAlpah []  = {0.4454423966781518,0.46303899247008184,0.27746597159487285};//{0.4799631840073914,0.4944929532552769,0.2857067163030646}; 
		double redStateAlpha       []  = {71.25302335371902,0.201959003605188,0.20299272598532336};//{71.25302335371902,0.201959003605188,0.20299272598532336}
		double redToGreenEdgeAlpha []  = {0.39728078185044136,0.3710506533397517,0.2490053255191496};
		//double redFollowingMixedAlpha   [] = {10.5200496933818615,0.1772646945273629,0.17979061871698085};//{4.5200496933818615,0.1772646945273629,0.17979061871698085};
		double junctionStateAlpha []   = {0.37553647891207315, 0.40882228834461276,0.3426240669996672};

		
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
		
		for(String key:statesAndDirPar.keySet() ){
			System.out.println("state is : " + key);
		}
		
		
		File serDir       =  new File(serDirName);
		for(File aSerFile: serDir.listFiles()){
			String oneFileName = aSerFile.getName();
			if(oneFileName.startsWith(".")){
				continue;
			}
			BLOCK block = new BLOCK(aSerFile,aliDirName,statesAndDirPar);
			double [][] observationProfiles  = block.observedProfiles;
			double [][] emissionProbability  = block.emissionMatrx;
			Sequence seq                     = block.blockSeq;
			Matrix2D m                       = block.blockMatrix;
			String blockId                   = block.blockId;
			System.out.println("block id is " + blockId);
			
			app.setObservationProfiles(observationProfiles);
			app.setEmissionProbability(emissionProbability);
			
			//viterbi
		
			ArrayList<String> viterbiPath   = block.viterbi(states, startProbability, observationProfiles, emissionProbability, transitionProbability).first;
			//ArrayList<Double> viebiScores   = block.viterbi(states, startProbability, observationProfiles, emissionProbability, transitionProbability).second;
			
			ArrayList<Double> viterbiPathMappedToAlignment = app.mapViterbiPathToAlignment(seq,m, viterbiPath);
			String outputFileName  = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/ViterbiPaths/"+blockId+".txt";
			app.printoutViterbiPath(viterbiPathMappedToAlignment, outputFileName);
			
			String chainDir = app.getChainDir();
			app.wirteChainsIntoAGffFile(viterbiPathMappedToAlignment, blockId, chainDir);
			
			//app.printoutViterbiPath(viterbiPathMappedToAlignment, "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/ViterbiPaths/ViterbiPath_withSevenStatesHashMap.txt");
			
			
			System.out.println("viterbiPath size is " + viterbiPath.size());
			System.out.println();
			
			//write resutls into a  file:
			String viterbiPathAndScoresOutputFileName =  "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/ViterbiPaths/viterbiPathAndScores"+"_"+blockId+".txt";
			
			boolean printoutViterbiScores = getPrintoutViterviScores();
			
			if(printoutViterbiScores){
				try{
					FileOutputStream fo = new FileOutputStream(viterbiPathAndScoresOutputFileName);
					PrintStream ps = new PrintStream(fo);
					for(int i = 0; i<viterbiPath.size(); i++){
						//String oneline = Integer.toString(i) + " " + viterbiPath.get(i) + " " + Double.toString(viebiScores.get(i)); 
						String oneline = Integer.toString(i) + " " + viterbiPath.get(i); 
						ps.println(oneline);
					}
					ps.close();
					fo.close();
				}
				catch(Exception e){
				}
			}

		}



	}/*main*/

}/*VITERBI*/
