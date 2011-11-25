/**
 * 
 */
package valuepostprocessing;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.BioException;
import org.biojava.bio.program.gff.GFFEntrySet;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.GFFTools;
import org.biojava.bio.program.gff.GFFWriter;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.biojava.utils.ParserException;
import org.bjv2.util.IntSet;
import org.bjv2.util.IntSets;

import net.derkholm.nmica.utils.CliTools;

/**
 * @author hk3 18-05-2011
 * This class is to compare composure predictions (chains) with DnaseI. For this:  it prints TP, TN, FP and FN into an output file (each line corresponds to one alignment block) 
 *
 */
public class EvaluateComposurePredictions {

	/**
	 * @param args
	 */
	
	private double chainThreshold;// =  0.7;
	private String compDir;//           = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/Eisenlab_Composure_Output_Model_20051212_With_Window_Length_350_temp";
	private String chainFile;//        =  "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/Chains_From_MarkovModel_ChainFinder_In_Real_Coordinates/ChainesInRealCoordinages.gff"; // note this chains is assumed to be in real coordiantes
	private String gsFile;//            =  "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/Unionof_DnaseI_CodingRegions_and_RedFlyEnhancers/UnionOfbdtnpDnaseAccAllStages.gff"; // note: this is the gold standard file (for example DnaseI), must be in fly release five coordinates (or the same release with chains coordinates).
	private String pathToOutputDir;//   =  "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/Eisenlab_Composure_Output_Model_20051212_With_Window_Length_350_Stats";
	private String comparisonTarget;// =  "DnaseI"; // Note: only "DnaseI", "CodingRegions" and "RedFlyEnhancers" are accepted.
	private double gsScoreThreshold;// =  0.0;// it is assumed that records in RedFlyEnhancers, CodingRegions and DnaseI all have scores greater than or equal to zero.
	private int   offset; // for a gs record r we accept rStart-offset and rEnd+offset.
	private int thresholdForBlockLengths;// = 2000; // some blocks are too short, this is to ignore those short blocks
	
	
	//SETTERS
	public void setOffset(int os){
		this.offset = os;
	}/*setOffset*/
	
	public void setCompDir(String cd){
		this.compDir = cd;
	}/*setCompDir*/
	
	public void setGsScoreThreshold(double score){
		this.gsScoreThreshold = score;
	}/*setGsScoreThreshold*/
	
	public void setChainFile(String cf){
		this.chainFile = cf;
	}/*setChainFile*/
	
	public void setGsFile(String gsf){
		this.gsFile = gsf;
	}/*setGsFile*/
	
	public void setPathToOutputDir(String pod){
		this.pathToOutputDir = pod;
	}/*setPathToOutputDir*/
	
	public void setChainThreshold(double ct){
		this.chainThreshold = ct;
	}/*setChainThreshold*/
	
	public void setComparisonTarget(String ct){
		this.comparisonTarget = ct;
	}/*setComparisonTarget*/
	
	public void setThresholdForBlockLengths(int l){
		this.thresholdForBlockLengths = l;
	}/*setThresholdForBlockLengths*/
	
	//GETTERS
	
	public int getThresholdForBlockLengths(){
		return thresholdForBlockLengths;
	}/*getThresholdForBlockLengths*/
	
	public int getOffset(){
		return offset;
	}/*getOffset*/
	
	public String getCompDir(){
		return compDir;
	}/*getCompDir*/
	public String getPathToOutputDir(){
		return pathToOutputDir;
	} /*getPathToOutputDir*/
	
	public String getChainFile(){
		return chainFile;
	}/*getChainFile*/
	
	public String getGsFile(){
		return gsFile;
	}/*getGsFile*/
	
	public  double getGsScoreThreshold(){
		return gsScoreThreshold;
	}/*getGsScoreThreshold*/
	
	public double getChainThreshold(){
		return chainThreshold;
	}/*getChainThreshold*/
	
	public String getComparisonTarget(){
		return comparisonTarget;
	}/*getComparisonTarget*/
	
	public static void main(String[] args) throws Exception{
		// TODO Auto-generated method stub
		EvaluateComposurePredictions app = new EvaluateComposurePredictions();
		args = CliTools.configureBean(app, args);
		app.run(args);
		

	}/*main*/
	
	private void run(String[] args) throws ParserException, BioException, IOException{
		
		ArrayList <String > listOfCompFiles = getListOfNamesOfFilesInADirectory(compDir);
		UnionOfFeaturesInAGffFile aGFFObject = new UnionOfFeaturesInAGffFile();
		
		double chainThr = getChainThreshold();
		double gsThr    = getGsScoreThreshold();
		String chainFileName = getChainFile(); 
		String gsFileName    = getGsFile();
		int lThreshold        = getThresholdForBlockLengths();
		String pathToStats    = getPathToOutputDir();
		
		int numberOfRepeats =0;
		String senOutputFile = pathToStats+"/sensitivity.txt";
		String specOutputFile = pathToStats + "/specificity.txt";
		FileWriter senWriter    = new FileWriter(senOutputFile);
		FileWriter specWriter   = new FileWriter(specOutputFile);

		while(numberOfRepeats<1){
			int OverallLengthOfChains = 0;
			int OverallLengthOFGs     = 0;
			int TP                    = 0;
			int FP                    = 0;
			int FN                    = 0;
			

			//int transformationLength = generator.nextInt(500);
			//randomWriter.append(Integer.toString(transformationLength));
			//randomWriter.append(",");
			for(int i=0; i<listOfCompFiles.size(); i++){
			
				String blockId = listOfCompFiles.get(i).replaceAll(".comp", "");
				int bLength = getBlockLength(blockId);
				if(bLength< lThreshold){
					continue;
				}
				
				Block  block = new Block();
				block.setBlockid(blockId);
				//System.out.println("block  " + i + " is " + blockId +  " ************************* " );
				//String chr = block.getBlockChromosome();
				//int start = block.getBlockStart();
				//int end   = block.getBlockEnd();
				//int oneBlockLength =  end- start;
		
			//make two temp files for comparison (if theres is any record from chains or gs in the range of this block, it is to print in these files to be compared)
				String chainTempFileName = pathToStats+"/chain_"+blockId+".gff" ;
				String gsTempFileName    = pathToStats+"/gs_"+blockId+".gff" ;
			
			
			
				setOffset(0);
				int chainOffset = getOffset();
				printRecordsWithinThisAlingmentBlockIntoATempFile(chainFileName, blockId, chainTempFileName,chainThr,chainOffset,false);
				setOffset(0);
				int gsOffset = getOffset();
				printRecordsWithinThisAlingmentBlockIntoATempFile(gsFileName,    blockId, gsTempFileName,gsThr,gsOffset, false );//normally boolean must be false, put it true only you want to shuffle locations.
			
				String outputFileName = getOutputFileName(comparisonTarget);
			 
				//getStatisticsVersion2(chainTempFileName,gsTempFileName,outputFileName,blockId,chainThr);
				
				BLOCKSTATS bs = new BLOCKSTATS(chainTempFileName,gsTempFileName,blockId,chainThr);
				int tp = bs.tp;
				int fp  = bs.fp;
				int fn  = bs.fn;
				int chainsLength = bs.totalLengthOfChains;
				int gsLength     = bs.totalLengthOfGoldStandardRecords;
				int blockLength = bs.blockLength;
				double probChains = bs.prob_of_chains;
				double probGS     = bs.prob_of_gs;
				double probOverlaps = bs.prov_of_overlaps;
				double sen = Double.parseDouble(Integer.toString(tp))/Double.parseDouble(Integer.toString(gsLength));
				double spec = Double.parseDouble(Integer.toString(tp))/Double.parseDouble(Integer.toString(chainsLength));
				OverallLengthOfChains = OverallLengthOfChains + chainsLength;
				OverallLengthOFGs     = OverallLengthOFGs + gsLength; 
				TP                    = TP + tp;
				FP                    = FP + fp;
				FN                    = FN + fn;
				
				//System.out.println(" tp = " + tp +  " totalLenghtOfChains = " + chainsLength + " totalLengthOfGoldStandardRecords= " + " blockLength = " +blockLength );
				//System.out.println("prov_of_chains = " +probChains + " prob_of_gs = " + probGS + " probOverlaps = " + probOverlaps );
				//System.out.println(" Sen = " + sen + " Spec = " + spec);
				
				
				
			
				aGFFObject.deleteFile(chainTempFileName);
				aGFFObject.deleteFile(gsTempFileName);
			//System.out.println("chainTempFile " + chainTempFileName + " gsTempFile = " +gsTempFileName +" output = " +  outputFileName + " bId = " +   blockId );
			}
			double SEN = Double.parseDouble(Integer.toString(TP))/Double.parseDouble(Integer.toString(OverallLengthOFGs));
			double SPEC = Double.parseDouble(Integer.toString(TP))/Double.parseDouble(Integer.toString(OverallLengthOfChains));
			double SEN_2 = Double.parseDouble(Integer.toString(TP))/Double.parseDouble(Integer.toString(TP+FN));
			double SPEC_2 = Double.parseDouble(Integer.toString(TP))/Double.parseDouble(Integer.toString(TP+FP));
			senWriter.append(Double.toString(SEN));
			senWriter.append(",");
			specWriter.append(Double.toString(SPEC));
			specWriter.append(",");
			System.out.println(" SEN = " + SEN + " SPEC = " + SPEC);
			System.out.println("TP = " + TP + " FP = "+ FP + " FN = " + FN );
			//System.out.println( " SEN_2 = " + SEN_2 + " SPEC_2 = " + SPEC_2);
			
			
			numberOfRepeats++;
			
		}
		//randomWriter.flush();
		senWriter.flush();
		specWriter.flush();
		//randomWriter.close();
		senWriter.close();
		specWriter.close();
		
		//System.out.println("total length of chain = " + OverallLengthOfChains);
	}/*run*/
	
	 public static int  getBlockLength(String blockId){// it is assumed block id is in  chr_start_end form
		String [] blockSpecifiations = blockId.split("_");
		int start = Integer.parseInt(blockSpecifiations[1]);
		int end  =  Integer.parseInt(blockSpecifiations[2]);
		int blockLength = end-start;
		return blockLength;
	}/*getBlockLength*/
	
	 
	 public static class BLOCKSTATS{
		 public int tp   = 0; 
		 public int fp   = 0;
		 public int fn   = 0;
		 public int tn   = 0;
		 public int blockLength = 0;
		 public int totalLengthOfChains = 0;
		 public int totalLengthOfGoldStandardRecords = 0;
		 public double prob_of_chains ;
		 public	double prob_of_gs ;
		 public	double prov_of_overlaps;
		 

		 
		 
		 
		 public BLOCKSTATS(String chainFileName, String gsFileName, String blockId, double chainThreshold) throws IOException, ParserException, BioException{
			 this.blockLength = getBlockLength(blockId);
			 FileInputStream fisChain  = new FileInputStream(new File(chainFileName)); 
			 FileInputStream fisGs     = new FileInputStream(new File(gsFileName));

			 int chainEmty             = fisChain.read();  //this is to check if file is empty fisChain.available();
			 int gsEmty                = fisGs.read();;//fisGs.available();
			 if(chainEmty == -1 && gsEmty == -1){//there is no chain no gs record
				 this.tn = blockLength;
				 this.prob_of_chains   = 0.0;
				 this.prob_of_gs       = 0.0;
				 this.prov_of_overlaps = 0.0;
			 }
			 else if(chainEmty == -1 &&  gsEmty != -1){//there is no chain
				 this.totalLengthOfGoldStandardRecords = getOverallLengthOfFeatures(gsFileName);
				 this.fn = getOverallLengthOfFeatures(gsFileName);
				 this.tn = blockLength-fn;
				 this.prob_of_chains =  0;
				 this.prov_of_overlaps = 0;
				 this.prob_of_gs = Double.parseDouble(Integer.toString(fn))/Double.parseDouble(Integer.toString(blockLength));
			 }
			 
			 else if(gsEmty == -1 && chainEmty != -1){//there is no gs
				 
				 this.totalLengthOfChains  = getOverallLengthOfFeatures(chainFileName);
				 this.fp = totalLengthOfChains;
				 this.tn = blockLength-fp;
				 this.prob_of_chains = Double.parseDouble(Integer.toString(fp))/Double.parseDouble(Integer.toString(blockLength));
				 this.prob_of_gs  = 0.0;
				 this.prov_of_overlaps = 0.0;
			 }
			 else if(gsEmty != -1 && chainEmty != -1){//both files exist
			 		
				 this.totalLengthOfChains              = getOverallLengthOfFeatures(chainFileName);
				 this.totalLengthOfGoldStandardRecords = getOverallLengthOfFeatures(gsFileName);
			
			
				 //IntSet UnionOfFeatures = IntSets.none(); //this union of all chains and gs records in within this block
				 IntSet UnionOfOverlaps = IntSets.none(); // this is union of overlaps (ie intesection of chain and gs records)
				 int    totalOverlapLength =0;
				 
			
				 File chainGffFile = new File(chainFileName);
				 GFFEntrySet chaingff = GFFTools.readGFF(chainGffFile);
			
				 File gsGffFile  = new File(gsFileName);
				 GFFEntrySet gsgff = GFFTools.readGFF(gsGffFile);
			
			
			
				 for(Iterator<?> chainGffI = chaingff.lineIterator(); chainGffI.hasNext() ;){
					 Object anObject = chainGffI.next();
					 if(anObject instanceof GFFRecord){
						 GFFRecord chainRecord = (GFFRecord) anObject;
						 int oneChainStart = chainRecord.getStart();
						 int oneChainEnd   = chainRecord.getEnd();
						 IntSet OneChainRange = IntSets.range(oneChainStart, oneChainEnd);
						 for(Iterator<?> gsGffI = gsgff.lineIterator(); gsGffI.hasNext();){
							 Object anotherObject = gsGffI.next();
							 if(anotherObject instanceof GFFRecord){
								 GFFRecord gsRecord = (GFFRecord) anotherObject;
								 int oneGsStart = gsRecord.getStart();
								 int oneGsEnd   = gsRecord.getEnd();
								 IntSet oneGsRange = IntSets.range(oneGsStart, oneGsEnd);
								 //IntSet oneUnion = getUnionOfTwoIntSets(OneChainRange,oneGsRange);
								 //UnionOfFeatures = getUnionOfTwoIntSets(UnionOfFeatures,oneUnion);
								 IntSet OneIntesection = IntSets.intersection(OneChainRange,oneGsRange);
								 if(!OneIntesection.isEmpty()){
									 int overlapStar = Math.max(oneChainStart,oneGsStart);
									 int overlapEnd  = Math.min(oneChainEnd, oneGsEnd);
									 int oneOverlapLength = overlapEnd- overlapStar;
									 totalOverlapLength   = totalOverlapLength+oneOverlapLength;
								 }
								 UnionOfOverlaps   = getUnionOfTwoIntSets(UnionOfOverlaps,OneIntesection);		
							//System.out.println("chainSt = " + oneChainStart + " chainEnd = "+ oneChainEnd +  " gsStart = " + oneGsStart + " gsEnd " + oneGsEnd);
							 }					
						 }/*end of gs iterator*/
					 }
				 }/*end of  chain iterator*/
				 	this.tp  = totalOverlapLength;//UnionOfOverlaps.size()-1;
				 	this.fp = totalLengthOfChains- tp;
				 	this.fn  = totalLengthOfGoldStandardRecords- tp;
				 	this.tn = blockLength -totalLengthOfChains-totalLengthOfGoldStandardRecords +tp ;
				 	this.prob_of_chains  = Double.parseDouble(Integer.toString(totalLengthOfChains))/Double.parseDouble(Integer.toString(blockLength));
				 	this.prob_of_gs      = Double.parseDouble(Integer.toString(totalLengthOfGoldStandardRecords))/Double.parseDouble(Integer.toString(blockLength));
				 	this.prov_of_overlaps = prob_of_chains*prob_of_gs;
			 		}
			 else{
				 System.err.println("classification");
			 }
			 
			 fisChain.close();
			 fisGs.close();

		 }/*BLOCKSTATS-constructor*/
		 
		 
	 }/*BLOCKSTATS*/
	 
	
	public void getStatisticsVersion2(String chainFileName, String gsFileName, String outputFileName, String blockId, double chainThreshold) 
	throws IOException, ParserException, BioException{
		int tp = 0;
		int fp = 0;
		int fn = 0;
		int tn = 0;
		int blockLength                      = getBlockLength(blockId);
		int totalLengthOfChains              = 0;
		int totalLengthOfGoldStandardRecords = 0;
		outputFileName                       =  modifyOutputFileName(outputFileName, chainThreshold);
		
		 FileInputStream fisChain  = new FileInputStream(new File(chainFileName)); 
		 FileInputStream fisGs     = new FileInputStream(new File(gsFileName));
		 int chainEmty             = fisChain.read();  //this is to check if file is empty fisChain.available();
		 int gsEmty                = fisGs.read();;//fisGs.available();
		 
		 if(chainEmty == -1 && gsEmty == -1){//there is no chain no gs record
			 tn = blockLength;
			 writeBlockStats(outputFileName,blockId,tp,fp,fn,tn);
		 }
		 else if(chainEmty == -1 &&  gsEmty != -1){//there is no chain
			 
			 fn = getOverallLengthOfFeatures(gsFileName);
			 tn = blockLength-fn;
			 writeBlockStats(outputFileName,blockId,tp,fp,fn,tn);	
			 System.out.println("length of gs = " + fn);
		 }
		 else if(gsEmty == -1 && chainEmty != -1){//there is no gs
			 fp = getOverallLengthOfFeatures(chainFileName);
			 tn = blockLength-fp;
			 writeBlockStats(outputFileName,blockId,tp,fp,fn,tn);
			 System.out.println("length of chains = " + fp);
		 }
		 else if(gsEmty != -1 && chainEmty != -1){//both files exist
		 		
			 totalLengthOfChains              = getOverallLengthOfFeatures(chainFileName);
			 totalLengthOfGoldStandardRecords = getOverallLengthOfFeatures(gsFileName);
		
		
			 //IntSet UnionOfFeatures = IntSets.none(); //this union of all chains and gs records in within this block
			 IntSet UnionOfOverlaps = IntSets.none(); // this is union of overlaps (ie intesection of chain and gs records)
			 int    totalOverlapLength =0;
			 
		
			 File chainGffFile = new File(chainFileName);
			 GFFEntrySet chaingff = GFFTools.readGFF(chainGffFile);
		
			 File gsGffFile  = new File(gsFileName);
			 GFFEntrySet gsgff = GFFTools.readGFF(gsGffFile);
		
		
		
			 for(Iterator<?> chainGffI = chaingff.lineIterator(); chainGffI.hasNext() ;){
				 Object anObject = chainGffI.next();
				 if(anObject instanceof GFFRecord){
					 GFFRecord chainRecord = (GFFRecord) anObject;
					 int oneChainStart = chainRecord.getStart();
					 int oneChainEnd   = chainRecord.getEnd();
					 IntSet OneChainRange = IntSets.range(oneChainStart, oneChainEnd);
					 for(Iterator<?> gsGffI = gsgff.lineIterator(); gsGffI.hasNext();){
						 Object anotherObject = gsGffI.next();
						 if(anotherObject instanceof GFFRecord){
							 GFFRecord gsRecord = (GFFRecord) anotherObject;
							 int oneGsStart = gsRecord.getStart();
							 int oneGsEnd   = gsRecord.getEnd();
							 IntSet oneGsRange = IntSets.range(oneGsStart, oneGsEnd);
							 //IntSet oneUnion = getUnionOfTwoIntSets(OneChainRange,oneGsRange);
							 //UnionOfFeatures = getUnionOfTwoIntSets(UnionOfFeatures,oneUnion);
							 IntSet OneIntesection = IntSets.intersection(OneChainRange,oneGsRange);
							 if(!OneIntesection.isEmpty()){
								 int overlapStar = Math.max(oneChainStart,oneGsStart);
								 int overlapEnd  = Math.min(oneChainEnd, oneGsEnd);
								 int oneOverlapLength = overlapEnd- overlapStar;
								 totalOverlapLength   = totalOverlapLength+oneOverlapLength;
							 }
							 UnionOfOverlaps   = getUnionOfTwoIntSets(UnionOfOverlaps,OneIntesection);		
						//System.out.println("chainSt = " + oneChainStart + " chainEnd = "+ oneChainEnd +  " gsStart = " + oneGsStart + " gsEnd " + oneGsEnd);
						 }					
					 }/*end of gs iterator*/
				 }
			 }/*end of  chain iterator*/
			 	tp  = totalOverlapLength;//UnionOfOverlaps.size()-1;
			 	fp = totalLengthOfChains- tp;
			 	fn  = totalLengthOfGoldStandardRecords- tp;
			 	tn = blockLength -totalLengthOfChains-totalLengthOfGoldStandardRecords +tp ;
			 	writeBlockStats(outputFileName,blockId,tp,fp,fn,tn);
			 	double probability_of_chain = Double.parseDouble(Integer.toString(totalLengthOfChains))/Double.parseDouble(Integer.toString(blockLength));
			 	double probability_of_gs    = Double.parseDouble(Integer.toString(totalLengthOfGoldStandardRecords))/Double.parseDouble(Integer.toString(blockLength));
			 	System.out.println("tp = " + tp + " totalLenghtOfChains = " + totalLengthOfChains + " totalLengthOfGoldStandardRecords = " + totalLengthOfGoldStandardRecords + " blockLength = " + blockLength);
			 	System.out.println("prob_of_chains = " + probability_of_chain + " prob_of_gs = " + probability_of_gs + " prob_of_overlap =" + probability_of_chain*probability_of_gs);
			 	
			 	//System.out.println("tp = " + tp + " fp = " + fp + " fn = " + fn + " tn = " + tn);
			 	
			 	double Sen = Double.parseDouble(Integer.toString(tp))/Double.parseDouble(Integer.toString(tp+fn));
			 	double Spec = Double.parseDouble(Integer.toString(tp))/Double.parseDouble(Integer.toString(tp+fp));
			 	System.out.println("Sen = " + Sen + " Spec = " + Spec);
			 	
		 		}
		 else{
			 System.err.println("classification");
		 }
		 fisChain.close();
		 fisGs.close();
		//System.out.println("blokc length is " + blockLength + " tp = " + tp + " fp =" + fp + " fn = " + fn + " tn = " + tn );
	}/*getStatisticsVersion2*/
	
	
	
	
	
	public String getOutputFileName(String comparsionTarget){
		String fileName = null;
		if(comparsionTarget.equals("DnaseI")){
			fileName = pathToOutputDir+"/DnaseIStats.txt";
		}
		else if(comparsionTarget.equals("CodingRegions")){
			fileName = pathToOutputDir+"/GenesStats.txt";
		}
		else if(comparsionTarget.equals("RedFlyEnhancers")){
			fileName = pathToOutputDir+"/EnhancersStats.txt";
		}
		else{
			System.err.println("Unkown comparison target! The accepted Strings are  only DnaseI or CodingRegions or RedFlyEnhancers ");
			System.exit(0);
		}
		return fileName;
	}/*getOutputFileName*/

	
	public void writeBlockStats(String outputFileName , String blockID, int truePositive, int falsePositive,int falseNegative, int trueNegative ) throws IOException{
		FileWriter fw     = new FileWriter(outputFileName,true);
        BufferedWriter bw = new BufferedWriter(fw);
        
 
		String line = blockID+"\t"+truePositive+"\t"+falsePositive+"\t"+falseNegative+"\t"+trueNegative;
		bw.write(line);
		bw.newLine();
		bw.close();
		
	}/*writeBlockStats*/

	
	public static int getOverallLengthOfFeatures(String inputFile) throws IOException, ParserException, BioException {
		int overallLength = 0;
		File gffFile = new File(inputFile);
		//System.out.println("From inside the function  input file is " +  inputFile);
		GFFEntrySet gff  = GFFTools.readGFF(gffFile);
		
		for(Iterator<?> gffi = gff.lineIterator(); gffi.hasNext();){
			Object o = gffi.next();
			if (o instanceof GFFRecord){
				GFFRecord gro = (GFFRecord) o;
				int oneStart = gro.getStart();
				int oneEnd   = gro.getEnd();
				int oneLength = oneEnd-oneStart;
				overallLength = overallLength + oneLength;
			}
		}
		

		return overallLength;
	}/*getOverallLengthOfFeatures*/

	
	public  static IntSet getUnionOfTwoIntSets(IntSet A, IntSet B){//there seems to me that the IntSets.union doesnt work properly. when a set contains the other, it doesnt return mathematical union of two sets. This is t get rid of that.
		
		IntSet  C = IntSets.none();//C = null;
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

	
	public static String  modifyOutputFileName(String outputFileName, double x){
		String [] fileNameSplitted = outputFileName.split("\\.");
		String outputFile = fileNameSplitted[0]+"_"+ x + "."+ fileNameSplitted[1]; 
		return outputFile;
	}/*modifyOutputFileName*/

	
	
	public ArrayList<String> getListOfNamesOfFilesInADirectory(String dirName ){
		ArrayList<String> listOfNamesOfFilesInCompDir = new ArrayList<String>();
		
		File dir = new File(dirName);
		String [] children = dir.list();
		if(children == null){
			System.err.println("Either directory is empty or is not a directory!");
		}
		else{
			for(String child:children){
				if(  !child.startsWith(".") ){
					listOfNamesOfFilesInCompDir.add(child);
				}
			}
		}
		
		
		return listOfNamesOfFilesInCompDir;
	} /*getListOfNamesOfFilesInADirectory*/
	
	

	
	public void printRecordsWithinThisAlingmentBlockIntoATempFile(String aGffFile, String blockId, String outputFileName, double scoreThreshold, int offset, boolean shufleLengths)//lengthOfTransromation, if required locations is transformed in order to get the significant of measurement in real data.
	throws ParserException, BioException, IOException{
		String [] idSpecifications = blockId.split("\\_");
		String chr = idSpecifications[0];
		int    blockStart = Integer.parseInt( idSpecifications[1]);
		int   blockEnd   =  Integer.parseInt(idSpecifications[2]);
		IntSet blockRange = IntSets.range(blockStart, blockEnd);
		Random     generateRandom = new Random();
			
		File    tempFile    = new File(outputFileName);
		PrintWriter pw      = new PrintWriter(tempFile);
		GFFWriter  gffw     = new GFFWriter(pw);
		
	
		File inputGffFile   = new File(aGffFile);
		GFFEntrySet  gff    = GFFTools.readGFF(inputGffFile);
		
		for(Iterator<?> gffi = gff.lineIterator(); gffi.hasNext();){
			Object o = gffi.next();
			if(o instanceof GFFRecord){
				GFFRecord gro = (GFFRecord) o;
				String oneChr = gro.getSeqName();
				double oneScore = gro.getScore();
				
				if( (oneChr.equals(chr) )  && (oneScore >=  scoreThreshold )){
					SimpleGFFRecord r = new SimpleGFFRecord(gro);
					int oneStart = gro.getStart()-offset;
					int oneEnd   = gro.getEnd()+offset;
					if(oneStart < 0 ){
						oneStart=0;
						System.out.println("Warining! one record modified!");
					}
					//System.out.println("ouside of shufle loop!");

					IntSet oneRecordRange = IntSets.range(oneStart, oneEnd);
					IntSet oneIntersection = IntSets.intersection(blockRange,oneRecordRange);
					if(!oneIntersection.isEmpty()){
						int overlapStart = oneIntersection.getMin();
						int overlapEnd   = oneIntersection.getMax();
						if(shufleLengths){
							int overlapLength = overlapEnd-overlapStart;
							int blockLength      = blockEnd-blockStart;
							int upperBoundForStart = blockLength-overlapLength;
							int randomStart = generateRandom.nextInt(upperBoundForStart)+blockStart;
							overlapStart = randomStart;
							overlapEnd   = overlapStart + overlapLength;
							//System.out.println("blockS " + blockStart + " rStart = " + overlapStart + " rEnd = " + overlapEnd + " bEnd = " + blockEnd);
							
						}
						r.setStart(overlapStart);
						r.setEnd(overlapEnd);
						r.setScore(oneScore);
						gffw.recordLine(r);
						//System.out.println("from gro = chr " + gro.getSeqName() + " st = " + gro.getStart() + " end = " + gro.getEnd() );
						//System.out.println("from r = chr " + r.getSeqName() + " st = " + r.getStart() + " end = " + r.getEnd() );
					}
				}
				
				
			}			
		}
		pw.flush();
		pw.close();
		
	}/*printRecordsWithinThisAlingmentBlockIntoATempFile*/
	
	
	private  class Block{
		private String blockid;
		
		public void setBlockid(String id){
			this.blockid = id;
		}/*setBlockid*/
		
		public String  getBlockid(){
			return blockid;
		}/*getBlockid*/
		
		public String getBlockChromosome(){
			String [] specifications = this.blockid.split("\\_");
			return specifications[0];
		}/*getBlockChromosome*/
		
		public int getBlockStart(){
			String [] specifications = this.blockid.split("\\_");
			return Integer.parseInt(specifications[1]);
		}/*getBlockStart*/
		
		public int getBlockEnd(){
			String [] specifications = this.blockid.split("\\_");
			return Integer.parseInt(specifications[2]);
		}
	}/*Block*/
		
}/*EvaluateComposurePredictions*/