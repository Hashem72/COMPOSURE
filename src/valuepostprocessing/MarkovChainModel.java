/**
 * 
 */
package valuepostprocessing;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;

import net.derkholm.nmica.matrix.Matrix2D;
import net.derkholm.nmica.matrix.SimpleMatrix2D;
import net.derkholm.nmica.utils.CliTools;
import net.derkholm.nmica.utils.ConfigurationException;

/**
 * @author hk3  27-05-2011
 * this will read a positive and negative models from the given files and it will make the beta matrix (see Richard Durbins book page 52) and then associate each sequence x with a score s(s)
 */
public class MarkovChainModel {

	/**
	 * @param args
	 */
	
	private String serDir                  = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/Eisenlab_Composure_Output_Model_20051212_With_Window_Length_350_SER";
	private String aliDir                  = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ALIGNMENTS_EISENLAB/2L";
	private String positiveModelFileName   = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/MarkovModels/EnhancersPositive.txt";
	private String negativeModelFileName   = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/MarkovModels/EnhancersNegative.txt";	
	private String pathToGffFiles          = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/Enhancers_Mapped_to_Alignments";
	private String outputFilePosSet        = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/Eisenlab_Composure_Output_Model_20051212_With_Window_Length_350_Stats/EnhancersPosMCScores.txt";
	private String outputFileNegSet        = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/Eisenlab_Composure_Output_Model_20051212_With_Window_Length_350_Stats/EnhancersNegMCScores.txt";
	private int minLength         = 100; // sequence shorter (after removing gaps) than this will not be accepted for training  the model.

	
	
	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		MarkovChainModel app  = new MarkovChainModel();
		args = CliTools.configureBean(app, args);
		app.run(args);
	}/*main*/
	
	private void run(String [] args) throws Exception{
		
		String posModelFileName       = getPositiveModelFileName();
		String negModelFileName       = getNegativeModelFileName();
		String pathToAliDir           = getAliDir();
		int    minLength              = getMinLength();
		Matrix2D posTransitionMatrix  = getTransitionMatrix(posModelFileName);
		Matrix2D negTransitionMatrix  = getTransitionMatrix(negModelFileName);
		Matrix2D logLikelihoodRatioMatrix = getLogLikelihoodRationsMatrix(posTransitionMatrix,negTransitionMatrix);
		for(int i = 0;i<3; i++){
			for(int j= 0; j<3; j++){
				System.out.print(logLikelihoodRatioMatrix.get(i, j) + "\t");
			}
			System.out.println();
		}
		
		PositiveModel aPosObj  = new PositiveModel();
		NegativeModel aNegObj  = new NegativeModel();
		ArrayList<Integer> listOfLengthOfPosStrings = aNegObj.getlistOfLengthsOfAllGoldStandardRecords(serDir, pathToGffFiles);
		ArrayList<String> posSetAllStrings = aPosObj.getAllStringOfSymbolsAssociatedToAllSerFiles(serDir, pathToGffFiles,minLength);
		ArrayList<String> negSetAllStrings = aNegObj.makeNegativeSetOfStrings(listOfLengthOfPosStrings, pathToGffFiles, serDir,pathToAliDir,minLength);
		printScoresOfAllSequenesIntoAFile(posSetAllStrings,getOutputFilePosSet(),logLikelihoodRatioMatrix, true);
		printScoresOfAllSequenesIntoAFile(negSetAllStrings,getOutputFileNegSet(),logLikelihoodRatioMatrix,true);
		
		//double seqScore = getSequenceScore(testStr, logLikelihoodRatioMatrix);
		
	}/*run*/
	
	
	/************************      METHODS ***************/
	public void setPositiveModelFileName(String pm){
		this.positiveModelFileName = pm;
	}/*setPositiveModelFileName*/
	
	public void setNegativeModelFileName(String nm){
		this.negativeModelFileName= nm;
	}/*setNegativeModelFileName*/
	
	public void setMinLength(int minLen){
		this.minLength = minLen;
	}/*setMinLength*/
	
	public void setAliDir(String ad){
		this.aliDir = ad;
	}/*setAliDir*/
	
	public String getAliDir(){
		return aliDir;
	}/*getAliDir*/
	public int getMinLength(){
		return minLength;
	}/*getMinLength*/
	
	
	public String getPositiveModelFileName(){
		return positiveModelFileName;
	}/*getPositiveModelFileName*/

	public String getNegativeModelFileName(){
		return negativeModelFileName;
	}/*getNegativeModelFileName*/
	
	public void setSerDir(String sd){
		this.serDir = sd;
	}/*setSerDir*/
	
	public void setPatToGffFiles(String path){
		this.pathToGffFiles = path;
	}/*setPatToGffFiles*/
	
	public void setOutputFilePosSet(String outputFile){
		this.outputFilePosSet = outputFile;
	}/*setOutputFilePosSet*/
	
	public void setOutputFileNegSet(String outputFileNeg){
		this.outputFileNegSet = outputFileNeg;
	}/*setOutputFileNegSet*/
	
	public String getOutputFilePosSet(){
		return outputFilePosSet;
	}/*getOutputFilePosSet*/
	
	public String  getOutputFileNegSet(){
		return outputFileNegSet;
	}/*getOutputFileNegSet*/
	
	
	public String getSerDir(){
		return serDir;
	}/*getSerDir*/
	
	
	
	public String getPathToGffFiles(){
		return pathToGffFiles;
	}/*getPathToGffFiles*/
	
	public void printScoresOfAllSequenesIntoAFile(ArrayList<String> listOfStrings, String outputFileName, Matrix2D logRatioMatrix, boolean normalize){//give list of all sequences ( either in positive set or negative set), calculate the scores and print them into a file
		ArrayList<Double> scores = new ArrayList<Double>();
		int counter = 0;
		for(String oneStr: listOfStrings){
			double oneScore = getSequenceScore(oneStr, logRatioMatrix, normalize);
			System.out.println("seq " + counter + " is in process");
			scores.add(oneScore);
			counter++;
		}
		
		try{
			//FileOutputStream fos = new FileOutputStream(outputFileName);
			//DataOutputStream dos = new DataOutputStream(fos);
			  FileWriter fw = new FileWriter(outputFileName);
			  BufferedWriter bw = new BufferedWriter(fw);

			for(double d: scores){
				bw.write(Double.toString(d)+",");
				//dos.writeChars(",");
			}
			bw.close();
			
		}
		catch(IOException e){
			System.out.println("IOException : " + e);
		}
	}/*printScoresOfAllSequenesIntoAFile*/
	
	public double getSequenceScore(String seq, Matrix2D logRatioScores, boolean normalize){
		double strScore = 0;
		int strLength = seq.length();
		if(strLength <= 2){
			System.out.println("seq has length = " + strLength + " and seq itself is " + seq);
		}
		for(int i = 0; i<strLength-1; i++){
			String subStr = seq.substring(i, i+2);
			double oneScore =0;
			if(subStr.equals("RR")){
				oneScore = logRatioScores.get(0, 0);			
			}
			else if(subStr.equals("RG")){
				oneScore = logRatioScores.get(0,1);
			}
			else if(subStr.equals("RB")){
				oneScore = logRatioScores.get(0,2);
			}
			else if(subStr.equals("GR")){
				oneScore = logRatioScores.get(1,0);
			}
			else if(subStr.equals("GG")){
				oneScore = logRatioScores.get(1,1);
			}
			else if(subStr.equals("GB")){
				oneScore = logRatioScores.get(1,2);
			}
			else if(subStr.equals("BR")){
				oneScore = logRatioScores.get(2,0);
			}
			else if(subStr.equals("BG")){
				oneScore = logRatioScores.get(2,2);
			}
			else if(subStr.equals("BB")){
				oneScore = logRatioScores.get(2,2);
			}
			else{
				System.err.println("Error! found unkonw charecter whereas expecting only R, G or B");
			}
			strScore = strScore + oneScore;
			//System.out.println("one score is : " + oneScore);
		}
		if (normalize){
			strScore = strScore/strLength;
		}
		
		return strScore;
	}/*getSequenceScore*/
	
	public Matrix2D getLogLikelihoodRationsMatrix(Matrix2D posTransMatrix, Matrix2D negTranMatrix){//Note: this is what is called Beta Matrix in R Durbin's book page 52
		Matrix2D betaMatrix = new SimpleMatrix2D(3,3);
		for (int i = 0; i<3; i++){
			for (int j = 0; j<3; j++){
				double onePosElement = posTransMatrix.get(i, j);
				double oneNegElement = negTranMatrix.get(i, j);
				double oneLogLikelihoodRatio = log2(onePosElement/oneNegElement);
				betaMatrix.set(i, j, oneLogLikelihoodRatio);
				
			}
		}		
		return betaMatrix;
		
	}/*getLogLikelihoodRationsMatrix*/
	
	public double log2(double x){
		return(Math.log(x)/Math.log(2));
	}/*log2*/
	
	public Matrix2D getTransitionMatrix(String posModelFileName){ //note that here alphabet is length is 3 therefore matrix is 3 by 3
		Matrix2D transitionMatrx = new SimpleMatrix2D(3,3);
		
		try{
			FileInputStream fstream  = new FileInputStream(posModelFileName);
			DataInputStream din      = new DataInputStream(fstream);
			BufferedReader br        = new BufferedReader(new InputStreamReader(din));
			String strLine;
			int row = 0;
			while((strLine = br.readLine()) != null ){
				
				String [] dataStr = strLine.split("\\t");
				for(int i = 0; i < dataStr.length; i++){
					transitionMatrx.set(row, i, Double.parseDouble(dataStr[i]));
					//System.out.print( "row = " + row + " data = " + dataStr[i]+ "--");
				}
				row++;
				if(row>3){
					System.err.println("Error: Expecting a 3 by 3 matrix but found a higher dimension matirx!");
					//System.out.println(" row = " + row);
				}
				//System.out.println();
				
			}
			
		}catch (Exception e){
			System.err.println("Error: " + e.getMessage());
		}
		
		return transitionMatrx;
		
	}/*getTransitionMatrix*/
		
}/*MarkovChainModel*/
