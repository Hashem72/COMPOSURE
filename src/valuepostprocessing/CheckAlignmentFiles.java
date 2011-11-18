/**
 * 
 */
package valuepostprocessing;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.NoSuchElementException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.io.SeqIOTools;

import net.derkholm.nmica.utils.CliTools;

/**
 * @author hk3
 *
 */
public class CheckAlignmentFiles {

	/**
	 * @param args,
	 * 02-06-2011 by Hashem Koohy
	 * 
	 * for a given alignment block coordinates (length of dmel sequences including gaps) and length of real alignment blocks (excluding gaps, which in this analysis is given as chr_min_max), max-min must
	 * be equal to length of alignment blocks - number of gaps. But I have found out this is not true! I dont know why, maybe something wrong with aligment. By the way this is to check, and if this
	 * criteria is not satisfied, delete comp file and ser file! this will prevent getting wrong information.
	 */
	
	private String aliDir= "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ALIGNMENTS_EISENLAB/2L";
	private String compDir =  "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/Eisenlab_Composure_Output_Model_20051212_With_Window_Length_350";
	private String serDir = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ANALYSIS_EISENLAB/2L/Eisenlab_Composure_Output_Model_20051212_With_Window_Length_350_SER";
	
	
	public void setAliDir (String ad){
		this.aliDir = ad;
	}/*setAliDir*/
	
	public void setCompDir(String cd){
		this.compDir = cd;
	}/*setCompDir*/
	
	public void setSerDir(String sd){
		this.serDir = sd;
	}/*setSerDir*/
	
	public String  getAliDir(){
		return aliDir;
	}/*getAliDir*/
	
	public String getCompDir(){
		return compDir;
	}/*getCompDir*/
	
	public String getSerDir(){
		return serDir;
	}/*getSerDir*/
	
	public static void main(String[] args)  throws Exception{
		// TODO Auto-generated method stub
		CheckAlignmentFiles app = new CheckAlignmentFiles();
		args = CliTools.configureBean(app, args);
		app.run(args);

	}/*main*/
	
	private void run(String [] args) throws NoSuchElementException, FileNotFoundException, BioException{
		
		Pattern fnPattern = Pattern.compile("(([^_]+)_([0-9]+)_([0-9]+)).fa");
		File aliDir = new File(getAliDir());
		
		for(File aliFile: aliDir.listFiles()){
			if(aliFile.getName().startsWith(".")){
				continue;
			}
			
            Matcher fnMatcher = fnPattern.matcher(aliFile.getName());
            if(!fnMatcher.matches()){
            	System.out.println("Bad File Name!");
            }
            String  baseName = fnMatcher.group(1);
            String chr = fnMatcher.group(2);
            int    min = Integer.parseInt(fnMatcher.group(3));
            int    max = Integer.parseInt(fnMatcher.group(4));
            int    lengthofBlock = max- min;
            Sequence seq = SeqIOTools.readFastaDNA(new BufferedReader(new FileReader(aliFile))).nextSequence();
            int aliLength = seq.length();
            if( aliLength < lengthofBlock){
            	//System.out.println("for block = " + baseName + " blokc length is "+ lengthofBlock + " which is longer than aliLength which is " + aliLength);
            	String compFileName = getCompDir()+"/"+ baseName+".comp";
            	File compFile = new File(compFileName);
            	String serFileName  = getSerDir()+"/"+baseName+".ser";
            	File serFile   = new File(serFileName);
            	if(compFile.exists()){
            		compFile.delete();
            		System.out.println(" file " + compFile.getName() + " deleted!");	
            	}
            	if(serFile.exists()){
            		serFile.delete();
            		System.out.println("file " + serFile.getName() + " deleted!");
            	}
            	
            	
            	//File compFile = new File()
            }
		}
		
	}/*run*/

}/*CheckAlignmentFiles*/
