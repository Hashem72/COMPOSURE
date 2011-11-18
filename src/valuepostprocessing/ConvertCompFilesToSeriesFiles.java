/**
 * 
 */
package valuepostprocessing;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import net.derkholm.nmica.matrix.Matrix2D;

import org.biojava.bio.program.gff.GFFEntrySet;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.GFFTools;
import org.biojava.bio.program.gff.GFFWriter;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.FastaFormat;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.Symbol;

import series.SeriesIO;


/**
 * @author hk3  28-03-2011
 *This is to convert each comp file in a directory  into a ser file and print them into another directory. Input arguments for this code  are dir of comp files and a dir that ser files are meant
 *to be written into.
 *ser files are a kind of binary files with a 6 byte 'Series' on top of the file. 
 */
public class ConvertCompFilesToSeriesFiles {

	/**
	 * @param args
	 */
	private String compDirName = null;
	private String serDirName = null;
	
	//setters:
	public void setCompDirName(String comDir){
		this.compDirName = comDir;
	}/**/
	
	public void setSerDirName(String serDir){
		this.serDirName = serDir;
	}/*setSerDirName*/
	
	//getters:
	public String getCompDirName(){
		return compDirName;
	}/*getCompDirName*/
	
	public String getSerDirName(){
		return serDirName;
	}/*getSerDirName*/
	
	public static void main(String[] args)  throws IOException, Exception{
		// TODO Auto-generated method stub
		
		String compDirName = args[0];
		String serDirName  = args[1];
		
		File inputFile = new File(compDirName);
		
		
		Pattern fileNamePattern = Pattern.compile(".comp");
		for(File aCompFile : inputFile.listFiles()){
			Matcher fileNameMatcher = fileNamePattern.matcher(aCompFile.getName());
			if(fileNameMatcher.find()){
				System.out.println("file "+ aCompFile.toString() + "  is being done!");
				
				//make output file name:
				String [] fileName = aCompFile.getName().split("\\.");
				String outputFile = serDirName+"/"+fileName[0]+".ser";
				
				File serFile = new File(outputFile);
				boolean serFileExists = serFile.exists();
				
				if(!serFileExists){
				
				//instantiate a comp object
				ConvertCompFileToSeriesFile comp = new ConvertCompFileToSeriesFile();
				comp.setCompFile(aCompFile.getName());
				String thisFile = comp.getCompFile();
				thisFile        =  compDirName+"/" +thisFile;
				
				
				//convert and write it into the appropreite file:
				comp.setBinFile(outputFile);
				//System.out.println("thisfile is " + thisFile);
				Matrix2D  thisData = comp.readMatrix(thisFile);
				File outPutFile = new File(outputFile);
				SeriesIO.writeSeries(thisData,outPutFile);	
				}
			}
		}
	}/*main*/

}/*ConvertCompFilesToSeriesFiles*/
