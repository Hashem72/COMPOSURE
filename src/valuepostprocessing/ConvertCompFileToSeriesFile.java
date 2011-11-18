/**
 * 
 */
package valuepostprocessing;



import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import series.SeriesIO;

import net.derkholm.nmica.matrix.Matrix2D;
import net.derkholm.nmica.matrix.SimpleMatrix2D;




/**
 * @author hk3 24-03-2011
 *
 */
public class ConvertCompFileToSeriesFile {

	/**
	 * @param args
	 */
	
	private String compFile;
	private String binFile;
	
	//setters
	public void setCompFile(String cf){
		this.compFile = cf;
	}/*setCompFile*/
	
	public void setBinFile(String bf){
		this.binFile = bf;
	}/*setBinFile*/
	
	//getters
	public String getCompFile(){
		return compFile;
	}/*getCompFile*/
	
	
	public String getBinFile(){
		return binFile;
	}/*getBinFile*/
	
	
	public int getFileNumberOfLines(String file) throws IOException{
	    InputStream is = new BufferedInputStream(new FileInputStream(file));
	    try {
	        byte[] c = new byte[1024];
	        int count = 0;
	        int readChars = 0;
	        while ((readChars = is.read(c)) != -1) {
	            for (int i = 0; i < readChars; ++i) {
	                if (c[i] == '\n')
	                    ++count;
	            }
	        }
	        return count;
	    } finally {
	        is.close();
	    }

	}/*getFileNumberOfLines*/
	
	public Matrix2D readMatrix(String file) throws Exception{
		int numberOfRows = getFileNumberOfLines(file);
		int numberOfCols = 3;
		Matrix2D matrix = new SimpleMatrix2D(numberOfRows,numberOfCols);
		try{
			FileInputStream fs = new FileInputStream(file);
			BufferedReader br =  new BufferedReader(new InputStreamReader(fs));
			for(int l=0; l<numberOfRows; l++){
				String thisline = br.readLine();
				String [] states = thisline.split("\t");
				
				states = thisline.split("\t");
				for(int k= 0;k<numberOfCols; k++){
					matrix.set(l, k, Double.parseDouble(states[k+1]));//note: first colm is the line number therefore is discarded
				}
			}

		}
		catch(FileNotFoundException fe){
			fe.fillInStackTrace();
		}
		return matrix;
		
	}/*readMatrix*/
		
		
	
	 
		
	public static void main(String[] args) throws IOException,Exception {
		// TODO Auto-generated method stub
		String compFile = args[0];
		String binFile  = args[1];
		
		ConvertCompFileToSeriesFile comp = new ConvertCompFileToSeriesFile();
		
		comp.setCompFile(compFile);
		comp.setBinFile(binFile);
		
		int numberOfLines = comp.getFileNumberOfLines(compFile);
		Matrix2D data = comp.readMatrix(compFile);
		File outPutFile = new File(binFile);
		SeriesIO.writeSeries(data, outPutFile);
	}/*main*/

}/*ConvertCompFileToSeriesFile*/
