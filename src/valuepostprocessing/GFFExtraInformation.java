/**
 * 
 */
package valuepostprocessing;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;

import org.biojava.bio.program.gff.GFFEntrySet;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.GFFTools;

import java.util.Collections;


/**
 * @author hk3
 * 11-05-2011
 * Given a Gff file, I want number of records, average length and overall length of features
 *
 */
public class GFFExtraInformation {

	/**
	 * @param args
	 */
	
	private String fileName = null;
	
	//SETTERS:
	public void setFileName(String fn){
		this.fileName = fn;
	}/*setFileName*/
	
	
	//GETTERS:
	public String getFileName(){
		return fileName;
	}/*getFileName*/
	
	public static void main(String[] args)  throws Exception{
		// TODO Auto-generated method stub
		String inputFileName = args[0];
		
		GFFExtraInformation app = new GFFExtraInformation();
		app.setFileName(inputFileName);
		
		System.out.println(app.getFileName());
		File gffFile = new File(app.getFileName());
		GFFEntrySet gff  = GFFTools.readGFF(gffFile);
		
		int lineCounter   = 0;
		int overallLength = 0;
		ArrayList<Integer> lengthOfRecords = new ArrayList<Integer>();
		double  average       = 0;
		for (Iterator<?> gffi = gff.lineIterator(); gffi.hasNext();){
			Object o = gffi.next();
			if(o instanceof GFFRecord){
				GFFRecord gro = (GFFRecord) o;
				int oneStrart = gro.getStart();
				int oneEnd    = gro.getEnd();
				int oneLength = oneEnd - oneStrart;
				lengthOfRecords.add(oneLength);
				lineCounter   = lineCounter +1;
				overallLength = overallLength + oneLength;
			}
		}
		average = overallLength/lineCounter;
		
		Object obj_max = Collections.max(lengthOfRecords);
		Object obj_min = Collections.min(lengthOfRecords);
		
		System.out.println("Number of record = " + lineCounter + " average length of records = " + average + " total length of records = " +  overallLength + " maximum legth = " +obj_max + " minimum length = " + obj_min );
		
		
	}/*main*/

}/*GFFExtraInformation*/





