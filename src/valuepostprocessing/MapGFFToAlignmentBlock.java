/**
 * 
 */
package valuepostprocessing;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.program.gff.GFFEntrySet;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.GFFTools;
import org.biojava.bio.program.gff.GFFWriter;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.Symbol;










/**
 * @author hk3
 * Note: this is to map some features such as real enhancers or genes or chromatin accessibile positions which are  in gff format into an alignment file (block).
 * Note:  the mapGffToAlignment method in this class is a modification of the  run method in  data.fly.MapGFFToAlignment written by Thomas Down.
 *
 */
public class MapGFFToAlignmentBlock {
	//instances
	private int blockNumber;
	private String mapFile ;
	private String enhancersInputGffFile;
	private String dnaseIInputGffFile;
	private String genesInputGfffile;
	private String enhancersOutputFile;
	private String dnaseIOutputFile;
	private String genesOutputFile;
	private String alignmentFile;
	private String targetName = "dmel";

	
	//geters
	public int getBlockNumber(){
		return blockNumber;
	}/*getBlockNumber*/
	
	public String getAlignmentFile(){
		return alignmentFile;
	}/*getAlignmentFile*/
	
	
	public String getMapFile(){
		return mapFile;
	}/*getMapFile*/
	
	public String getEnhancersInputGffFile(){
		return enhancersInputGffFile;
	}/*getEnhancersInputGffFile*/
	
	public String getDnaseIInputGffFile(){
		return dnaseIInputGffFile;
	}/*getDnaseIInputGffFile*/

	public String getGenesIputGffFile(){
		return genesInputGfffile;
	}/*getGenesIputGffFile*/
		
	public String getEnhancersOutputFiel(){
			return enhancersOutputFile;
		}/*getEnhancersOutputFiel*/
		
	public String getGenesOutputFile(){
			return genesOutputFile;
		}/*getGenesOutputFile*/
	public String getDnaseIOutputFile(){
		return dnaseIOutputFile;
	}/*getDnaseIOutputFile*/
	
	//setters
	public void setMapFile (String mf){
		this.mapFile = mf;
	}/*setMapFile*/
	
	public void setBlockNumber(int b){
		this.blockNumber = b;
	}/*setBlockNumber*/
	
	public void setAlignmentFile(String af){
		this.alignmentFile = af;
	}/*setAlignmentFile*/
	
	
	public void setEnhancersInputGffFile(String eif){
		this.enhancersInputGffFile = eif;
	}/*setEnhancersInputGffFile*/
	
	public void setEnhancersOutputFile(String eof){
		this.enhancersOutputFile = eof;
	}/*setEnhancersOutputFile*/
	
	
	
	public void setDnaseIInputGffFile(String dif){
		this.dnaseIInputGffFile = dif;
	}/*setDnaseIInputGffFile*/
	
	public void setDnaseIOutputFile(String dof){
		this.dnaseIOutputFile = dof;
	}/*setDnaseIOutputFile*/
	
	public void setGenesInputGfffile(String gif){
		this.genesInputGfffile = gif;
	}/*setGenesInputGfffile*/
	public void setGenesOutputFile(String gof){
		this.genesOutputFile = gof;
	}/*setGenesOutputFile*/
	
	//delet file
	public void deleteFile(String file){
		File f = new File(file);
		boolean exists = f.exists();
		if(exists){
			boolean success = f.delete();
			if(!success){
				System.out.println("Deletion failed.");
				System.exit(0);
			}
		}
	}/*deleteFile*/

	
	//get block specifications from the map file. Note that the block number is the line number in the map file.
	public String [] getSpecifications(String file, int lineNumber) throws IOException{
		String [] specifications = null;
		try{
			FileInputStream fs = new FileInputStream(file);
			BufferedReader br = new BufferedReader(new InputStreamReader(fs));
			String lineOfInterest = null;
			for(int i=0; i<lineNumber; i++){
				lineOfInterest = br.readLine();
			}
			specifications = lineOfInterest.split("\t");
		}
		catch(FileNotFoundException fe){
			fe.fillInStackTrace();	
	}/*getSpecifications*/
		return specifications;
}/*getSpecifications*/
	
	
    
    public void mapGffToAlignment(String [] blockSpecifications, String inputGffFile,  String alignmentFile, String outputGffFile) throws IOException,Exception{
    	
    	String chr = blockSpecifications[1];
    	int min    = Integer.parseInt(blockSpecifications[2]);
    	int max    = Integer.parseInt(blockSpecifications[3]);
    	File gffFile        = new File(inputGffFile);
    	File aliFile        = new File(alignmentFile);
    	File gffoutputFile  = new File(outputGffFile);
    	
        GFFEntrySet gff = GFFTools.readGFF(new BufferedReader(new FileReader(gffFile)));
        
        Map<GFFRecord,Mapping> localGff = new HashMap<GFFRecord,Mapping>();
        Set<Integer> interesting = new HashSet<Integer>();
        for (Iterator<?> gffi = gff.lineIterator(); gffi.hasNext(); ) {
            Object o = gffi.next();
            if (o instanceof GFFRecord) {
                GFFRecord gro = (GFFRecord) o;
                     //if ((chr == null || gro.getSeqName().equals(chr)) && gro.getEnd() >= min && gro.getStart() <= max) {
                    //Mapping m = new Mapping(Math.max(min, gro.getStart()), Math.min(max, gro.getEnd())); //WARNING:this two lines changed by hashem in order to get m.aMax being updated for the min point! hasehm 19-04-2011
                if ((chr == null || gro.getSeqName().equals(chr)) && gro.getEnd() > min && gro.getStart() < max) {
                	Mapping m = new Mapping(Math.max(min+1, gro.getStart()), Math.min(max, gro.getEnd())); 
                	interesting.add(m.rMax);
                	interesting.add(m.rMin);
                    localGff.put(gro, m);
                }
            }
        }
        
        if (localGff.size() > 0) {
            Sequence seq = null;
            SequenceIterator seqI = SeqIOTools.readFastaDNA(new BufferedReader(new FileReader(aliFile)));
            while (seqI.hasNext()) { 
            	Sequence ss = seqI.nextSequence(); 
            	if (ss.getName().equals(targetName)) {
            		seq = ss;
            		break;
            	}
            }
            if (seq == null) {
            	System.err.println("Couldn't find " + seq);
            	return;
            }
            Symbol gap = seq.getAlphabet().getGapSymbol();
            // int gPos = min - 1;
            int gPos = min;     // WARNING: Quick workaround for an untraced off-by-one error.
  
           for (int aPos = 1; aPos <= seq.length(); ++aPos) {
                if (seq.symbolAt(aPos) != gap) {
                    ++gPos;
                    if (interesting.contains(gPos)) {
	                    for (Mapping m : localGff.values()) {
	                        if (m.rMin == gPos) {
	                            m.aMin = aPos;
	                        }
	                        if (m.rMax == gPos) {
	                            m.aMax = aPos;
	                        }
	                    }
                    }
                }
            }
        }

        PrintWriter pw = new PrintWriter(new FileWriter(gffoutputFile));
        GFFWriter gffw = new GFFWriter(pw);
        
        
        //temp commented allowed to print out for debugging purposes! must me commented out
         //pw = new PrintWriter(new OutputStreamWriter(System.out));
         //gffw = new GFFWriter(pw);


        for (Map.Entry<GFFRecord,Mapping> me : localGff.entrySet()) {
            SimpleGFFRecord r2 = new SimpleGFFRecord(me.getKey());
            r2.setSeqName(targetName);
            r2.setStart(me.getValue().aMin);
            r2.setEnd(me.getValue().aMax);
            if(r2.getStart() < 0 || r2.getEnd() <0){
            	System.out.println("seqName = " + me.getKey().getSeqName()+ " start = " + me.getKey().getStart() + " end = "+ me.getKey().getEnd());
            	System.out.println("WARNING: SOME GFF RECORDS ARE NOT MAPPED PROPERLY!");
            }
            // r2.setScore(100.0);
             gffw.recordLine(r2);
            
           
            //this  if has been put here for debugging purposes Hashem - 15-04-2011
           // if(r2.getStart() < 0){
            //	System.out.println("WARNING: In MapGFFToAlignmentBlock start of a feature became negative and therefore ignored");
            	//System.out.println("name " + r2.getFeature()+ " start " + r2.getStart()+ " end " + r2.getEnd());           	
            //}
            
        }
        pw.flush();

    }/*mapGffToAlignment*/

	

    private static class Mapping {
        public final int rMin, rMax;
        public int aMin = -1, aMax = -1;
        
        public Mapping(int min, int max)
        {
            this.rMin = min;
            this.rMax = max;
        }/*Mapping-Constructor*/
    }/*Mapping-class*/

    
	
}/*MapGFFToAlignmentBlock*/
