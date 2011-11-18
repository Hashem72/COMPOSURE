/**
 * 
 */
package valuepostprocessing;

import java.io.File;
import java.io.FileWriter;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.biojava.bio.program.gff.GFFEntrySet;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.GFFTools;
import org.biojava.bio.program.gff.GFFWriter;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.bjv2.util.IntSet;
import org.bjv2.util.IntSets;
import org.bjv2.util.SmallMap;

//import composure2.ChainFinder2;

import series.SeriesIO;

import net.derkholm.nmica.matrix.Matrix2D;
import net.derkholm.nmica.matrix.SimpleMatrix2D;
import net.derkholm.nmica.utils.CliTools;


/**
 * @author hk3 25-03-2011
 * This is a modification of the ChainFider2 by Thomas Down. The Key Difference here is that it would be able to write output into a gff file if it is asked to do so.
 * Note: the input for this Class must be a series file (a binary file with 6 bytes "Series" on top of that).
 *
 */
public class ChainFinder3 {

	private File comp; 
	private File seed;
	private String output= null;
	private double maxC2 = 0.5;//originally was 0.4
	private double minC0 = 0.9;
	private double minC1 = 0.7; //orginally was 0.6
	private int coverageC1 = 1;//originally was 1
	private boolean blocks = false;
	private String chainPrefix = "chain";
	private int chainLengthThreshold = 100; // this is to ignore chains shorter than this threshold
	
	public void setChainLengthThreshold(int clt){
		this.chainLengthThreshold = clt;
	}
	
	public void setBlocks(boolean b) {
		this.blocks = b;
	}

	public void setChainPrefix(String s) {
		this.chainPrefix = s;
	}
	
	public void setCoverageC1(int i) {
		this.coverageC1 = i;
	}
	
	public void setComp(File comp) {
		this.comp = comp;
	}

	public void setMaxC2(double maxC2) {
		this.maxC2 = maxC2;
	}

	public void setMinC0(double minC0) {
		this.minC0 = minC0;
	}

	public void setMinC1(double minC1) {
		this.minC1 = minC1;
	}

	public void setSeed(File seed) {
		this.seed = seed;
	}
	
	public void setOutput(String op){
		this.output = op;
	}/*setOutput*/
	
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

	

	/**
	 * @param args
	 */
	public static void main(String[] args) 
		throws Exception
	{
		ChainFinder3 app = new ChainFinder3();
		args = CliTools.configureBean(app, args);
		app.run(args);
	}
	
	private void run(String[] args)
		throws Exception
	{
		
		
		Matrix2D comp = SeriesIO.readSeries(this.comp);
		
		
		GFFEntrySet seeds;
		if (this.seed != null) {
			seeds = GFFTools.readGFF(this.seed);
		} else {
			int offset = 1;
			Matrix2D seedm = new SimpleMatrix2D(comp.rows(), 1);
			for (int i = 0; i < comp.rows() - offset; ++i) {
				seedm.set(i, 0, Math.max(0, (comp.get(i, 0) - comp.get(i + offset, 0)) * (comp.get(i + offset, 1) - comp.get(i, 1))));
			}
			double seedThreshold =0.001; // originally was:  0.001; 
		    seeds = new GFFEntrySet();
		    int blockStart = -1;
		    for (int i = 0; i < seedm.rows(); ++i) {
		    	if (seedm.get(i, 0) >= seedThreshold) { 
		    		if (blockStart < 0) {
		    			blockStart = i;
		    		}
		    	} else {
		    		if (blockStart > 0) {
		    			SimpleGFFRecord r = new SimpleGFFRecord();
		    			r.setStart(blockStart);
		    			r.setEnd(i - 1);
		    			seeds.add(r);
		    			blockStart = -1; 
		    		}
		    	}/*else*/
		    	
		    }
		}
	    IntSet covered = IntSets.none();
	    
	    PrintWriter pw;	    
	    if(output != null){ //given the output file name, the data must be printed into it.
	    	// but first delete the old file if exists
	    	deleteFile(this.output);

	    	File gffoutputFile  = new File(this.output);	    	
	    	 pw = new PrintWriter(new FileWriter(gffoutputFile));
	    }
	    else{	    
	     pw = new PrintWriter(new OutputStreamWriter(System.out));
	    }
	    
	    GFFWriter gffw = new GFFWriter(pw);
	    
	    int chainIdSeed = 0; 
	    for (Iterator<?> i = seeds.lineIterator(); i.hasNext(); ) {
	    	Object o = i.next();
	    	if (o instanceof GFFRecord) {
	    		GFFRecord seedRecord = (GFFRecord) o;
	    		IntSet seed = IntSets.range(seedRecord.getStart(), seedRecord.getEnd());
	    		boolean overlaps = false;
	    		for (int ii : seed) {
	    			if (covered.contains(ii)) {
	    				overlaps=true;
	    			}
	    		}
	    		if (!overlaps) {
		    		int chainMin = seed.getMax();
		    		while (comp.get(chainMin, 2) < maxC2 && chainMin > 1) {
		    			--chainMin;
		    		}
		    		++chainMin;
		    		int chainMax = seed.getMin();
		    		while (comp.get(chainMax, 2) < maxC2 && chainMax < comp.rows() - 1) {
		    			++chainMax;
		    		}
		    		--chainMax;
		    		if (chainMax > chainMin) {
		    			IntSet chain = IntSets.range(chainMin, chainMax);
		    			if (chanMax(comp, chain, 0) >= minC0 && chanOverThreshold(comp, chain, 1, minC1) >= coverageC1) { // NOTE: when coverageC1 is considered as 1, then the condition one and 
		    				// two in this if clause are equivalent, apart from the fact that the first is on chan 0 and the second is on chan 1; 
		    				String chainId = String.format("%s.%05d", chainPrefix, ++chainIdSeed);
		    				double score = chanMax(comp, chain, 1);
		    				
		    				SimpleGFFRecord r = new SimpleGFFRecord();
		    				r.setSource("composure");
		    				r.setFeature("chain");
		    				r.setStart(chain.getMin());
		    				r.setEnd(chain.getMax());
		    				r.setScore(score);
		    				Map<String,List<String>> gaga = new SmallMap<String,List<String>>();
		    				gaga.put("chain.id", Collections.singletonList(chainId));
		    				r.setGroupAttributes(gaga);
		    				
		    				//ignore very short chains
		    				int chainStart = r.getStart();
		    				int chainEnd   = r.getEnd();
		    				int chainlength = chainEnd-chainStart;
		    				if(chainlength > chainLengthThreshold){
		    				gffw.recordLine(r);
		    				}
		    				
		    				covered = IntSets.union(covered, chain);

		    				
		    				if (blocks) {
			    				int blockIdSeed = 0;
			    				int blockStart = -1;
			    				for (int bpos = chainMin; bpos <= chainMax; ++bpos) {
			    					if (comp.get(bpos, 0) >= minC1) {
			    						if (blockStart < 0) {
			    							blockStart = bpos;
			    						}
			    					} else {
			    						if (blockStart >= 0) {
			    							emitBlockRecord(gffw, chainId, String.format("%s.%d", chainId, ++blockIdSeed), blockStart, bpos - 1, score);
			    							blockStart = -1;
			    						}
			    					}
			    				}
			    				if (blockStart >= 0) {
			    					emitBlockRecord(gffw, chainId, String.format("%s.%d", chainId, ++blockIdSeed), blockStart, chainMax, score);
			    				}
		    				}
		    			}
		    		}
	    		}/*if !overlaps*/
	    	}/*o instanceof*/
	    }/*Iterator i*/
	    
	    pw.flush();
	}
	
	

	private void emitBlockRecord(GFFWriter gffw, String chainId, String blockId, int blockStart, int blockEnd, double score) 
		throws Exception
	{
		SimpleGFFRecord r = new SimpleGFFRecord();
		r.setSource("composure");
		r.setFeature("block");
		r.setStart(blockStart);
		r.setEnd(blockEnd);
		r.setScore(score);
		Map<String,List<String>> gaga = new SmallMap<String,List<String>>();
		gaga.put("chain.id", Collections.singletonList(chainId));
		gaga.put("block.id", Collections.singletonList(blockId));
		r.setGroupAttributes(gaga);
		gffw.recordLine(r);
	}

	private double chanMax(Matrix2D m, IntSet target, int channel)
	{
		double max = Double.NEGATIVE_INFINITY;
		for (int i : target) {
			max = Math.max(max, m.get(i, channel));
		}
		return max;
	}
	
	private int chanOverThreshold(Matrix2D m, IntSet target, int channel, double threshold) {
		int cnt = 0;
		for (int i : target) {
			if (m.get(i, channel) >= threshold) {
				++cnt;
			}
		}
		return cnt;
	}



}/*ChainFinder3*/
