/**
 * 
 */
package valuepostprocessing;

/**
 * @author hk3
 *13-04-2011
 */
public class Pair {
	private double firstComponent  = 0.0;
	private double secondComponent = 0.0;
	
	//SETTERS
	public void setFirstComponent(double fc){
		this.firstComponent = fc;
	}/*setFirstComponent*/
	
	public void setSecondComponent(double sc){
		this.secondComponent = sc;
	}/*setSecondComponent*/
	
	//GETTERS
	public double getFirstComponent(){
		return firstComponent;
	}/*getFirstComponent*/
	
	public double getSecondComponent(){
		return secondComponent;
	}/*getSecondComponent*/
	

}/*Pair*/
