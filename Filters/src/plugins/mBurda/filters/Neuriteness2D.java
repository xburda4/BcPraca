package plugins.mBurda.filters;

import java.awt.image.BufferedImage;
import icy.image.IcyBufferedImage;
import icy.image.colormodel.IcyColorModel;
import icy.type.DataType;

public class Neuriteness2D extends Filter{
	public IcyBufferedImage ret;
	private double[][] neuriteness2D;
	private double alpha;
	
	/**@param src grayscale blurred image with one channel
	 * @param alpha float
	 * */
	public Neuriteness2D(BufferedImage src,double alpha){
		Filter.source = IcyBufferedImage.createFrom(src);
		this.alpha = alpha;
		this.neuriteness2D = new double[source.getWidth()][source.getHeight()];
	}
	
	/** Computes neuriteness for 2D grayscale images.Saves it to a matrix neuriteness2D
	 * */
	private void computeNeuriteness2D(){
		
		double minEig = 0,curEig = 0;
		for(int y=0;y<source.getHeight();y++){
			for(int x=0;x<source.getWidth();x++){
				double[] eigens = getEigenValues(x,y);
				
				double psi1 = eigens[0]+alpha*eigens[1];
				double psi2 = eigens[1]+alpha*eigens[0];
				
				//Max absolute value
				if(Math.abs(psi1)>Math.abs(psi2)) curEig = psi1;
				else curEig = psi2;
				//Get minimal eigenvalue through the image
				if(curEig < minEig || (x == 0 && y == 0)) minEig = curEig;
				//Save current max eigenvalue at x,y coordinates
				neuriteness2D[x][y] = curEig;
				}
			}
		for(int x=0;x<neuriteness2D.length;x++){
			for(int y=0;y<neuriteness2D[x].length;y++){
				curEig = neuriteness2D[x][y];
				//Counting of neuriteness
				if(curEig < 0) neuriteness2D[x][y] = curEig/minEig;
				else neuriteness2D[x][y] = 0;
			}
		}
	}
	
	
	private void computeNeuritenessWithPhase2D(){
		
		double minEig = 0,curEig = 0;
		for(int y=0;y<source.getHeight();y++){
			for(int x=0;x<source.getWidth();x++){
				double[] eigens = getEigenValuesWithPhase(x,y);
				
				double psi1 = eigens[0]+alpha*eigens[1];
				double psi2 = eigens[1]+alpha*eigens[0];
				
				//Max absolute value
				if(Math.abs(psi1)>Math.abs(psi2)) curEig = psi1;
				else curEig = psi2;
				//Get minimal eigenvalue through the image
				if(curEig < minEig || (x == 0 && y == 0)) minEig = curEig;
				//Save current max eigenvalue at x,y coordinates
				neuriteness2D[x][y] = curEig;
				}
			}
		for(int x=0;x<neuriteness2D.length;x++){
			for(int y=0;y<neuriteness2D[x].length;y++){
				curEig = neuriteness2D[x][y];
				//Counting of neuriteness
				if(curEig < 0) neuriteness2D[x][y] = curEig/minEig;
				else neuriteness2D[x][y] = 0;
			}
		}
	}
	
	/** Rewrites data from neuriteness2D matrix to image ret.
	 * */
	public void makeImage2D(){
		computeNeuriteness2D();
		ret = new IcyBufferedImage(source.getWidth(),source.getHeight(),IcyColorModel.createInstance(1, DataType.DOUBLE));
		ret.beginUpdate();
		for(int y=0;y<source.getHeight();y++){
			for(int x=0;x<source.getWidth();x++){
				ret.setData(x, y, 0, neuriteness2D[x][y]);
			}
		}
		ret.endUpdate();
	 }
	
	public void makeImageWithPhase2D(){
		computeNeuritenessWithPhase2D();
		ret = new IcyBufferedImage(source.getWidth(),source.getHeight(),IcyColorModel.createInstance(1, DataType.DOUBLE));
		ret.beginUpdate();
		for(int y=0;y<source.getHeight();y++){
			for(int x=0;x<source.getWidth();x++){
				ret.setData(x, y, 0, neuriteness2D[x][y]);
			}
		}
		ret.endUpdate();
	 }
	
}
