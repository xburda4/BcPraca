package plugins.mBurda.filters;

import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;
import Jama.Matrix;
import icy.image.IcyBufferedImage;
import icy.image.colormodel.IcyColorModel;
import icy.type.DataType;

public class Neuriteness2D{
	private BufferedImage source;
	public IcyBufferedImage ret;
	private double[][] neuriteness2D;
	private double alpha;
	/**@param src grayscale blurred image with one channel
	 * @param alpha float
	 * */
	public Neuriteness2D(BufferedImage src,double alpha){
		this.source = src;
		this.alpha = alpha;
	}
	
	/** Computes neuriteness for 2D grayscale images.Saves it to a matrix neuriteness2D
	 * */
	private void computeNeuriteness2D(){
		Matrix hessian = new Matrix(2,2);
		
		this.neuriteness2D = new double[source.getWidth()][source.getHeight()];
		WritableRaster raster = source.getRaster();
		double minEig = 0,curEig = 0;
		for(int y=0;y<source.getHeight();y++){
			for(int x=0;x<source.getWidth();x++){
				if(x+1<source.getWidth() && x-1>=0){
					//Horizontal approximation
					double a = (raster.getSample(x+1, y, 0)+raster.getSample(x-1, y, 0)-2*raster.getSample(x, y, 0));
					hessian.set(0, 0, a); 	
				} else hessian.set(0, 0, 0);
				if(y+1 < source.getHeight() && y-1>=0)
				{
					//Vertical approximation
					double a = (raster.getSample(x, y+1, 0)+raster.getSample(x, y-1, 0)-2*raster.getSample(x, y, 0));
						hessian.set(1, 1, a);
				} else hessian.set(1, 1, 0);
				if(x+1<source.getWidth() && x-1>=0 && y+1 < source.getHeight() && y-1>=0)
				{
					//Diagonal approximation
					double a = (raster.getSample(x+1, y+1, 0)+raster.getSample(x-1, y-1, 0)-raster.getSample(x+1, y-1, 0)-raster.getSample(x-1, y+1, 0))/4;
						hessian.set(0,1,a);
				} else hessian.set(0, 1, 0);
				hessian.set(1, 0, hessian.get(0,1));
				
				
				double[] eigens = hessian.eig().getRealEigenvalues();
				
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
	
//	public void steerFilter(int sigma){
//		for(int y=0;y<source.getHeight();y++){
//			for(int x=0;x<source.getWidth();x++){
//				for(int xx = 0;xx<2*sigma+1;xx++){
//					double g = Math.exp(-(xx*xx)/(2*sigma*sigma));
//					double gp =  -(xx/sigma)*Math.exp(-(xx*xx)/(2*sigma*sigma));
//					
//					imageX[x][y] = ;
//					imageY[x][y] = ;
//		}
//			}
//			}
//		
//	}
}
