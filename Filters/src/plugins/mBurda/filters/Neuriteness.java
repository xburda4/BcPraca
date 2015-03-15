package plugins.mBurda.filters;

import java.awt.Graphics;
import java.awt.image.BufferedImage;

import org.jdesktop.swingx.image.GaussianBlurFilter;

import Jama.Matrix;

import icy.image.IcyBufferedImage;
import icy.image.colormodel.IcyColorModel;
import icy.type.DataType;

public class Neuriteness {
	private IcyBufferedImage source;
	private int scale;
	public IcyBufferedImage ret;
	private double[][] neuriteness2D;
	private double[][][] neuriteness3D;
	private double alpha = 4;
	
	public Neuriteness(IcyBufferedImage src,int scale){
		this.scale = scale;
		this.source = src;
	}
	
	private IcyBufferedImage getGrayScale(BufferedImage inputImage){
	    BufferedImage img = new BufferedImage(inputImage.getWidth(), inputImage.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
	    Graphics g = img.getGraphics();
	    g.drawImage(inputImage, 0, 0, null);
	    g.dispose();
	    
	    return IcyBufferedImage.createFrom(img);
	}
	
	private void computeNeuriteness2D(){
		Matrix hessian = new Matrix(2,2);
		
		GaussianBlurFilter gauss = new GaussianBlurFilter(scale);
		IcyBufferedImage blurred = getGrayScale(IcyBufferedImage.createFrom(gauss.filter(source, null)));
		this.neuriteness2D = new double[source.getWidth()][source.getHeight()];
		
		double minEig = 0,curEig = 0;
		for(int y=0;y<blurred.getHeight();y++){
			for(int x=0;x<blurred.getWidth();x++){
				if(x+1<blurred.getWidth() && x-1>=0){
					double a = (blurred.getData(x+1, y, 0)+blurred.getData(x-1, y, 0)-2*blurred.getData(x, y, 0));
					hessian.set(0, 0, a); 	
				} else hessian.set(0, 0, 0);
				if(y+1 < blurred.getHeight() && y-1>=0)
				{
					double a = (blurred.getData(x, y+1, 0)+blurred.getData(x, y-1, 0)-2*blurred.getData(x, y, 0));
						hessian.set(1, 1, a);
				} else hessian.set(1, 1, 0);
				if(x+1<blurred.getWidth() && x-1>=0 && y+1 < blurred.getHeight() && y-1>=0)
				{
					double a = (blurred.getData(x+1, y+1, 0)+blurred.getData(x-1, y-1, 0)-blurred.getData(x+1, y-1, 0)-blurred.getData(x-1, y+1, 0))/4;
						hessian.set(0,1,a);
				} else hessian.set(0, 1, 0);
				hessian.set(1, 0, hessian.get(0,1));
				
				double[] eigens = hessian.eig().getRealEigenvalues();
				curEig = Math.max(eigens[0]+alpha*eigens[1], eigens[1]+alpha*eigens[0]);
				if(x == 0 && y == 0) minEig = curEig;
				else if(curEig < minEig) minEig = curEig;
				neuriteness2D[x][y] = curEig;
				}
			}
		curEig = 0;
		for(int x=0;x<neuriteness2D.length;x++){
			for(int y=0;y<neuriteness2D[x].length;y++){
				curEig = neuriteness2D[x][y];
				if(curEig < 0) neuriteness2D[x][y] = curEig/minEig;
				else neuriteness2D[x][y] = 0;
			}
		}
	}
	
	public void makeImage2D(){
		computeNeuriteness2D();
		ret = new IcyBufferedImage(source.getWidth(),source.getHeight(),IcyColorModel.createInstance(1, DataType.DOUBLE));
		for(int y=0;y<source.getHeight();y++){
			for(int x=0;x<source.getWidth();x++){
				ret.setData(x, y, 0, neuriteness2D[x][y]);
			}
		}
	 }
}
