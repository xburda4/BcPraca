package plugins.mBurda.filters;

import java.awt.Graphics;
import java.awt.image.BufferedImage;

import org.jdesktop.swingx.image.GaussianBlurFilter;

import Jama.Matrix;
import icy.image.IcyBufferedImage;
import icy.image.colormodel.IcyColorModel;
import icy.type.DataType;


public class Vesselness {
	public IcyBufferedImage source;
	public int scale;
	public IcyBufferedImage ret;
	public double[][] vesselness2D;
	public double[][][] vesselness3D;
	public double beta=0.5,brightThresh=7.5;
	public int ampli = 1,spreadX = 2;
	
	public Vesselness(IcyBufferedImage source,int scale){
		this.source = source;
		this.scale = scale;
	}
	
	
	/*IcyBufferedImage blur2D(){
		
		double[] gauss2D = (new GaussianKernel()).computeGaussian(scale);
		ret.copyData(source);
			for(int y=0;y<ret.getHeight();y++){
				for(int x=0;x<ret.getWidth();x++){
					ret.setData(x, y, 0, 0);
						for(int kX=-scale;kX<=scale;kX++){
							if(x-kX >= 0 && x-kX < source.getWidth()){
							ret.setData(x, y,0, ret.getData(x, y,0)+(source.getData(x-kX, y, 0)*gauss2D[scale+kX]));
							}
					}
				}
		}
			for(int y=0;y<ret.getHeight();y++){
			for(int x=0;x<ret.getWidth();x++){
				ret.setData(x, y, 0, 0);
					for(int kX=-scale;kX<=scale;kX++){
						if(y-kX >= 0 && y-kX < source.getHeight()){
						ret.setData(x, y,0, ret.getData(x, y,0)+(source.getData(x, y-kX, 0)*gauss2D[scale+kX]));
						}
				}
			}
	}
			
		return ret;
	}*/
	
	public IcyBufferedImage getGrayScale(BufferedImage inputImage){
	    BufferedImage img = new BufferedImage(inputImage.getWidth(), inputImage.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
	    Graphics g = img.getGraphics();
	    g.drawImage(inputImage, 0, 0, null);
	    g.dispose();
	    
	    return IcyBufferedImage.createFrom(img);
	}
	
	void computeVesselness2D(){
		Matrix hessian = new Matrix(2,2);
		//IcyBufferedImage blurred = blur2D();
		
		GaussianBlurFilter gauss = new GaussianBlurFilter(scale);
		IcyBufferedImage blurred = getGrayScale(IcyBufferedImage.createFrom(gauss.filter(source, null)));
		this.vesselness2D = new double[source.getWidth()][source.getHeight()];
		
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
				double eig1;
				double eig2;
				if(Math.abs(eigens[0]) <= Math.abs(eigens[1])){
					eig1= eigens[0];
					eig2= eigens[1];
				} else {
					eig1= eigens[1];
					eig2= eigens[0];
				}
				double a = Math.exp(-((eig1*eig1)/(2*beta*beta*eig2*eig2)));
				double b = (1 - Math.exp(-((eig1*eig1+eig2*eig2)/2*brightThresh*brightThresh)));
				this.vesselness2D[x][y] =  a * b;
				}
			}
	}
	
public void makeImage(){
		computeVesselness2D();
		ret = new IcyBufferedImage(source.getWidth(),source.getHeight(),IcyColorModel.createInstance(1, DataType.DOUBLE));
		for(int y=0;y<source.getHeight();y++){
			for(int x=0;x<source.getWidth();x++){
				ret.setData(x, y, 0, vesselness2D[x][y]);
			}
		}
	 }
}
