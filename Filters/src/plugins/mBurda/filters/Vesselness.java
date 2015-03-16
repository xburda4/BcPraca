package plugins.mBurda.filters;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;

import org.jdesktop.swingx.image.GaussianBlurFilter;

import Jama.Matrix;
import icy.image.IcyBufferedImage;
import icy.image.colormodel.IcyColorModel;
import icy.type.DataType;


public class Vesselness{
	public IcyBufferedImage source;
	private int scale = 0;
	public IcyBufferedImage ret;
	public double[][] vesselness2D;
	public double[][][] vesselness3D;
	public double beta,brightThresh;
	//public int ampli = 1,spreadX = 2;
	
	public Vesselness(IcyBufferedImage source,int scale,double beta,double bright){
		this.source = source;
		this.scale = scale;
		this.beta = beta;
		this.brightThresh = bright;
	}
	public Vesselness(IcyBufferedImage source){
		this.source = source;
	}
	
	private BufferedImage getGrayScale(BufferedImage original){
		BufferedImage image = new BufferedImage(original.getWidth(), original.getHeight(),  
			    BufferedImage.TYPE_BYTE_GRAY);  
			Graphics g = image.getGraphics();  
			g.drawImage(original, 0, 0, null);  
			g.dispose(); 
			return image;
	}
	
	private void computeVesselness2D(int scale){
		Matrix hessian = new Matrix(2,2);
		//IcyBufferedImage blurred = blur2D();
		
//		GaussianBlurFilter gauss = new GaussianBlurFilter(scale);
//		BufferedImage blurred = getGrayScale(gauss.filter(source, null));
		IcyBufferedImage blurred = source;
		this.vesselness2D = new double[source.getWidth()][source.getHeight()];
		WritableRaster raster = blurred.getRaster();
		for(int y=0;y<blurred.getHeight();y++){
			for(int x=0;x<blurred.getWidth();x++){
				if(x+1<blurred.getWidth() && x-1>=0){
					
					double a = (raster.getSample(x+1, y, 0)+raster.getSample(x-1, y, 0)-2*raster.getSample(x, y, 0));
					hessian.set(0, 0, a); 	
				} else hessian.set(0, 0, 0);
				if(y+1 < blurred.getHeight() && y-1>=0)
				{
					double a = (raster.getSample(x, y+1, 0)+raster.getSample(x, y-1, 0)-2*raster.getSample(x, y, 0));
						hessian.set(1, 1, a);
				} else hessian.set(1, 1, 0);
				if(x+1<blurred.getWidth() && x-1>=0 && y+1 < blurred.getHeight() && y-1>=0)
				{
					double a = (raster.getSample(x+1, y+1, 0)+raster.getSample(x-1, y-1, 0)-raster.getSample(x+1, y-1, 0)-raster.getSample(x-1, y+1, 0))/4;
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
	
public void makeImage2D(){
		computeVesselness2D(scale);
		ret = new IcyBufferedImage(source.getWidth(),source.getHeight(),IcyColorModel.createInstance(1, DataType.DOUBLE));
		for(int y=0;y<source.getHeight();y++){
			for(int x=0;x<source.getWidth();x++){
				ret.setData(x, y, 0, vesselness2D[x][y]);
			}
		}
	 }

public void setScale(int scale){
	this.scale = scale;
}
}
