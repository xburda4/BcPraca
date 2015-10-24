package plugins.mBurda.filters;


import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;
import Jama.Matrix;
import icy.image.IcyBufferedImage;
import icy.image.colormodel.IcyColorModel;
import icy.type.DataType;


public class Vesselness2D{
	public BufferedImage source;
	public IcyBufferedImage ret;
	public double[][] vesselness2D;
	public double beta,brightThresh;
	public double[][] phaseCong; //prepísať na static + triedu z ktorej budú vesselness a neuriteness dediť
	
	/**@param source blurred grayscale image with one channel
	 * @param beta,bright thresholds controling sensitivity of line measurement
	 * */
	public Vesselness2D(BufferedImage source,double beta,double bright){
		this.source = source;
		this.beta = beta;
		this.brightThresh = bright;
	}
	public Vesselness2D(BufferedImage source){
		this.source = source;
		beta = 1.5;
		brightThresh = 0.08;
	}
	
	/**Computes vesselness for 2D grayscale images. Saves it to a matrix vesselness2D
	 * */
	private void computeVesselness(){
		Matrix hessian = new Matrix(2,2);
		
		
		WritableRaster raster = source.getRaster();
		this.vesselness2D = new double[source.getWidth()][source.getHeight()];
		
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
				double eig1;
				double eig2;
				if(Math.abs(eigens[0]) <= Math.abs(eigens[1])){
					eig1= eigens[0];
					eig2= eigens[1];
				} else {
					eig1= eigens[1];
					eig2= eigens[0];
				}
				double disparsity = Math.exp(-((eig1*eig1)/(2*beta*beta*eig2*eig2)));
				double relativeBrightness = (1 - Math.exp(-((eig1*eig1+eig2*eig2)/2*brightThresh*brightThresh)));
				if(eig2 < 0) this.vesselness2D[x][y] = 0;
				else this.vesselness2D[x][y] =  disparsity * relativeBrightness;
				}
			}
	}

	private void computeVesselnessWithPhase(){
		Matrix hessian = new Matrix(2,2);
		
		
		WritableRaster raster = source.getRaster();
		this.vesselness2D = new double[source.getWidth()][source.getHeight()];
		
		for(int y=0;y<source.getHeight();y++){
			for(int x=0;x<source.getWidth();x++){
				if(x+1<source.getWidth() && x-1>=0){
					//Horizontal approximation
					double a = (raster.getSample(x+1, y, 0)*phaseCong[x+1][y]+raster.getSample(x-1, y, 0)*phaseCong[x-1][y]-2*raster.getSample(x, y, 0)*phaseCong[x][y]);
					hessian.set(0, 0, a/**phaseCong[x][y]*/); 	
				} else hessian.set(0, 0, 0);
				if(y+1 < source.getHeight() && y-1>=0)
				{
					//Vertical approximation
					double a = (raster.getSample(x, y+1, 0)*phaseCong[x][y+1]+raster.getSample(x, y-1, 0)*phaseCong[x][y-1]-2*raster.getSample(x, y, 0)*phaseCong[x][y]);
						hessian.set(1, 1, a/**phaseCong[x][y]*/);
				} else hessian.set(1, 1, 0);
				if(x+1<source.getWidth() && x-1>=0 && y+1 < source.getHeight() && y-1>=0)
				{
					//Diagonal approximation
					double a = (raster.getSample(x+1, y+1, 0)*phaseCong[x+1][y+1]+raster.getSample(x-1, y-1, 0)*phaseCong[x-1][y-1]-raster.getSample(x+1, y-1, 0)*phaseCong[x+1][y-1]-raster.getSample(x-1, y+1, 0)*phaseCong[x-1][y+1])/4;
						hessian.set(0,1,a/**phaseCong[x][y]*/);
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
				double disparsity = Math.exp(-((eig1*eig1)/(2*beta*beta*eig2*eig2)));
				double relativeBrightness = (1 - Math.exp(-((eig1*eig1+eig2*eig2)/2*brightThresh*brightThresh)));
				if(eig2 < 0) this.vesselness2D[x][y] = 0;
				else this.vesselness2D[x][y] =  disparsity * relativeBrightness/* * phaseCong[y][x]*/;
				}
			}
	}

	
	/** Rewrites data from vesselness2D matrix to image ret.
	 * */
	/*program sa zaobíde bez funkcie,dá sa prepísať aj na koniec computeVesselness2D funkcie
	 * */
public void makeImageWithPhase2D(){
		computeVesselnessWithPhase();
		ret = new IcyBufferedImage(source.getWidth(),source.getHeight(),IcyColorModel.createInstance(1, DataType.DOUBLE));
		ret.beginUpdate();
		for(int y=0;y<source.getHeight();y++){
			for(int x=0;x<source.getWidth();x++){
				ret.setData(x, y, 0, vesselness2D[x][y]);
			}
		}
		ret.endUpdate();
	 }

public void makeImage2D(){
	computeVesselness();
	ret = new IcyBufferedImage(source.getWidth(),source.getHeight(),IcyColorModel.createInstance(1, DataType.DOUBLE));
	ret.beginUpdate();
	for(int y=0;y<source.getHeight();y++){
		for(int x=0;x<source.getWidth();x++){
			ret.setData(x, y, 0, vesselness2D[x][y]);
		}
	}
	ret.endUpdate();
 }
}
