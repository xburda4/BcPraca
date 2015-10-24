package plugins.mBurda.filters;


import java.awt.image.WritableRaster;

import Jama.Matrix;
import icy.image.IcyBufferedImage;
import icy.image.colormodel.IcyColorModel;
import icy.sequence.Sequence;
import icy.sequence.SequenceUtil;
import icy.type.DataType;


public class Vesselness3D{
	public Sequence source;
	public Sequence ret;
	public double[][][] vesselness3D;
	public double beta,brightThresh,gamma;
	public double[][][] phaseCong;
	
	/**@param source blurred grayscale image with one channel
	 * @param beta,bright,gamma thresholds controling sensitivity of line measurement
	 * */
	public Vesselness3D(Sequence source,double beta,double bright,double gamma){
		this.source = SequenceUtil.getCopy(source);
		this.beta = beta;
		this.brightThresh = bright;
		this.gamma = gamma;
	}
	
	public Vesselness3D(Sequence source){
		this.source = SequenceUtil.getCopy(source);
	}
	
	
	/**Computes vesselness for 3D grayscale images. Saves it to a matrix vesselness2D
	 * */
	private void computeVesselness(){
		Matrix hessian = new Matrix(3,3);
		
		
		WritableRaster raster;
		this.vesselness3D = new double[source.getSizeT()][source.getWidth()][source.getHeight()];
		WritableRaster rasMinus=null,rasPlus=null;
		for(int z=0;z<source.getSizeT();z++){
			if(z+1<source.getSizeT() && z-1>=0){
				rasMinus = source.getImage(z-1, 0).getRaster();
				rasPlus = source.getImage(z+1, 0).getRaster();
			}
			raster = source.getImage(z, 0).getRaster();
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
					
					if(z+1<source.getSizeT() && z-1>=0 && y+1 < source.getHeight() && y-1>=0){
						double a = (rasPlus.getSample(x, y+1, 0)+rasMinus.getSample(x, y-1, 0)-rasPlus.getSample(x, y-1, 0)-rasMinus.getSample(x, y+1, 0))/4;
						hessian.set(1,2,a);
						
					}else hessian.set(1, 2, 0);
					hessian.set(2, 1, hessian.get(1,2));
					
					if(z+1<source.getSizeT() && z-1>=0){
						
						double a = (rasPlus.getSample(x, y, 0)+rasMinus.getSample(x, y, 0)-2*raster.getSample(x, y, 0));
						hessian.set(2, 2, a); 	
					} else hessian.set(2, 2, 0);
					
					if(z+1<source.getSizeT() && z-1>=0 && x+1 < source.getWidth() && x-1>=0){
						double a = (rasPlus.getSample(x+1, y, 0)+rasMinus.getSample(x-1, y, 0)-rasPlus.getSample(x-1, y, 0)-rasMinus.getSample(x+1, y, 0))/4;
						hessian.set(0,2,a);
						
					}else hessian.set(0, 2, 0);
					hessian.set(2, 0, hessian.get(0,2));
					
					
					
				double[] eigens = hessian.eig().getRealEigenvalues();
				double eig1,eig2,eig3;
				
				eig3 = Math.max(Math.max(Math.abs(eigens[2]),Math.abs(eigens[1])), Math.abs(eigens[0]));
				eig1 = Math.min(Math.min(Math.abs(eigens[2]), Math.abs(eigens[1])), Math.abs(eigens[0]));
				eig2 = eigens[2]+eigens[1]+eigens[0]-eig1-eig3;
				double a = eig2/eig3,b = eig1/Math.sqrt(eig2*eig3);
				double disparsity = Math.exp(-((b*b)/(2*beta*beta)));
				double relativeBrightness = (1 - Math.exp(-(Math.sqrt(eig1*eig1+eig2*eig2+eig3*eig3)/2*brightThresh*brightThresh)));
				//change name
				double third = (1-Math.exp((a*a)/(2*gamma*gamma)));
				if(eig2 < 0 && eig3 < 0) this.vesselness3D[z][x][y] = 0;
				else this.vesselness3D[z][x][y] =  disparsity * relativeBrightness * third;
				}
			}
		}
	}

	
	/** Rewrites data from vesselness2D matrix to image ret.
	 * */
	/*program sa zaobíde bez funkcie,dá sa prepísať aj na koniec computeVesselness2D funkcie
	 * */
public void makeImage3D(){
		computeVesselness();
		ret = new Sequence("Vesselness3D");
		IcyBufferedImage img = null;
		WritableRaster ras = null;
		ret.beginUpdate();
		for(int z=0;z<source.getSizeT();z++){
			img = new IcyBufferedImage(source.getWidth(),source.getHeight(),IcyColorModel.createInstance(1, DataType.DOUBLE));
			img.beginUpdate();
			ras = img.getRaster();
				for(int y=0;y<source.getHeight();y++){
					for(int x=0;x<source.getWidth();x++){
				ras.setSample(x, y, 0, vesselness3D[z][x][y]);
					}
			}
			ret.addImage(img);
			
			
			}
		img.endUpdate();
		ret.endUpdate();
		}
}
