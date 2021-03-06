package plugins.mBurda.filters;


import java.awt.image.WritableRaster;

import Jama.Matrix;
import icy.image.IcyBufferedImage;
import icy.image.colormodel.IcyColorModel;
import icy.sequence.Sequence;
import icy.sequence.SequenceUtil;
import icy.type.DataType;


public class Neuriteness3D{
	private Sequence source;
	public Sequence ret;
	public double[][][] neuriteness3D;
	public double alpha,gamma;
	public double[][][] phaseCong;
	
	/**@param source blurred grayscale image with one channel
	 * @param beta,bright,gamma thresholds controling sensitivity of line measurement
	 * */
	public Neuriteness3D(Sequence source,double alpha,double gamma){
		this.source = SequenceUtil.getCopy(source);
		this.alpha = alpha;
		this.gamma = gamma;
	}
	
	public Neuriteness3D(Sequence source){
		this.source = SequenceUtil.getCopy(source);
		gamma = 0.33;
		alpha = 0;
	}
	
	
	/**Computes neuriteness for 3D grayscale images. Saves it to a matrix vesselness3D
	 * */
	private void computeNeuriteness(){
		Matrix hessian = new Matrix(3,3);
		
		
		WritableRaster raster;
		WritableRaster rasMinus=null,rasPlus=null;
		this.neuriteness3D = new double[source.getSizeT()][source.getWidth()][source.getHeight()];
		
		double minEig=0;
		for(int z=0;z<source.getSizeT();z++){
			raster = source.getImage(z, 0).getRaster();
			if(z+1<source.getSizeT() && z-1>=0){
				rasMinus = source.getImage(z-1, 0).getRaster();
				rasPlus = source.getImage(z+1, 0).getRaster();
			}
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
					
				double temp0 = hessian.get(0, 0),temp1 = hessian.get(1, 1),temp2 = hessian.get(2, 2);
					
				hessian.set(0, 1, (1-gamma)*hessian.get(0, 1));
				hessian.set(0, 2, (1-gamma)*hessian.get(0, 2));
				hessian.set(1, 2, (1-gamma)*hessian.get(1, 2));
				hessian.set(0, 0, gamma/2*temp1+gamma/2*temp2+temp0);
				hessian.set(1, 1, gamma/2*temp0+gamma/2*temp2+temp1);
				hessian.set(2, 2, gamma/2*temp1+gamma/2*temp0+temp2);
				hessian.set(1, 0, hessian.get(0, 1));
				hessian.set(2, 0, hessian.get(0, 2));
				hessian.set(2, 1, hessian.get(1, 2));
				
				double[] eigens = hessian.eig().getRealEigenvalues();
				if(z==0 && x == 0 && y == 0) minEig = Math.min(Math.min(Math.abs(eigens[1]),Math.abs(eigens[0])),Math.abs(eigens[2]));
				else minEig = Math.min(minEig,Math.min(Math.min(Math.abs(eigens[1]),Math.abs(eigens[0])),Math.abs(eigens[2])));
				
				neuriteness3D[z][x][y] = Math.max(Math.max(Math.abs(eigens[1]),Math.abs(eigens[0])),Math.abs(eigens[2]));
				}
			}
		}
		for(int z=0;z<source.getSizeT();z++){
			for(int y=0;y<source.getHeight();y++){
				for(int x=0;x<source.getWidth();x++){
					neuriteness3D[z][x][y] /= minEig;
					}
				}
			}
	}

	
	/** Rewrites data from neuriteness3D matrix to image ret.
	 * */
	/*program sa zaobíde bez funkcie,dá sa prepísať aj na koniec predošlej funkcie
	 * */
public void makeImage3D(){
		computeNeuriteness();
		ret = new Sequence("Neuriteness3D");
		IcyBufferedImage img = null;
		WritableRaster ras = null;
		ret.beginUpdate();
		for(int z=0;z<source.getSizeT();z++){
			img = new IcyBufferedImage(source.getWidth(),source.getHeight(),IcyColorModel.createInstance(1, DataType.DOUBLE));
			img.beginUpdate();
			ras = img.getRaster();
				for(int y=0;y<source.getHeight();y++){
					for(int x=0;x<source.getWidth();x++){
				ras.setSample(x, y, 0, neuriteness3D[z][x][y]);
					}
			}
			ret.addImage(img);
			
			
			}
		img.endUpdate();
		ret.endUpdate();
	}
}
