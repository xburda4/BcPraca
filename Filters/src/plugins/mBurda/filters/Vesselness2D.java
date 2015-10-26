package plugins.mBurda.filters;


import java.awt.image.BufferedImage;
import icy.image.IcyBufferedImage;
import icy.image.colormodel.IcyColorModel;
import icy.type.DataType;


public class Vesselness2D extends Filter{
	//public BufferedImage source;
	public IcyBufferedImage ret;
	public double[][] vesselness2D;
	public double beta,brightThresh;
	//public double[][] phaseCong; //prepísať na static + triedu z ktorej budú vesselness a neuriteness dediť
	
	/**@param source blurred grayscale image with one channel
	 * @param beta,bright thresholds controling sensitivity of line measurement
	 * */
	public Vesselness2D(BufferedImage source,double beta,double bright){
		Filter.source = IcyBufferedImage.createFrom(source);
		this.beta = beta;
		this.brightThresh = bright;
		vesselness2D = new double[source.getWidth()][source.getHeight()];
	}
	public Vesselness2D(BufferedImage source){
		Filter.source = IcyBufferedImage.createFrom(source);
		beta = 1.5;
		brightThresh = 0.08;
		vesselness2D = new double[source.getWidth()][source.getHeight()];
	}
	
	/**Computes vesselness for 2D grayscale images. Saves it to a matrix vesselness2D
	 * */
	private void computeVesselness2D(){
		
		for(int y=0;y<source.getHeight();y++){
			for(int x=0;x<source.getWidth();x++){
				double[] eigens = getEigenValues(x,y);
				
				if(eigens[1] < 0) this.vesselness2D[x][y] = 0;
				else this.vesselness2D[x][y] =  computeVesselness(eigens);
				}
			}
	}

	private void computeVesselnessWithPhase(){
		
		for(int y=0;y<source.getHeight();y++){
			for(int x=0;x<source.getWidth();x++){
				double[] eigens = getEigenValuesWithPhase(x,y);
				
				if(eigens[1] < 0) 
					this.vesselness2D[x][y] = 0;
				else {
					this.vesselness2D[x][y] =  computeVesselness(eigens);}
				}
			}
	}

	private double computeVesselness(double[] eigens){
		double disparsity = Math.exp(-((eigens[0]*eigens[0])/(2*beta*beta*eigens[1]*eigens[1])));
		double relativeBrightness = (1 - Math.exp(-((eigens[0]*eigens[0]+eigens[1]*eigens[2])/2*brightThresh*brightThresh)));
		
		return disparsity * relativeBrightness;
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
	computeVesselness2D();
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
