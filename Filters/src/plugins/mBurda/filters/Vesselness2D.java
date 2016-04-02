package plugins.mBurda.filters;


import java.awt.image.BufferedImage;
import icy.image.IcyBufferedImage;
import icy.image.colormodel.IcyColorModel;
import icy.type.DataType;


public class Vesselness2D extends Filter{
	public double beta,brightThresh;
	
	/**@param source blurred grayscale image with one channel
	 * @param beta,bright thresholds controling sensitivity of line measurement
	 * */
	public Vesselness2D(BufferedImage source,double beta,double bright){
		Filter.source = IcyBufferedImage.createFrom(source);
		this.beta = beta;
		this.brightThresh = bright;
	}
	
	public Vesselness2D(BufferedImage source){
		Filter.source = IcyBufferedImage.createFrom(source);
		beta = 1.5;
		brightThresh = 0.008;
	}
	
	/**
	 * Computes vesselness for point in image
	 * @param x x-coordinate
	 * @param y y-coordinate
	 * @return vesselness at point (x,y)
	 */
	private double computeVesselnessForPoint(int x,int y){
		double[] eigens = getEigenValues(x,y);
		
		/*if(eigens[1] < 0) return 0;
		else */return computeVesselness(eigens);
	}
	
	/**
	 * Computes vesselness for point in image with phase congruency
	 * @param x x-coordinate
	 * @param y y-coordinate
	 * @return vesselness at point (x,y)
	 */
	private double computeVesselnessWithPhaseForPoint(int x,int y){
		double[] eigens = getEigenValuesWithPhase(x,y);
		return computeVesselness(eigens);//eigens[1];
		/*
		if(eigens[1] < 0) return 0;
		else return computeVesselness(eigens);
	*/}
	
//	private void computeVesselnessWithPhase(){
//		
//		for(int y=0;y<source.getHeight();y++){
//			for(int x=0;x<source.getWidth();x++){
//				double[] eigens = getEigenValuesWithPhase(x,y);
//				
//				if(eigens[1] < 0) 
//					this.vesselness2D[x][y] = 0;
//				else {
//					this.vesselness2D[x][y] =  computeVesselness(eigens);}
//				}
//			}
//	}

	/**
	 * 
	 * @param eigens array of eigen-values
	 * @return vesselness
	 */
	private double computeVesselness(double[] eigens){
		double disparsity = Math.exp(-((eigens[0]*eigens[0])/(2*beta*beta*eigens[1]*eigens[1])));
		double relativeBrightness = (1 - Math.exp(-((eigens[0]*eigens[0]+eigens[1]*eigens[1])/2*brightThresh*brightThresh)));
		
		return (disparsity * relativeBrightness);
	}
	
	/** Rewrites data from vesselness2D matrix to image ret.
	 * */
	/*program sa zaobíde bez funkcie,dá sa prepísať aj na koniec computeVesselness2D funkcie
	 * */
public IcyBufferedImage makeImageWithPhase2D(){
		IcyBufferedImage ret = new IcyBufferedImage(source.getWidth(),source.getHeight(),IcyColorModel.createInstance(1, DataType.DOUBLE));
		ret.beginUpdate();
		for(int y=0;y<source.getHeight();y++){
			for(int x=0;x<source.getWidth();x++){
				ret.setData(x, y, 0, computeVesselnessWithPhaseForPoint(x, y));
			}
		}
		ret.endUpdate();
		return ret;
	 }

public IcyBufferedImage makeImage2D(){
	IcyBufferedImage ret = new IcyBufferedImage(source.getWidth(),source.getHeight(),IcyColorModel.createInstance(1, DataType.DOUBLE));
	ret.beginUpdate();
	
	for(int y=0;y<source.getHeight();y++){
		for(int x=0;x<source.getWidth();x++){
			
			ret.setData(x, y, 0, computeVesselnessForPoint(x, y));
		}
	}
	ret.endUpdate();
	return ret;
 }
}
