package plugins.mBurda.filters;

import icy.image.IcyBufferedImage;
import icy.image.colormodel.IcyColorModel;
import icy.plugin.abstract_.PluginActionable;
import icy.sequence.Sequence;
import icy.type.DataType;

public class Filters extends PluginActionable {
	@Override
	public void run() {
		
		
		Sequence seq = getActiveSequence();
		
		Vesselness3D ves = new Vesselness3D(seq,0.2,0,3);
		ves.makeImage3D();
		addSequence(ves.ret);
		
		Neuriteness3D neu = new Neuriteness3D(seq,50,0.33);
		neu.makeImage3D();
		addSequence(neu.ret);
		
		//rozmazanie obrazu s okolim 3
//		GaussianBlurFilter gauss = new GaussianBlurFilter(3);
//		BufferedImage image = getGrayScale(gauss.filter(getActiveSequence().getImage(0, 0), null));
//		
		
//		BufferedImage img = getGrayScale(getActiveSequence().getImage(0, 0));
//		
//		
//		double[][] image = new double[img.getRaster().getHeight()][img.getRaster().getWidth()];
//		for(int y = 0;y<img.getRaster().getHeight();y++){
//			for(int x=0;x<img.getRaster().getWidth();x++){
//				image[y][x] = img.getRaster().getSampleDouble(x, y, 0);
//			}
//		}
		
		
		/**
		 * 2D results on eye image
		 * */
//	{
//		ComplexMatrix matrix = Computations.FourierTransform2D(image,false);
//		
//		double[] scs = {0.0015,0.002,0.005,0.009,0.013};
//		double[] ors = {0,90,180,270};
//		
//		Vesselness2D ves = new Vesselness2D(image,5.5,0.007);
//		
//		Neuriteness2D neu = new Neuriteness2D(image,-5);
//		
//		ves.phaseCong = Computations.getPhaseCong(matrix, scs, ors, -2, -55 ,8);
//		neu.phaseCong = Computations.getPhaseCong(matrix, scs, ors, -2, -55 ,8);
//		ves.makeImage2D();
//		addSequence(new Sequence("Vesselness",ves.ret));
//		ves.makeImageWithPhase2D();
//		addSequence(new Sequence("Vesselness with phase",ves.ret));
//		neu.makeImage2D();
//		addSequence(new Sequence("Neuriteness",neu.ret));
//		neu.makeImageWithPhase2D();
//		addSequence(new Sequence("Neuriteness with phase",neu.ret));
//	}
	}

	public static IcyBufferedImage makeImage2D(double[][] source){
		long start = System.nanoTime();
		IcyBufferedImage ret = new IcyBufferedImage(source[0].length,source.length,IcyColorModel.createInstance(1, DataType.DOUBLE));
		
		//Without this,it would take forever
		ret.beginUpdate();
		
		for(int y=0;y< source.length;y++){
			for(int x=0;x<source[y].length;x++){
				ret.setData(x, y, 0, source[y][x]);
			}
		}
		ret.endUpdate();
		
		System.out.println("Time of making an image is " + (System.nanoTime() - start));
		
		return ret;
	 }
}


