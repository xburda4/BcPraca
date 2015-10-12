package plugins.mBurda.filters;


import java.awt.Graphics;
import java.awt.image.BufferedImage;

import org.jdesktop.swingx.image.GaussianBlurFilter;

import flanagan.complex.ComplexMatrix;

import icy.gui.dialog.MessageDialog;
import icy.image.IcyBufferedImage;
import icy.image.colormodel.IcyColorModel;
import icy.plugin.abstract_.PluginActionable;
import icy.sequence.Sequence;
import icy.type.DataType;

public class Filters extends PluginActionable {
	@Override
	public void run() {
		System.out.println("Depth is " + getActiveSequence().getSizeT() );
		
		//rozmazanie obrazu s okolim 
//		GaussianBlurFilter gauss = new GaussianBlurFilter(5);
//		BufferedImage img = getGrayScale(gauss.filter(getActiveSequence().getImage(0, 0), null));
//
//		Vesselness ves = new Vesselness(blurred,1.5,0.08);
//		ves.makeImage2D();
//		addSequence(new Sequence("Vesselness",ves.ret));
//
//		Neuriteness neu = new Neuriteness(blurred,-5);
//		neu.makeImage2D();
//		addSequence(new Sequence("Neuriteness",neu.ret));
		
		
//		BufferedImage img = getGrayScale(getActiveSequence().getImage(0, 0));
//		
//		
//		double[][] image = new double[img.getRaster().getHeight()][img.getRaster().getWidth()];
//		for(int y = 0;y<img.getRaster().getHeight();y++){
//			for(int x=0;x<img.getRaster().getWidth();x++){
//				image[y][x] = img.getRaster().getSampleDouble(x, y, 0);
//			}
//		}
//		ComplexMatrix matrix = Computations.FourierTransform2D(image,false);
//		double[] scs = {0.0015,0.02,0.05};
//		double[] ors = {0,90,180,270};
//		
//		
//		
//		Vesselness ves = new Vesselness(img,1.5,0.08);
//		ves.phaseCong = Computations.getPhaseCong(matrix, scs, ors, 2, 5, 4);
//		
//		addSequence(new Sequence("Final Phase Cong",makeImage2D(ves.phaseCong)));
//		addSequence(new Sequence("invFTT",makeImage2D(Computations.InverseFourierTransform2D(matrix, true))));
//		
//		addSequence(new Sequence("Magnitude",makeImage2D(Computations.fttToDoubleArr(matrix,true))));
//		addSequence(new Sequence("Phase",makeImage2D(Computations.fttToDoubleArr(matrix, false))));
		
//		ves.makeImage2D();
//		addSequence(new Sequence("Final product - Vesselness",ves.ret));

		
		MessageDialog.showDialog("Filt is done !");
	}

	public IcyBufferedImage makeImage2D(double[][] source){
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
	
	private BufferedImage getGrayScale(BufferedImage original){
		BufferedImage image = new BufferedImage(original.getWidth(), original.getHeight(),
				BufferedImage.TYPE_BYTE_GRAY);  
		Graphics g = image.getGraphics();
		g.drawImage(original, 0, 0, null);
		g.dispose(); 
		return image;
	}
}


