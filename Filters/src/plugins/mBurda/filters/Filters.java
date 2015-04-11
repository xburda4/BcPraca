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
		//rozmazanie obrazu s okolim 5
//		GaussianBlurFilter gauss = new GaussianBlurFilter(3);
//		BufferedImage blurred = getGrayScale(gauss.filter(getActiveSequence().getImage(0, 0), null));
//
//		Vesselness ves = new Vesselness(blurred,1.5,0.08);
//		ves.makeImage2D();
//		addSequence(new Sequence("Vesselness",ves.ret));
//
//		Neuriteness neu = new Neuriteness(blurred,-5);
//		neu.makeImage2D();
//		addSequence(new Sequence("Neuriteness",neu.ret));
		
		
		BufferedImage img = getGrayScale(getActiveSequence().getImage(0, 0));
		double[][] image = new double[img.getRaster().getHeight()][img.getRaster().getWidth()];
		for(int y = 0;y<img.getRaster().getHeight();y++){
			for(int x=0;x<img.getRaster().getWidth();x++){
				image[y][x] = img.getRaster().getSampleDouble(x, y, 0);
			}
		}
		ComplexMatrix matrix= Computations.FourierTransform2D(image);
		image = Computations.InverseFourierTransform2D(matrix);
		addSequence(new Sequence("FFT",makeImage2D(image)));
		
		MessageDialog.showDialog("Filt is done !");
	}

	public IcyBufferedImage makeImage2D(double[][] source){
		IcyBufferedImage ret = new IcyBufferedImage(source[0].length,source.length,IcyColorModel.createInstance(1, DataType.DOUBLE));
		for(int y=0;y<source.length;y++){
			for(int x=0;x<source[y].length;x++){
				ret.setData(x, y, 0, source[y][x]);
			}
		}
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


