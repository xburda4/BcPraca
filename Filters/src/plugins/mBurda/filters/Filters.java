package plugins.mBurda.filters;


import java.awt.Graphics;
import java.awt.image.BufferedImage;
import org.jdesktop.swingx.image.GaussianBlurFilter;
import icy.plugin.abstract_.PluginActionable;
import icy.sequence.Sequence;

public class Filters extends PluginActionable {
	@Override
	public void run() {
					//rozmazanie obrazu s okolim 9
					GaussianBlurFilter gauss = new GaussianBlurFilter(9);
					BufferedImage blurred = getGrayScale(gauss.filter(getActiveSequence().getImage(0, 0), null));
					
					Vesselness ves = new Vesselness(blurred,5,0.85);
					ves.makeImage2D();
					addSequence(new Sequence(ves.ret));

					Neuriteness neu = new Neuriteness(blurred,6);
					neu.makeImage2D();
					addSequence(new Sequence(neu.ret));

		//MessageDialog.showDialog("Filt is done !");
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
