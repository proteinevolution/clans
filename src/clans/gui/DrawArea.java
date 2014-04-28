package clans.gui;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.JPanel;

import clans.UnusedWeirdClassPreviouslyCalledSelectclass;
import clans.model.ClusterData;
import clans.model.SequenceGroup;
import clans.model.proteins.MinimalAttractionValue;

public class DrawArea extends JPanel implements java.awt.print.Printable, ComponentListener {

	private static final long serialVersionUID = -9084230540557021638L;
	
	private ProgramWindow parent;
	private ClusterData data;
	
	/** 
	 * The current zoom.
	 */
	private float zoomFactor = 1f;
	/**
	 * Pad the graph with this many pixels on all four sides when in 100% zoom. 
	 */
	final private float PADDING_PIXEL = 20;
	/**
	 * Pre-zoom is used to achieve an approximation of {@code PADDING_PIXEL} pixels as padding irrespective of the
	 * window size.
	 */
	private float preZoomFactor;

	int drawWidth;
	int drawHeight;
	
	/**
	 * Movement of the graph inside the DrawArea. Used to move the graph with or without active zoom.
	 */
	int xTranslate = 0;
	int yTranslate = 0;
	
	Color mainColor = Color.black;
	Color innerAreaBackgroundColor = Color.white;
	
	private final Color selectionRectangleColor = Color.orange;
	
	Color selectedSequenceCircleColor = Color.red;
	private final Color blastQueryCircleColor = Color.magenta;
	Color blastHitCircleColor = Color.green;
	Color blastHitNumberColor = new Color(0.4f, 0.4f, 0.4f);
	
	protected java.awt.Font myfont = getFont();
	private java.awt.Font smallfont = new java.awt.Font("Monospace", 0, 9);
	private int fontsize = myfont.getSize();
	private int fontwidth = (int) (fontsize / 2);
	
	// mostly unsorted fields follow
	private HashMap<String, MinimalAttractionValue> frustration;// color code for the frustration of connections
	
	private int[] originpos = new int[2];
	private int colornum = 10;

	int stereoangle = 4;
	private double[][] box;
	private double[][] boxdata;
	double xscale = 1;
	double yscale = 1;
	boolean rotatez = false;
	
	private double minX = 1;
	private double minY = 1;
	private double maxX = 1;
	private double maxY = 1;
	
	private double xMove = 0;
	private double yMove = 0;
	private double zMove = 0;

	private int[] draw = new int[4];// the coordinates of the field selected by the mouse
	
	private String refdbstring = "";

	int[] xposarr;
	int[] yposarr;
	int posnum;

	/**
	 * the current mouse rotation
	 */
	double[][] tmprotmtx = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };	
	
	/**
	 * Creates a new DrawArea.
	 * 
	 * @param parent
	 *            The parent window.
	 * @param data
	 *            The cluster data.
	 */
	public DrawArea(ProgramWindow parent, ClusterData data) {
		this.parent = parent;
		this.data = data;
				
		this.setOpaque(true);// should be set by default, but can't hurt to re-set
		this.setDoubleBuffered(true);// should be set by default, but can't hurt to re-set
		
		this.setBackground(innerAreaBackgroundColor);
		
		this.addComponentListener(this);
	}
	

	void init() {
		data.posarr = new float[0][0];
		
		data.colorarr = makecolors(colornum);
		
		data.colorcutoffs = new float[colornum];
		if (data.usescval) {
			for (int i = 0; i < colornum; i++) {
				data.colorcutoffs[i] = (float) (((float) i / (float) (colornum)) * data.p2attfactor);
			}
		} else {
			for (int i = 0; i < colornum; i++) {
				data.colorcutoffs[i] = ((float) i / (float) (colornum));
			}
		}
		
		refdbstring = "";
		for (int i = 0; i < data.referencedb.length; i++) {
			refdbstring += data.referencedb[i] + "; ";
		}
		
		// set the drawbox coordinates (a box of size 1 in x,y,z space)
		box = new double[8][3];
		box[0][0] = 0;
		box[0][1] = 0;
		box[0][2] = 0;

		box[1][0] = 1;
		box[1][1] = 0;
		box[1][2] = 0;

		box[2][0] = 1;
		box[2][1] = 1;
		box[2][2] = 0;

		box[3][0] = 0;
		box[3][1] = 1;
		box[3][2] = 0;

		box[4][0] = 0;
		box[4][1] = 0;
		box[4][2] = 1;

		box[5][0] = 1;
		box[5][1] = 0;
		box[5][2] = 1;

		box[6][0] = 1;
		box[6][1] = 1;
		box[6][2] = 1;

		box[7][0] = 0;
		box[7][1] = 1;
		box[7][2] = 1;

		boxdata = new double[8][3];
		originpos[0] = 0;
		originpos[1] = 0;
	}

	/**
	 * Creates a grayscale color gradient.
	 * 
	 * @param count
	 *            The number of colors to create.
	 * @return The colors.
	 */
	Color[] makecolors(int count) {
		Color[] retarr = new Color[count];
		
		int worstred = 230;
		int worstgreen = 230;
		int worstblue = 230;
		int bestred = 0;
		int bestgreen = 0;
		int bestblue = 0;
		
		float redstep = ((float) (bestred - worstred)) / ((float) count);
		float greenstep = ((float) (bestgreen - worstgreen)) / ((float) count);
		float bluestep = ((float) (bestblue - worstblue)) / ((float) count);

		// now compute the stepwise gradient
		for (int i = 0; i < count; i++) {
			retarr[i] = new Color(
					(int) (worstred + (i * redstep)),
					(int) (worstgreen + (i * greenstep)),
					(int) (worstblue + (i * bluestep)));
		}

		return retarr;
	}

	// ----------------------------------------------------------------------

	public int print(Graphics g, java.awt.print.PageFormat pf, int pi)
			throws java.awt.print.PrinterException {
		double xscalel = pf.getImageableWidth() / ((double) getWidth());
		double yscalel = pf.getImageableHeight() / ((double) getHeight());
		if (pi >= 1) {
			return java.awt.print.Printable.NO_SUCH_PAGE;
		}
		Graphics2D g2d = (Graphics2D) g;
		g2d.translate(pf.getImageableX(), pf.getImageableY());
		g2d.scale(xscalel, yscalel);
		setDoubleBuffered(false);
		paintComponent(g2d);
		setDoubleBuffered(true);
		return java.awt.print.Printable.PAGE_EXISTS;
	}

	/**
	 * 
	 */
	public void paintComponent(Graphics g) {

		// the overlay has to be done first as otherwise the dataLock might prevent it from showing
		if (parent.messageOverlayActive) {
			parent.messageOverlay.activate_overlay();
		}

		g.setFont(myfont);
		
		g.setColor(innerAreaBackgroundColor);
		g.fillRect(0, 0, getWidth(), getHeight());
		
		drawImage(g, true);
		
		if (parent.stereocheckboxmenuitem.isSelected()) {
			// now rotate the graph by the stereoangle along the x-axis
			parent.mousemove[0] += ((float) stereoangle / 360f) * getWidth() / 2;
			
			drawImage(g, false);
			
			// now rotate back to what it was before
			parent.mousemove[0] -= ((float) stereoangle / 360f) * getWidth() / 2;
		}
		
		parent.recalc = true;

		if (parent.showInfoOnDrawAreaCheckbox.isSelected()) {
			drawInfoMessages(g);
		}
		
		if (parent.showorigcheckbox.isSelected()) {
			g.setColor(Color.red);
			g.drawString("X(0,0,0)", originpos[0], originpos[1]);
			
			drawCoordinateSystemn(g);
		}
	}
	
	private void drawImage(Graphics g, boolean is_first_stereo_image) {
		if (!parent.contains_data()) {
			return;
		}
		
		drawWidth = getWidth();
		if (parent.stereocheckboxmenuitem.isSelected()) {
			drawWidth /= 2;
		}
		
		drawHeight = getHeight();
		
		if (is_first_stereo_image) {
			g.translate(xTranslate, yTranslate);	
		} else {
			g.translate(xTranslate + drawWidth, yTranslate);
		}
		
		drawWidth *= zoomFactor;
		drawHeight *= zoomFactor;
		
		g.setColor(mainColor);

		getangles();
		
		drawdata(g);
	}

	private void drawInfoMessages(Graphics g) {
		g.setColor(mainColor);
		g.drawString(" INFO: blast=" + data.blastpath + " refdb=" + refdbstring, (colornum + 9) * fontwidth
				- xTranslate, fontsize - yTranslate);
		
		g.setFont(smallfont);
		g.drawString("x:" + minX + " to " + maxX, fontwidth - xTranslate, getHeight() - fontsize - 3 - yTranslate);
		g.drawString("y:" + minY + " to " + maxY, fontwidth - xTranslate, getHeight() - yTranslate - 3);
		
		g.drawString((float) ((int) (data.myrotmtx[0][0] * 100)) / 100 + ","
				+ (float) ((int) (data.myrotmtx[0][1] * 100)) / 100 + ","
				+ (float) ((int) (data.myrotmtx[0][2] * 100))/ 100,
				xTranslate, fontsize - yTranslate);
		g.drawString((float) ((int) (data.myrotmtx[1][0] * 100)) / 100 + ","
				+ (float) ((int) (data.myrotmtx[1][1] * 100)) / 100 + ","
				+ (float) ((int) (data.myrotmtx[1][2] * 100)) / 100,
				xTranslate, 2 * fontsize - yTranslate);
		g.drawString((float) ((int) (data.myrotmtx[2][0] * 100)) / 100 + ","
				+ (float) ((int) (data.myrotmtx[2][1] * 100)) / 100 + ","
				+ (float) ((int) (data.myrotmtx[2][2] * 100)) / 100,
				xTranslate, 3 * fontsize - yTranslate);
	}
	
	void getangles() {
		// set alphaxy alphaxz and alphayz and their sin and cos (the rotation angles)
		// first compute the angle xy by taking the last rotaion angles and the last mouse movement

		// add minimal value to avoid division by 0 errors
		xMove = java.lang.Math.PI / 2 * ((double) parent.mousemove[0] / (drawWidth / 2)) + Double.MIN_VALUE; 
		yMove = java.lang.Math.PI / 2 * ((double) parent.mousemove[1] / (drawWidth / 2)) + Double.MIN_VALUE;
		
		// if rotatez the calculate the rotation angle for z axis
		if (rotatez || (parent.cluster2dbutton.isSelected())) {
			zMove = yMove + xMove;
			tmprotmtx[0][0] = java.lang.Math.cos(zMove);
			tmprotmtx[0][1] = -java.lang.Math.sin(zMove);
			tmprotmtx[0][2] = 0;
			tmprotmtx[1][0] = java.lang.Math.sin(zMove);
			tmprotmtx[1][1] = java.lang.Math.cos(zMove);
			tmprotmtx[1][2] = 0;
			tmprotmtx[2][0] = 0;
			tmprotmtx[2][1] = 0;
			tmprotmtx[2][2] = 1;
		} else {// get the parameters for the xy rotation
			if (xMove == 0 && yMove == 0) {
				return;
			}
			double rotz = java.lang.Math.atan(yMove / xMove) - java.lang.Math.PI / 2;// the amount I have to rotate
																						// in xy
			double rotx = java.lang.Math.sqrt(xMove * xMove + yMove * yMove);// the amount to rotate by
			if ((xMove < 0)) {
				rotx = -rotx;
			}
			double cosrotz = java.lang.Math.cos(rotz);
			double sinrotz = java.lang.Math.sin(rotz);
			double cosrotx = java.lang.Math.cos(rotx);
			double sinrotx = java.lang.Math.sin(rotx);
			tmprotmtx[0][0] = ((cosrotz * cosrotz) + (-sinrotz * cosrotx * -sinrotz));
			tmprotmtx[0][1] = ((cosrotz * sinrotz) + (-sinrotz * cosrotx * cosrotz));
			tmprotmtx[0][2] = (-sinrotz * -sinrotx);
			tmprotmtx[1][0] = ((sinrotz * cosrotz) + (cosrotz * cosrotx * -sinrotz));
			tmprotmtx[1][1] = ((sinrotz * sinrotz) + (cosrotz * cosrotx * cosrotz));
			tmprotmtx[1][2] = (cosrotz * -sinrotx);
			tmprotmtx[2][0] = (sinrotx * -sinrotz);
			tmprotmtx[2][1] = (sinrotx * cosrotz);
			tmprotmtx[2][2] = cosrotx;
		}
		xMove = 0;
		yMove = 0;
		zMove = 0;
		// now multiply tmprotmtx to rotmtx and save in myrotmtx
		data.myrotmtx[0][0] = (tmprotmtx[0][0] * data.rotmtx[0][0]) + (tmprotmtx[0][1] * data.rotmtx[1][0])
				+ (tmprotmtx[0][2] * data.rotmtx[2][0]);
		data.myrotmtx[0][1] = (tmprotmtx[0][0] * data.rotmtx[0][1]) + (tmprotmtx[0][1] * data.rotmtx[1][1])
				+ (tmprotmtx[0][2] * data.rotmtx[2][1]);
		data.myrotmtx[0][2] = (tmprotmtx[0][0] * data.rotmtx[0][2]) + (tmprotmtx[0][1] * data.rotmtx[1][2])
				+ (tmprotmtx[0][2] * data.rotmtx[2][2]);
		data.myrotmtx[1][0] = (tmprotmtx[1][0] * data.rotmtx[0][0]) + (tmprotmtx[1][1] * data.rotmtx[1][0])
				+ (tmprotmtx[1][2] * data.rotmtx[2][0]);
		data.myrotmtx[1][1] = (tmprotmtx[1][0] * data.rotmtx[0][1]) + (tmprotmtx[1][1] * data.rotmtx[1][1])
				+ (tmprotmtx[1][2] * data.rotmtx[2][1]);
		data.myrotmtx[1][2] = (tmprotmtx[1][0] * data.rotmtx[0][2]) + (tmprotmtx[1][1] * data.rotmtx[1][2])
				+ (tmprotmtx[1][2] * data.rotmtx[2][2]);
		data.myrotmtx[2][0] = (tmprotmtx[2][0] * data.rotmtx[0][0]) + (tmprotmtx[2][1] * data.rotmtx[1][0])
				+ (tmprotmtx[2][2] * data.rotmtx[2][0]);
		data.myrotmtx[2][1] = (tmprotmtx[2][0] * data.rotmtx[0][1]) + (tmprotmtx[2][1] * data.rotmtx[1][1])
				+ (tmprotmtx[2][2] * data.rotmtx[2][1]);
		data.myrotmtx[2][2] = (tmprotmtx[2][0] * data.rotmtx[0][2]) + (tmprotmtx[2][1] * data.rotmtx[1][2])
				+ (tmprotmtx[2][2] * data.rotmtx[2][2]);
	}

	/**
	 * Show the coordinate system orientation	
	 * @param g
	 */
	void drawCoordinateSystemn(Graphics g) {
		// do the xz rotation
		double minxl = data.myrotmtx[0][0] * box[0][0] + data.myrotmtx[0][1] * box[0][1] + data.myrotmtx[0][2]
				* box[0][2];
		double maxxl = minxl;
		double minyl = data.myrotmtx[1][0] * box[0][0] + data.myrotmtx[1][1] * box[0][1] + data.myrotmtx[1][2]
				* box[0][2];
		double maxyl = minyl;
		boxdata[3][0] = data.myrotmtx[0][0] * box[3][0] + data.myrotmtx[0][1] * box[3][1] + data.myrotmtx[0][2]
				* box[3][2];
		boxdata[3][1] = data.myrotmtx[1][0] * box[3][0] + data.myrotmtx[1][1] * box[3][1] + data.myrotmtx[1][2]
				* box[3][2];
		boxdata[3][2] = data.myrotmtx[2][0] * box[3][0] + data.myrotmtx[2][1] * box[3][1] + data.myrotmtx[2][2]
				* box[3][2];
		if (boxdata[3][0] < minxl) {
			minxl = boxdata[3][0];
		} else if (boxdata[3][0] > maxxl) {
			maxxl = boxdata[3][0];
		}
		if (boxdata[3][1] < minyl) {
			minyl = boxdata[3][1];
		} else if (boxdata[3][1] > maxyl) {
			maxyl = boxdata[3][1];
		}
		boxdata[2][0] = data.myrotmtx[0][0] * box[2][0] + data.myrotmtx[0][1] * box[2][1] + data.myrotmtx[0][2]
				* box[2][2];
		boxdata[2][1] = data.myrotmtx[1][0] * box[2][0] + data.myrotmtx[1][1] * box[2][1] + data.myrotmtx[1][2]
				* box[2][2];
		boxdata[2][2] = data.myrotmtx[2][0] * box[2][0] + data.myrotmtx[2][1] * box[2][1] + data.myrotmtx[2][2]
				* box[2][2];
		if (boxdata[2][0] < minxl) {
			minxl = boxdata[2][0];
		} else if (boxdata[2][0] > maxxl) {
			maxxl = boxdata[2][0];
		}
		if (boxdata[2][1] < minyl) {
			minyl = boxdata[2][1];
		} else if (boxdata[2][1] > maxyl) {
			maxyl = boxdata[2][1];
		}
		boxdata[7][0] = data.myrotmtx[0][0] * box[7][0] + data.myrotmtx[0][1] * box[7][1] + data.myrotmtx[0][2]
				* box[7][2];
		boxdata[7][1] = data.myrotmtx[1][0] * box[7][0] + data.myrotmtx[1][1] * box[7][1] + data.myrotmtx[1][2]
				* box[7][2];
		boxdata[7][2] = data.myrotmtx[2][0] * box[7][0] + data.myrotmtx[2][1] * box[7][1] + data.myrotmtx[2][2]
				* box[7][2];
		if (boxdata[7][0] < minxl) {
			minxl = boxdata[7][0];
		} else if (boxdata[7][0] > maxxl) {
			maxxl = boxdata[7][0];
		}
		if (boxdata[7][1] < minyl) {
			minyl = boxdata[7][1];
		} else if (boxdata[7][1] > maxyl) {
			maxyl = boxdata[7][1];
		}
		double xoffset = -minxl;
		double yoffset = -minyl;
		double xfac = 50;
		double yfac = 50;

		g.setColor(mainColor);
		// now I have the coordinates; now draw the box
		g.drawLine((int) ((boxdata[0][0] + xoffset) * xfac), (int) ((boxdata[0][1] + yoffset) * yfac),
				(int) ((boxdata[3][0] + xoffset) * xfac), (int) ((boxdata[3][1] + yoffset) * yfac));
		g.drawLine((int) ((boxdata[2][0] + xoffset) * xfac), (int) ((boxdata[2][1] + yoffset) * yfac),
				(int) ((boxdata[3][0] + xoffset) * xfac), (int) ((boxdata[3][1] + yoffset) * yfac));
		g.drawLine((int) ((boxdata[3][0] + xoffset) * xfac), (int) ((boxdata[3][1] + yoffset) * yfac),
				(int) ((boxdata[7][0] + xoffset) * xfac), (int) ((boxdata[7][1] + yoffset) * yfac));
		
		g.setColor(Color.red);
		g.drawString("X", (int) ((boxdata[2][0] + xoffset) * xfac), (int) ((boxdata[2][1] + yoffset) * yfac));
		g.drawString("Y", (int) ((boxdata[0][0] + xoffset) * xfac), (int) ((boxdata[0][1] + yoffset) * yfac));
		g.drawString("Z", (int) ((boxdata[7][0] + xoffset) * xfac), (int) ((boxdata[7][1] + yoffset) * yfac));
		g.drawString("0,0,0", (int) ((boxdata[3][0] + xoffset) * xfac), (int) ((boxdata[3][1] + yoffset) * yfac));
	}


	/**
	 * Precomputes some aspects of the data for drawing.
	 */
	void prepareData() {
		// use the data from posarr and the rotation angles xz and yz to
		double xfac = 1;
		double yfac = 1;
		double xoffset;
		double yoffset;
		float[][] tposarrtmp = data.posarrtmp;
		float[][] posarr = data.posarr;

		// compute the new positions and save them to posarrtmp
		for (int i = data.elements; --i >= 0;) {
			// same as for boxdata
			tposarrtmp[i][0] = (float) (data.myrotmtx[0][0] * posarr[i][0] + data.myrotmtx[0][1] * posarr[i][1] + data.myrotmtx[0][2]
					* posarr[i][2]);
			tposarrtmp[i][1] = (float) (data.myrotmtx[1][0] * posarr[i][0] + data.myrotmtx[1][1] * posarr[i][1] + data.myrotmtx[1][2]
					* posarr[i][2]);
			data.drawarrtmp[i][0] = (int) tposarrtmp[i][0];
			data.drawarrtmp[i][1] = (int) tposarrtmp[i][1];
			// I don't need to calculate the 3rd dimention
			// tposarrtmp[i][2]=(float)(myrotmtx[2][0]*posarr[i][0]+myrotmtx[2][1]*posarr[i][1]+myrotmtx[2][2]*posarr[i][2]);
			// now add perspective using posarrtmp[i][2]
			// maxdiff is the maximum value for any coordinate in the system
		}
		
		// we can use the same code for zoom-on-selected and normal display if we set considered_sequences to the range
		// 0..<data.elements> in the normal case
		int[] considered_sequences;
		if (parent.isZoomingOnSelectedSequences() && (parent.getNumberOfSelectedSequences() > 1)) {
			considered_sequences = data.selectedSequencesIndices;

		} else {
			considered_sequences = new int[data.elements];
			for (int i = 0; i < considered_sequences.length; i++) {
				considered_sequences[i] = i;
			}
		}
		
		// determine the X and Y range we need to consider
		maxX = tposarrtmp[considered_sequences[0]][0];
		maxY = tposarrtmp[considered_sequences[0]][1];
		minX = maxX;
		minY = maxY;
		for (int i = considered_sequences.length; --i >= 0;) {
			if (maxX < tposarrtmp[considered_sequences[i]][0]) {
				maxX = tposarrtmp[considered_sequences[i]][0];
			} else if (minX > tposarrtmp[considered_sequences[i]][0]) {
				minX = tposarrtmp[considered_sequences[i]][0];
			}
			if (maxY < tposarrtmp[considered_sequences[i]][1]) {
				maxY = tposarrtmp[considered_sequences[i]][1];
			} else if (minY > tposarrtmp[considered_sequences[i]][1]) {
				minY = tposarrtmp[considered_sequences[i]][1];
			}
		}
		
		xfac = drawWidth / (maxX - minX);
		yfac = drawHeight / (maxY - minY);
		xoffset = (-minX);
		yoffset = (-minY);
		if (maxX - minX == 0) {
			System.out.println("isZero (in Zoom)!");
		}
		
		// even though we might be zoomed in on selected, unselected sequences must be shown, too
		for (int i = 0; i < data.elements; i++) {
			tposarrtmp[i][0] = (float) (((tposarrtmp[i][0] + xoffset) * xfac));
			tposarrtmp[i][1] = (float) (((tposarrtmp[i][1] + yoffset) * yfac));
			data.drawarrtmp[i][0] = (int) tposarrtmp[i][0];
			data.drawarrtmp[i][1] = (int) tposarrtmp[i][1];
			// our projection ignores the 3rd dimension
		}
		originpos[0] = (int) ((xoffset * xfac));
		originpos[1] = (int) ((yoffset * yfac));
		
		xscale = xfac;
		yscale = yfac;
	}

	// ---------------------------------------

	void drawdata(Graphics g1d) {
		synchronized (parent.dataLock) {
			drawdata_synchronized(g1d);
		}
	}

	void drawdata_synchronized(Graphics g1d) {

		Graphics2D g = (Graphics2D) g1d;
		
		if (parent.antialiasingcheckboxmenuitem.isSelected()) {
			g.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING, java.awt.RenderingHints.VALUE_ANTIALIAS_ON);
		} else {
			g.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING,
					java.awt.RenderingHints.VALUE_ANTIALIAS_OFF);
		}
		
		
		if (data.elements > 0) {
			prepareData();

			if (parent.showConnectionsCheckbox.isSelected()) {
				// here i am looking for a speedup by eliminating the multiple "for" loops
				// to this effect the drawing order is computed in an earlier step. here I just have to loop through
				// the 2d array and draw the elements in the according color

				if (parent.showConnectionsFrustrationColoredCheckbox.isSelected()) {
					drawConnectionsFrustrationColored(g);

				} else {
					drawConnections(g);
				}
			}
			
			if (parent.showSequenceNamesCheckbox.isSelected()) {
				drawSequenceNames(g);
				
			} else if (parent.showSequenceNumbersCheckbox.isSelected()) {
				drawSequenceNumbers(g);
			}
			
			// now draw the sequence dots
			if (parent.dotsfirst) {
				if (parent.colorDotsBySequenceLengthCheckbox.isSelected()) {
					drawDotsColoredBySequenceLength(g);
				
				} else {
					drawDots(g);
				}
			}
			
			// draw the sequences from older selection steps
			// now draw the elements in the selectvec
			drawOldSelection(g);

			// draw the sequence groups
			if (data.showSequenceGroups) {
				drawSequenceGroups(g);
			}

			if (parent.groupseqs != null) {// draw the sequences from the currently selected group
				drawSequencesOfCurrentlySelectedGroup(g);
			}
			
			// draw the current selection
			if (parent.clusterconf == null) {
				drawCurrentlySelectedSequences(g);

			} else {
				drawCurrentlySelectedSequenceClusterconf(g);
			}
						
			if (parent.blastselectseqs.length > 0) {
				drawBlastHits(g);
			}
			
			if (!parent.dotsfirst) {
				if (parent.colorDotsBySequenceLengthCheckbox.isSelected()) {
					drawDotsColoredBySequenceLength(g);
				
				} else {
					drawDots(g);
				}
			}
		}
		
		if (parent.inSelectionMode() && parent.mouse_is_pressed) {
			drawSelectionRectangle(g);
		}
	}
	
	/**
	 * Draw connections with the setup color gradient.
	 * 
	 * @param g
	 */
	private void drawConnections(Graphics g) {
		// color the lines by their blast P-value
		int j;
		int[] currdraw;
		int vecsize;
		if (data.draworder == null || data.draworder.size() < 1) {
			data.draworder = getLineDrawOrder(data.attractionValues, colornum);
		}
		for (int i = 0; i < colornum; i++) {
			vecsize = data.draworder.get(i).size();
			g.setColor(data.colorarr[i]);
			g.drawString("-", (5 + i) * fontwidth - xTranslate, fontsize - yTranslate);
			// draw all the vector elements
			for (j = vecsize; --j >= 0;) {
				currdraw = data.draworder.get(i).get(j);
				g.drawLine(data.drawarrtmp[currdraw[0]][0], data.drawarrtmp[currdraw[0]][1],
						data.drawarrtmp[currdraw[1]][0], data.drawarrtmp[currdraw[1]][1]);
			}
		}
		g.setColor(mainColor);
		g.drawString("Worst", -xTranslate, fontsize - yTranslate);
		g.drawString("Best", ((colornum + 5) * fontwidth) - xTranslate, fontsize - yTranslate);
	}
	
	/**
	 * Draw connections colored by their "frustration" (i.e. too long or too short for P-value)
	 * 
	 * @param g
	 */
	private void drawConnectionsFrustrationColored(Graphics g) {
		String key;
		int[] tmparr;
		MinimalAttractionValue currfrust;
		
		frustration = getfrustration(data.attractionValues, data.posarr);
		
		if (data.draworder.size() < 1) {
			// draw the lines in the same order as before
			data.draworder = getLineDrawOrder(data.attractionValues, colornum);
		}
		
		for (int i = 0; i < colornum; i++) {

			// draw all the vector elements
			for (int j = 0; j < data.draworder.get(i).size(); j++) {
				tmparr = (int[]) data.draworder.get(i).get(j);
				key = tmparr[0] + "_" + tmparr[1];

				if (frustration.containsKey(key)) {
					currfrust = (MinimalAttractionValue) frustration.get(key);
			
					if (currfrust.att > 0) {
						// then the coloring should be in blue(i.e. too short for attval)
						g.setColor(new Color((1 - currfrust.att), (1 - currfrust.att), 1));
					} else {
						// then the coloring should be in red(i.e. too long for attval)
						g.setColor(new Color(1, (1 + currfrust.att), (1 + currfrust.att)));
					}
					
					g.drawLine(data.drawarrtmp[currfrust.query][0], data.drawarrtmp[currfrust.query][1],
							data.drawarrtmp[currfrust.hit][0], data.drawarrtmp[currfrust.hit][1]);
				
				} else {
					System.err.println("no value for frustration key '" + key + "'");
				}
			}
		}
		
		g.setColor(mainColor);
		g.drawString("Worst", -xTranslate, fontsize - yTranslate);
		g.drawString("Best", ((colornum + 5) * fontwidth) - xTranslate, fontsize - yTranslate);
	}

	/**
	 * Draws the name of each sequence next to its dot.
	 * 
	 * @param g
	 */
	private void drawSequenceNames(Graphics g) {
		String[] namearr = data.sequence_names;
		if (parent.getNumberOfSelectedSequences() == 0) {
			for (int i = 0; i < data.elements; i++) {
				g.drawString(String.valueOf(i) + "-" + namearr[i], (int) data.drawarrtmp[i][0],
						(int) data.drawarrtmp[i][1]);
			}
		} else {
			for (int i = 0; i < parent.getNumberOfSelectedSequences(); i++) {
				g.drawString(String.valueOf(i) + "-" + namearr[data.selectedSequencesIndices[i]],
						(int) data.drawarrtmp[data.selectedSequencesIndices[i]][0],
						(int) data.drawarrtmp[data.selectedSequencesIndices[i]][1]);
			}
		}
	}
	
	/**
	 * Draws the number of each sequence next to its dot.
	 * 
	 * @param g
	 */
	private void drawSequenceNumbers(Graphics g) {
		for (int i = 0; i < data.elements; i++) {
			g.drawString(String.valueOf(i), (int) data.drawarrtmp[i][0], (int) data.drawarrtmp[i][1]);
		}
	}
	
	/**
	 * Draws the dots for each sequence in the foreground color.
	 * 
	 * @param g
	 */
	private void drawDots(Graphics g) {

		float halfdot = (float) (data.dotsize / 2);

		g.setColor(mainColor);

		for (int i = 0; i < data.elements; i++) {
			g.fillOval((int) (data.drawarrtmp[i][0] - halfdot), (int) (data.drawarrtmp[i][1] - halfdot), data.dotsize,
					data.dotsize);
		}
	}

	/**
	 * Draws the dots for each sequence colored by sequence length
	 * 
	 * @param g
	 */
	private void drawDotsColoredBySequenceLength(Graphics g) {

		float halfdot = (float) (data.dotsize / 2);

		for (int i = 0; i < data.elements; i++) {
			g.setColor(new Color(1 - data.seqlengths[i], 1 - data.seqlengths[i], data.seqlengths[i]));
			g.fillOval((int) (data.drawarrtmp[i][0] - halfdot), (int) (data.drawarrtmp[i][1] - halfdot), data.dotsize,
					data.dotsize);
		}
	}
	
	/**
	 * TODO: find out what this method is good for.
	 * 
	 * @param g
	 */
	private void drawOldSelection(Graphics g) {

		float halfoval = (float) (data.ovalsize / 2);

		UnusedWeirdClassPreviouslyCalledSelectclass tmpselectclass;

		for (int i = parent.selectvec.size() - 1; i >= 0; i--) {
			tmpselectclass = parent.selectvec.elementAt(i);
			g.setColor(tmpselectclass.color);
			int tmpelements = tmpselectclass.selectednames.length;
			for (int j = 0; j < tmpelements; j++) {
				g.fillOval((int) (data.drawarrtmp[tmpselectclass.selectednames[j]][0] - halfoval),
						(int) (data.drawarrtmp[tmpselectclass.selectednames[j]][1] - halfoval), data.ovalsize,
						data.ovalsize);
			}
		}
	}

	/**
	 * Draws the sequence groups.
	 * 
	 * @param g
	 */
	private void drawSequenceGroups(Graphics g) {
		SequenceGroup group;

		// draw groups from bottom to top, so that the topmost ones get painted over the lower ones
		for (int i = data.seqgroupsvec.size() - 1; i >= 0; i--) {
			group = data.seqgroupsvec.elementAt(i);

			if (group.hide) {
				continue;
			}

			g.setColor(group.color);

			if (group.type == 0) {
				
				drawShapeOval(g, group);

			} else {
				drawShapePolygon(g, group);
			}
		}
	}
	
	/**
	 * Draws the group's oval shape.
	 * 
	 * @param g
	 * @param group
	 *            The sequence group.
	 */
	private void drawShapeOval(Graphics g, SequenceGroup group) {

		float halfgroup = (float) group.size / 2;

		for (int j = group.sequences.length - 1; j >= 0; j--) {

			if (group.sequences[j] < data.elements) {
				g.fillOval((int) (data.drawarrtmp[group.sequences[j]][0] - halfgroup),
						(int) (data.drawarrtmp[group.sequences[j]][1] - halfgroup), group.size, group.size);

			} else {
				System.err.println("sequence number " + group.sequences[j]
						+ " is not present in file; removing entry from group");
				group.remove(j);
			}
		}
	}
	
	/**
	 * Draws the group's polygon shape.
	 * 
	 * @param g
	 * @param group
	 *            The sequence group.
	 */
	private void drawShapePolygon(Graphics g, SequenceGroup group) {

		posnum = group.polygon[2][0];// the number of points
		xposarr = new int[posnum];
		yposarr = new int[posnum];

		for (int j = group.sequences.length - 1; j >= 0; j--) {

			if (group.sequences[j] < data.elements) {

				for (int k = 0; k < posnum; k++) {
					xposarr[k] = (int) (group.polygon[0][k] + data.drawarrtmp[group.sequences[j]][0]);
					yposarr[k] = (int) (group.polygon[1][k] + data.drawarrtmp[group.sequences[j]][1]);
				}

				g.fillPolygon(xposarr, yposarr, posnum);

			} else {
				System.err.println("sequence number " + group.sequences[j]
						+ " is not present in file; removing entry from group");
				group.remove(j);
			}
		}
	}
	
	/**
	 * Draws the sequence circles for the sequences in the currently selected groups.
	 * 
	 * @param g
	 */
	private void drawSequencesOfCurrentlySelectedGroup(Graphics g) {

		float halfoval = (float) (data.ovalsize / 2);

		g.setColor(parent.groupseqscolor);
		for (int i = parent.groupseqs.length; --i >= 0;) {
			g.fillOval((int) (data.drawarrtmp[parent.groupseqs[i]][0] - halfoval),
					(int) (data.drawarrtmp[parent.groupseqs[i]][1] - halfoval), data.ovalsize, data.ovalsize);
		}
	}
	
	/**
	 * Draws the currently selected sequences.
	 * 
	 * @param g
	 */
	private void drawCurrentlySelectedSequences(Graphics g) {

		float halfoval = (float) (data.ovalsize / 2);
		
		for (int i = parent.getNumberOfSelectedSequences() - 1; i >= 0; i--) {
			g.setColor(selectedSequenceCircleColor);
			g.fillOval((int) (data.drawarrtmp[data.selectedSequencesIndices[i]][0] - halfoval),
					(int) (data.drawarrtmp[data.selectedSequencesIndices[i]][1] - halfoval), data.ovalsize,
					data.ovalsize);
			
			g.setColor(mainColor);
			g.drawOval((int) (data.drawarrtmp[data.selectedSequencesIndices[i]][0] - halfoval),
					(int) (data.drawarrtmp[data.selectedSequencesIndices[i]][1] - halfoval), data.ovalsize,
					data.ovalsize);
		}
	}
	
	private void drawCurrentlySelectedSequenceClusterconf(Graphics g) {

		float halfoval = (float) (data.ovalsize / 2);

		for (int i = parent.getNumberOfSelectedSequences() - 1; i >= 0; i--) {
			if (parent.clusterconf[i] < 0) {
				parent.clusterconf[i] = 0;
			}
			g.setColor(new Color((int) (selectedSequenceCircleColor.getRed() * parent.clusterconf[i]),
					(int) (selectedSequenceCircleColor.getGreen() * parent.clusterconf[i]),
					(int) (selectedSequenceCircleColor.getBlue() * parent.clusterconf[i])));
			g.fillOval((int) (data.drawarrtmp[data.selectedSequencesIndices[i]][0] - halfoval),
					(int) (data.drawarrtmp[data.selectedSequencesIndices[i]][1] - halfoval), data.ovalsize,
					data.ovalsize);
			
			g.setColor(mainColor);
			g.drawOval((int) (data.drawarrtmp[data.selectedSequencesIndices[i]][0] - halfoval),
					(int) (data.drawarrtmp[data.selectedSequencesIndices[i]][1] - halfoval), data.ovalsize,
					data.ovalsize);
			
		}
	}
	
	/**
	 * Draws Blast Hits
	 * @param g
	 */
	private void drawBlastHits(Graphics g) {
		float halfoval = (float) (data.ovalsize / 2);

		if (parent.showBlastHitNumberCheckbox.isSelected()) {
			g.setColor(blastHitNumberColor);
			// write the names to screen
			for (int i = parent.blastselectseqs.length - 1; i >= 0; i--) {
				g.drawString(String.valueOf(parent.blastselectseqs[i]),
						(int) data.drawarrtmp[parent.blastselectseqs[i]][0],
						(int) data.drawarrtmp[parent.blastselectseqs[i]][1]);
			}
		}

		// draw the blast query (entry 0) in one color, then the hits in another one and outline them
		g.setColor(blastQueryCircleColor);
		for (int i = 0; i < parent.blastselectseqs.length; i++) {
			if (i > 0) {
				g.setColor(blastHitCircleColor);
			}
			
			g.fillOval((int) (data.drawarrtmp[parent.blastselectseqs[i]][0] - halfoval),
					(int) (data.drawarrtmp[parent.blastselectseqs[i]][1] - halfoval), data.ovalsize, data.ovalsize);
			
			g.setColor(mainColor);
			g.drawOval((int) (data.drawarrtmp[parent.blastselectseqs[i]][0] - halfoval),
					(int) (data.drawarrtmp[parent.blastselectseqs[i]][1] - halfoval), data.ovalsize, data.ovalsize);
		}
	}
	
	/**
	 * Draws the selection rectangle while making a selection.
	 * 
	 * @param g
	 */
	private void drawSelectionRectangle(Graphics g) {

		g.setColor(selectionRectangleColor);

		if (parent.currmousepos[0] < parent.selectstart[0]) {
			draw[0] = parent.currmousepos[0];
			draw[2] = parent.selectstart[0];

		} else {
			draw[0] = parent.selectstart[0];
			draw[2] = parent.currmousepos[0];
		}

		if (parent.currmousepos[1] < parent.selectstart[1]) {
			draw[1] = parent.currmousepos[1];
			draw[3] = parent.selectstart[1];

		} else {
			draw[1] = parent.selectstart[1];
			draw[3] = parent.currmousepos[1];
		}

		g.drawRect(draw[0], draw[1], (draw[2] - draw[0]), (draw[3] - draw[1]));
	}

	/**
	 * bin the lines connecting datapoints by the color they are assigned (makes subsequent drawing quicker)
	 * 
	 * @param myattvals
	 * @param colornum
	 * @return
	 */
	private ArrayList<ArrayList<int[]>> getLineDrawOrder(MinimalAttractionValue[] myattvals, int colornum) {

		ArrayList<ArrayList<int[]>> retarr = new ArrayList<ArrayList<int[]>>();
		for (int i = 0; i < colornum; i++) {
			retarr.add(new ArrayList<int[]>());
		}
		int[] mydraw = new int[2];// connection between sequences i & j

		float[] colorcutoffs = data.colorcutoffs;
		for (int i = myattvals.length; --i >= 0;) {

			if ((myattvals[i].att == 0) || (myattvals[i].att <= colorcutoffs[0])) {// if I am below the lowest to
																					// draw
				continue;
			}

			for (int j = 0; j < 10; j++) {
				if (j < 9) {
					if (myattvals[i].att <= colorcutoffs[j + 1]) { // from j to j + 1
						mydraw = new int[2];
						mydraw[0] = myattvals[i].query;
						mydraw[1] = myattvals[i].hit;
						retarr.get(j).add(mydraw);
						break;
					}
				} else { // j == 9
					if (myattvals[i].att > colorcutoffs[j]) { // better than value at j
						mydraw = new int[2];
						mydraw[0] = myattvals[i].query;
						mydraw[1] = myattvals[i].hit;
						retarr.get(j).add(mydraw);
						break;
					}
				}
			}
		}

		return retarr;
	}

	/**
	 * Gets wether each connection is longer or shorter than expected based on the attraction value.
	 * 
	 * @param attvals
	 * @param posarr
	 * @return
	 */
	private HashMap<String, MinimalAttractionValue> getfrustration(MinimalAttractionValue[] attvals, float[][] posarr) {
		HashMap<String, MinimalAttractionValue> retarr = new HashMap<String, MinimalAttractionValue>(
				(int) (attvals.length / 0.8) + 1, 0.8f);
		MinimalAttractionValue[] minattarr = new MinimalAttractionValue[attvals.length];
		int query_index, hit_index;
		float attraction;
		
		float avgatt = 0;
		float avgdist = 0;
		float maxval = 0;
		int counter = 0;
		
		String key;
		float tmpfloat;
		float tmpx, tmpy, tmpz;
		float tmpval;
		
		for (int i = 0; i < attvals.length; i++) {
			attraction = attvals[i].att;
			query_index = attvals[i].query;
			hit_index = attvals[i].hit;
			key = query_index + "_" + hit_index;
			if (attraction > 0) {// else don't do the calculations as I won't be drawing a line
				counter++;
				tmpx = posarr[query_index][0] - posarr[hit_index][0];
				tmpy = posarr[query_index][1] - posarr[hit_index][1];
				tmpz = posarr[query_index][2] - posarr[hit_index][2];
				tmpval = (float) java.lang.Math.sqrt(tmpx * tmpx + tmpy * tmpy + tmpz * tmpz);
				minattarr[i] = new MinimalAttractionValue(query_index, hit_index, tmpval);
				retarr.put(key, minattarr[i]);
				avgdist += tmpval;
				avgatt += attraction;
			}
		}
		
		avgdist /= counter;
		avgatt /= counter;
		tmpfloat = avgatt * avgdist;

		// note: the higher the attval, the smaller the dist!
		for (int i = 0; i < attvals.length; i++) {
			attraction = attvals[i].att;
			query_index = attvals[i].query;
			hit_index = attvals[i].hit;
			if (attraction > 0) {
				minattarr[i].att = (tmpfloat - (attraction * minattarr[i].att));
				if (maxval < minattarr[i].att) {
					maxval = minattarr[i].att;
				}
				if (maxval < -minattarr[i].att) {
					maxval = -minattarr[i].att;
				}
			}
		}

		for (int i = 0; i < attvals.length; i++) {
			if (attvals[i].att > 0) {
				minattarr[i].att /= maxval;
			}
		}
		
		return retarr;
	}

	/**
	 * @return the current zoom factor without pre-zoom, which is handles completely inside the class.
	 */
	public float getZoomFactor() {
		return zoomFactor / preZoomFactor;
	}
	
	/**
	 * Updates the zoom of the view to the currently set zoom value while keeping the same area centered.
	 * 
	 * @param old_zoom
	 *            The previous zoom value.
	 */
	protected void setZoom(float new_zoom_factor) throws IllegalStateException {
		if (new_zoom_factor <= 0) {
			throw new IllegalStateException("the zoom factor must be a float value >= 0, was "
					+ String.format("%s", new_zoom_factor));
		}
		
		// the old zoom factor always remains with the preZoomFactor current when the old zoom was set. 
		float old_zoom_factor = zoomFactor;
		zoomFactor = new_zoom_factor * preZoomFactor;
		
		int image_center_x = (int) ((getWidth() / 2 - xTranslate) / old_zoom_factor);
		int image_center_y = (int) ((getHeight() / 2 - yTranslate) / old_zoom_factor);

		xTranslate = (int) (getWidth() / 2 - (image_center_x * zoomFactor));
		yTranslate = (int) (getHeight() / 2 -(image_center_y * zoomFactor));

		repaint();
	}
	
	/**
	 * Resets the zoom to none.
	 */
	protected void resetZoom() {
		setZoom(1);
		centerGraph();
	}

	/**
	 * Centers the graph while maintaining the current zoom level.
	 */
	protected void centerGraph() {
		xTranslate = (int) (getWidth() / 2 * (1 - zoomFactor));
		yTranslate = (int) (getHeight() / 2 * (1 - zoomFactor));

		repaint();
	}
	
	/**
	 * Recomputes the pre-zoom factor to adjust to changed window size. Pre-zoom is used to achieve a fixed pixel
	 * padding.
	 */
	private void updatePreZoomFactor() {
		float old_zoom_factor = getZoomFactor();
		preZoomFactor = (getHeight() - 2 * PADDING_PIXEL) / getHeight();
		setZoom(old_zoom_factor);
	}

	/**
	 * Whenever the DrawArea is resized, we need to recompute the pre-zoom factor.
	 */
	@Override
	public void componentResized(ComponentEvent e) {
		updatePreZoomFactor();
		
		// try to repair the graph centering whenever there is no zoom.
		// TODO: find problem in setZoom(int) then remove this. 
		if (getZoomFactor() == 1) {
			centerGraph();
		}
	}

	@Override
	public void componentMoved(ComponentEvent e) {
	}

	@Override
	public void componentShown(ComponentEvent e) {
	}

	@Override
	public void componentHidden(ComponentEvent e) {
	}
}