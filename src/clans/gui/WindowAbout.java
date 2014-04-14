package clans.gui;

public class WindowAbout extends javax.swing.JDialog {

	/**
	 * 
	 */
	private static final long serialVersionUID = -7085934144448877843L;

	private javax.swing.JTextArea textarea;
	private javax.swing.JTextField topfield;

	public WindowAbout(java.awt.Frame parent, boolean modal) {
		super(parent, modal);
		initComponents();
	}

	private void initComponents() {

		textarea = new javax.swing.JTextArea();
		topfield = new javax.swing.JTextField();

		setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);

		textarea.setBackground(new java.awt.Color(204, 204, 204));
		textarea.setEditable(false);
		textarea.setLineWrap(true);
		textarea.setPreferredSize(new java.awt.Dimension(600, 400));

		textarea.setText("About CLANS.: Clans was developed as part of my PhD in biology at the  Max-Planck-Institute for Developmental Biology in Tuebingen, Germany. I frequently use this program and continuously extend it with features I  feel are liable to make my life easier.  Hope you enjoy using it.\nThe program has now been extended to allow analysis of microarray expression data and overlay functional annotation on the maps.\nMany thanks also go to the NCBI team who constantly keep changing the blast output format (which gives me the opportunity to change the blast parser over and over).\nTancred Frickey");

		getContentPane().add(textarea, java.awt.BorderLayout.CENTER);

		topfield.setEditable(false);
		topfield.setText("CLANS (2014-04-07)");

		getContentPane().add(topfield, java.awt.BorderLayout.NORTH);

		pack();
	}

}
