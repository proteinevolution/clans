package clans.gui;

/**
 * 
 * @author tancred
 */
public class WindowShowCopyPasteableSequences extends javax.swing.JDialog {

	/**
	 * 
	 */
	private static final long serialVersionUID = 8972158369206381300L;

	private javax.swing.JScrollPane jScrollPane1;
	private javax.swing.JTextArea textarea;

	public WindowShowCopyPasteableSequences(java.awt.Frame parent, StringBuffer text) {
		super(parent, false);
		this.setTitle("Sequence window");
		initComponents();
		textarea.setText(text.toString());
		pack();
	}

	private void initComponents() {
		jScrollPane1 = new javax.swing.JScrollPane();
		textarea = new javax.swing.JTextArea();

		getContentPane().setLayout(
				new javax.swing.BoxLayout(getContentPane(),
						javax.swing.BoxLayout.X_AXIS));

		setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
		jScrollPane1.setPreferredSize(new java.awt.Dimension(200, 400));
		jScrollPane1.setViewportView(textarea);

		getContentPane().add(jScrollPane1);

		pack();
	}
}