/*
 * aboutwindow.java
 *
 * Created on July 8, 2004, 3:13 PM
 */
package clans;
/**
 *
 * @author  tancred
 */
public class aboutwindow extends javax.swing.JDialog {
    
    /** Creates new form aboutwindow */
    public aboutwindow(java.awt.Frame parent, boolean modal) {
        super(parent, modal);
        initComponents();
    }
    
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        textarea = new javax.swing.JTextArea();
        topfield = new javax.swing.JTextField();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);

        textarea.setBackground(new java.awt.Color(204, 204, 204));
        textarea.setEditable(false);
        textarea.setLineWrap(true);
        textarea.setText("About CLANS.: Clans was developed as part of my PhD in biology at the  Max-Planck-Institute for Developmental Biology in Tuebingen, Germany. I frequently use this program and continuously extend it with features I  feel are liable to make my life easier.  Hope you enjoy using it.\nThe program has now been extended to allow analysis of microarray expression data and overlay functional annotation on the maps.\nMany thanks also go to the NCBI team who constantly keep changing the blast output format (which gives me the opportunity to change the blast parser over and over).\nTancred Frickey");
        getContentPane().add(textarea, java.awt.BorderLayout.CENTER);

        topfield.setEditable(false);
        topfield.setText("CLANS (27.08.2012)");
        getContentPane().add(topfield, java.awt.BorderLayout.NORTH);

        pack();
    }// </editor-fold>//GEN-END:initComponents
    
    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        new aboutwindow(new javax.swing.JFrame(), true).setVisible(true);
    }
    
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JTextArea textarea;
    private javax.swing.JTextField topfield;
    // End of variables declaration//GEN-END:variables
    
}
