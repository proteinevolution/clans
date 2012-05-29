/*
 * querydialog.java
 *
 * Created on July 7, 2003, 10:53 AM
 */
package clans;
/**
 *
 * @author  tancred
 */
public class querydialog extends javax.swing.JDialog {
    
    /** Creates new form querydialog */
    public querydialog(java.awt.Frame parent, boolean modal) {
        super(parent, modal);
        initComponents();
    }
    
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    private void initComponents() {//GEN-BEGIN:initComponents
        Label = new javax.swing.JLabel();
        queryfield = new javax.swing.JTextField();
        okbutton = new javax.swing.JButton();

        getContentPane().setLayout(new java.awt.GridLayout(3, 1));

        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosing(java.awt.event.WindowEvent evt) {
                closeDialog(evt);
            }
        });

        Label.setText("Enter a query:");
        getContentPane().add(Label);

        queryfield.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                queryfieldActionPerformed(evt);
            }
        });

        getContentPane().add(queryfield);

        okbutton.setText("OK");
        okbutton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                okbuttonActionPerformed(evt);
            }
        });

        getContentPane().add(okbutton);

        pack();
    }//GEN-END:initComponents
    
    private void okbuttonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_okbuttonActionPerformed
        querystring=queryfield.getText();
        setVisible(false);
        dispose();
    }//GEN-LAST:event_okbuttonActionPerformed
    
    private void queryfieldActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_queryfieldActionPerformed
        querystring=queryfield.getText();
        setVisible(false);
        dispose();
    }//GEN-LAST:event_queryfieldActionPerformed
    
    /** Closes the dialog */
    private void closeDialog(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_closeDialog
        setVisible(false);
        dispose();
    }//GEN-LAST:event_closeDialog
    
    /**
     * @param args the command line arguments
     */
    //public static void main(String args[]) {
    //   new querydialog(new javax.swing.JFrame(), true).show();
    //}
    
    
    static String querystring="";
    
    public static String getquery(){
        new querydialog(new javax.swing.JFrame(), true).setVisible(true);
        return querystring;
    }
    
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JLabel Label;
    private javax.swing.JButton okbutton;
    private javax.swing.JTextField queryfield;
    // End of variables declaration//GEN-END:variables
    
}
