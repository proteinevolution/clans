/*
 * errwindow.java
 *
 * Created on July 9, 2004, 9:26 AM
 */
package clans;
/**
 *
 * @author  tancred
 */
public class errwindow extends javax.swing.JDialog {
    
    /** Creates new form errwindow */
    public errwindow(java.awt.Frame parent, boolean modal,String text) {
        super(parent, modal);
        this.text=text;
        initComponents();
    }
    
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    private void initComponents() {//GEN-BEGIN:initComponents
        jPanel1 = new javax.swing.JPanel();
        jScrollPane1 = new javax.swing.JScrollPane();
        errarea = new javax.swing.JTextArea();
        errarea.setText(text);
        jPanel2 = new javax.swing.JPanel();
        jTextField1 = new javax.swing.JTextField();
        jPanel3 = new javax.swing.JPanel();
        okbutton = new javax.swing.JButton();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        jPanel1.setLayout(new javax.swing.BoxLayout(jPanel1, javax.swing.BoxLayout.X_AXIS));

        jPanel1.setPreferredSize(new java.awt.Dimension(600, 400));
        errarea.setEditable(false);
        jScrollPane1.setViewportView(errarea);

        jPanel1.add(jScrollPane1);

        getContentPane().add(jPanel1, java.awt.BorderLayout.CENTER);

        jPanel2.setLayout(new javax.swing.BoxLayout(jPanel2, javax.swing.BoxLayout.X_AXIS));

        jTextField1.setEditable(false);
        jTextField1.setText("Errors occurred so far (check command line)");
        jPanel2.add(jTextField1);

        getContentPane().add(jPanel2, java.awt.BorderLayout.NORTH);

        jPanel3.setLayout(new java.awt.GridLayout(1, 0));

        okbutton.setText("OK (Close)");
        okbutton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                okbuttonActionPerformed(evt);
            }
        });

        jPanel3.add(okbutton);

        getContentPane().add(jPanel3, java.awt.BorderLayout.SOUTH);

        pack();
    }//GEN-END:initComponents

    private void okbuttonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_okbuttonActionPerformed
        this.setVisible(false);
        dispose();
    }//GEN-LAST:event_okbuttonActionPerformed
    
    /**
     * @param args the command line arguments
     */
    //public static void main(String args[]) {
    //    new errwindow(new javax.swing.JFrame(), true).show();
    //}
    
    String text;
    
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JTextArea errarea;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JPanel jPanel3;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JTextField jTextField1;
    private javax.swing.JButton okbutton;
    // End of variables declaration//GEN-END:variables
    
}
