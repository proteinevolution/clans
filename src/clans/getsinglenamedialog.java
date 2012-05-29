/*
 * shownamedialog.java
 *
 * Created on December 17, 2002, 3:04 PM
 */
package clans;
import java.util.*;
/**
 *
 * @author  tancred
 */
public class getsinglenamedialog extends javax.swing.JDialog {
    
    /** Creates new form shownamedialog */
    public getsinglenamedialog(String[] namesarr, java.awt.Frame parent) {
        super(parent,true);
        this.namesarr=numberarr(namesarr);
        initComponents();
        this.jScrollPane2.getVerticalScrollBar().setUnitIncrement((int)(jPanel1.getHeight()/java.lang.reflect.Array.getLength(namesarr)));
    }
    
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    private void initComponents() {//GEN-BEGIN:initComponents
        jScrollPane2 = new javax.swing.JScrollPane();
        jPanel1 = new javax.swing.JPanel();
        seqnamelist = new javax.swing.JList(namesarr);
        buttonpanel = new javax.swing.JPanel();
        okbutton = new javax.swing.JButton();
        clearbutton = new javax.swing.JButton();
        searchbutton = new javax.swing.JButton();

        setTitle("Selected Sequences (3D)");
        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosing(java.awt.event.WindowEvent evt) {
                closeDialog(evt);
            }
        });

        jScrollPane2.setPreferredSize(new java.awt.Dimension(100, 500));
        jPanel1.setLayout(new java.awt.BorderLayout());

        seqnamelist.setFont(new java.awt.Font("Monospaced", 0, 10));
        seqnamelist.setSelectionMode(javax.swing.ListSelectionModel.SINGLE_SELECTION);
        jPanel1.add(seqnamelist, java.awt.BorderLayout.CENTER);

        jScrollPane2.setViewportView(jPanel1);

        getContentPane().add(jScrollPane2, java.awt.BorderLayout.CENTER);

        okbutton.setText("OK");
        okbutton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                okbuttonActionPerformed(evt);
            }
        });

        buttonpanel.add(okbutton);

        clearbutton.setText("Clear");
        clearbutton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                clearbuttonActionPerformed(evt);
            }
        });

        buttonpanel.add(clearbutton);

        searchbutton.setText("Search");
        searchbutton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                searchbuttonActionPerformed(evt);
            }
        });

        buttonpanel.add(searchbutton);

        getContentPane().add(buttonpanel, java.awt.BorderLayout.SOUTH);

        pack();
    }//GEN-END:initComponents

    private void searchbuttonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_searchbuttonActionPerformed
        // Add your handling code here:
        String query=getquery();
        int[] selectedseqs=getmatches(query,namesarr);
        if(java.lang.reflect.Array.getLength(selectedseqs)==0){
            javax.swing.JOptionPane.showMessageDialog(this,"No sequences found for '"+query+"'","ERROR",javax.swing.JOptionPane.ERROR_MESSAGE);
            return;
        }
        seqnamelist.setSelectedIndices(selectedseqs);
    }//GEN-LAST:event_searchbuttonActionPerformed

    private void clearbuttonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_clearbuttonActionPerformed
        seqnamelist.setSelectedIndices(new int[0]);
    }//GEN-LAST:event_clearbuttonActionPerformed

    private void okbuttonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_okbuttonActionPerformed
        refseq=seqnamelist.getSelectedIndex();
        setVisible(false);
        dispose();
    }//GEN-LAST:event_okbuttonActionPerformed
    
    /** Closes the dialog */
    private void closeDialog(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_closeDialog
        setVisible(false);
        dispose();
    }//GEN-LAST:event_closeDialog
    
    public static int getrefseq(String[] seqnames) {
        new getsinglenamedialog(seqnames,new javax.swing.JFrame()).setVisible(true);
        return refseq;
    }
    
    
    String[] namesarr;
    static int refseq=-1;    
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JPanel buttonpanel;
    private javax.swing.JButton searchbutton;
    private javax.swing.JButton clearbutton;
    public javax.swing.JList seqnamelist;
    private javax.swing.JButton okbutton;
    private javax.swing.JScrollPane jScrollPane2;
    private javax.swing.JPanel jPanel1;
    // End of variables declaration//GEN-END:variables
    
    String getquery(){
        //get the regexp query string 
        String retstring="";
        retstring="(?i).*"+querydialog.getquery()+".*";
        return retstring;
    }//end getquery
    
    //--------------------------------------------------------------------------
    
    int[] getmatches(String query, String[] namesarr){
        //get the index of all names that match query
        int namesnum=java.lang.reflect.Array.getLength(namesarr);
        Vector tmpvec=new Vector();
        for(int i=0;i<namesnum;i++){
            if(namesarr[i].matches(query)){
                tmpvec.addElement(new Integer(i));
            }
        }//end for i
        int[] retarr=new int[tmpvec.size()];
        for(int i=0;i<tmpvec.size();i++){
            retarr[i]=((Integer)tmpvec.elementAt(i)).intValue();
        }
        return retarr;
    }//end getmatches
    
    //--------------------------------------------------------------------------
    
    String[] numberarr(String[] innames){
        int elements=java.lang.reflect.Array.getLength(innames);
        String[] retarr=new String[elements];
        int numlength=(String.valueOf(elements).length())+1;
        StringBuffer tmpstrbuff=new StringBuffer();
        for(int i=0;i<elements;i++){
            tmpstrbuff.setLength(0);
            tmpstrbuff.append(i);
            for(int j=tmpstrbuff.length();j<numlength;j++){
                tmpstrbuff.append(" ");
            }
            retarr[i]=tmpstrbuff+innames[i];
        }// end for i
        return retarr;
    }// end printnames

    void setselected(int[] selectednames){
        seqnamelist.setSelectedIndices(selectednames);
    }
    
}// end class
