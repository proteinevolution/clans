/*
 * affyselectwindow.java
 *
 * Created on January 16, 2006, 12:26 PM
 */
package clans;
import javax.swing.*;
/**
 *
 * @author  tancred
 */
public class affyselectwindow extends javax.swing.JFrame {
    
    /** Creates new form affyselectwindow */
    public affyselectwindow(affyplotdialog parent) {
        this.parent=parent;
        myelements=new listelem[parent.draw1.datnum];
        for(int i=parent.draw1.datnum;--i>=0;){
            myelements[i]=new listelem();
            myelements[i].name=((datapoint[])parent.datlist.get(i))[0].name;
            myelements[i].value=parent.draw1.drawdat[i];
        }//end for i
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
        jList1 = new javax.swing.JList(myelements);
        jList1.setCellRenderer(new listcellrenderer());
        jPanel2 = new javax.swing.JPanel();
        selectbutton = new javax.swing.JButton();
        deselectbutton = new javax.swing.JButton();
        changecolorbutton = new javax.swing.JButton();

        setTitle("Select Identifiers");
        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosing(java.awt.event.WindowEvent evt) {
                exitForm(evt);
            }
        });

        jPanel1.setLayout(new javax.swing.BoxLayout(jPanel1, javax.swing.BoxLayout.X_AXIS));

        jList1.addListSelectionListener(new javax.swing.event.ListSelectionListener() {
            public void valueChanged(javax.swing.event.ListSelectionEvent evt) {
                jList1ValueChanged(evt);
            }
        });

        jScrollPane1.setViewportView(jList1);

        jPanel1.add(jScrollPane1);

        getContentPane().add(jPanel1, java.awt.BorderLayout.CENTER);

        jPanel2.setLayout(new java.awt.GridLayout(1, 0));

        selectbutton.setText("Show");
        selectbutton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                selectbuttonActionPerformed(evt);
            }
        });

        jPanel2.add(selectbutton);

        deselectbutton.setText("Hide");
        deselectbutton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                deselectbuttonActionPerformed(evt);
            }
        });

        jPanel2.add(deselectbutton);

        changecolorbutton.setText("Change Color");
        changecolorbutton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                changecolorbuttonActionPerformed(evt);
            }
        });

        jPanel2.add(changecolorbutton);

        getContentPane().add(jPanel2, java.awt.BorderLayout.SOUTH);

        pack();
    }//GEN-END:initComponents

    private void jList1ValueChanged(javax.swing.event.ListSelectionEvent evt) {//GEN-FIRST:event_jList1ValueChanged
        // update the selection change in parent
        int[] selecteds=jList1.getSelectedIndices();
        int num=java.lang.reflect.Array.getLength(selecteds);
        if(num==0){
            return;
        }
        for(int i=parent.draw1.datnum;--i>=0;){
            parent.draw1.selecteds[i]=0;
        }//end for i
        for(int i=0;i<num;i++){
            parent.draw1.selecteds[selecteds[i]]=1;
        }//end for i
        parent.repaint();
    }//GEN-LAST:event_jList1ValueChanged

    private void changecolorbuttonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_changecolorbuttonActionPerformed
        // change the color of this element
        int[] selecteds=jList1.getSelectedIndices();
        int num=java.lang.reflect.Array.getLength(selecteds);
        if(num==0){
            return;
        }
        java.awt.Color newcolor=JColorChooser.showDialog(this,"Select new color",parent.draw1.colorarr[selecteds[0]]);
        for(int i=0;i<num;i++){
            parent.draw1.colorarr[selecteds[i]]=newcolor;
        }//end for i
        parent.repaint();
        repaint();
    }//GEN-LAST:event_changecolorbuttonActionPerformed
    
    private void deselectbuttonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_deselectbuttonActionPerformed
        int[] selecteds=jList1.getSelectedIndices();
        int num=java.lang.reflect.Array.getLength(selecteds);
        if(num==0){
            return;
        }
        for(int i=0;i<num;i++){
            parent.draw1.drawdat[selecteds[i]]=0;
            myelements[selecteds[i]].value=0;
        }//end for i
        parent.repaint();
        repaint();
    }//GEN-LAST:event_deselectbuttonActionPerformed
    
    private void selectbuttonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_selectbuttonActionPerformed
        int[] selecteds=jList1.getSelectedIndices();
        int num=java.lang.reflect.Array.getLength(selecteds);
        if(num==0){
            return;
        }
        for(int i=0;i<num;i++){
            parent.draw1.drawdat[selecteds[i]]=1;
            myelements[selecteds[i]].value=1;
        }//end for i
        parent.repaint();
        repaint();
    }//GEN-LAST:event_selectbuttonActionPerformed
    
    /** Exit the Application */
    private void exitForm(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_exitForm
        this.setVisible(false);
        this.dispose();
    }//GEN-LAST:event_exitForm
    
    /**
     * @param args the command line arguments
     */
    
    affyplotdialog parent;
    listelem[] myelements=null;
    
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton changecolorbutton;
    private javax.swing.JButton deselectbutton;
    public javax.swing.JList jList1;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JButton selectbutton;
    // End of variables declaration//GEN-END:variables
    
    class listelem{
        public listelem(){}
        
        String name="none";
        int value=1;
        
    }
        
    class listcellrenderer extends javax.swing.JLabel implements javax.swing.ListCellRenderer {
        // This is the only method defined by ListCellRenderer.
        // We just reconfigure the JLabel each time we're called.
        
        public java.awt.Component getListCellRendererComponent(
        JList list,
        Object value,            // value to display
        int index,               // cell index
        boolean isSelected,      // is the cell selected
        boolean cellHasFocus)    // the list and the cell have the focus
        {
            String s = ((listelem)value).name;
            if(((listelem)value).value!=1){
                s="(hidden) "+s;
            }
            setText(s);
            setForeground(parent.draw1.colorarr[index]);
            if (isSelected){
                setBackground(list.getSelectionBackground());
                //setForeground(list.getSelectionForeground());
            }else{
                setBackground(list.getBackground());
                //setForeground(list.getForeground());
            }
            setEnabled(list.isEnabled());
            setFont(list.getFont());
            setOpaque(true);
            return this;
        }
        
        
    }
    
}//end class

