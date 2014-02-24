package clans;

public class DialogChangeConnectionColors extends javax.swing.JDialog {

    /**
     * 
     */
    private static final long serialVersionUID = 1895643871437384612L;

    public DialogChangeConnectionColors(ClusteringWithGui parent, boolean modal) {
        super(parent, modal);
        this.parent = parent;
        this.colornum = colorarr.length;
        draw1 = new drawpanel(colorarr);
        initComponents();
        draw1.drawwidth = mainpanel.getWidth();
        draw1.drawheight = mainpanel.getHeight();
        draw1.elementwidth = ((float) draw1.drawwidth) / ((float) colornum);
        //System.out.println("test");
        if (parent.data.blasthits != null) {
            data2d = false;
            if (parent.data.usescval) {
                field1.setText(String.valueOf((((((parent.data.colorcutoffs[0]) * parent.data.p2attfactor) + parent.data.p2attoffset)))));
                field2.setText(String.valueOf((((((parent.data.colorcutoffs[1]) * parent.data.p2attfactor) + parent.data.p2attoffset)))));
                field3.setText(String.valueOf((((((parent.data.colorcutoffs[2]) * parent.data.p2attfactor) + parent.data.p2attoffset)))));
                field4.setText(String.valueOf((((((parent.data.colorcutoffs[3]) * parent.data.p2attfactor) + parent.data.p2attoffset)))));
                field5.setText(String.valueOf((((((parent.data.colorcutoffs[4]) * parent.data.p2attfactor) + parent.data.p2attoffset)))));
                field6.setText(String.valueOf((((((parent.data.colorcutoffs[5]) * parent.data.p2attfactor) + parent.data.p2attoffset)))));
                field7.setText(String.valueOf((((((parent.data.colorcutoffs[6]) * parent.data.p2attfactor) + parent.data.p2attoffset)))));
                field8.setText(String.valueOf((((((parent.data.colorcutoffs[7]) * parent.data.p2attfactor) + parent.data.p2attoffset)))));
                field9.setText(String.valueOf((((((parent.data.colorcutoffs[8]) * parent.data.p2attfactor) + parent.data.p2attoffset)))));
                field10.setText(String.valueOf((((((parent.data.colorcutoffs[9]) * parent.data.p2attfactor) + parent.data.p2attoffset)))));
            } else {
                field1.setText(String.valueOf((int) (((((parent.data.colorcutoffs[0]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10) + 0.5)));
                field2.setText(String.valueOf((int) (((((parent.data.colorcutoffs[1]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10) + 0.5)));
                field3.setText(String.valueOf((int) (((((parent.data.colorcutoffs[2]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10) + 0.5)));
                field4.setText(String.valueOf((int) (((((parent.data.colorcutoffs[3]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10) + 0.5)));
                field5.setText(String.valueOf((int) (((((parent.data.colorcutoffs[4]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10) + 0.5)));
                field6.setText(String.valueOf((int) (((((parent.data.colorcutoffs[5]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10) + 0.5)));
                field7.setText(String.valueOf((int) (((((parent.data.colorcutoffs[6]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10) + 0.5)));
                field8.setText(String.valueOf((int) (((((parent.data.colorcutoffs[7]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10) + 0.5)));
                field9.setText(String.valueOf((int) (((((parent.data.colorcutoffs[8]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10) + 0.5)));
                field10.setText(String.valueOf((int) (((((parent.data.colorcutoffs[9]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10) + 0.5)));
            }
        } else {//if I have only attraction values
            data2d = true;
            field1.setText(String.valueOf(parent.data.colorcutoffs[0]));//*parent.p2attfactor));
            field2.setText(String.valueOf(parent.data.colorcutoffs[1]));//*parent.p2attfactor));
            field3.setText(String.valueOf(parent.data.colorcutoffs[2]));//*parent.p2attfactor));
            field4.setText(String.valueOf(parent.data.colorcutoffs[3]));//*parent.p2attfactor));
            field5.setText(String.valueOf(parent.data.colorcutoffs[4]));//*parent.p2attfactor));
            field6.setText(String.valueOf(parent.data.colorcutoffs[5]));//*parent.p2attfactor));
            field7.setText(String.valueOf(parent.data.colorcutoffs[6]));//*parent.p2attfactor));
            field8.setText(String.valueOf(parent.data.colorcutoffs[7]));//*parent.p2attfactor));
            field9.setText(String.valueOf(parent.data.colorcutoffs[8]));//*parent.p2attfactor));
            field10.setText(String.valueOf(parent.data.colorcutoffs[9]));//*parent.p2attfactor));
        }
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    private void initComponents() {//GEN-BEGIN:initComponents
        jButton1 = new javax.swing.JButton();
        mainpanel = new javax.swing.JPanel();
        buttonpanel = new javax.swing.JPanel();
        worstlabel = new javax.swing.JLabel();
        updatebutton = new javax.swing.JButton();
        gradientbutton = new javax.swing.JButton();
        valgradientbutton = new javax.swing.JButton();
        closebutton = new javax.swing.JButton();
        bestlabel = new javax.swing.JLabel();
        jPanel1 = new javax.swing.JPanel();
        jPanel2 = new javax.swing.JPanel();
        infotextfield = new javax.swing.JTextField();
        textfieldpanel = new javax.swing.JPanel();
        field1 = new javax.swing.JTextField();
        field2 = new javax.swing.JTextField();
        field3 = new javax.swing.JTextField();
        field4 = new javax.swing.JTextField();
        field5 = new javax.swing.JTextField();
        field6 = new javax.swing.JTextField();
        field7 = new javax.swing.JTextField();
        field8 = new javax.swing.JTextField();
        field9 = new javax.swing.JTextField();
        field10 = new javax.swing.JTextField();

        jButton1.setText("jButton1");

        setTitle("Change Colors");
        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosing(java.awt.event.WindowEvent evt) {
                closeDialog(evt);
            }
        });

        mainpanel.setLayout(new javax.swing.BoxLayout(mainpanel, javax.swing.BoxLayout.X_AXIS));

        mainpanel.setBorder(new javax.swing.border.LineBorder(new java.awt.Color(0, 0, 0)));
        mainpanel.setPreferredSize(new java.awt.Dimension(400, 20));
        mainpanel.add(draw1);
        mainpanel.addHierarchyBoundsListener(new java.awt.event.HierarchyBoundsListener() {
            public void ancestorMoved(java.awt.event.HierarchyEvent evt) {
            }
            public void ancestorResized(java.awt.event.HierarchyEvent evt) {
                mainpanelAncestorResized(evt);
            }
        });
        mainpanel.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseReleased(java.awt.event.MouseEvent evt) {
                mainpanelMouseReleased(evt);
            }
        });

        getContentPane().add(mainpanel, java.awt.BorderLayout.CENTER);

        buttonpanel.setLayout(new java.awt.GridLayout(1, 0));

        worstlabel.setText("Worst");
        buttonpanel.add(worstlabel);

        updatebutton.setText("UPDATE");
        updatebutton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                updatebuttonActionPerformed(evt);
            }
        });

        buttonpanel.add(updatebutton);

        gradientbutton.setText("Color gradient");
        gradientbutton.setToolTipText("select the worst and best color, rest is computed");
        gradientbutton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                gradientbuttonActionPerformed(evt);
            }
        });

        buttonpanel.add(gradientbutton);

        valgradientbutton.setText("Value gradient");
        valgradientbutton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                valgradientbuttonActionPerformed(evt);
            }
        });

        buttonpanel.add(valgradientbutton);

        closebutton.setText("CLOSE");
        closebutton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                closebuttonActionPerformed(evt);
            }
        });

        buttonpanel.add(closebutton);

        bestlabel.setHorizontalAlignment(javax.swing.SwingConstants.RIGHT);
        bestlabel.setText("Best");
        bestlabel.setHorizontalTextPosition(javax.swing.SwingConstants.LEADING);
        buttonpanel.add(bestlabel);

        getContentPane().add(buttonpanel, java.awt.BorderLayout.SOUTH);

        jPanel1.setLayout(new java.awt.GridLayout(2, 0));

        jPanel2.setLayout(new javax.swing.BoxLayout(jPanel2, javax.swing.BoxLayout.X_AXIS));

        infotextfield.setEditable(false);
        infotextfield.setText("Connections with P-value better than 1E-X are drawn in the corresponding color");
        jPanel2.add(infotextfield);

        jPanel1.add(jPanel2);

        textfieldpanel.setLayout(new java.awt.GridLayout(1, 0));

        textfieldpanel.add(field1);

        textfieldpanel.add(field2);

        textfieldpanel.add(field3);

        textfieldpanel.add(field4);

        textfieldpanel.add(field5);

        textfieldpanel.add(field6);

        textfieldpanel.add(field7);

        textfieldpanel.add(field8);

        textfieldpanel.add(field9);

        textfieldpanel.add(field10);

        jPanel1.add(textfieldpanel);

        getContentPane().add(jPanel1, java.awt.BorderLayout.NORTH);

        pack();
    }//GEN-END:initComponents

    private void valgradientbuttonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_valgradientbuttonActionPerformed
        // calculate a gradient of the values in the textfields
        String tmpstr = "";
        try {
            tmpstr = field1.getText();
            float worst = Float.parseFloat(tmpstr);
            tmpstr = field10.getText();
            float best = Float.parseFloat(tmpstr);
            float interval = (best - worst) / 9;
            field2.setText(String.valueOf(worst + interval));
            field3.setText(String.valueOf(worst + 2 * interval));
            field4.setText(String.valueOf(worst + 3 * interval));
            field5.setText(String.valueOf(worst + 4 * interval));
            field6.setText(String.valueOf(worst + 5 * interval));
            field7.setText(String.valueOf(worst + 6 * interval));
            field8.setText(String.valueOf(worst + 7 * interval));
            field9.setText(String.valueOf(worst + 8 * interval));
        } catch (NumberFormatException ne) {
            javax.swing.JOptionPane.showMessageDialog(this, "ERROR, unable to parse float from '" + tmpstr + "'");
        }
        repaint();
        parent.repaint();
    }//GEN-LAST:event_valgradientbuttonActionPerformed

    private void gradientbuttonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_gradientbuttonActionPerformed
        // calculate a color gradient
        int worstred = colorarr[0].getRed();
        int worstgreen = colorarr[0].getGreen();
        int worstblue = colorarr[0].getBlue();
        int bestred = colorarr[colornum - 1].getRed();
        int bestgreen = colorarr[colornum - 1].getGreen();
        int bestblue = colorarr[colornum - 1].getBlue();
        //now compute the stepwise gradient
        float redstep = ((float) (bestred - worstred)) / ((float) colornum);
        float greenstep = ((float) (bestgreen - worstgreen)) / ((float) colornum);
        float bluestep = ((float) (bestblue - worstblue)) / ((float) colornum);
        for (int i = 1; i < colornum - 1; i++) {
            colorarr[i] = new java.awt.Color((int) (worstred + (i * redstep)), (int) (worstgreen + (i * greenstep)), (int) (worstblue + (i * bluestep)));
        }
        repaint();
        parent.repaint();
    }//GEN-LAST:event_gradientbuttonActionPerformed

    private void closebuttonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_closebuttonActionPerformed
        // Add your handling code here:
        setVisible(false);
        dispose();
    }//GEN-LAST:event_closebuttonActionPerformed

    private void updatebuttonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_updatebuttonActionPerformed
        // Add your handling code here:
        //read all the data from the textfields
        if (data2d == false) {
            if (parent.data.usescval == false) {//if the data is in P-values
                String tmpstr = field1.getText();
                try {
                    parent.data.colorcutoffs[0] = (float) ((Float.parseFloat(tmpstr) * ln10 - parent.data.p2attoffset) / parent.data.p2attfactor);
                    //System.out.println("cut0="+parent.draw1.colorcutoffs[0]);
                } catch (NumberFormatException ne) {
                    System.err.println("unable to parse float from " + tmpstr);
                    field1.setText(String.valueOf((int) ((((parent.data.colorcutoffs[0]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10)));
                    return;
                }

                tmpstr = field2.getText();
                try {
                    parent.data.colorcutoffs[1] = (float) ((Float.parseFloat(tmpstr) * ln10 - parent.data.p2attoffset) / parent.data.p2attfactor);
                    //System.out.println("cut1="+parent.draw1.colorcutoffs[1]);
                } catch (NumberFormatException ne) {
                    System.err.println("unable to parse float from " + tmpstr);
                    field2.setText(String.valueOf((int) ((((parent.data.colorcutoffs[1]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10)));
                    return;
                }

                tmpstr = field3.getText();
                try {
                    parent.data.colorcutoffs[2] = (float) ((Float.parseFloat(tmpstr) * ln10 - parent.data.p2attoffset) / parent.data.p2attfactor);
                    //System.out.println("cut2="+parent.draw1.colorcutoffs[2]);
                } catch (NumberFormatException ne) {
                    System.err.println("unable to parse float from " + tmpstr);
                    field3.setText(String.valueOf((int) ((((parent.data.colorcutoffs[2]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10)));
                    return;
                }

                tmpstr = field4.getText();
                try {
                    parent.data.colorcutoffs[3] = (float) ((Float.parseFloat(tmpstr) * ln10 - parent.data.p2attoffset) / parent.data.p2attfactor);
                    //System.out.println("cut3="+parent.draw1.colorcutoffs[3]);
                } catch (NumberFormatException ne) {
                    System.err.println("unable to parse float from " + tmpstr);
                    field4.setText(String.valueOf((int) ((((parent.data.colorcutoffs[3]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10)));
                    return;
                }

                tmpstr = field5.getText();
                try {
                    parent.data.colorcutoffs[4] = (float) ((Float.parseFloat(tmpstr) * ln10 - parent.data.p2attoffset) / parent.data.p2attfactor);
                    //System.out.println("cut4="+parent.draw1.colorcutoffs[4]);
                } catch (NumberFormatException ne) {
                    System.err.println("unable to parse float from " + tmpstr);
                    field5.setText(String.valueOf((int) ((((parent.data.colorcutoffs[4]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10)));
                    return;
                }

                tmpstr = field6.getText();
                try {
                    parent.data.colorcutoffs[5] = (float) ((Float.parseFloat(tmpstr) * ln10 - parent.data.p2attoffset) / parent.data.p2attfactor);
                    //System.out.println("cut5="+parent.draw1.colorcutoffs[5]);
                } catch (NumberFormatException ne) {
                    System.err.println("unable to parse float from " + tmpstr);
                    field6.setText(String.valueOf((int) ((((parent.data.colorcutoffs[5]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10)));
                    return;
                }

                tmpstr = field7.getText();
                try {
                    parent.data.colorcutoffs[6] = (float) ((Float.parseFloat(tmpstr) * ln10 - parent.data.p2attoffset) / parent.data.p2attfactor);
                    //System.out.println("cut6="+parent.draw1.colorcutoffs[6]);
                } catch (NumberFormatException ne) {
                    System.err.println("unable to parse float from " + tmpstr);
                    field7.setText(String.valueOf((int) ((((parent.data.colorcutoffs[6]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10)));
                    return;
                }

                tmpstr = field8.getText();
                try {
                    parent.data.colorcutoffs[7] = (float) ((Float.parseFloat(tmpstr) * ln10 - parent.data.p2attoffset) / parent.data.p2attfactor);
                    //System.out.println("cut7="+parent.draw1.colorcutoffs[7]);
                } catch (NumberFormatException ne) {
                    System.err.println("unable to parse float from " + tmpstr);
                    field8.setText(String.valueOf((int) ((((parent.data.colorcutoffs[7]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10)));
                    return;
                }

                tmpstr = field9.getText();
                try {
                    parent.data.colorcutoffs[8] = (float) ((Float.parseFloat(tmpstr) * ln10 - parent.data.p2attoffset) / parent.data.p2attfactor);
                    //System.out.println("cut8="+parent.draw1.colorcutoffs[8]);
                } catch (NumberFormatException ne) {
                    System.err.println("unable to parse float from " + tmpstr);
                    field9.setText(String.valueOf((int) ((((parent.data.colorcutoffs[8]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10)));
                    return;
                }

                tmpstr = field10.getText();
                try {
                    parent.data.colorcutoffs[9] = (float) ((Float.parseFloat(tmpstr) * ln10 - parent.data.p2attoffset) / parent.data.p2attfactor);
                    //System.out.println("cut9="+parent.draw1.colorcutoffs[9]);
                } catch (NumberFormatException ne) {
                    System.err.println("unable to parse float from " + tmpstr);
                    field10.setText(String.valueOf((int) ((((parent.data.colorcutoffs[9]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10)));
                    return;
                }
            } else {//if data is in scval
                String tmpstr = field1.getText();
                try {
                    parent.data.colorcutoffs[0] = (float) ((Float.parseFloat(tmpstr) - parent.data.p2attoffset) / parent.data.p2attfactor);
                    System.out.println("cut0="+parent.data.colorcutoffs[0]);
                } catch (NumberFormatException ne) {
                    System.err.println("unable to parse float from " + tmpstr);
                    field1.setText(String.valueOf((int) ((((parent.data.colorcutoffs[0]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10)));
                    return;
                }

                tmpstr = field2.getText();
                try {
                    parent.data.colorcutoffs[1] = (float) ((Float.parseFloat(tmpstr)- parent.data.p2attoffset) / parent.data.p2attfactor);
                    System.out.println("cut1="+parent.data.colorcutoffs[1]);
                } catch (NumberFormatException ne) {
                    System.err.println("unable to parse float from " + tmpstr);
                    field2.setText(String.valueOf((int) ((((parent.data.colorcutoffs[1]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10)));
                    return;
                }

                tmpstr = field3.getText();
                try {
                    parent.data.colorcutoffs[2] = (float) ((Float.parseFloat(tmpstr) - parent.data.p2attoffset) / parent.data.p2attfactor);
                    System.out.println("cut2="+parent.data.colorcutoffs[2]);
                } catch (NumberFormatException ne) {
                    System.err.println("unable to parse float from " + tmpstr);
                    field3.setText(String.valueOf((int) ((((parent.data.colorcutoffs[2]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10)));
                    return;
                }

                tmpstr = field4.getText();
                try {
                    parent.data.colorcutoffs[3] = (float) ((Float.parseFloat(tmpstr) - parent.data.p2attoffset) / parent.data.p2attfactor);
                    System.out.println("cut3="+parent.data.colorcutoffs[3]);
                } catch (NumberFormatException ne) {
                    System.err.println("unable to parse float from " + tmpstr);
                    field4.setText(String.valueOf((int) ((((parent.data.colorcutoffs[3]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10)));
                    return;
                }

                tmpstr = field5.getText();
                try {
                    parent.data.colorcutoffs[4] = (float) ((Float.parseFloat(tmpstr) - parent.data.p2attoffset) / parent.data.p2attfactor);
                    System.out.println("cut4="+parent.data.colorcutoffs[4]);
                } catch (NumberFormatException ne) {
                    System.err.println("unable to parse float from " + tmpstr);
                    field5.setText(String.valueOf((int) ((((parent.data.colorcutoffs[4]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10)));
                    return;
                }

                tmpstr = field6.getText();
                try {
                    parent.data.colorcutoffs[5] = (float) ((Float.parseFloat(tmpstr) - parent.data.p2attoffset) / parent.data.p2attfactor);
                    System.out.println("cut5="+parent.data.colorcutoffs[5]);
                } catch (NumberFormatException ne) {
                    System.err.println("unable to parse float from " + tmpstr);
                    field6.setText(String.valueOf((int) ((((parent.data.colorcutoffs[5]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10)));
                    return;
                }

                tmpstr = field7.getText();
                try {
                    parent.data.colorcutoffs[6] = (float) ((Float.parseFloat(tmpstr) - parent.data.p2attoffset) / parent.data.p2attfactor);
                    System.out.println("cut6="+parent.data.colorcutoffs[6]);
                } catch (NumberFormatException ne) {
                    System.err.println("unable to parse float from " + tmpstr);
                    field7.setText(String.valueOf((int) ((((parent.data.colorcutoffs[6]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10)));
                    return;
                }

                tmpstr = field8.getText();
                try {
                    parent.data.colorcutoffs[7] = (float) ((Float.parseFloat(tmpstr) - parent.data.p2attoffset) / parent.data.p2attfactor);
                    System.out.println("cut7="+parent.data.colorcutoffs[7]);
                } catch (NumberFormatException ne) {
                    System.err.println("unable to parse float from " + tmpstr);
                    field8.setText(String.valueOf((int) ((((parent.data.colorcutoffs[7]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10)));
                    return;
                }

                tmpstr = field9.getText();
                try {
                    parent.data.colorcutoffs[8] = (float) ((Float.parseFloat(tmpstr) - parent.data.p2attoffset) / parent.data.p2attfactor);
                    System.out.println("cut8="+parent.data.colorcutoffs[8]);
                } catch (NumberFormatException ne) {
                    System.err.println("unable to parse float from " + tmpstr);
                    field9.setText(String.valueOf((int) ((((parent.data.colorcutoffs[8]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10)));
                    return;
                }

                tmpstr = field10.getText();
                try {
                    parent.data.colorcutoffs[9] = (float) ((Float.parseFloat(tmpstr) - parent.data.p2attoffset) / parent.data.p2attfactor);
                    System.out.println("cut9="+parent.data.colorcutoffs[9]);
                } catch (NumberFormatException ne) {
                    System.err.println("unable to parse float from " + tmpstr);
                    field10.setText(String.valueOf((int) ((((parent.data.colorcutoffs[9]) * parent.data.p2attfactor) + parent.data.p2attoffset) / ln10)));
                    return;
                }
            }
        } else {//if parent is in attraction value mode
            String tmpstr = field1.getText();
            try {
                parent.data.colorcutoffs[0] = (float) (Float.parseFloat(tmpstr) / parent.data.p2attfactor);
                //System.out.println("cut0="+parent..colorcutoffs[0]);
            } catch (NumberFormatException ne) {
                System.err.println("unable to parse float from " + tmpstr);
                field1.setText(String.valueOf(parent.data.colorcutoffs[0] * parent.data.p2attfactor));
                return;
            }

            tmpstr = field2.getText();
            try {
                parent.data.colorcutoffs[1] = (float) (Float.parseFloat(tmpstr) / parent.data.p2attfactor);
                //System.out.println("cut1="+parent..colorcutoffs[1]);
            } catch (NumberFormatException ne) {
                System.err.println("unable to parse float from " + tmpstr);
                field2.setText(String.valueOf(parent.data.colorcutoffs[1] * parent.data.p2attfactor));
                return;
            }

            tmpstr = field3.getText();
            try {
                parent.data.colorcutoffs[2] = (float) (Float.parseFloat(tmpstr) / parent.data.p2attfactor);
                //System.out.println("cut2="+parent..colorcutoffs[2]);
            } catch (NumberFormatException ne) {
                System.err.println("unable to parse float from " + tmpstr);
                field3.setText(String.valueOf(parent.data.colorcutoffs[2] * parent.data.p2attfactor));
                return;
            }

            tmpstr = field4.getText();
            try {
                parent.data.colorcutoffs[3] = (float) (Float.parseFloat(tmpstr) / parent.data.p2attfactor);
                //System.out.println("cut3="+parent..colorcutoffs[3]);
            } catch (NumberFormatException ne) {
                System.err.println("unable to parse float from " + tmpstr);
                field4.setText(String.valueOf(parent.data.colorcutoffs[3] * parent.data.p2attfactor));
                return;
            }

            tmpstr = field5.getText();
            try {
                parent.data.colorcutoffs[4] = (float) (Float.parseFloat(tmpstr) / parent.data.p2attfactor);
                //System.out.println("cut4="+parent..colorcutoffs[4]);
            } catch (NumberFormatException ne) {
                System.err.println("unable to parse float from " + tmpstr);
                field5.setText(String.valueOf(parent.data.colorcutoffs[4] * parent.data.p2attfactor));
                return;
            }

            tmpstr = field6.getText();
            try {
                parent.data.colorcutoffs[5] = (float) (Float.parseFloat(tmpstr) / parent.data.p2attfactor);
                //System.out.println("cut5="+parent..colorcutoffs[5]);
            } catch (NumberFormatException ne) {
                System.err.println("unable to parse float from " + tmpstr);
                field6.setText(String.valueOf(parent.data.colorcutoffs[5] * parent.data.p2attfactor));
                return;
            }

            tmpstr = field7.getText();
            try {
                parent.data.colorcutoffs[6] = (float) (Float.parseFloat(tmpstr) / parent.data.p2attfactor);
                //System.out.println("cut6="+parent..colorcutoffs[6]);
            } catch (NumberFormatException ne) {
                System.err.println("unable to parse float from " + tmpstr);
                field7.setText(String.valueOf(parent.data.colorcutoffs[6] * parent.data.p2attfactor));
                return;
            }

            tmpstr = field8.getText();
            try {
                parent.data.colorcutoffs[7] = (float) (Float.parseFloat(tmpstr) / parent.data.p2attfactor);
                //System.out.println("cut7="+parent..colorcutoffs[7]);
            } catch (NumberFormatException ne) {
                System.err.println("unable to parse float from " + tmpstr);
                field8.setText(String.valueOf(parent.data.colorcutoffs[7] * parent.data.p2attfactor));
                return;
            }

            tmpstr = field9.getText();
            try {
                parent.data.colorcutoffs[8] = (float) (Float.parseFloat(tmpstr) / parent.data.p2attfactor);
                //System.out.println("cut8="+parent..colorcutoffs[8]);
            } catch (NumberFormatException ne) {
                System.err.println("unable to parse float from " + tmpstr);
                field9.setText(String.valueOf(parent.data.colorcutoffs[8] * parent.data.p2attfactor));
                return;
            }

            tmpstr = field10.getText();
            try {
                parent.data.colorcutoffs[9] = (float) (Float.parseFloat(tmpstr) / parent.data.p2attfactor);
                //System.out.println("cut9="+parent.draw1.colorcutoffs[9]);
            } catch (NumberFormatException ne) {
                System.err.println("unable to parse float from " + tmpstr);
                field10.setText(String.valueOf(parent.data.colorcutoffs[9] * parent.data.p2attfactor));
                return;
            }
        }

        parent.data.resetDrawOrder();
        parent.repaint();
    }//GEN-LAST:event_updatebuttonActionPerformed

	/**
	 * open and handle color chooser dialog for a clicked element of the color bar
	 * 
	 * @param evt
	 */
	private void mainpanelMouseReleased(java.awt.event.MouseEvent evt) {

		int xval = evt.getX();
		int colorelement = (int) (xval / draw1.elementwidth); // determine clicked element

		colorarr[colorelement] = parent.safe_change_color_dialog("Select New Color", colorarr[colorelement]);
		
		repaint();
		parent.repaint();
	}

    private void mainpanelAncestorResized(java.awt.event.HierarchyEvent evt) {//GEN-FIRST:event_mainpanelAncestorResized
        // Add your handling code here:
        draw1.drawwidth = mainpanel.getWidth();
        draw1.drawheight = mainpanel.getHeight();
        draw1.elementwidth = ((float) draw1.drawwidth) / ((float) colornum);
        repaint();
    }//GEN-LAST:event_mainpanelAncestorResized

    /** Closes the dialog */
    private void closeDialog(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_closeDialog
        setVisible(false);
        dispose();
    }//GEN-LAST:event_closeDialog
    /**
     * @param args the command line arguments
     */
    static java.awt.Color[] colorarr = new java.awt.Color[2];
    drawpanel draw1;
    javax.swing.JColorChooser colorchooser = new javax.swing.JColorChooser();
    ClusteringWithGui parent;
    int colornum;
    static double ln10 = java.lang.Math.log(10);
    boolean data2d = false;//has parent loaded in 2d

    public static void changecolor(ClusteringWithGui parent, java.awt.Color[] incolorarr) {
        //this should replace the colors in colorarr with new ones
        colorarr = incolorarr;
        new DialogChangeConnectionColors(parent, false).setVisible(true);
    }//en getnewcolors
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JLabel bestlabel;
    private javax.swing.JPanel buttonpanel;
    private javax.swing.JButton closebutton;
    private javax.swing.JTextField field1;
    private javax.swing.JTextField field10;
    private javax.swing.JTextField field2;
    private javax.swing.JTextField field3;
    private javax.swing.JTextField field4;
    private javax.swing.JTextField field5;
    private javax.swing.JTextField field6;
    private javax.swing.JTextField field7;
    private javax.swing.JTextField field8;
    private javax.swing.JTextField field9;
    private javax.swing.JButton gradientbutton;
    private javax.swing.JTextField infotextfield;
    private javax.swing.JButton jButton1;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JPanel mainpanel;
    private javax.swing.JPanel textfieldpanel;
    private javax.swing.JButton updatebutton;
    private javax.swing.JButton valgradientbutton;
    private javax.swing.JLabel worstlabel;
    // End of variables declaration//GEN-END:variables

    class drawpanel extends javax.swing.JPanel {

        /**
         * 
         */
        private static final long serialVersionUID = 5838939531675395908L;

        public drawpanel(java.awt.Color[] colorarr) {
            this.colorarr = colorarr;
        }
        java.awt.Color[] colorarr;
        float elementwidth;
        int drawwidth;
        int drawheight;

        public void paintComponent(java.awt.Graphics g) {
            for (int i = 0; i < colornum; i++) {
                g.setColor(colorarr[i]);
                g.fillRect((int) (i * elementwidth), 0, (int) elementwidth, drawheight);
                g.setColor(java.awt.Color.black);
                g.drawRect((int) (i * elementwidth), 0, (int) elementwidth, drawheight);
            }//end for i
        }//end paintcomponent
    }//end class drawpanel
}//end class

