/*
 * clustertest.java
 *
 * Created on December 15, 2002, 3:47 PM
 */
package clans;
import java.util.*;
import javax.swing.*;

import java.awt.event.KeyEvent;
import java.io.*;
/**
 *
 * @author  tancred
 */
public class clustermain_graphics extends javax.swing.JFrame {
    
    /** Creates new form clustertest */
    public clustermain_graphics() {
        initComponents();
    }
    
    public clustermain_graphics(clusterdata data){//minhsp[] blasthits,aaseq[] inaln, String[] namearr,HashMap nameshash,double eval,double pval,float scval,int verbose,int cpu,boolean savepos, String cmd, String blastpath,boolean addblastvbparam, String formatdbpath,String[] referencedb,StringBuffer errbuff,String loadsaved) {
    	draw1=new drawpanel();
        //System.out.println("started clustertest");
        initComponents();
        System.out.println("done init");
        this.data=data;
        data.nographics = false;
        //this.blasthits=blasthits;
        //this.namearr=namearr;
        //this.nameshash=nameshash;
        data.mineval=data.eval;
        data.minpval=data.pval;
        if(data.scval>=0){//in that case use a score cutoff
            data.minpval=data.scval;
            data.usescval=true;
            button_cutoff_value.setText("Use SC-vals better than");
            textfield_cutoff_value.setText("0");
            evalueitem.setText("SC-value plot");
            attvalcompcheckbox.setSelected(false);
        }else{
            data.usescval=false;
        }
        //this.verbose=verbose;
        //this.inaln=rmgaps(inaln);
        //this.cpu=cpu;
        //this.savepos=savepos;
        //this.cmd=cmd;
        //this.blastpath=blastpath;
        //this.addblastvbparam=addblastvbparam;
        //this.formatdbpath=formatdbpath;
        //System.out.println("formatdbpath="+formatdbpath);
        //this.referencedb=referencedb;
        //if(loadsaved!=null){
        //    this.loadsaved=loadsaved;
        //}
        //calculate the ralive sequence lengths
        //int seqnum=java.lang.reflect.Array.getLength(namearr);
        System.out.println("seqnum="+data.seqnum);
        data.seqlengths=new float[data.seqnum];
        float maxlength=0;
        for(int i=0;i<data.seqnum;i++){
            data.seqlengths[i]=data.sequences[i].seq.length();
            if(data.seqlengths[i]>maxlength){
                maxlength=data.seqlengths[i];
            }
        }//end for i
        for(int i=0;i<data.seqnum;i++){
            data.seqlengths[i]/=maxlength;
        }//end for i
        //now all my sequences have a value assigned between 0 and 1 reflecting their length
        //data.movethreads=new getmovethread[cpu];
        textfield_cutoff_value.setText(String.valueOf(data.minpval));
        textfield_info_min_blast_evalue.setText(String.valueOf(data.minpval));
        System.out.println("initializing, please wait");
        //now initialize the stuff
        mousestart[0]=0;
        mousestart[1]=0;
        mousemove[0]=0;
        mousemove[1]=0;
        draw1.init();
        //and now initialize a run
        initgraph();
        /*currcool=1;
        if(myoptionswindow!=null){
            String tmpstr="";
            try{
                tmpstr=myoptionswindow.attfield.getText();
                attfactor=java.lang.Float.parseFloat(tmpstr);
                tmpstr=myoptionswindow.repfield.getText();
                repfactor=java.lang.Float.parseFloat(tmpstr);
                tmpstr=myoptionswindow.dampfield.getText();
                dampening=java.lang.Float.parseFloat(tmpstr);
                tmpstr=myoptionswindow.coolfield.getText();
                cooling=java.lang.Double.parseDouble(tmpstr);
                tmpstr=myoptionswindow.minattfield.getText();
                minattract=java.lang.Double.parseDouble(tmpstr);
                tmpstr=myoptionswindow.maxmovefield.getText();
                maxmove=java.lang.Float.parseFloat(tmpstr);
                tmpstr=myoptionswindow.attvalpowtextfield.getText();
                attvalpow=java.lang.Integer.parseInt(tmpstr);
                tmpstr=myoptionswindow.repvalpowtextfield.getText();
                repvalpow=java.lang.Integer.parseInt(tmpstr);
            }catch(NumberFormatException e){
                System.err.println("ERROR "+tmpstr+" is not a number");
            }
        }
        if(hidebelowold!=hidebelow){
            //draw1.draworder=new Vector[0];
            draw1.draworder=new ArrayList[0];
            hidebelowold=hidebelow;
        }
        mythread.stop=true;
        //System.out.println("starting cluster3d");
        myposarr=cluster3d(blasthits,3);
        draw1.posarr=myposarr;
        if(myoptionswindow!=null){
            myoptionswindow.currcoolfield.setText(String.valueOf(currcool));
        }
        stopbutton.setText("Start run");
        rounds=0;
        mousemove[0]=0;
        mousemove[1]=0;
        */
        button_select_all_or_clear.setText("Select All");
        //System.out.println("done all; displaying now");
        if(new File("positionfile.dat").canRead()){
            System.out.println("reading former data");
            data.myposarr=readsave.readpos();
            data.posarr=data.myposarr;
        }
        System.out.println("done init.");
        if(data.errbuff.length()>0){
            //If I have had errors up to this point
            new errwindow(this,true,data.errbuff.toString()).setVisible(true);
        }
        if(data.input_filename!=null){
            System.out.println("loading data from "+data.input_filename);
            loaddata(data.input_filename);
        }
    }// end init
    
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    private void initComponents() {//GEN-BEGIN:initComponents

        graphpanel = new javax.swing.JPanel();
        buttonpanel = new javax.swing.JPanel();
        drawbuttonpanel = new javax.swing.JPanel();
        button_initialize = new javax.swing.JButton();
        button_start_stop_resume = new javax.swing.JButton();
        button_show_selected = new javax.swing.JButton();
        button_select_move = new javax.swing.JToggleButton();
        button_cutoff_value = new javax.swing.JButton();
        textfield_cutoff_value = new javax.swing.JTextField();
        textfield_info_min_blast_evalue = new javax.swing.JTextField();
        button_select_all_or_clear = new javax.swing.JButton();
        checkbox_show_names = new javax.swing.JCheckBox();
        checkbox_show_numbers = new javax.swing.JCheckBox();
        checkbox_show_connections = new javax.swing.JCheckBox();
        button_zoom_on_selected = new javax.swing.JButton();
        jMenuBar1 = new javax.swing.JMenuBar();
        menu_file = new javax.swing.JMenu();
        loadmenuitem = new javax.swing.JMenuItem();
        savemenuitem = new javax.swing.JMenuItem();
        saveattvalsmenuitem = new javax.swing.JMenuItem();
        addseqsmenuitem = new javax.swing.JMenuItem();
        savemtxmenuitem = new javax.swing.JMenuItem();
        save2dmenuitem = new javax.swing.JMenuItem();
        printmenuitem = new javax.swing.JMenuItem();
        loadalternatemenuitem = new javax.swing.JMenuItem();
        loadtabsmenuitem = new javax.swing.JMenuItem();
        loadgroupsmenuitem = new javax.swing.JMenuItem();
        menu_misc = new javax.swing.JMenu();
        getseqsmenuitem = new javax.swing.JMenuItem();
        hidesingletonsmenuitem = new javax.swing.JMenuItem();
        getchildmenuitem = new javax.swing.JMenuItem();
        getparentmenuitem = new javax.swing.JMenuItem();
        setrotmenuitem = new javax.swing.JMenuItem();
        attvalcompcheckbox = new javax.swing.JCheckBoxMenuItem();
        moveselectedonly = new javax.swing.JCheckBoxMenuItem();
        cluster2dbutton = new javax.swing.JCheckBoxMenuItem();
        rescalepvaluescheckbox = new javax.swing.JCheckBoxMenuItem();
        skipdrawingrounds = new javax.swing.JMenuItem();
        menu_draw = new javax.swing.JMenu();
        changefontmenuitem = new javax.swing.JMenuItem();
        getdotsizemenuitem = new javax.swing.JMenuItem();
        getovalsizemenuitem = new javax.swing.JMenuItem();
        changecolormenuitem = new javax.swing.JMenuItem();
        changefgcolormenuitem = new javax.swing.JMenuItem();
        changebgcolormenuitem = new javax.swing.JMenuItem();
        changeselectcolormenuitem = new javax.swing.JMenuItem();
        changenumbercolor = new javax.swing.JMenuItem();
        changeblastcolor = new javax.swing.JMenuItem();
        lengthcolormenuitem = new javax.swing.JCheckBoxMenuItem();
        colorfrustrationcheckbox = new javax.swing.JCheckBoxMenuItem();
        showorigcheckbox = new javax.swing.JCheckBoxMenuItem();
        showinfocheckbox = new javax.swing.JCheckBoxMenuItem();
        shownamesselectcheckbox = new javax.swing.JCheckBoxMenuItem();
        showblasthitnamescheckbox = new javax.swing.JCheckBoxMenuItem();
        zoommenuitem = new javax.swing.JMenuItem();
        centermenuitem = new javax.swing.JMenuItem();
        antialiasingcheckboxmenuitem = new javax.swing.JCheckBoxMenuItem();
        stereocheckboxmenuitem = new javax.swing.JCheckBoxMenuItem();
        stereoanglemenuitem = new javax.swing.JMenuItem();
        menu_windows = new javax.swing.JMenu();
        showoptionsmenuitem = new javax.swing.JMenuItem();
        sequencesitem = new javax.swing.JMenuItem();
        evalueitem = new javax.swing.JMenuItem();
        getblasthitsmenuitem = new javax.swing.JMenuItem();
        clustermenuitem = new javax.swing.JMenuItem();
        getseqsforselectedhits = new javax.swing.JMenuItem();
        seqscoloring = new javax.swing.JMenuItem();
        showseqsmenuitem = new javax.swing.JMenuItem();
        rotationmenuitem = new javax.swing.JMenuItem();
        affymenuitem = new javax.swing.JMenuItem();
        mapmanmenuitem = new javax.swing.JMenuItem();
        taxonomymenuitem = new javax.swing.JMenuItem();
        menu_help = new javax.swing.JMenu();
        aboutmenuitem = new javax.swing.JMenuItem();
        helpmenuitem = new javax.swing.JMenuItem();

        setTitle("3D-View");
        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosing(java.awt.event.WindowEvent evt) {
                exitForm(evt);
            }
        });

        graphpanel.setPreferredSize(new java.awt.Dimension(640, 480));
        graphpanel.add(draw1);
        graphpanel.addMouseMotionListener(new java.awt.event.MouseMotionAdapter() {
            public void mouseDragged(java.awt.event.MouseEvent evt) {
                graphpanelMouseDragged(evt);
            }
        });
        graphpanel.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mousePressed(java.awt.event.MouseEvent evt) {
                graphpanelMousePressed(evt);
            }
            public void mouseReleased(java.awt.event.MouseEvent evt) {
                graphpanelMouseReleased(evt);
            }
        });
        graphpanel.addHierarchyBoundsListener(new java.awt.event.HierarchyBoundsListener() {
            public void ancestorMoved(java.awt.event.HierarchyEvent evt) {
            }
            public void ancestorResized(java.awt.event.HierarchyEvent evt) {
                graphpanelAncestorResized(evt);
            }
        });
        graphpanel.addMouseWheelListener(new java.awt.event.MouseWheelListener() {
            public void mouseWheelMoved(java.awt.event.MouseWheelEvent evt) {
                graphpanelMouseWheelMoved(evt);
            }
        });
        graphpanel.setLayout(new javax.swing.BoxLayout(graphpanel, javax.swing.BoxLayout.LINE_AXIS));
        getContentPane().add(graphpanel, java.awt.BorderLayout.CENTER);

        buttonpanel.setLayout(new java.awt.GridLayout(1, 0));

        drawbuttonpanel.setLayout(new java.awt.GridLayout(0, 4));

        button_initialize.setText("Initialize");
        button_initialize.setToolTipText("initialize a new run");
        button_initialize.setMnemonic(KeyEvent.VK_I);
        button_initialize.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                button_initializeActionPerformed(evt);
            }
        });
        drawbuttonpanel.add(button_initialize);

        button_start_stop_resume.setText("Stop");
        button_start_stop_resume.setToolTipText("start/resume/stop the current run");
        button_start_stop_resume.setMnemonic(KeyEvent.VK_S);
        button_start_stop_resume.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                button_start_stop_resumeActionPerformed(evt);
            }
        });
        drawbuttonpanel.add(button_start_stop_resume);

        button_show_selected.setText("Show selected");
        button_show_selected.setMnemonic(KeyEvent.VK_O);
        button_show_selected.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                button_show_selectedActionPerformed(evt);
            }
        });
        drawbuttonpanel.add(button_show_selected);

        button_select_move.setText("select/MOVE");
        button_select_move.setToolTipText("Toggle between moving the world and selecting sequences");
        button_select_move.setMnemonic(KeyEvent.VK_V);
        button_select_move.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                button_select_moveActionPerformed(evt);
            }
        });
        drawbuttonpanel.add(button_select_move);

        button_cutoff_value.setText("Use P-values better than:");
        button_cutoff_value.setMnemonic(KeyEvent.VK_B);
        button_cutoff_value.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                button_cutoff_valueActionPerformed(evt);
            }
        });
        drawbuttonpanel.add(button_cutoff_value);

        textfield_cutoff_value.setText("1");
        textfield_cutoff_value.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                textfield_cutoff_valueActionPerformed(evt);
            }
        });
        drawbuttonpanel.add(textfield_cutoff_value);

        textfield_info_min_blast_evalue.setEditable(false);
        textfield_info_min_blast_evalue.setText("1");
        textfield_info_min_blast_evalue.setToolTipText("evalue limit used for blast");
        drawbuttonpanel.add(textfield_info_min_blast_evalue);

        button_select_all_or_clear.setText("Clear Selection");
        button_select_all_or_clear.setMnemonic(KeyEvent.VK_A);
        button_select_all_or_clear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                button_clear_selectionActionPerformed(evt);
            }
        });
        drawbuttonpanel.add(button_select_all_or_clear);

        checkbox_show_names.setText("show names");
        checkbox_show_names.setToolTipText("show sequence names");
        checkbox_show_names.setMnemonic(KeyEvent.VK_N);
        checkbox_show_names.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                checkbox_show_namesItemStateChanged(evt);
            }
        });
        drawbuttonpanel.add(checkbox_show_names);

        checkbox_show_numbers.setText("show numbers");
        checkbox_show_numbers.setMnemonic(KeyEvent.VK_U);
        checkbox_show_numbers.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                checkbox_show_numbersActionPerformed(evt);
            }
        });
        drawbuttonpanel.add(checkbox_show_numbers);

        checkbox_show_connections.setText("show connections");
        checkbox_show_connections.setToolTipText("draw lines for all connections better than the selected cutoff");
        checkbox_show_connections.setMnemonic(KeyEvent.VK_T);
        checkbox_show_connections.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                checkbox_show_connectionsItemStateChanged(evt);
            }
        });
        drawbuttonpanel.add(checkbox_show_connections);

        button_zoom_on_selected.setText("Zoom on selected");
        button_zoom_on_selected.setMnemonic(KeyEvent.VK_Z);
        button_zoom_on_selected.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                button_zoom_on_selectedActionPerformed(evt);
            }
        });
        drawbuttonpanel.add(button_zoom_on_selected);

        buttonpanel.add(drawbuttonpanel);

        getContentPane().add(buttonpanel, java.awt.BorderLayout.SOUTH);

        menu_file.setText("File");
        menu_file.setMnemonic(KeyEvent.VK_F);

        loadmenuitem.setText("Load Run");
        loadmenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                loadmenuitemActionPerformed(evt);
            }
        });
        menu_file.add(loadmenuitem);

        savemenuitem.setText("Save Run");
        savemenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                savemenuitemActionPerformed(evt);
            }
        });
        menu_file.add(savemenuitem);

        saveattvalsmenuitem.setText("Save attraction values to file");
        saveattvalsmenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                saveattvalsmenuitemActionPerformed(evt);
            }
        });
        menu_file.add(saveattvalsmenuitem);

        addseqsmenuitem.setText("Add Sequences");
        addseqsmenuitem.setToolTipText("You have to do that from the command line");
        addseqsmenuitem.setEnabled(false);
        addseqsmenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                addseqsmenuitemActionPerformed(evt);
            }
        });
        menu_file.add(addseqsmenuitem);

        savemtxmenuitem.setText("Save blast matrix pP-values");
        savemtxmenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                savemtxmenuitemActionPerformed(evt);
            }
        });
        menu_file.add(savemtxmenuitem);

        save2dmenuitem.setText("Save 2d graph data");
        save2dmenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                save2dmenuitemActionPerformed(evt);
            }
        });
        menu_file.add(save2dmenuitem);

        printmenuitem.setText("Print view");
        printmenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                printmenuitemActionPerformed(evt);
            }
        });
        menu_file.add(printmenuitem);

        loadalternatemenuitem.setText("Load data in matrix format");
        loadalternatemenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                loadalternatemenuitemActionPerformed(evt);
            }
        });
        menu_file.add(loadalternatemenuitem);

        loadtabsmenuitem.setText("Load tabular data");
        loadtabsmenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                loadtabsmenuitemActionPerformed(evt);
            }
        });
        menu_file.add(loadtabsmenuitem);

        loadgroupsmenuitem.setText("Append sequence groups from file");
        loadgroupsmenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                loadgroupsmenuitemActionPerformed(evt);
            }
        });
        menu_file.add(loadgroupsmenuitem);

        
        jMenuBar1.add(menu_file);

        
        menu_misc.setText("Misc");
        menu_misc.setMnemonic(KeyEvent.VK_M);

        getseqsmenuitem.setText("Extract selected sequences");
        getseqsmenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                getseqsmenuitemActionPerformed(evt);
            }
        });
        menu_misc.add(getseqsmenuitem);

        hidesingletonsmenuitem.setText("Hide singletons");
        hidesingletonsmenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                hidesingletonsmenuitemActionPerformed(evt);
            }
        });
        menu_misc.add(hidesingletonsmenuitem);

        getchildmenuitem.setText("Use selected subset");
        getchildmenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                getchildmenuitemActionPerformed(evt);
            }
        });
        menu_misc.add(getchildmenuitem);

        getparentmenuitem.setText("Use parent group (0)");
        getparentmenuitem.setEnabled(false);
        getparentmenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                getparentmenuitemActionPerformed(evt);
            }
        });
        menu_misc.add(getparentmenuitem);

        setrotmenuitem.setText("Set rotation values");
        setrotmenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                setrotmenuitemActionPerformed(evt);
            }
        });
        menu_misc.add(setrotmenuitem);

        attvalcompcheckbox.setSelected(true);
        attvalcompcheckbox.setText("Complex attraction");
        attvalcompcheckbox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                attvalcompcheckboxActionPerformed(evt);
            }
        });
        menu_misc.add(attvalcompcheckbox);

        moveselectedonly.setText("Optimize only selected sequences");
        moveselectedonly.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                moveselectedonlyActionPerformed(evt);
            }
        });
        menu_misc.add(moveselectedonly);

        cluster2dbutton.setText("Cluster in 2D");
        cluster2dbutton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                cluster2dbuttonActionPerformed(evt);
            }
        });
        menu_misc.add(cluster2dbutton);

        rescalepvaluescheckbox.setText("Rescale attraction values");
        rescalepvaluescheckbox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                rescalepvaluescheckboxActionPerformed(evt);
            }
        });
        menu_misc.add(rescalepvaluescheckbox);

        skipdrawingrounds.setText("Only draw every Nth round (speedup)");
        skipdrawingrounds.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                skipdrawingroundsActionPerformed(evt);
            }
        });
        menu_misc.add(skipdrawingrounds);

        jMenuBar1.add(menu_misc);

        
        menu_draw.setText("Draw");
        menu_draw.setMnemonic(KeyEvent.VK_D);

        changefontmenuitem.setText("Change Font");
        changefontmenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                changefontmenuitemActionPerformed(evt);
            }
        });
        menu_draw.add(changefontmenuitem);

        getdotsizemenuitem.setText("Set dot size");
        getdotsizemenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                getdotsizemenuitemActionPerformed(evt);
            }
        });
        menu_draw.add(getdotsizemenuitem);

        getovalsizemenuitem.setText("Set selected circle size");
        getovalsizemenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                getovalsizemenuitemActionPerformed(evt);
            }
        });
        menu_draw.add(getovalsizemenuitem);

        changecolormenuitem.setText("Change color (dot connections)");
        changecolormenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                changecolormenuitemActionPerformed(evt);
            }
        });
        menu_draw.add(changecolormenuitem);

        changefgcolormenuitem.setText("Change color (Foreground)");
        changefgcolormenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                changefgcolormenuitemActionPerformed(evt);
            }
        });
        menu_draw.add(changefgcolormenuitem);

        changebgcolormenuitem.setText("Change color (Background)");
        changebgcolormenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                changebgcolormenuitemActionPerformed(evt);
            }
        });
        menu_draw.add(changebgcolormenuitem);

        changeselectcolormenuitem.setText("Change color (Selecteds)");
        changeselectcolormenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                changeselectcolormenuitemActionPerformed(evt);
            }
        });
        menu_draw.add(changeselectcolormenuitem);

        changenumbercolor.setText("Change color (BLAST hit numbers)");
        changenumbercolor.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                changenumbercolorActionPerformed(evt);
            }
        });
        menu_draw.add(changenumbercolor);

        changeblastcolor.setText("Change color (BLAST hit circles)");
        changeblastcolor.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                changeblastcolorActionPerformed(evt);
            }
        });
        menu_draw.add(changeblastcolor);

        lengthcolormenuitem.setText("Color dots by sequence length (yellow=short, blue=long)");
        menu_draw.add(lengthcolormenuitem);

        colorfrustrationcheckbox.setText("Color by edge \"frustration\" (red=too long, blue=too short)");
        colorfrustrationcheckbox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                colorfrustrationcheckboxActionPerformed(evt);
            }
        });
        menu_draw.add(colorfrustrationcheckbox);

        showorigcheckbox.setText("Show origin");
        menu_draw.add(showorigcheckbox);

        showinfocheckbox.setSelected(true);
        showinfocheckbox.setText("Show info");
        showinfocheckbox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                showinfocheckboxActionPerformed(evt);
            }
        });
        menu_draw.add(showinfocheckbox);

        shownamesselectcheckbox.setText("Show names while selecting");
        menu_draw.add(shownamesselectcheckbox);

        showblasthitnamescheckbox.setText("Show hsp sequence numbers");
        menu_draw.add(showblasthitnamescheckbox);

        zoommenuitem.setText("Zoom");
        zoommenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                zoommenuitemActionPerformed(evt);
            }
        });
        menu_draw.add(zoommenuitem);

        centermenuitem.setText("Center graph");
        centermenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                centermenuitemActionPerformed(evt);
            }
        });
        menu_draw.add(centermenuitem);

        antialiasingcheckboxmenuitem.setText("Antialiasing (slow !)");
        antialiasingcheckboxmenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                antialiasingcheckboxmenuitemActionPerformed(evt);
            }
        });
        menu_draw.add(antialiasingcheckboxmenuitem);

        stereocheckboxmenuitem.setText("Stereo");
        stereocheckboxmenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                stereocheckboxmenuitemActionPerformed(evt);
            }
        });
        menu_draw.add(stereocheckboxmenuitem);

        stereoanglemenuitem.setText("Change stereo angle (0-360)");
        stereoanglemenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                stereoanglemenuitemActionPerformed(evt);
            }
        });
        menu_draw.add(stereoanglemenuitem);

        jMenuBar1.add(menu_draw);

        
        menu_windows.setText("Windows");
        menu_windows.setMnemonic(KeyEvent.VK_W);

        showoptionsmenuitem.setText("Show options window");
        showoptionsmenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                showoptionsmenuitemActionPerformed(evt);
            }
        });
        menu_windows.add(showoptionsmenuitem);

        sequencesitem.setText("Selecteds");
        sequencesitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                sequencesitemActionPerformed(evt);
            }
        });
        menu_windows.add(sequencesitem);

        evalueitem.setText("P-value plot");
        evalueitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                evalueitemActionPerformed(evt);
            }
        });
        menu_windows.add(evalueitem);

        getblasthitsmenuitem.setText("Show blast hits for sequence:");
        getblasthitsmenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                getblasthitsmenuitemActionPerformed(evt);
            }
        });
        menu_windows.add(getblasthitsmenuitem);

        clustermenuitem.setText("find clusters");
        clustermenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                clustermenuitemActionPerformed(evt);
            }
        });
        menu_windows.add(clustermenuitem);

        getseqsforselectedhits.setText("Get sequence with hits from/to selected");
        getseqsforselectedhits.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                getseqsforselectedhitsActionPerformed(evt);
            }
        });
        menu_windows.add(getseqsforselectedhits);

        seqscoloring.setText("Edit Groups");
        seqscoloring.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                seqscoloringActionPerformed(evt);
            }
        });
        menu_windows.add(seqscoloring);

        showseqsmenuitem.setText("Show selected sequences as text (copy/pastable)");
        showseqsmenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                showseqsmenuitemActionPerformed(evt);
            }
        });
        menu_windows.add(showseqsmenuitem);

        rotationmenuitem.setText("Rotation");
        rotationmenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                rotationmenuitemActionPerformed(evt);
            }
        });
        menu_windows.add(rotationmenuitem);

        affymenuitem.setText("Microarray_data");
        affymenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                affymenuitemActionPerformed(evt);
            }
        });
        menu_windows.add(affymenuitem);

        mapmanmenuitem.setText("Functional mapping");
        mapmanmenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                mapmanmenuitemActionPerformed(evt);
            }
        });
        menu_windows.add(mapmanmenuitem);

        taxonomymenuitem.setText("Taxonomy");
        taxonomymenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                taxonomymenuitemActionPerformed(evt);
            }
        });
        menu_windows.add(taxonomymenuitem);

        jMenuBar1.add(menu_windows);

        
        menu_help.setText("Help");
        menu_help.setMnemonic(KeyEvent.VK_H);

        aboutmenuitem.setText("About");
        aboutmenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                aboutmenuitemActionPerformed(evt);
            }
        });
        menu_help.add(aboutmenuitem);

        helpmenuitem.setText("Help");
        helpmenuitem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                helpmenuitemActionPerformed(evt);
            }
        });
        menu_help.add(helpmenuitem);

        jMenuBar1.add(menu_help);

        setJMenuBar(jMenuBar1);

        pack();
    }//GEN-END:initComponents
    
    private void loadtabsmenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_loadtabsmenuitemActionPerformed
        if(checkstop()==false){//check for stop of run
            return;
        }
        groupseqs=null;
        int returnVal = fc.showOpenDialog(this);
        if(returnVal == JFileChooser.APPROVE_OPTION) {
            //update loadtabdata(fc.getSelectedFile().getAbsolutePath());
            if(myseqgroupwindow!=null){
                myseqgroupwindow.setVisible(false);
                myseqgroupwindow.dispose();
            }
            if(mymapfunctiondialog!=null){
                mymapfunctiondialog.setVisible(false);
                mymapfunctiondialog.dispose();
            }
            repaint();
        }
    }//GEN-LAST:event_loadtabsmenuitemActionPerformed
    
    private void stereocheckboxmenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_stereocheckboxmenuitemActionPerformed
        repaint();
    }//GEN-LAST:event_stereocheckboxmenuitemActionPerformed
    
    private void stereoanglemenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_stereoanglemenuitemActionPerformed
        // change the angle of stereo vision
        String tmpstr="";
        try{
            tmpstr=JOptionPane.showInputDialog(this,"Enter the new angle (int):",String.valueOf(draw1.stereoangle));
            if(tmpstr!=null){
                draw1.stereoangle=(int)(Float.parseFloat(tmpstr));//someone might enter a float and it's not time critical
            }
        }catch (NumberFormatException ne){
            javax.swing.JOptionPane.showMessageDialog(this,"ERROR, unable to parse integer from '"+tmpstr+"'");
        }
        repaint();
    }//GEN-LAST:event_stereoanglemenuitemActionPerformed
    
    private void mapmanmenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_mapmanmenuitemActionPerformed
        if(mymapfunctiondialog!=null){
            mymapfunctiondialog.setVisible(false);
            mymapfunctiondialog.dispose();
        }
        mymapfunctiondialog=new mapfunctiondialog_tab(this);
        mymapfunctiondialog.setVisible(true);
    }//GEN-LAST:event_mapmanmenuitemActionPerformed
    
    private void graphpanelMouseWheelMoved(java.awt.event.MouseWheelEvent evt) {//GEN-FIRST:event_graphpanelMouseWheelMoved
        //add zooming ability
        if(button_select_move.isSelected()==false){
            if(evt.isShiftDown()){
                if(evt.isControlDown()){
                    data.zoomfactor+=((float)evt.getWheelRotation())/10;
                }else{
                    data.zoomfactor+=((float)evt.getWheelRotation())/100;
                }
                repaint();
            }
        }
    }//GEN-LAST:event_graphpanelMouseWheelMoved
    
    private void affymenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_affymenuitemActionPerformed
        if(myaffydialog!=null){
            myaffydialog.setVisible(false);
            myaffydialog.dispose();
        }
        myaffydialog=new affydialog(this);
        myaffydialog.setVisible(true);
    }//GEN-LAST:event_affymenuitemActionPerformed
    
    private void rotationmenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_rotationmenuitemActionPerformed
        // open a window making specific rotating possible
        if(myrotationdialog!=null){
            myrotationdialog.setVisible(false);
            myrotationdialog.dispose();
        }
        myrotationdialog=new rotationdialog(this);
        myrotationdialog.setVisible(true);
    }//GEN-LAST:event_rotationmenuitemActionPerformed
    
    private void loadgroupsmenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_loadgroupsmenuitemActionPerformed
        // load additional sequence groups from a file
        int returnVal = fc.showOpenDialog(this);
        if(returnVal == JFileChooser.APPROVE_OPTION) {
            data.seqgroupsvec=customutils.loadgroups(data.seqgroupsvec,data.namearr,fc.getSelectedFile());
            if(myseqgroupwindow!=null){
                myseqgroupwindow.setVisible(false);
                myseqgroupwindow.dispose();
            }
            repaint();
        }
    }//GEN-LAST:event_loadgroupsmenuitemActionPerformed
    
    private void addseqsmenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_addseqsmenuitemActionPerformed
        javax.swing.JOptionPane.showMessageDialog(this,"You currently have to do that from the command line");
    }//GEN-LAST:event_addseqsmenuitemActionPerformed
    
    private void showoptionsmenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_showoptionsmenuitemActionPerformed
        // show the options window
        if(myoptionswindow!=null){
            myoptionswindow.setVisible(false);
            myoptionswindow.dispose();
        }
        myoptionswindow=new optionsdialog(this);
        myoptionswindow.setVisible(true);
    }//GEN-LAST:event_showoptionsmenuitemActionPerformed
    
    private void skipdrawingroundsActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_skipdrawingroundsActionPerformed
        // only draw every Nth round
        String tmpstr;
        tmpstr=javax.swing.JOptionPane.showInputDialog(this,"Draw each Nth round. N=",String.valueOf(skiprounds));
        try{
            if(tmpstr!=null){
                skiprounds=Integer.parseInt(tmpstr);
            }
        }catch (NumberFormatException ne){
            javax.swing.JOptionPane.showMessageDialog(this, "ERROR, unable to parse integer from '"+tmpstr+"'");
        }
    }//GEN-LAST:event_skipdrawingroundsActionPerformed
    
    private void colorfrustrationcheckboxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_colorfrustrationcheckboxActionPerformed
        repaint();//rest is done in paintComponent function
    }//GEN-LAST:event_colorfrustrationcheckboxActionPerformed
    
    private void helpmenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_helpmenuitemActionPerformed
        // "help" was clicked
        new helpwindow(this,true).setVisible(true);
    }//GEN-LAST:event_helpmenuitemActionPerformed
    
    private void aboutmenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_aboutmenuitemActionPerformed
        // the about clans menu was selected
        new aboutwindow(this,true).setVisible(true);
    }//GEN-LAST:event_aboutmenuitemActionPerformed
    
    private void antialiasingcheckboxmenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_antialiasingcheckboxmenuitemActionPerformed
        repaint();
    }//GEN-LAST:event_antialiasingcheckboxmenuitemActionPerformed
    
    private void showseqsmenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_showseqsmenuitemActionPerformed
        // open a text field with the selected sequences
        int seqnum=java.lang.reflect.Array.getLength(data.selectednames);
        if(seqnum<1){
            javax.swing.JOptionPane.showMessageDialog(this,"Please select some sequences");
        }else{
            StringBuffer outbuff=new StringBuffer();
            for(int i=0;i<seqnum;i++){
                outbuff.append(">"+data.namearr[data.selectednames[i]]+" "+data.selectednames[i]+"\n");
                outbuff.append(data.sequences[data.selectednames[i]].seq+"\n");
            }//end for i
            new showsequences(new javax.swing.JFrame(),outbuff).setVisible(true);
        }
        
    }//GEN-LAST:event_showseqsmenuitemActionPerformed
    
    private void centermenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_centermenuitemActionPerformed
        draw1.xtranslate=(int)(((graphpanel.getWidth()-2*draw1.xadd)-draw1.drawwidth)/2);
        draw1.ytranslate=(int)(((graphpanel.getHeight()-2*draw1.yadd)-draw1.drawheight)/2);
        repaint();
    }//GEN-LAST:event_centermenuitemActionPerformed
    
    private void zoommenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_zoommenuitemActionPerformed
        String tmpstr="";
        float oldzoom=data.zoomfactor;
        try{
            tmpstr=javax.swing.JOptionPane.showInputDialog(this,"New zoom factor (in percent)",String.valueOf((int)(data.zoomfactor*100)));
            if(tmpstr!=null){
                data.zoomfactor=((float)(Integer.parseInt(tmpstr)))/100;
            }
        }catch (NumberFormatException ne){
            javax.swing.JOptionPane.showMessageDialog(this,"Error; cannot parse float from '"+tmpstr+"'");
            //System.err.println("ERROR parsing number from "+tmpstr);
            return;
        }
        //now keep the center of the screen in focus
        int panelwidth=graphpanel.getWidth()-2*draw1.xadd;
        int panelheight=graphpanel.getHeight()-2*draw1.yadd;
        int imagecenterx=(int)((panelwidth/2-draw1.xtranslate)/oldzoom);
        int imagecentery=(int)((panelheight/2-draw1.ytranslate)/oldzoom);
        draw1.xtranslate=(int)(-(imagecenterx*data.zoomfactor)+panelwidth/2);
        draw1.ytranslate=(int)(-(imagecentery*data.zoomfactor)+panelheight/2);
        repaint();
    }//GEN-LAST:event_zoommenuitemActionPerformed
    
    private void changefontmenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_changefontmenuitemActionPerformed
        // change the current font
        draw1.myfont=fontchooserdialog.getfont("Select Font",draw1.myfont);
        repaint();
    }//GEN-LAST:event_changefontmenuitemActionPerformed
    
    private void seqscoloringActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_seqscoloringActionPerformed
        // open a window showing the vector containing the selected groups of sequences
        if(myseqgroupwindow!=null){
            myseqgroupwindow.setVisible(false);
            myseqgroupwindow.dispose();
        }
        myseqgroupwindow=new seqgroupwindow(this);
        myseqgroupwindow.setVisible(true);
    }//GEN-LAST:event_seqscoloringActionPerformed
    
    private void getseqsforselectedhitsActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_getseqsforselectedhitsActionPerformed
        // get the sequences with hits to or from the selected set of sequences
        // output in a "sequences" window the currently selected and in a second frame the
        // sequences with hits to the selected.
        int[] blasthitsarr=showblasthitsforselected.getblasthits(data.myattvals,data.selectednames,data.namearr);//get the blast hits
        new showblasthitsforselected(this,blasthitsarr,data.selectednames).setVisible(true);
    }//GEN-LAST:event_getseqsforselectedhitsActionPerformed
    
    private void clustermenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_clustermenuitemActionPerformed
        // detect clusters in this dataset
        //use minpval and the current attraction values to get the clustering
        //don't define optionsvec as it contains both strings and numbers in a defined order
        Vector optionsvec=new Vector();
        new clusteroptionsdialog(this,optionsvec).setVisible(true);
        if(optionsvec.size()==0){//if I canceled
            return;
        }
        boolean didbootstrap=false;
        String tmpstr=(String)optionsvec.remove(0);
        if(tmpstr.equals("convex")){
            tmpstr=(String)optionsvec.remove(0);
            try{
                float sigmafac=Float.parseFloat(tmpstr);
                tmpstr=(String)optionsvec.remove(0);
                int minseqnum=Integer.parseInt(tmpstr);
                System.out.println("searching for convex clusters");
                Vector <cluster>clustervec=findclusters.getconvex(data.myattvals,sigmafac,minseqnum,data.elements);
                System.out.println("done searching for clusters; opening window");
                if(((String)optionsvec.remove(0)).equalsIgnoreCase("true")){//if do bootstrap
                    didbootstrap=true;
                    tmpstr=(String)optionsvec.remove(0);
                    int replicates=Integer.parseInt(tmpstr);
                    tmpstr=(String)optionsvec.remove(0);
                    float remove=Float.parseFloat(tmpstr);
                    if(remove>1){
                        remove/=100;
                    }
                    if(bootstrapcluster.bootstrapconvex(data.myattvals,clustervec,"convex",replicates,remove,sigmafac,minseqnum,data.elements)==false){
                        javax.swing.JOptionPane.showMessageDialog(this,"ERROR while bootstrapping");
                        return;
                    }
                }
                new clusterwindow(this,clustervec,"convex: "+minseqnum+";"+sigmafac,didbootstrap).setVisible(true);
            }catch (NumberFormatException ne){
                javax.swing.JOptionPane.showMessageDialog(this,"Unable to parse float from "+tmpstr);
            }
        }else if(tmpstr.equals("linkage")){
            tmpstr=(String)optionsvec.remove(0);
            try{
                int minlinkage=Integer.parseInt(tmpstr);
                tmpstr=(String)optionsvec.remove(0);
                int minseqnum=Integer.parseInt(tmpstr);
                System.out.println("searching for linkage clusters");
                Vector <cluster>clustervec=findclusters.getlinkage(data.myattvals,minlinkage,minseqnum,data.elements);
                System.out.println("done searching for clusters; opening window");
                if(((String)optionsvec.remove(0)).equalsIgnoreCase("true")){//if do bootstrap
                    didbootstrap=true;
                    tmpstr=(String)optionsvec.remove(0);
                    int replicates=Integer.parseInt(tmpstr);
                    tmpstr=(String)optionsvec.remove(0);
                    float remove=Float.parseFloat(tmpstr);
                    if(remove>1){
                        remove/=100;
                    }
                    if(bootstrapcluster.bootstraplinkage(data.myattvals,clustervec,"linkage",replicates,remove,minlinkage,minseqnum,data.elements)==false){
                        javax.swing.JOptionPane.showMessageDialog(this,"Error while bootstrapping");
                        return;
                    }
                }
                new clusterwindow(this,clustervec,"linkage: "+minseqnum+";"+minlinkage,didbootstrap).setVisible(true);
            }catch (NumberFormatException ne){
                javax.swing.JOptionPane.showMessageDialog(this,"Unable to parse int from "+tmpstr);
            }
        }else if(tmpstr.equals("network")){
            tmpstr=(String)optionsvec.remove(0);
            try{
                int minseqnum=Integer.parseInt(tmpstr);
                boolean dooffset=false;
                boolean globalaverage=false;
                tmpstr=(String)optionsvec.remove(0);
                if(tmpstr.equalsIgnoreCase("true")){
                    dooffset=true;
                }
                tmpstr=(String)optionsvec.remove(0);
                if(tmpstr.equalsIgnoreCase("true")){
                    globalaverage=true;
                }
                int maxrounds=((Integer)optionsvec.remove(optionsvec.size()-1)).intValue();
                System.out.println("searching for network clusters, maxrounds="+maxrounds);
                Vector <cluster>clustervec=findclusters.getnetwork(data.myattvals,minseqnum,dooffset,globalaverage,data.elements,maxrounds);
                System.out.println("done searching for clusters; opening window");
                if(((String)optionsvec.remove(0)).equalsIgnoreCase("true")){//if do bootstrap
                    didbootstrap=true;
                    tmpstr=(String)optionsvec.remove(0);
                    int replicates=Integer.parseInt(tmpstr);
                    tmpstr=(String)optionsvec.remove(0);
                    float remove=Float.parseFloat(tmpstr);
                    if(remove>1){
                        remove/=100;
                    }
                    if(bootstrapcluster.bootstrapnetwork(data.myattvals,clustervec,"network",replicates,remove,minseqnum,dooffset,globalaverage,data.elements,maxrounds)==false){
                        javax.swing.JOptionPane.showMessageDialog(this,"Error while bootstrapping");
                        return;
                    }
                }
                new clusterwindow(this,clustervec,"network:"+minseqnum+";"+dooffset+";"+globalaverage,didbootstrap).setVisible(true);
            }catch(NumberFormatException ne){
                javax.swing.JOptionPane.showMessageDialog(this,"Unable to parse int from "+tmpstr);
            }
        }else{
            javax.swing.JOptionPane.showMessageDialog(this,"Error in selecting clustering method: "+tmpstr);
        }
    }//GEN-LAST:event_clustermenuitemActionPerformed
    
    private void loadalternatemenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_loadalternatemenuitemActionPerformed
        // load data from a file with matrix info (oly one value per pair)
        if(checkstop()==false){//check for stop of run
            return;
        }
        groupseqs=null;
        int returnVal = fc.showOpenDialog(this);
        if(returnVal == JFileChooser.APPROVE_OPTION) {
            saverunobject saveddata=customutils.loadrunalternate(fc.getSelectedFile());
            if(saveddata.file!=null){//if the data was read all right
                repaint="Error Loading Data";
                data.sequences=ClusterMethods.remove_gaps_from_sequences(saveddata.inaln);
                int seqs=java.lang.reflect.Array.getLength(data.sequences);
                data.nameshash=new HashMap((int)(seqs/0.75)+1,(float)0.75);//holds info about which name is which array number
                data.namearr=new String[seqs];
                data.myposarr=new float[seqs][3];
                Random rand=ClusterMethods.rand;
                for(int i=0;i<seqs;i++){
                    data.namearr[i]=data.sequences[i].name.trim();
                    data.sequences[i].name=new String("sequence"+i);
                    data.nameshash.put(data.sequences[i].name,new Integer(i));
                    data.myposarr[i][0]=rand.nextFloat();
                    data.myposarr[i][1]=rand.nextFloat();
                    data.myposarr[i][2]=rand.nextFloat();
                }
                data.blasthits=null;
                data.orgattvals=null;
                data.myattvals=saveddata.attvals;
                data.myattvals=saveddata.attvals;
                button_cutoff_value.setText("Use Attraction values better than");
                textfield_cutoff_value.setText("0");
                data.elements=java.lang.reflect.Array.getLength(data.namearr);
                //now symmetrize and normalize the attvals to range from -1 to +1
                float minval=0;
                float maxval=0;
                for(int i=java.lang.reflect.Array.getLength(data.myattvals)-1;i>=0;i--){
                    if(data.myattvals[i].att>maxval){
                        maxval=data.myattvals[i].att;
                    }else if(data.myattvals[i].att<minval){
                        minval=data.myattvals[i].att;
                    }
                }//end for i
                if(-minval>maxval){//decide wether to divide by maxval or -minval
                    maxval=-minval;
                }
                for(int i=java.lang.reflect.Array.getLength(data.myattvals)-1;i>=0;i--){
                    //now normalize the values
                    data.myattvals[i].att/=maxval;
                }//end for i
                data.p2attfactor=maxval;
                data.selectednames=new int[0];
                data.seqgroupsvec=saveddata.seqgroupsvec;
                data.posarr=data.myposarr;
                data.lastmovearr=new float[data.elements][data.dimensions];
                data.mymovearr=new float[data.elements][data.dimensions];
                data.posarrtmp=new float[data.elements][data.dimensions];
                data.drawarrtmp=new int[data.elements][data.dimensions];
                data.draworder=new ArrayList[0];
                data.attvalsimple=true;
                repaint=null;
                data.minpval=1;
                textfield_cutoff_value.setText("1");
                textfield_info_min_blast_evalue.setText("1");
            }else{//if the data had errors
                JOptionPane.showMessageDialog(this,"Error reading data","Error reading",JOptionPane.ERROR_MESSAGE);
                return;
            }
        }
        int seqnum=java.lang.reflect.Array.getLength(data.namearr);
        System.out.println("seqnum="+seqnum);
        data.seqlengths=new float[seqnum];
        float maxlength=0;
        for(int i=0;i<seqnum;i++){
            data.seqlengths[i]=data.sequences[i].seq.replaceAll("-","").length();
            if(data.seqlengths[i]>maxlength){
                maxlength=data.seqlengths[i];
            }
        }//end for i
        if(maxlength>0){
            for(int i=0;i<seqnum;i++){
                data.seqlengths[i]/=maxlength;
            }//end for i
        }
        //now all my sequences have a value assigned between 0 and 1 reflecting their length
        repaint();
        
    }//GEN-LAST:event_loadalternatemenuitemActionPerformed
    
    private void cluster2dbuttonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cluster2dbuttonActionPerformed
        // get rid of all z-axis information and set the rotation matrices to 1,1,1 (x,y,z)
        //set the rotation matrices to 1,1,1
        data.rotmtx[0][0]=1;
        data.rotmtx[0][1]=0;
        data.rotmtx[0][2]=0;
        data.rotmtx[1][0]=0;
        data.rotmtx[1][1]=1;
        data.rotmtx[1][2]=0;
        data.rotmtx[2][0]=0;
        data.rotmtx[2][1]=0;
        data.rotmtx[2][2]=1;
        draw1.tmprotmtx[0][0]=1;
        draw1.tmprotmtx[0][1]=0;
        draw1.tmprotmtx[0][2]=0;
        draw1.tmprotmtx[1][0]=0;
        draw1.tmprotmtx[1][1]=1;
        draw1.tmprotmtx[1][2]=0;
        draw1.tmprotmtx[2][0]=0;
        draw1.tmprotmtx[2][1]=0;
        draw1.tmprotmtx[2][2]=1;
        data.myrotmtx[0][0]=1;
        data.myrotmtx[0][1]=0;
        data.myrotmtx[0][2]=0;
        data.myrotmtx[1][0]=0;
        data.myrotmtx[1][1]=1;
        data.myrotmtx[1][2]=0;
        data.myrotmtx[2][0]=0;
        data.myrotmtx[2][1]=0;
        data.myrotmtx[2][2]=1;
        //now set all the z-values to "0"
        for(int i=java.lang.reflect.Array.getLength(data.myposarr);--i>=0;){
            data.myposarr[i][2]=0;
        }//end for i
        if(cluster2dbutton.isSelected()){
            data.cluster2d=true;
        }else{
            data.cluster2d=false;
        }
        repaint();
    }//GEN-LAST:event_cluster2dbuttonActionPerformed
    
    private void printmenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_printmenuitemActionPerformed
        // print the current view
        if(checkstop()==false){//check for stop of run
            return;
        }
        java.awt.print.PrinterJob printJob = java.awt.print.PrinterJob.getPrinterJob();
        printJob.setPrintable(draw1);
        if (printJob.printDialog()) {
            try{
                printJob.print();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
    }//GEN-LAST:event_printmenuitemActionPerformed
    
    private void changeselectcolormenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_changeselectcolormenuitemActionPerformed
        // change color for the selected sequence circles
        draw1.selectedcolor=JColorChooser.showDialog(this,"Choose a new foreground color", draw1.selectedcolor);
        repaint();
    }//GEN-LAST:event_changeselectcolormenuitemActionPerformed
    
    private void changeblastcolorActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_changeblastcolorActionPerformed
        // change the color of the blast hit circles
        draw1.blastcirclecolor=JColorChooser.showDialog(this,"Choose a new foreground color", draw1.blastcirclecolor);
        repaint();
    }//GEN-LAST:event_changeblastcolorActionPerformed
    
    private void button_show_selectedActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_showselectbuttonActionPerformed
        // open the sequences dialog and display the names for only the selected sequences
    	if (!this.contains_data(true)) {
			return;
		}
    	if(shownames!=null){
            shownames.setVisible(false);
            shownames.dispose();
        }
        shownames=new shownamedialog(data.namearr,this);
        shownames.setVisible(true);
    }//GEN-LAST:event_showselectbuttonActionPerformed
    
    private void changefgcolormenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_changefgcolormenuitemActionPerformed
        // change the foreground color
        draw1.fgcolor=JColorChooser.showDialog(this,"Choose a new foreground color", draw1.fgcolor);
        repaint();
    }//GEN-LAST:event_changefgcolormenuitemActionPerformed
    
    private void changebgcolormenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_changebgcolormenuitemActionPerformed
        // change the background color
        draw1.bgcolor=JColorChooser.showDialog(this,"Choose a new background color", draw1.bgcolor);
        repaint();
    }//GEN-LAST:event_changebgcolormenuitemActionPerformed
    
    private void save2dmenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_save2dmenuitemActionPerformed
        // save the 2d coordinates of the graph points to a file
        if(checkstop()==false){//check for stop of run
            return;
        }
        int returnVal = fc.showSaveDialog(this);
        if(returnVal == JFileChooser.APPROVE_OPTION) {
            File savefile=fc.getSelectedFile();
            try{
                PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter(savefile)));
                //first write the sequence names
                int namenum=java.lang.reflect.Array.getLength(data.namearr);
                outwrite.println("ID\tNAME\tX\tY");
                for(int i=0;i<namenum;i++){
                    outwrite.println(i+"\t"+data.namearr[i]+"\t"+(data.posarrtmp[i][0]-draw1.xadd)/draw1.drawwidth+"\t"+(data.posarrtmp[i][1]-draw1.yadd)/draw1.drawheight);
                }//end for i names
                outwrite.close();
            }catch (IOException ioe){
                javax.swing.JOptionPane.showMessageDialog(this,"IOERROR writing to '"+savefile.getName()+"'");
            }
        }
    }//GEN-LAST:event_save2dmenuitemActionPerformed
    
    private void changenumbercolorActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_changenumbercolorActionPerformed
        try{
            draw1.blasthitcolor=JColorChooser.showDialog(this, "Select New Color",draw1.blasthitcolor);
        }catch (java.awt.HeadlessException e){
            javax.swing.JOptionPane.showMessageDialog(this,"HeadlessException");
            System.err.println("HeadlessException!");
        }
    }//GEN-LAST:event_changenumbercolorActionPerformed
    
    private void button_cutoff_valueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_minpvalbuttonActionPerformed
        // same as pressing return in the corresponding text field
    	if (!this.contains_data(true)) {
			return;
		}
    	if(checkstop()==false){//check for stop of run
            return;
        }
        try{
            data.minpval=java.lang.Double.parseDouble(textfield_cutoff_value.getText());
        }catch (NumberFormatException e){
            javax.swing.JOptionPane.showMessageDialog(this,"ERROR; unable to parse double from '"+textfield_cutoff_value.getText()+"'");
            return;
        }
        //if I have a valid number update the attvals data
        if(data.blasthits!=null){
            synchronized(data.myattvals){
                data.myattvals=ClusterMethods.compute_attraction_values(data.blasthits,data.minpval,data);
            }
        }else if(data.myattvals!=null){//remove all attvals below the specified value
            if(data.orgattvals==null){
                data.orgattvals=data.myattvals;
            }
            data.myattvals=ClusterMethods.filter_attraction_values(data.orgattvals,data.minpval);
        }
        data.draworder=new ArrayList[0];
        repaint();
    }//GEN-LAST:event_minpvalbuttonActionPerformed
    
    private void savemtxmenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_savemtxmenuitemActionPerformed
        // save the blast matrix to file
        if(checkstop()==false){//check for stop of run
            return;
        }
        int returnVal = fc.showSaveDialog(this);
        if(returnVal == JFileChooser.APPROVE_OPTION) {
            File savefile=fc.getSelectedFile();
            int elementsl=java.lang.reflect.Array.getLength(data.namearr);
            try{
                PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter(savefile)));
                outwrite.println("sequences="+elementsl);
                outwrite.println("<seqs>");
                for(int i=0;i<elementsl;i++){
                    outwrite.println(">"+data.namearr[i]);
                }//end for i
                outwrite.println();
                outwrite.println("</seqs>");
                if(data.seqgroupsvec.size()>0){
                    outwrite.println("#user defined sequence groups");
                    outwrite.println("<seqgroups>");
                    seqgroup mygroup;
                    for(int i=data.seqgroupsvec.size()-1;i>=0;i--){
                        mygroup=(seqgroup)data.seqgroupsvec.elementAt(i);
                        outwrite.println("name="+mygroup.name);
                        outwrite.println("color="+mygroup.color.getRed()+";"+mygroup.color.getGreen()+";"+mygroup.color.getBlue());
                        outwrite.print("numbers=");
                        for(int j=java.lang.reflect.Array.getLength(mygroup.sequences)-1;j>=0;j--){
                            outwrite.print(mygroup.sequences[j]+";");
                        }//end for j
                        outwrite.println();
                    }//end for i
                    outwrite.println("</seqgroups>");
                }//end if seqgroupsvec.size>0
                outwrite.println("<att>");
                minattvals[] myattvals=data.myattvals;
                int datnum=java.lang.reflect.Array.getLength(myattvals);
                for(int i=0;i<datnum;i++){
                    outwrite.println(myattvals[i].query+" "+myattvals[i].hit+" :"+myattvals[i].att);
                }//end for i
                outwrite.println("</att>");
                outwrite.close();
                System.out.println("done");
            }catch (IOException ioe){
                javax.swing.JOptionPane.showMessageDialog(this,"IOERROR writing to '"+savefile.getAbsolutePath()+"'");
            }
        }
        
    }//GEN-LAST:event_savemtxmenuitemActionPerformed
    
    private void getdotsizemenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_getdotsizemenuitemActionPerformed
        String tmpstr="";
        try{
            tmpstr=JOptionPane.showInputDialog(this,"Enter the new size (int):",String.valueOf(data.dotsize));
            if(tmpstr!=null){
                data.dotsize=(int)(Float.parseFloat(tmpstr));//someone might enter a float and it's not time critical
            }
        }catch (NumberFormatException ne){
            javax.swing.JOptionPane.showMessageDialog(this,"ERROR, unable to parse integer from '"+tmpstr+"'");
        }
        repaint();
    }//GEN-LAST:event_getdotsizemenuitemActionPerformed
    
    private void showinfocheckboxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_showinfocheckboxActionPerformed
        repaint();
    }//GEN-LAST:event_showinfocheckboxActionPerformed
    
    private void setrotmenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_setrotmenuitemActionPerformed
        // specify the rotation values
        String tmpstr="";
        try{
            tmpstr=JOptionPane.showInputDialog(this,"Enter the new rotation values: x , y , z (9 values total)");
            if(tmpstr==null){
                return;
            }
            String[] tmparr=tmpstr.split(",");
            if(java.lang.reflect.Array.getLength(tmparr)!=9){
                JOptionPane.showMessageDialog(this,"You have to enter nine values separated by commas ','","Error",JOptionPane.ERROR_MESSAGE);
                return;
            }
            for(int i=0;i<3;i++){
                for(int j=0;j<3;j++){
                    data.myrotmtx[i][j]=Double.parseDouble(tmparr[i*3+j]);
                }
            }//end for i
            data.rotmtx[0][0]=data.myrotmtx[0][0];
            data.rotmtx[0][1]=data.myrotmtx[0][1];
            data.rotmtx[0][2]=data.myrotmtx[0][2];
            data.rotmtx[1][0]=data.myrotmtx[1][0];
            data.rotmtx[1][1]=data.myrotmtx[1][1];
            data.rotmtx[1][2]=data.myrotmtx[1][2];
            data.rotmtx[2][0]=data.myrotmtx[2][0];
            data.rotmtx[2][1]=data.myrotmtx[2][1];
            data.rotmtx[2][2]=data.myrotmtx[2][2];
        }catch (NumberFormatException ne){
            javax.swing.JOptionPane.showMessageDialog(this,"ERROR, unable to parse double from '"+tmpstr+"'");
        }
        repaint();
    }//GEN-LAST:event_setrotmenuitemActionPerformed
    
    private void getovalsizemenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_getovalsizemenuitemActionPerformed
        // set the size for the circles for selected sequences
        String tmpstr="";
        try{
            tmpstr=JOptionPane.showInputDialog(this,"Enter the new size(int):",String.valueOf(data.ovalsize));
            if(tmpstr!=null){
                data.ovalsize=Integer.parseInt(tmpstr);
            }
        }catch (NumberFormatException ne){
            javax.swing.JOptionPane.showMessageDialog(this,"ERROR, unable to parse int from '"+tmpstr+"'");
        }
        repaint();
    }//GEN-LAST:event_getovalsizemenuitemActionPerformed
    
    private void getblasthitsmenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_getblasthitsmenuitemActionPerformed
        // this should pop up a sequence selection menu, do a blast run for this sequence against
        // the present database, extract the hsp's (sequence regions and evalues)
        // and map it on to the selected sequence.
        //i.e. region 1-200 hits cluster A, region 210-300 cluster b, ergo 2 domains.
        int referenceseqnum=getsinglenamedialog.getrefseq(data.namearr);
        if(referenceseqnum==-1){
            javax.swing.JOptionPane.showMessageDialog(this,"Please select a sequence");
            return;
        }
        //get the blast hits to this sequence
        hsp[] thishsp=(new viewblasthitsutils()).gethsps(referenceseqnum,data.sequences,data.cmd,data.formatdbpath,data.blastpath,data.addblastvbparam,data.referencedb,data.mineval,data.minpval);
        //plot these blast hits on to the sequence
        viewblasthits myview=new viewblasthits(this,thishsp,referenceseqnum,data.namearr,data.sequences[referenceseqnum],data.nameshash);
        viewblasthitsvec.addElement(myview);
        myview.setVisible(true);
    }//GEN-LAST:event_getblasthitsmenuitemActionPerformed
    
    private void savemenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_savemenuitemActionPerformed
        // save the current run to file
        if(checkstop()==false){//check for stop of run
            return;
        }
        int returnVal = fc.showSaveDialog(this);
        if(returnVal == JFileChooser.APPROVE_OPTION) {
            this.setTitle("Clustering of "+fc.getSelectedFile().getName());
            savetofile(fc.getSelectedFile());
            //saverunobject myrun=new saverunobject();
            //myrun.file=fc.getSelectedFile();
            //myrun.inaln=inaln;
            //myrun.blasthits=blasthits;
            //myrun.attvals=myattvals;
            //myrun.posarr=myposarr;
            //myrun.maxmove=maxmove;
            //myrun.pval=minpval;
            //myrun.usescval=usescval;
            //if(attvalcompcheckbox.isSelected()){
            //    myrun.complexatt=true;
            //}else{
            //    myrun.complexatt=false;
            //}
            //myrun.rotmtx=draw1.rotmtx;
            //myrun.seqgroupsvec=seqgroupsvec;
            //myrun.cooling=cooling;
            //myrun.currcool=currcool;
            //myrun.attfactor=attfactor;
            //myrun.attvalpow=attvalpow;
            //myrun.repfactor=repfactor;
            //myrun.repvalpow=repvalpow;
            //myrun.dampening=dampening;
            //myrun.minattract=minattract;
            //myrun.blastpath=blastpath;
            //myrun.formatdbpath=formatdbpath;
            //myrun.dotsize=draw1.dotsize;
            //myrun.ovalsize=draw1.ovalsize;
            //myrun.groupsize=draw1.groupsize;
            //myrun.mapfiles=mapfiles;
            //myrun.lookupfiles=lookupfiles;
            //myrun.usefoldchange=usefoldchange;
            //myrun.avgfoldchange=avgfoldchange;
            //myrun.affyfiles=affyfiles;
            //myrun.namesdmp_file=namesdmp_file;
            //myrun.nodesdmp_file=nodesdmp_file;
            //if(cluster2dbutton.isSelected()){
            //    myrun.cluster2d=true;
            //}else{
            //    myrun.cluster2d=false;
            //}
            //myrun.colorarr=draw1.colorarr;
            //myrun.colorcutoffs=draw1.colorcutoffs;
            //customutils.saverun(myrun,namearr);
            //myrun=null;
        }
    }//GEN-LAST:event_savemenuitemActionPerformed
    
    private void loadmenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_loadmenuitemActionPerformed
        // load my custom format for a run from disk
        if(checkstop()==false){//check for stop of run
            return;
        }
        groupseqs=null;
        int returnVal = fc.showOpenDialog(this);
        if(returnVal == JFileChooser.APPROVE_OPTION) {
            loaddata(fc.getSelectedFile().getAbsolutePath());
            if(myseqgroupwindow!=null){
                myseqgroupwindow.setVisible(false);
                myseqgroupwindow.dispose();
            }
            if(mymapfunctiondialog!=null){
                mymapfunctiondialog.setVisible(false);
                mymapfunctiondialog.dispose();
            }
            repaint();
        }
    }//GEN-LAST:event_loadmenuitemActionPerformed
    
    private void hidesingletonsmenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_hidesingletonsmenuitemActionPerformed
        // hide sequences that have no hit whatsoever to others in the dataset
        if(checkstop()==false){//check for stop of run
            return;
        }
        groupseqs=null;
        int sequences=java.lang.reflect.Array.getLength(data.namearr);
        int j;
        //look through all sequences to see if they have an attraction value to some other sequence
        boolean[] keepseqs=new boolean[sequences];
        for(int i=0;i<sequences;i++){
            keepseqs[i]=false;
        }
        minattvals[] myattvals=data.myattvals;
        int attnum=java.lang.reflect.Array.getLength(myattvals);
        for(int i=0;i<attnum;i++){
            if(myattvals[i].att!=0){
                keepseqs[myattvals[i].query]=true;
                keepseqs[myattvals[i].hit]=true;
            }
        }//end for i
        int counter=0;
        for(int i=0;i<sequences;i++){
            if(keepseqs[i]==true){
                counter++;
            }
        }
        data.selectednames=new int[counter];
        counter=0;
        for(int i=0;i<sequences;i++){
            if(keepseqs[i]==true){
                data.selectednames[counter]=i;
                counter++;
            }
        }
        int selectednum=java.lang.reflect.Array.getLength(data.selectednames);
        if(selectednum<2){
            return;
        }
        if(shownames!=null){
            shownames.setVisible(false);
            shownames.dispose();
            shownames=null;
        }
        if(myseqgroupwindow!=null){
            myseqgroupwindow.setVisible(false);
            myseqgroupwindow.dispose();
            myseqgroupwindow=null;
        }
        if(viewblasthitsvec.size()>0){
            for(int i=0;i<viewblasthitsvec.size();i++){
                blastselectseqs=new int[0];
                ((viewblasthits)viewblasthitsvec.elementAt(i)).setVisible(false);
                ((viewblasthits)viewblasthitsvec.elementAt(i)).dispose();
            }
            viewblasthitsvec.setSize(0);
        }
        level++;
        getparentmenuitem.setEnabled(true);
        getparentmenuitem.setText("use parent group ("+level+")");
        parentnameshash.addElement(data.nameshash);
        data.nameshash=new HashMap((int)(selectednum/0.75)+1,(float)0.75);//holds info about which name is which array number
        for(int i=0;i<selectednum;i++){
            data.nameshash.put(data.sequences[data.selectednames[i]].name,new Integer(i));
        }
        if(data.blasthits!=null){
            parentblasthits.addElement(data.blasthits);
            data.blasthits=zoomdata.getblasthitsubset(data.blasthits,data.selectednames);
        }else{
            if(data.orgattvals==null){
                data.orgattvals=data.myattvals;
            }
            parentblasthits.addElement(data.orgattvals);
            data.orgattvals=zoomdata.getmyattvalssubset(data.orgattvals,data.selectednames);
            data.myattvals=ClusterMethods.filter_attraction_values(data.orgattvals,data.minpval);
        }
        parentmovearr.addElement(data.mymovearr);
        data.mymovearr=zoomdata.getmymovearrsubset(data.mymovearr,data.selectednames);
        if(java.lang.reflect.Array.getLength(data.mymovearr)>0){
            data.lastmovearr=new float[java.lang.reflect.Array.getLength(data.mymovearr)][java.lang.reflect.Array.getLength(data.mymovearr[0])];
        }
        parentposarr.addElement(data.myposarr);
        data.myposarr=zoomdata.getmyposarrsubset(data.myposarr,data.selectednames);
        parentaln.addElement(data.sequences);
        data.sequences=zoomdata.getinalnsubset(data.sequences,data.selectednames);
        parentnamearr.addElement(data.namearr);
        data.namearr=zoomdata.getnamearrsubset(data.namearr,data.selectednames);
        data.selectednames=new int[0];
        if(data.blasthits!=null){
            synchronized(data.myattvals){
                data.myattvals=ClusterMethods.compute_attraction_values(data.blasthits,data.minpval,data);
            }
        }
        data.elements=java.lang.reflect.Array.getLength(data.namearr);
        data.posarr=data.myposarr;
        data.posarrtmp=new float[data.elements][data.dimensions];
        data.drawarrtmp=new int[data.elements][data.dimensions];
        data.draworder=new ArrayList[0];
        repaint();
    }//GEN-LAST:event_hidesingletonsmenuitemActionPerformed
    
    private void getparentmenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_getparentmenuitemActionPerformed
        // zoom out one level (use all sequences from the level before for further computation
        if(checkstop()==false){//check for stop of run
            return;
        }
        groupseqs=null;
        if(shownames!=null){
            shownames.setVisible(false);
            shownames.dispose();
            shownames=null;
        }
        if(myseqgroupwindow!=null){
            myseqgroupwindow.setVisible(false);
            myseqgroupwindow.dispose();
            myseqgroupwindow=null;
        }
        if(viewblasthitsvec.size()>0){
            for(int i=0;i<viewblasthitsvec.size();i++){
                blastselectseqs=new int[0];
                ((viewblasthits)viewblasthitsvec.elementAt(i)).setVisible(false);
                ((viewblasthits)viewblasthitsvec.elementAt(i)).dispose();
            }
            viewblasthitsvec.setSize(0);
        }
        level--;
        if(level==0){
            getparentmenuitem.setEnabled(false);
        }
        getparentmenuitem.setText("use parent group ("+level+")");
        if(data.blasthits!=null){
            data.blasthits=(minhsp[])parentblasthits.elementAt(level);
            parentblasthits.removeElementAt(level);
        }else{
            data.orgattvals=(minattvals[]) parentblasthits.elementAt(level);
            parentblasthits.removeElementAt(level);
            data.myattvals=ClusterMethods.filter_attraction_values(data.orgattvals,data.minpval);
        }
        data.mymovearr=(float[][])parentmovearr.elementAt(level);
        if(java.lang.reflect.Array.getLength(data.mymovearr)>0){
            data.lastmovearr=new float[java.lang.reflect.Array.getLength(data.mymovearr)][java.lang.reflect.Array.getLength(data.mymovearr[0])];
        }
        parentmovearr.removeElementAt(level);
        data.myposarr=(float[][])parentposarr.elementAt(level);
        parentposarr.removeElementAt(level);
        data.sequences=(AminoAcidSequence[])parentaln.elementAt(level);
        parentaln.removeElementAt(level);
        data.namearr=(String[])parentnamearr.elementAt(level);
        parentnamearr.removeElementAt(level);
        data.nameshash=(HashMap)parentnameshash.elementAt(level);
        parentnameshash.removeElementAt(level);
        data.weights=(float[])parentweights.elementAt(level);
        parentweights.removeElementAt(level);
        data.selectednames=new int[0];
        data.elements=java.lang.reflect.Array.getLength(data.namearr);
        if(data.blasthits!=null){
            synchronized(data.myattvals){
                data.myattvals=ClusterMethods.compute_attraction_values(data.blasthits,data.minpval,data);
            }
        }
        if(myseqgroupwindow!=null){
            myseqgroupwindow.setVisible(false);
            myseqgroupwindow.dispose();
        }
        if(parentseqgroups.size()>level){
            data.seqgroupsvec=(Vector)parentseqgroups.remove(level);
        }
        if(mymapfunctiondialog!=null){
            mymapfunctiondialog.makenameshash();
        }
        data.posarr=data.myposarr;
        data.posarrtmp=new float[data.elements][data.dimensions];
        data.drawarrtmp=new int[data.elements][data.dimensions];
        data.draworder=new ArrayList[0];
        repaint();
    }//GEN-LAST:event_getparentmenuitemActionPerformed
    
    private void getchildmenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_getchildmenuitemActionPerformed
        // for further computation use only the selected sequences.(zoom in one level)
        if(checkstop()==false){//check for stop of run
            return;
        }
        groupseqs=null;
        int selectednum=java.lang.reflect.Array.getLength(data.selectednames);
        if(selectednum<2){
            return;
        }
        if(shownames!=null){
            shownames.setVisible(false);
            shownames.dispose();
            shownames=null;
        }
        if(myseqgroupwindow!=null){
            myseqgroupwindow.setVisible(false);
            myseqgroupwindow.dispose();
            myseqgroupwindow=null;
        }
        if(myseqgroupwindow!=null){
            myseqgroupwindow.setVisible(false);
            myseqgroupwindow.dispose();
            data.showseqgroups=false;
        }
        selectvec.clear();
        if(viewblasthitsvec.size()>0){
            for(int i=0;i<viewblasthitsvec.size();i++){
                blastselectseqs=new int[0];
                ((viewblasthits)viewblasthitsvec.elementAt(i)).setVisible(false);
                ((viewblasthits)viewblasthitsvec.elementAt(i)).dispose();
            }
            viewblasthitsvec.setSize(0);
        }
        level++;
        getparentmenuitem.setEnabled(true);
        getparentmenuitem.setText("use parent group ("+level+")");
        parentnameshash.addElement(data.nameshash);
        data.nameshash=new HashMap((int)(selectednum/0.75)+1,(float)0.75);//holds info about which name is which array number
        for(int i=0;i<selectednum;i++){
            data.nameshash.put(data.sequences[data.selectednames[i]].name,new Integer(i));
        }
        if(data.blasthits!=null){
            parentblasthits.addElement(data.blasthits);
            data.blasthits=zoomdata.getblasthitsubset(data.blasthits,data.selectednames);
        }
        parentmovearr.addElement(data.mymovearr);
        data.mymovearr=zoomdata.getmymovearrsubset(data.mymovearr,data.selectednames);
        parentposarr.addElement(data.myposarr);
        data.myposarr=zoomdata.getmyposarrsubset(data.myposarr,data.selectednames);
        parentaln.addElement(data.sequences);
        data.sequences=zoomdata.getinalnsubset(data.sequences,data.selectednames);
        parentnamearr.addElement(data.namearr);
        data.namearr=zoomdata.getnamearrsubset(data.namearr,data.selectednames);
        parentweights.addElement(data.weights);
        data.weights=zoomdata.getweightssubset(data.weights,data.selectednames);
        data.elements=java.lang.reflect.Array.getLength(data.namearr);
        if(data.blasthits==null){
            if(data.orgattvals==null){
                data.orgattvals=data.myattvals;
            }
            parentblasthits.addElement(data.orgattvals);
            data.orgattvals=zoomdata.getmyattvalssubset(data.orgattvals,data.selectednames);
            data.myattvals=ClusterMethods.filter_attraction_values(data.orgattvals,data.minpval);
        }
        parentseqgroups.addElement(data.seqgroupsvec);
        data.seqgroupsvec=new Vector();
        data.selectednames=new int[0];
        synchronized(data.myattvals){
            data.myattvals=ClusterMethods.compute_attraction_values(data.blasthits,data.minpval,data);
        }
        data.posarr=data.myposarr;
        data.posarrtmp=new float[data.elements][data.dimensions];
        data.drawarrtmp=new int[data.elements][data.dimensions];
        data.draworder=new ArrayList[0];
        if(mymapfunctiondialog!=null){
            mymapfunctiondialog.makenameshash();
        }
        repaint();
    }//GEN-LAST:event_getchildmenuitemActionPerformed
    
    private void evalueitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_evalueitemActionPerformed
        // Add your handling code here:
        if(data.blasthits!=null){
            if(data.usescval){//remember I'm in score mode
                eplotdialog eplot=new eplotdialog(data.blasthits,data.minpval,true);
                eplot.setVisible(true);
            }else{//i'm in P-value mode
                eplotdialog eplot=new eplotdialog(data.blasthits,data.minpval,false);
                eplot.setVisible(true);
            }
        }else{
            attplotdialog attplot=new attplotdialog(data.myattvals,data.minpval);
            attplot.setVisible(true);
        }
    }//GEN-LAST:event_evalueitemActionPerformed
    
    private void sequencesitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_sequencesitemActionPerformed
        // Add your handling code here:
        if(shownames!=null){
            shownames.setVisible(false);
            shownames.dispose();
        }
        shownames=new shownamedialog(data.namearr,this);
        shownames.setVisible(true);
    }//GEN-LAST:event_sequencesitemActionPerformed
    
    private void getseqsmenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_getseqsmenuitemActionPerformed
        // Add your handling code here:
        //save the currently selected sequences to file
        if(checkstop()==false){//check for stop of run
            return;
        }
        int returnVal = fc.showSaveDialog(this);
        if(returnVal == JFileChooser.APPROVE_OPTION) {
            AminoAcidSequence[] selectedseqs=getselectedseqs();
            if(java.lang.reflect.Array.getLength(selectedseqs)==0){
                javax.swing.JOptionPane.showMessageDialog(this,"Please select some sequences");
                return;
            }
            printout.printfasta(selectedseqs,fc.getSelectedFile());
        }
    }//GEN-LAST:event_getseqsmenuitemActionPerformed
    
    private void changecolormenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_changecolormenuitemActionPerformed
        // change the colors
        changecolordialog.changecolor(this,data.colorarr);
        repaint();
    }//GEN-LAST:event_changecolormenuitemActionPerformed
    
    private void checkbox_show_numbersActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_shownumberscheckboxActionPerformed
        // Add your handling code here:
        mousemove[0]=0;
        mousemove[1]=0;
        repaint();
    }//GEN-LAST:event_shownumberscheckboxActionPerformed
    
    private void button_clear_selectionActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_clearselectbuttonActionPerformed
    	if (!this.contains_data(true)) {
			return;
		}
        // Add your handling code here:
        if(java.lang.reflect.Array.getLength(data.selectednames)>0){//if sth. is selected clear the selection
            data.selectednames=new int[0];
            if(shownames!=null){
                shownames.seqnamelist.setSelectedIndices(data.selectednames);//clear all selected
            }
            if(zoom==true){
                zoom=false;
                button_zoom_on_selected.setText("Zoom on selected");
            }
            button_select_all_or_clear.setText("Select All");
        }else{//if nothing is selected, select all
            int alnseqs=java.lang.reflect.Array.getLength(data.sequences);
            data.selectednames=new int[alnseqs];
            for(int i=0;i<alnseqs;i++){
                data.selectednames[i]=i;
            }
            if(shownames!=null){
                shownames.seqnamelist.setSelectedIndices(data.selectednames);
            }
            button_select_all_or_clear.setText("Clear Selection");
        }
        repaint();
    }//GEN-LAST:event_clearselectbuttonActionPerformed
    
    private void button_zoom_on_selectedActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_zoombuttonActionPerformed
        //see if I want to zoom in on selected sequences
    	if (!this.contains_data(true)) {
			return;
		}
    	if(zoom==true){
            button_zoom_on_selected.setText("Zoom on selected");
            zoom=false;
        }else{
            button_zoom_on_selected.setText("Show all");
            zoom=true;
            if(java.lang.reflect.Array.getLength(data.selectednames)<1){
                zoom=false;
                JOptionPane.showMessageDialog(null,"Please select some sequences","Message",JOptionPane.INFORMATION_MESSAGE);
                button_zoom_on_selected.setText("Zoom on selected");
            }
        }
        data.zoomfactor=1;
        draw1.xtranslate=0;
        draw1.ytranslate=0;
        repaint();
    }//GEN-LAST:event_zoombuttonActionPerformed
    
    private void textfield_cutoff_valueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_minpvaltextfieldActionPerformed
    	if (!this.contains_data(true)) {
			return;
		}
        if(checkstop()==false){//check for stop of run
            return;
        }
        try{
            data.minpval=java.lang.Double.parseDouble(textfield_cutoff_value.getText());
        }catch (NumberFormatException e){
            javax.swing.JOptionPane.showMessageDialog(this,"ERROR, unable to parse double from '"+textfield_cutoff_value.getText()+"'");
            return;
        }
        //if I have a valid number update the attvals data
        if(data.blasthits!=null){
            synchronized(data.myattvals){
                data.myattvals=ClusterMethods.compute_attraction_values(data.blasthits,data.minpval,data);
            }
        }else if(data.myattvals!=null){//remove all attvals below the specified value
            if(data.orgattvals==null){
                System.out.println("setting orgattval=myattvals");
                data.orgattvals=data.myattvals;
            }
            data.myattvals=ClusterMethods.filter_attraction_values(data.orgattvals,data.minpval);
        }
        data.draworder=new ArrayList[0];
        repaint();
    }//GEN-LAST:event_minpvaltextfieldActionPerformed
    
    private void button_select_moveActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_selecttbuttonActionPerformed
        // Add your handling code here
        if(button_select_move.isSelected()){
            button_select_move.setText("SELECT/move");
        }else{
            button_select_move.setText("select/MOVE");
        }
    }//GEN-LAST:event_selecttbuttonActionPerformed
    
    private void graphpanelMouseReleased(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_graphpanelMouseReleased
        draw1.drawbox=false;
        mousemove[0]=0;
        mousemove[1]=0;
        mouse_is_pressed=false;
        if(moveseqs==true){
            moveseqs=false;
            int movex=evt.getX()-mousestart[0];
            int movey=evt.getY()-mousestart[1];
            moveselected(movex,movey);
        }else{
            if(button_select_move.isSelected()){
                int[] tmpreg=new int[4];
                tmpreg[0]=selectstart[0];
                tmpreg[1]=selectstart[1];
                tmpreg[2]=evt.getX()-draw1.xtranslate;
                tmpreg[3]=evt.getY()-draw1.ytranslate;
                if(evt.isAltDown()||evt.isControlDown()||evt.isShiftDown()){
                    //deselect all sequences in the selected square
                    updateselected(tmpreg,true);
                }else{
                    updateselected(tmpreg,false);
                }
                clusterconf=null;
            }
            //deepcopy the myrotmtx to rotmtx (do NOT assign reference)
            data.rotmtx[0][0]=data.myrotmtx[0][0];
            data.rotmtx[0][1]=data.myrotmtx[0][1];
            data.rotmtx[0][2]=data.myrotmtx[0][2];
            data.rotmtx[1][0]=data.myrotmtx[1][0];
            data.rotmtx[1][1]=data.myrotmtx[1][1];
            data.rotmtx[1][2]=data.myrotmtx[1][2];
            data.rotmtx[2][0]=data.myrotmtx[2][0];
            data.rotmtx[2][1]=data.myrotmtx[2][1];
            data.rotmtx[2][2]=data.myrotmtx[2][2];
            draw1.tmprotmtx=new double[3][3];
        }
        repaint();
    }//GEN-LAST:event_graphpanelMouseReleased
    
    private void graphpanelMousePressed(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_graphpanelMousePressed
        mouse_is_pressed=true;
        if(button_select_move.isSelected()==false){
            if(evt.isAltDown()||evt.isControlDown()||evt.isMetaDown()){//if I want to drag a sequence in 2d
                //move all selected sequences a certain amount
                moveseqs=true;
                mousestart[0]=evt.getX();
                mousestart[1]=evt.getY();
            }else if(evt.isShiftDown()){
                mousestart[0]=evt.getX();
                mousestart[1]=evt.getY();
                translate[0]=draw1.xtranslate;
                translate[1]=draw1.ytranslate;
            }else{
                draw1.drawbox=true;
                mousestart[0]=evt.getX();
                mousestart[1]=evt.getY();
                //if mouse is inside a certain area rotate x,y; outside rotate z
                if(stereocheckboxmenuitem.isSelected()==false){
                    if((mousestart[0]<draw1.xadd)||(mousestart[0]>draw1.drawwidth+draw1.xadd)||(mousestart[1]<draw1.yadd)||(mousestart[1]>draw1.drawheight+draw1.yadd)){
                        draw1.rotatez=true;
                    }else{
                        draw1.rotatez=false;
                    }
                }else{//if in stereo mode
                    if((mousestart[0]<draw1.xadd)||(mousestart[0]>draw1.drawwidth*2+draw1.xadd)||(mousestart[1]<draw1.yadd)||(mousestart[1]>draw1.drawheight+draw1.yadd)){
                        draw1.rotatez=true;
                    }else{
                        draw1.rotatez=false;
                    }
                }
                mousemove[0]=0;
                mousemove[1]=0;
            }
        }else{
            mousemove[0]=0;
            mousemove[1]=0;
            selectstart[0]=evt.getX()-draw1.xtranslate;
            selectstart[1]=evt.getY()-draw1.ytranslate;
            currmousepos[0]=selectstart[0];
            currmousepos[1]=selectstart[1];
        }
        repaint();
    }//GEN-LAST:event_graphpanelMousePressed
    
    private void graphpanelMouseDragged(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_graphpanelMouseDragged
        if(button_select_move.isSelected()==false){
            if(evt.isShiftDown()){
                draw1.xtranslate=evt.getX()-mousestart[0]+translate[0];
                draw1.ytranslate=evt.getY()-mousestart[1]+translate[1];
            }else if(moveseqs){
                
            }else{
                mousemove[0]=evt.getX()-mousestart[0];
                mousemove[1]=evt.getY()-mousestart[1];
            }
        }else{
            if(shownamesselectcheckbox.isSelected()){
                if(shownames==null){
                    shownames=new shownamedialog(data.namearr,this);
                    shownames.setVisible(true);
                }
                int[] tmpreg=new int[4];
                tmpreg[0]=selectstart[0];
                tmpreg[1]=selectstart[1];
                tmpreg[2]=evt.getX()-draw1.xtranslate;
                tmpreg[3]=evt.getY()-draw1.ytranslate;
                updatetmpselected(tmpreg);
            }
            mousemove[0]=0;
            mousemove[1]=0;
            currmousepos[0]=evt.getX()-draw1.xtranslate;
            currmousepos[1]=evt.getY()-draw1.ytranslate;
        }
        repaint();
    }//GEN-LAST:event_graphpanelMouseDragged
    
    private void checkbox_show_namesItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_shownamescheckboxItemStateChanged
        mousemove[0]=0;
        mousemove[1]=0;
        repaint();
    }//GEN-LAST:event_shownamescheckboxItemStateChanged
    
    private void button_start_stop_resumeActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_stopbuttonActionPerformed
    	if (!this.contains_data(true)) {
			return;
		}
    	startstopthread();
        /*//does the iteration and the stop for a thread
        mousemove[0]=0;
        mousemove[1]=0;
        if(mythread.stop==true){//if this thread was stopped
            String tmpstr="";
            if(mymovearr==null){
                javax.swing.JOptionPane.showMessageDialog(this,"WARNING: No data currently available; please load some data!");
                return;
            }
            if(myoptionswindow!=null){
                try{
                    tmpstr=myoptionswindow.coolfield.getText();
                    cooling=java.lang.Double.parseDouble(tmpstr);
                    tmpstr=myoptionswindow.attfield.getText();
                    attfactor=java.lang.Float.parseFloat(tmpstr);
                    tmpstr=myoptionswindow.repfield.getText();
                    repfactor=java.lang.Float.parseFloat(tmpstr);
                    tmpstr=myoptionswindow.dampfield.getText();
                    dampening=java.lang.Float.parseFloat(tmpstr);
                    tmpstr=myoptionswindow.minattfield.getText();
                    minattract=java.lang.Double.parseDouble(tmpstr);
                    tmpstr=myoptionswindow.maxmovefield.getText();
                    maxmove=java.lang.Float.parseFloat(tmpstr);
                    tmpstr=myoptionswindow.attvalpowtextfield.getText();
                    attvalpow=java.lang.Integer.parseInt(tmpstr);
                    tmpstr=myoptionswindow.repvalpowtextfield.getText();
                    repvalpow=java.lang.Integer.parseInt(tmpstr);
                }catch (NumberFormatException e){
                    javax.swing.JOptionPane.showMessageDialog(this,"ERROR, unable to parse number from '"+tmpstr+"'");
                }
            }
            if(hidebelowold!=hidebelow){
                draw1.draworder=new ArrayList[0];
                hidebelowold=hidebelow;
            }
            mythread.stop=true;
            mythread=new computethread();
            mythread.start();
            stopbutton.setText("Stop");
        }else{//if this thread is running
            mythread.stop=true;
            stopbutton.setEnabled(false);
            //is unavailable until the thread has stopped running
            //the thread then sets the text to "resume" and re-enables the button.
        }
        */
        repaint();
    }//GEN-LAST:event_stopbuttonActionPerformed
    
    private void checkbox_show_connectionsItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_showblasthitscheckboxItemStateChanged
        mousemove[0]=0;
        mousemove[1]=0;
        repaint();
    }//GEN-LAST:event_showblasthitscheckboxItemStateChanged
    
    private void graphpanelAncestorResized(java.awt.event.HierarchyEvent evt) {//GEN-FIRST:event_graphpanelAncestorResized
        mousemove[0]=0;
        mousemove[1]=0;
        repaint();
    }//GEN-LAST:event_graphpanelAncestorResized
    
    private void button_initializeActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_startbuttonActionPerformed
    	if (!this.contains_data(true)) {
			return;
		}
    	initgraph();
        /*currcool=1;
        String tmpstr="";
        if(myoptionswindow!=null){
            try{
                tmpstr=myoptionswindow.attfield.getText();
                attfactor=java.lang.Float.parseFloat(tmpstr);
                tmpstr=myoptionswindow.repfield.getText();
                repfactor=java.lang.Float.parseFloat(tmpstr);
                tmpstr=myoptionswindow.dampfield.getText();
                dampening=java.lang.Float.parseFloat(tmpstr);
                tmpstr=myoptionswindow.coolfield.getText();
                cooling=java.lang.Double.parseDouble(tmpstr);
                tmpstr=myoptionswindow.currcoolfield.getText();
                currcool=java.lang.Double.parseDouble(tmpstr);
                tmpstr=myoptionswindow.minattfield.getText();
                minattract=java.lang.Double.parseDouble(tmpstr);
                tmpstr=myoptionswindow.maxmovefield.getText();
                maxmove=java.lang.Float.parseFloat(tmpstr);
                tmpstr=myoptionswindow.attvalpowtextfield.getText();
                attvalpow=java.lang.Integer.parseInt(tmpstr);
                tmpstr=myoptionswindow.repvalpowtextfield.getText();
                repvalpow=java.lang.Integer.parseInt(tmpstr);
            }catch(NumberFormatException e){
                javax.swing.JOptionPane.showMessageDialog(this,"ERROR, unable to parse number from '"+tmpstr+"'");
            }
        }
        if(hidebelowold!=hidebelow){
            draw1.draworder=new ArrayList[0];
            hidebelowold=hidebelow;
        }
        changedvals=false;
        if(checkstop(false)==true){//if thread is stopped
            updatevals();
            repaint();
        }
        mythread.stop=true;
        myposarr=null;
        mymovearr=null;
        lastmovearr=null;
        myposarr=cluster3d(blasthits,3);
        draw1.posarr=myposarr;
        if(myoptionswindow!=null){
            myoptionswindow.currcoolfield.setText(String.valueOf(currcool));
        }
        stopbutton.setText("Start run");
        stopbutton.setEnabled(true);
        rounds=0;
        mousemove[0]=0;
        mousemove[1]=0;
        */
        repaint();
    }//GEN-LAST:event_startbuttonActionPerformed
    
    /** Exit the Application */
    private void exitForm(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_exitForm
        System.exit(0);
    }//GEN-LAST:event_exitForm

    private void taxonomymenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_taxonomymenuitemActionPerformed
        if(taxonomydialog!=null){
            taxonomydialog.setVisible(false);
            taxonomydialog.dispose();
        }
        taxonomydialog=new ncbitaxonomydialog(this);
        taxonomydialog.setVisible(true);
    }//GEN-LAST:event_taxonomymenuitemActionPerformed

    private void attvalcompcheckboxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_attvalcompcheckboxActionPerformed
        if(attvalcompcheckbox.isSelected()){
            data.attvalsimple=false;
        }else{
            data.attvalsimple=true;
        }
    }//GEN-LAST:event_attvalcompcheckboxActionPerformed

    private void rescalepvaluescheckboxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_rescalepvaluescheckboxActionPerformed
        if(rescalepvaluescheckbox.isSelected()){
            data.rescalepvalues=true;
        }else{
            data.rescalepvalues=false;
        }
    }//GEN-LAST:event_rescalepvaluescheckboxActionPerformed

    private void moveselectedonlyActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_moveselectedonlyActionPerformed
        if(moveselectedonly.isSelected()){
            data.moveselectedonly=true;
        }else{
            data.moveselectedonly=false;
        }
    }//GEN-LAST:event_moveselectedonlyActionPerformed

    private void saveattvalsmenuitemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_saveattvalsmenuitemActionPerformed
        // save the similarity matrix data to a file
        if(checkstop()==false){//check for stop of run
            return;
        }
        int returnVal = fc.showSaveDialog(this);
        if(returnVal == JFileChooser.APPROVE_OPTION) {
            this.setTitle("Clustering of "+fc.getSelectedFile().getName());
            customutils.saveattvalstofile(fc.getSelectedFile(),data.myattvals,data.nographics);
        }
    }//GEN-LAST:event_saveattvalsmenuitemActionPerformed
    
    /**
     * @param args the command line arguments
     */
    //used for backup of blasthit data when working with selection levels
    int level=0;
    Vector <float[][]> parentmovearr=new Vector<float[][]>();
    Vector <float[][]> parentposarr=new Vector<float[][]>();
    Vector <AminoAcidSequence[]> parentaln=new Vector<AminoAcidSequence[]>();
    Vector <String[]> parentnamearr=new Vector<String[]>();
    Vector <HashMap> parentnameshash=new Vector<HashMap>();
    //don't define parentblasthits as once I add minhsp's and once I add minattvals
    Vector parentblasthits=new Vector();
    Vector <Vector> parentseqgroups=new Vector<Vector>();
    Vector <float[]> parentweights=new Vector<float[]>();
    //end of zoom section variables
    
    //minhsp[] blasthits=null;
    //float[] weights=null;
    Vector <viewblasthits> viewblasthitsvec=new Vector<viewblasthits>();
    Vector <selectclass> selectvec=new Vector<selectclass>();
    //String[] namearr;
    //String loadsaved=null;// do I want to load a file on startup?
    //aaseq[] inaln;
    //HashMap nameshash;
    boolean zoom=false;
    boolean mouse_is_pressed=false;
    boolean moveseqs=false;//flag to set if I want to draw sequences in 3d-space
    //boolean changedvals=false;//if the input values (used in iteration) were changed
    //boolean attvalsimple=true;//how do I want to compute the attraction values
    boolean savepos=false;
    //boolean usescval=false;//am I usign P-values or scvalues?
    //String docalc="true";//to synchronize drawing and calculating
    boolean recalc=true;//synchronize drawing and calculating
    drawpanel draw1;
    shownamedialog shownames;
    seqgroupwindow myseqgroupwindow;
    //groupslegenddialog groupslegend=null;
    //boolean showseqgroups=false;
    //Vector <seqgroup> seqgroupsvec=new Vector<seqgroup>();
    computethread mythread=new computethread(this);
    //String cmd="";
    //String blastpath="blastall -p blastp";
    //boolean addblastvbparam=true;
    //String formatdbpath="formatdb -p T";
    //String[] referencedb=new String[0];
    //static float[][] myposarr;
    //static float[][] mymovearr;
    //static float[][] lastmovearr;
    //minattvals[] orgattvals=null;//the loaded attraction values
    //minattvals[] myattvals=null;
    //float[] seqlengths;
    //double minattract=1;//minimal attraction, centers the whole thing around origin
    //float maxmove=0.1f;//maximum movement per round for a point
    //double cooling=1;//how much the temp is reduced each round (multiplier)
    //double currcool=1;//the current cooling factor (movement multiplier)
    //double mineval=1;//maximum evalue to use
    //double minpval=1;//maximum p-value to use
    //double maxvalfound=0;//the maximum Pvalue found
    //float repfactor=(float)5;//repulsion, decreases in strength with distance**repfactor
    //float attfactor=(float)10;//attraction, increases in strength with distance**attfactor
    //int attvalpow=1;
    //int repvalpow=1;
    //float hidebelow=0;//hide the connections with values below
    //float hidebelowold=0;//to see if I need to recompute the drawdata
    //float dampening=(float)0.2;
    //int elements=0;
    //static int dimentions=3;
    //int rounds=0;
    //int roundslimit=-1;
    //int roundsdone=0;
    //int verbose=0;
    //int cpu=1;
    int skiprounds=1;
    int[] mousemove=new int[2];
    int[] translate=new int[2];
    int[] mousestart=new int[2];
    int[] selectstart=new int[2];
    int[] currmousepos=new int[2];
    //int[] selectednames=new int[0];//this gets accessed from the shownamesdialog
    int[] groupseqs=null;//used by seqgroupsdialog to show which group is currently selected
    java.awt.Color groupseqscolor=new java.awt.Color(0,0,0);
    int[] blastselectseqs=new int[0];
    float[] clusterconf=null;//confidence values for the selected sequences
    getmovethread[] movethreads;
    //static final Random rand=new Random(System.currentTimeMillis());
    static final JFileChooser fc=new JFileChooser(new File("."));
    String repaint=null;
    //static double ln10=java.lang.Math.log(10);//used to convert log_base_e to log_base_10
    //double p2attfactor=1;
    //double p2attoffset=0;
    optionsdialog myoptionswindow=null;//clans options
    rotationdialog myrotationdialog=null;//clans rotation values
    affydialog myaffydialog=null;//loads/shows affymetrix data
    mapfunctiondialog_tab mymapfunctiondialog=null;//loads/shows metabolic/functional mapping
    ncbitaxonomydialog taxonomydialog=null;//show the NCBI taxonomic keywords
    String namesdmp_file=null;//"names.dmp";
    String nodesdmp_file=null;//"nodes.dmp";
    ArrayList <java.io.File>mapfiles=new ArrayList<java.io.File>();
    ArrayList <java.io.File>lookupfiles=new ArrayList<java.io.File>();
    Vector affyfiles=null;
    boolean usefoldchange=false;
    boolean avgfoldchange=false;
    boolean dotsfirst=false;
    //boolean nographics=false;

    clusterdata data=null;
    
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JMenu menu_misc;
    private javax.swing.JMenuItem aboutmenuitem;
    private javax.swing.JMenuItem addseqsmenuitem;
    private javax.swing.JMenuItem affymenuitem;
    private javax.swing.JCheckBoxMenuItem antialiasingcheckboxmenuitem;
    private javax.swing.JCheckBoxMenuItem attvalcompcheckbox;
    private javax.swing.JPanel buttonpanel;
    private javax.swing.JMenuItem centermenuitem;
    private javax.swing.JMenuItem changebgcolormenuitem;
    private javax.swing.JMenuItem changeblastcolor;
    private javax.swing.JMenuItem changecolormenuitem;
    private javax.swing.JMenuItem changefgcolormenuitem;
    private javax.swing.JMenuItem changefontmenuitem;
    private javax.swing.JMenuItem changenumbercolor;
    private javax.swing.JMenuItem changeselectcolormenuitem;
    public javax.swing.JButton button_select_all_or_clear;
    private javax.swing.JCheckBoxMenuItem cluster2dbutton;
    private javax.swing.JMenuItem clustermenuitem;
    private javax.swing.JCheckBoxMenuItem colorfrustrationcheckbox;
    private javax.swing.JPanel drawbuttonpanel;
    private javax.swing.JMenu menu_draw;
    private javax.swing.JMenuItem evalueitem;
    private javax.swing.JMenu menu_file;
    private javax.swing.JMenuItem getblasthitsmenuitem;
    private javax.swing.JMenuItem getchildmenuitem;
    private javax.swing.JMenuItem getdotsizemenuitem;
    private javax.swing.JMenuItem getovalsizemenuitem;
    private javax.swing.JMenuItem getparentmenuitem;
    private javax.swing.JMenuItem getseqsforselectedhits;
    private javax.swing.JMenuItem getseqsmenuitem;
    private javax.swing.JPanel graphpanel;
    private javax.swing.JMenu menu_help;
    private javax.swing.JMenuItem helpmenuitem;
    private javax.swing.JMenuItem hidesingletonsmenuitem;
    private javax.swing.JMenuBar jMenuBar1;
    private javax.swing.JCheckBoxMenuItem lengthcolormenuitem;
    private javax.swing.JMenuItem loadalternatemenuitem;
    private javax.swing.JMenuItem loadgroupsmenuitem;
    private javax.swing.JMenuItem loadmenuitem;
    private javax.swing.JMenuItem loadtabsmenuitem;
    private javax.swing.JMenuItem mapmanmenuitem;
    private javax.swing.JButton button_cutoff_value;
    private javax.swing.JTextField textfield_cutoff_value;
    private javax.swing.JTextField textfield_info_min_blast_evalue;
    private javax.swing.JCheckBoxMenuItem moveselectedonly;
    private javax.swing.JMenuItem printmenuitem;
    private javax.swing.JCheckBoxMenuItem rescalepvaluescheckbox;
    private javax.swing.JMenuItem rotationmenuitem;
    private javax.swing.JMenuItem save2dmenuitem;
    private javax.swing.JMenuItem saveattvalsmenuitem;
    private javax.swing.JMenuItem savemenuitem;
    private javax.swing.JMenuItem savemtxmenuitem;
    private javax.swing.JToggleButton button_select_move;
    private javax.swing.JMenuItem seqscoloring;
    private javax.swing.JMenuItem sequencesitem;
    private javax.swing.JMenuItem setrotmenuitem;
    private javax.swing.JCheckBoxMenuItem showblasthitnamescheckbox;
    private javax.swing.JCheckBox checkbox_show_connections;
    private javax.swing.JCheckBoxMenuItem showinfocheckbox;
    private javax.swing.JCheckBox checkbox_show_names;
    private javax.swing.JCheckBoxMenuItem shownamesselectcheckbox;
    private javax.swing.JCheckBox checkbox_show_numbers;
    private javax.swing.JMenuItem showoptionsmenuitem;
    private javax.swing.JCheckBoxMenuItem showorigcheckbox;
    private javax.swing.JButton button_show_selected;
    private javax.swing.JMenuItem showseqsmenuitem;
    private javax.swing.JMenuItem skipdrawingrounds;
    private javax.swing.JButton button_initialize;
    private javax.swing.JMenuItem stereoanglemenuitem;
    private javax.swing.JCheckBoxMenuItem stereocheckboxmenuitem;
    private javax.swing.JButton button_start_stop_resume;
    private javax.swing.JMenuItem taxonomymenuitem;
    private javax.swing.JMenu menu_windows;
    private javax.swing.JButton button_zoom_on_selected;
    private javax.swing.JMenuItem zoommenuitem;
    // End of variables declaration//GEN-END:variables

    void savetofile(File savetofile){
        ClusterMethods.savetofile(savetofile,data);
        /*
        saverunobject myrun=new saverunobject();
        myrun.file=savetofile;
        myrun.inaln=inaln;
        myrun.blasthits=blasthits;
        myrun.attvals=myattvals;
        myrun.posarr=myposarr;
        myrun.maxmove=maxmove;
        myrun.pval=minpval;
        myrun.usescval=usescval;
        if(attvalcompcheckbox.isSelected()){
            myrun.complexatt=true;
        }else{
            myrun.complexatt=false;
        }
        myrun.rotmtx=draw1.rotmtx;
        myrun.seqgroupsvec=seqgroupsvec;
        myrun.cooling=cooling;
        myrun.currcool=currcool;
        myrun.attfactor=attfactor;
        myrun.attvalpow=attvalpow;
        myrun.repfactor=repfactor;
        myrun.repvalpow=repvalpow;
        myrun.dampening=dampening;
        myrun.minattract=minattract;
        myrun.blastpath=blastpath;
        myrun.formatdbpath=formatdbpath;
        myrun.dotsize=draw1.dotsize;
        myrun.ovalsize=draw1.ovalsize;
        myrun.groupsize=draw1.groupsize;
        myrun.mapfiles=mapfiles;
        myrun.lookupfiles=lookupfiles;
        myrun.usefoldchange=usefoldchange;
        myrun.avgfoldchange=avgfoldchange;
        myrun.affyfiles=affyfiles;
        myrun.namesdmp_file=namesdmp_file;
        myrun.nodesdmp_file=nodesdmp_file;
        if(cluster2dbutton.isSelected()){
            myrun.cluster2d=true;
        }else{
            myrun.cluster2d=false;
        }
        myrun.colorarr=draw1.colorarr;
        myrun.colorcutoffs=draw1.colorcutoffs;
        customutils.saverun(myrun,namearr);
        myrun=null;
         */
    }

    //--------------------------------------------------------------------------

    void startstopthread(){
        //does the iteration and the stop for a thread
        mousemove[0]=0;
        mousemove[1]=0;
        if(mythread.didrun==false || mythread.stop==true){//if this thread was never started or stopped
            String tmpstr="";
            if(data.mymovearr==null){
                javax.swing.JOptionPane.showMessageDialog(this,"WARNING: No data currently available; please load some data!");
                return;
            }
            if(myoptionswindow!=null){
                try{
                    tmpstr=myoptionswindow.coolfield.getText();
                    data.cooling=java.lang.Double.parseDouble(tmpstr);
                    tmpstr=myoptionswindow.attfield.getText();
                    data.attfactor=java.lang.Float.parseFloat(tmpstr);
                    tmpstr=myoptionswindow.repfield.getText();
                    data.repfactor=java.lang.Float.parseFloat(tmpstr);
                    tmpstr=myoptionswindow.dampfield.getText();
                    data.dampening=java.lang.Float.parseFloat(tmpstr);
                    tmpstr=myoptionswindow.minattfield.getText();
                    data.minattract=java.lang.Double.parseDouble(tmpstr);
                    tmpstr=myoptionswindow.maxmovefield.getText();
                    data.maxmove=java.lang.Float.parseFloat(tmpstr);
                    tmpstr=myoptionswindow.attvalpowtextfield.getText();
                    data.attvalpow=java.lang.Integer.parseInt(tmpstr);
                    tmpstr=myoptionswindow.repvalpowtextfield.getText();
                    data.repvalpow=java.lang.Integer.parseInt(tmpstr);
                }catch (NumberFormatException e){
                    javax.swing.JOptionPane.showMessageDialog(this,"ERROR, unable to parse number from '"+tmpstr+"'");
                }
            }
            if(data.hidebelowold!=data.hidebelow){
                data.draworder=new ArrayList[0];
                data.hidebelowold=data.hidebelow;
            }
            mythread.stop=true;
            mythread=new computethread(this);
            mythread.start();
            button_start_stop_resume.setText("Stop");
        }else{//if this thread is running
            mythread.stop=true;
            button_start_stop_resume.setEnabled(false);
            //is unavailable until the thread has stopped running
            //the thread then sets the text to "resume" and re-enables the button.
        }
    }//end startstopthread

    //--------------------------------------------------------------------------

    void initgraph(){
        //first some stuff specific to the graphical interface
        String tmpstr="";
        if(myoptionswindow!=null){
            try{
                tmpstr=myoptionswindow.attfield.getText();
                data.attfactor=java.lang.Float.parseFloat(tmpstr);
                tmpstr=myoptionswindow.repfield.getText();
                data.repfactor=java.lang.Float.parseFloat(tmpstr);
                tmpstr=myoptionswindow.dampfield.getText();
                data.dampening=java.lang.Float.parseFloat(tmpstr);
                tmpstr=myoptionswindow.coolfield.getText();
                data.cooling=java.lang.Double.parseDouble(tmpstr);
                tmpstr=myoptionswindow.currcoolfield.getText();
                data.currcool=java.lang.Double.parseDouble(tmpstr);
                tmpstr=myoptionswindow.minattfield.getText();
                data.minattract=java.lang.Double.parseDouble(tmpstr);
                tmpstr=myoptionswindow.maxmovefield.getText();
                data.maxmove=java.lang.Float.parseFloat(tmpstr);
                tmpstr=myoptionswindow.attvalpowtextfield.getText();
                data.attvalpow=java.lang.Integer.parseInt(tmpstr);
                tmpstr=myoptionswindow.repvalpowtextfield.getText();
                data.repvalpow=java.lang.Integer.parseInt(tmpstr);
            }catch(NumberFormatException e){
                javax.swing.JOptionPane.showMessageDialog(this,"ERROR, unable to parse number from '"+tmpstr+"'");
            }
        }
        if(data.hidebelowold!=data.hidebelow){
            data.draworder=new ArrayList[0];
            data.hidebelowold=data.hidebelow;
        }
        data.changedvals=false;
        if(checkstop(false)==true){//if thread is stopped
            updatevals();
            repaint();
        }
        ClusterMethods.setup_attraction_values_and_initialize(data);
        //currcool=1;
        //mythread.stop=true;
        //myposarr=null;
        //mymovearr=null;
        //lastmovearr=null;
        //myposarr=cluster3d(blasthits,3);
        //draw1.posarr=myposarr;
        if(myoptionswindow!=null){
            myoptionswindow.currcoolfield.setText(String.valueOf(data.currcool));
        }
        button_start_stop_resume.setText("Start run");
        button_start_stop_resume.setMnemonic(KeyEvent.VK_S);
        button_start_stop_resume.setEnabled(true);
        data.rounds=0;
        mousemove[0]=0;
        mousemove[1]=0;
    }//end initgraph

    //--------------------------------------------------------------------------
    
    void loaddata(String inname){
        if(mythread!=null && mythread.didrun==true && mythread.stop!=true){
            System.err.println("Warning, you should stop the clustering thread before loading another file; stopping thread now");
            mythread.stop=true;
        }//else everything is OK
        data.input_filename=inname;
        ClusterMethods.loaddata(data);
        textfield_info_min_blast_evalue.setText(String.valueOf(data.maxvalfound));
        if(data.blasthits==null){
            button_cutoff_value.setText("Use Attraction values better than");
            textfield_cutoff_value.setText("0");
            savemtxmenuitem.setText("Save Attraction values as matrix");
            evalueitem.setText("Attraction value plot");
        }else{
            if(data.usescval){
                button_cutoff_value.setText("Use SC-vals better than");
                textfield_cutoff_value.setText("0");
                evalueitem.setText("SC-value plot");
            }else{
                button_cutoff_value.setText("Use P-values better than");
                textfield_cutoff_value.setText("1");
                evalueitem.setText("P-value plot");
            }
        }

        if(data.complexatt==true){
            attvalcompcheckbox.setSelected(true);
        }else{
            attvalcompcheckbox.setSelected(false);
        }
        if(data.cluster2d){
            cluster2dbutton.setSelected(true);
        }else{
            cluster2dbutton.setSelected(false);
        }
        if(data.showinfo){
            showinfocheckbox.setSelected(true);
        }else{
            showinfocheckbox.setSelected(false);
        }
        textfield_cutoff_value.setText(String.valueOf(data.minpval));
        if(myoptionswindow!=null){
            myoptionswindow.coolfield.setText(String.valueOf(data.cooling));
            myoptionswindow.currcoolfield.setText(String.valueOf(data.currcool));
            myoptionswindow.attfield.setText(String.valueOf(data.attfactor));
            myoptionswindow.attvalpowtextfield.setText(String.valueOf(data.attvalpow));
            myoptionswindow.repfield.setText(String.valueOf(data.repfactor));
            myoptionswindow.repvalpowtextfield.setText(String.valueOf(data.repvalpow));
            myoptionswindow.dampfield.setText(String.valueOf(data.dampening));
            myoptionswindow.minattfield.setText(String.valueOf(data.minattract));
            myoptionswindow.maxmovefield.setText(String.valueOf(data.maxmove));
        }
        textfield_info_min_blast_evalue.setText(String.valueOf(data.maxvalfound));
        this.setTitle("Clustering of "+inname);
    }//end loaddata
    
    //--------------------------------------------------------------------------
    
    /*void loadtabdata(String inname){
        saverunobject saveddata=customutils.loadrunbiolayout(new File(inname));
        if(saveddata.file!=null){//if the data was read all right
            System.out.println("File loaded:"+inname);
            repaint="Error loading data;";
            inaln=rmgaps(saveddata.inaln);
            myposarr=saveddata.posarr;
            blasthits=saveddata.blasthits;
            usescval=saveddata.usescval;
            if(blasthits==null){
                //synchronized(myattvals){//first time I load myattvals; cannot be anything else; don't need to sync
                myattvals=saveddata.attvals;
                //}
                minpvalbutton.setText("Use Attraction values better than");
                minpvaltextfield.setText("0");
                savemtxmenuitem.setText("Save Attraction values as matrix");
                evalueitem.setText("Attraction value plot");
            }else{
                if(usescval){
                    minpvalbutton.setText("Use SC-vals better than");
                    minpvaltextfield.setText("0");
                    evalueitem.setText("SC-value plot");
                }else{
                    minpvalbutton.setText("Use P-values better than");
                    minpvaltextfield.setText("1");
                    evalueitem.setText("P-value plot");
                }
            }
            if(saveddata.complexatt==true){
                attvalcompcheckbox.setSelected(true);
            }else{
                attvalcompcheckbox.setSelected(false);
            }
            //now set the graph values
            maxmove=saveddata.maxmove;
            minpval=saveddata.pval;
            cooling=saveddata.cooling;
            currcool=saveddata.currcool;
            attfactor=saveddata.attfactor;
            repfactor=saveddata.repfactor;
            attvalpow=saveddata.attvalpow;
            repvalpow=saveddata.repvalpow;
            dampening=saveddata.dampening;
            minattract=saveddata.minattract;
            weights=saveddata.weights;
            mapfiles=saveddata.mapfiles;
            lookupfiles=saveddata.lookupfiles;
            affyfiles=saveddata.affyfiles;
            //be careful not to overwrite any blastpath and formatdbpath setting passed via command line
            if(blastpath.equals("blastall -p blastp")){//if it was not changed via command line
                blastpath=saveddata.blastpath;
            }
            if(formatdbpath.equals("formatdb -p T")){//if it was not changed via command line
                formatdbpath=saveddata.formatdbpath;
            }
            draw1.zoomfactor=saveddata.zoom;
            if(saveddata.cluster2d){
                cluster2dbutton.setSelected(true);
            }else{
                cluster2dbutton.setSelected(false);
            }
            if(saveddata.showinfo){
                showinfocheckbox.setSelected(true);
            }else{
                showinfocheckbox.setSelected(false);
            }
            minpvaltextfield.setText(String.valueOf(minpval));
            if(myoptionswindow!=null){
                myoptionswindow.coolfield.setText(String.valueOf(cooling));
                myoptionswindow.currcoolfield.setText(String.valueOf(currcool));
                myoptionswindow.attfield.setText(String.valueOf(attfactor));
                myoptionswindow.attvalpowtextfield.setText(String.valueOf(attvalpow));
                myoptionswindow.repfield.setText(String.valueOf(repfactor));
                myoptionswindow.repvalpowtextfield.setText(String.valueOf(repvalpow));
                myoptionswindow.dampfield.setText(String.valueOf(dampening));
                myoptionswindow.minattfield.setText(String.valueOf(minattract));
                myoptionswindow.maxmovefield.setText(String.valueOf(maxmove));
            }
            int seqs=java.lang.reflect.Array.getLength(inaln);
            nameshash=new HashMap((int)(seqs/0.75)+1,(float)0.75);//holds info about which name is which array number
            namearr=new String[seqs];
            for(int i=0;i<seqs;i++){
                namearr[i]=inaln[i].name.trim();
                inaln[i].name=new String("sequence"+i);
                nameshash.put(inaln[i].name,new Integer(i));
            }
            elements=java.lang.reflect.Array.getLength(namearr);
            selectednames=new int[0];
            draw1.posarr=myposarr;
            lastmovearr=new float[elements][dimentions];
            mymovearr=new float[elements][dimentions];
            draw1.posarrtmp=new float[elements][dimentions];
            draw1.drawarrtmp=new int[elements][dimentions];
            draw1.draworder=new ArrayList[0];
            draw1.myrotmtx=saveddata.rotmtx;
            draw1.rotmtx[0][0]=draw1.myrotmtx[0][0];
            draw1.rotmtx[0][1]=draw1.myrotmtx[0][1];
            draw1.rotmtx[0][2]=draw1.myrotmtx[0][2];
            draw1.rotmtx[1][0]=draw1.myrotmtx[1][0];
            draw1.rotmtx[1][1]=draw1.myrotmtx[1][1];
            draw1.rotmtx[1][2]=draw1.myrotmtx[1][2];
            draw1.rotmtx[2][0]=draw1.myrotmtx[2][0];
            draw1.rotmtx[2][1]=draw1.myrotmtx[2][1];
            draw1.rotmtx[2][2]=draw1.myrotmtx[2][2];
            orgattvals=null;
            //synchronized (myattvals){//first time I load myattvals; don't need to sync as nothing else can be using this yet
            myattvals=clustermethods.getattvals(data.blasthits,data.minpval,data);
            //}
            minsetpvaltextfield.setText(String.valueOf(maxvalfound));
            draw1.dotsize=saveddata.dotsize;
            draw1.ovalsize=saveddata.ovalsize;
            draw1.groupsize=saveddata.groupsize;
            draw1.polygons=makepolygons.get(draw1.groupsize);
            seqgroupsvec=saveddata.seqgroupsvec;
            if(seqgroupsvec.size()>0){
                showseqgroups=true;
            }
            changedvals=true;
            if(checkstop(false)==true){//if thread is stopped
                updatevals();
                repaint();
            }
            repaint=null;
        }else{//if the data had errors
            String retval=JOptionPane.showInputDialog(this,"Error reading data; Do you want to display anyway? (yes/NO)","Error reading",JOptionPane.ERROR_MESSAGE);
            if(retval.startsWith("y")||retval.startsWith("Y")){
                inaln=rmgaps(saveddata.inaln);
                myposarr=saveddata.posarr;
                blasthits=saveddata.blasthits;
                maxmove=saveddata.maxmove;
                minpval=saveddata.pval;
                if(saveddata.cluster2d){
                    cluster2dbutton.setSelected(true);
                }else{
                    cluster2dbutton.setSelected(false);
                }
                minpvaltextfield.setText(String.valueOf(minpval));
                int seqs=java.lang.reflect.Array.getLength(inaln);
                nameshash=new HashMap((int)(seqs/0.75)+1,(float)0.75);//holds info about which name is which array number
                namearr=new String[seqs];
                for(int i=0;i<seqs;i++){
                    namearr[i]=inaln[i].name.trim();
                    inaln[i].name=new String("sequence"+i);
                    nameshash.put(inaln[i].name,new Integer(i));
                }
                elements=java.lang.reflect.Array.getLength(namearr);
                selectednames=new int[0];
                draw1.posarr=myposarr;
                lastmovearr=new float[elements][dimentions];
                mymovearr=new float[elements][dimentions];
                draw1.posarrtmp=new float[elements][dimentions];
                draw1.drawarrtmp=new int[elements][dimentions];
                draw1.draworder=new ArrayList[0];
                draw1.myrotmtx=saveddata.rotmtx;
                draw1.rotmtx[0][0]=draw1.myrotmtx[0][0];
                draw1.rotmtx[0][1]=draw1.myrotmtx[0][1];
                draw1.rotmtx[0][2]=draw1.myrotmtx[0][2];
                draw1.rotmtx[1][0]=draw1.myrotmtx[1][0];
                draw1.rotmtx[1][1]=draw1.myrotmtx[1][1];
                draw1.rotmtx[1][2]=draw1.myrotmtx[1][2];
                draw1.rotmtx[2][0]=draw1.myrotmtx[2][0];
                draw1.rotmtx[2][1]=draw1.myrotmtx[2][1];
                draw1.rotmtx[2][2]=draw1.myrotmtx[2][2];
                orgattvals=null;
                //synchronized(myattvals){
                myattvals=clustermethods.getattvals(data.blasthits,data.minpval,data);
                //}
                minsetpvaltextfield.setText(String.valueOf(maxvalfound));
                seqgroupsvec=saveddata.seqgroupsvec;
                if(seqgroupsvec.size()>0){
                    showseqgroups=true;
                }
                draw1.dotsize=saveddata.dotsize;
                draw1.ovalsize=saveddata.ovalsize;
                draw1.groupsize=saveddata.groupsize;
                changedvals=true;
                if(checkstop(false)==true){//if thread is stopped
                    updatevals();
                    repaint();
                }
                repaint=null;
            }else{
                //do nothing
            }
        }
        if(saveddata.colorarr!=null){
            System.out.println("setting colorarr");
            draw1.colorarr=saveddata.colorarr;
        }
        if(saveddata.colorcutoffs!=null){
            System.out.println("setting colorcutoffs");
            draw1.colorcutoffs=saveddata.colorcutoffs;
        }
        int seqnum;
        seqnum=java.lang.reflect.Array.getLength(inaln);
        System.out.println("seqnum="+seqnum);
        seqlengths=new float[seqnum];
        float maxlength=0;
        for(int i=0;i<seqnum;i++){
            seqlengths[i]=inaln[i].seq.replaceAll("-","").length();
            if(seqlengths[i]>maxlength){
                maxlength=seqlengths[i];
            }
        }//end for i
        for(int i=0;i<seqnum;i++){
            seqlengths[i]/=maxlength;
        }//end for i
    }//end loadtabdata
    */
    //--------------------------------------------------------------------------
	boolean contains_data() {
		return this.data.mymovearr != null;
	}
	
	boolean contains_data(boolean showinfo) {
		if (this.contains_data()){
			return true;
		}
		
		if (showinfo) {
			javax.swing.JOptionPane.showMessageDialog(this,
					"you have to load data first", "no data loaded",
					javax.swing.JOptionPane.ERROR_MESSAGE);
		}
		return false;
	}

    boolean checkstop(boolean showinfo){
        //check to see if the clustering is running
        //if so, output an error message
        if(mythread.didrun==true && mythread.stop==false){//if the thread is running
            if(showinfo){
                javax.swing.JOptionPane.showMessageDialog(this,"you have to stop the clustering first","Stop first",javax.swing.JOptionPane.ERROR_MESSAGE);
            }
            return false;
        }
        return true;
    }
    
    boolean checkstop(){
        return checkstop(true);
    }//end checkstop
    
    //--------------------------------------------------------------------------
    /*
    aaseq[] rmgaps(aaseq[] inseqs){
        int seqnum=java.lang.reflect.Array.getLength(inseqs);
        for(int i=0;i<seqnum;i++){
            inseqs[i].seq=inseqs[i].seq.replaceAll("-","");
        }//end for i
        return inseqs;
    }//end rmgaps
    */
    //--------------------------------------------------------------------------
    
    void initaddedseqs(minhsp[] blastvec,AminoAcidSequence[]allaln,String[]allnamearr,HashMap allnameshash,int[]newnumarr,float[][]allposarr,float maxmove,double pval,boolean useselectedonly){
        //initialize the necessary variable for the case where new sequences are added to an existing run.
        data.sequences=allaln;
        data.myposarr=allposarr;
        data.blasthits=blastvec;
        data.maxmove=maxmove;
        data.minpval=pval;
        textfield_cutoff_value.setText(String.valueOf(data.minpval));
        textfield_info_min_blast_evalue.setText(String.valueOf(data.minpval));
        data.selectednames=newnumarr;
        int seqs=java.lang.reflect.Array.getLength(data.sequences);
        data.nameshash=allnameshash;
        data.namearr=allnamearr;
        data.elements=java.lang.reflect.Array.getLength(data.namearr);
        data.posarr=data.myposarr;
        data.lastmovearr=new float[data.elements][data.dimensions];
        data.mymovearr=new float[data.elements][data.dimensions];
        data.posarrtmp=new float[data.elements][data.dimensions];
        data.drawarrtmp=new int[data.elements][data.dimensions];
        //draw1.draworder=new Vector[0];
        data.draworder=new ArrayList[0];
        data.myattvals=ClusterMethods.compute_attraction_values(data.blasthits,data.minpval,data);
        moveselectedonly.setSelected(useselectedonly);
        int seqnum=java.lang.reflect.Array.getLength(data.sequences);
        System.out.println("seqnum="+seqnum);
        data.seqlengths=new float[seqnum];
        float maxlength=0;
        for(int i=0;i<seqnum;i++){
            data.seqlengths[i]=data.sequences[i].seq.replaceAll("-","").length();
            if(data.seqlengths[i]>maxlength){
                maxlength=data.seqlengths[i];
            }
        }//end for i
        for(int i=0;i<seqnum;i++){
            data.seqlengths[i]/=maxlength;
        }//end for i
        //now all my sequences have a value assigned between 0 and 1 reflecting their length
        
    }//end initaddedseqs
    
    //--------------------------------------------------------------------------
    
    AminoAcidSequence[] getselectedseqs(){
        //get the currently selected sequences
        ArrayList <AminoAcidSequence>tmpvec=new ArrayList<AminoAcidSequence>();
        int sequences=java.lang.reflect.Array.getLength(data.selectednames);
        AminoAcidSequence curraaseq;
        for(int i=0;i<sequences;i++){
            curraaseq=new AminoAcidSequence();
            curraaseq.name=data.namearr[data.selectednames[i]];
            curraaseq.seq=data.sequences[data.selectednames[i]].seq;
            tmpvec.add(curraaseq);
        }
        AminoAcidSequence[] retarr=(AminoAcidSequence[]) tmpvec.toArray(new AminoAcidSequence[0]);
        return retarr;
    }//end getselectedseqs
    
    //--------------------------------------------------------------------------
    /*
    float[][] cluster3d(minhsp[] indata){
        //take the hsp objects from indata and compute "attraction" values for all sequence pairs
        //once you have those try to cluster the data in 2d by "energy minimization" approach.
        //iterative approch, might want to specify maximum number of iterations
        return cluster3d(indata,10);//this runs for max 10 iterations
        //return value is an array with x and y coords for all sequences
    }//end cluster2d
    
    float[][] cluster3d(minhsp[] indata,int maxiter){
        //take the hsp objects from indata and compute "attraction" values for all sequence pairs
        //once you have those try to cluster the data in 2d by "energy minimization" approach.
        //iterative approch, might want to specify maximum number of iterations
        elements=java.lang.reflect.Array.getLength(namearr);
        if(elements==0){
            return new float[0][0];
        }
        dimentions=3;
        myposarr=new float[elements][dimentions];
        draw1.posarrtmp=new float[elements][dimentions];
        draw1.drawarrtmp=new int[elements][dimentions];
        mymovearr=new float[elements][dimentions];
        lastmovearr=new float[elements][dimentions];
        for(int i=elements;--i>=0;){
            lastmovearr[i][0]=0;
            lastmovearr[i][1]=0;
            lastmovearr[i][2]=0;
        }
        //compute the "attraction values
        if(indata!=null){
            if(myattvals==null){
                //synchronized(myattvals){//myattvals is null here; cannot sync on it
                myattvals=clustermethods.getattvals(indata,data.minpval,data);
                //}
            }
        }
        //now i have the matrix with the attraction values for each seqpair
        //next: seed the 2d environment randomly with the starting points for the sequences
        myposarr=seedstart(myposarr);
        //then iterate
        mymovearr=new float[elements][dimentions];
        for(int i=elements;--i>=0;){
            mymovearr[i][0]=0;
            mymovearr[i][1]=0;
            mymovearr[i][2]=0;
        }
        return myposarr;
    }// end cluster3d
    */
    void updatevals(){
        //System.out.println("udating values");
        String tmpstr="";
        if(myoptionswindow!=null){
            try{
                tmpstr=myoptionswindow.attfield.getText();
                data.attfactor=java.lang.Float.parseFloat(tmpstr);
                tmpstr=myoptionswindow.repfield.getText();
                data.repfactor=java.lang.Float.parseFloat(tmpstr);
                tmpstr=myoptionswindow.dampfield.getText();
                data.dampening=java.lang.Float.parseFloat(tmpstr);
                tmpstr=myoptionswindow.coolfield.getText();
                data.cooling=java.lang.Double.parseDouble(tmpstr);
                tmpstr=myoptionswindow.currcoolfield.getText();
                data.currcool=java.lang.Double.parseDouble(tmpstr);
                tmpstr=myoptionswindow.minattfield.getText();
                data.minattract=java.lang.Double.parseDouble(tmpstr);
                tmpstr=myoptionswindow.maxmovefield.getText();
                data.maxmove=java.lang.Float.parseFloat(tmpstr);
                tmpstr=myoptionswindow.attvalpowtextfield.getText();
                data.attvalpow=java.lang.Integer.parseInt(tmpstr);
                tmpstr=myoptionswindow.repvalpowtextfield.getText();
                data.repvalpow=java.lang.Integer.parseInt(tmpstr);
                tmpstr=myoptionswindow.roundstextfield.getText();
                data.roundslimit=java.lang.Integer.parseInt(tmpstr);
            }catch(NumberFormatException e){
                javax.swing.JOptionPane.showMessageDialog(this,"ERROR, unable to parse number from '"+tmpstr+"'");
                //System.err.println("ERROR "+tmpstr+" is not a number");
            }
        }
    }//end updatevals


    /*
    float[][] recluster3d(){
        if(data.changedvals){
            updatevals();
            if(hidebelowold!=hidebelow){
                //draw1.draworder=new Vector[0];
                draw1.draworder=new ArrayList[0];
                hidebelowold=hidebelow;
            }
            changedvals=false;
        }//end if changedvals
        clustermethods.recluster3d(data);

        //take the hsp objects from indata and compute "attraction" values for all sequence pairs
        //once you have those try to cluster the data in 2d by "energy minimization" approach.
        //iterative approch, might want to specify maximum number of iterations
        //use the positions array and the attracion/repulsion values to
        //compute movement vectors for each object
        //long time=System.currentTimeMillis();
        if(changedvals){
            updatevals();
            if(hidebelowold!=hidebelow){
                //draw1.draworder=new Vector[0];
                draw1.draworder=new ArrayList[0];
                hidebelowold=hidebelow;
            }
            changedvals=false;
        }//end if changedvals
        System.arraycopy(mymovearr, 0, lastmovearr, 0, elements);//move all values from mymovearr to lastmovearr
        for(int i=elements;--i>=0;){
            //lastmovearr[i][0]=mymovearr[i][0];
            //lastmovearr[i][1]=mymovearr[i][1];
            //lastmovearr[i][2]=mymovearr[i][2];
            mymovearr[i][0]=0;
            mymovearr[i][1]=0;
            mymovearr[i][2]=0;
        }
        int selnamenum=java.lang.reflect.Array.getLength(selectednames);
        int[] selectnames=new int[selnamenum];//(int[])selectednames.clone();//a copy of the array that won't change during this round of calculations
        System.arraycopy(selectednames,0, selectnames,0, selnamenum);
        String syncme="syncme";
        synchronized(myattvals){
            if(cpu==1){
                //mymovearr=getmovement(myposarr,myattvals,mymovearr);
                getmovement(myposarr,myattvals,mymovearr,selectnames);
                domove(myposarr,mymovearr,selectnames);
            }else{
                //compute the attraction values using multiple threads
                int selectnamesnum=java.lang.reflect.Array.getLength(selectednames);
                if((moveselectedonly.isSelected())&&(selectnamesnum>0)&&(selectnamesnum<elements)){
                    selectnamesnum=java.lang.reflect.Array.getLength(selectnames);
                    HashMap <String,Integer>tmphash=new HashMap<String,Integer>();
                    int threadnum=selectnamesnum/cpu;
                    int tmpval;
                    for(int i=selectnamesnum-1;i>=0;i--){
                        tmpval=(int)i/threadnum;
                        if(tmpval>cpu-1){
                            tmpval=cpu-1;
                        }
                        tmphash.put(String.valueOf(selectnames[i]),new Integer(tmpval));
                    }//end for i
                    for(int i=0;i<cpu;i++){
                        //movethreads[i]=new getmovethread(myposarr,myattvals,mymovearr,i,cpu,tmphash,selectnames,this);//start cpu threads to write the attraction values to mymovearr
                        movethreads[i]=new getmovethread(myposarr,myattvals,mymovearr,i,cpu,tmphash,selectnames,syncme,this);//start cpu threads to write the attraction values to mymovearr
                        movethreads[i].start();
                    }
                }else{
                    for(int i=0;i<cpu;i++){
                        //movethreads[i]=new getmovethread(myposarr,myattvals,mymovearr,i,cpu,this);//start cpu threads to write the attraction values to mymovearr
                        movethreads[i]=new getmovethread(myposarr,myattvals,mymovearr,i,cpu,syncme,this);//start cpu threads to write the attraction values to mymovearr
                        movethreads[i].start();
                    }
                }
                //now wait for all threads to finish
                boolean alldone=false;
                try{
                    //synchronized(this){
                    synchronized(syncme){
                        while (alldone==false){
                            alldone=true;
                            for(int i=0;i<cpu;i++){
                                if(movethreads[i].done==false){
                                    alldone=false;
                                    break;
                                }
                            }
                            if(alldone==false){//if I still have to wait for a thread
                                //this.wait();//release lock on this and wait
                                syncme.wait();//release lock on this and wait
                            }
                        }
                    }//end while alldone
                }catch (InterruptedException e){
                    System.err.println("Interrupted wait in searchblast");
                    e.printStackTrace();
                }
                domove(myposarr,mymovearr,selectnames);
            }//end else cpu>1
        }//end sync myattvals
        //then move each according to movevec
        //System.out.println("calctime="+(System.currentTimeMillis()-time));
        return myposarr;
    }// end recluster3d
    */
    //--------------------------------------------------------------------------
    /*
    void domove(float[][] posarr, float[][] movearr,int[] selectnames){
        //move all objects according to their movement vector.
        currcool=currcool*cooling;
        //float totalpos;
        float multdamp=1-dampening;
        int selectnamesnum=java.lang.reflect.Array.getLength(selectnames);
        if((moveselectedonly.isSelected())&&(selectnamesnum>0)&&(selectnamesnum!=elements)){
            //System.out.println("domoveselected");
            if(cluster2dbutton.isSelected()){
                for(int i=selectnamesnum;--i>=0;){
                    posarr[selectnames[i]][0]+=(currcool*((lastmovearr[selectnames[i]][0]*multdamp)+(movearr[selectnames[i]][0])));
                    posarr[selectnames[i]][1]+=(currcool*((lastmovearr[selectnames[i]][1]*multdamp)+(movearr[selectnames[i]][1])));
                    //posarr[selectednames[i]][2]+=(currcool*((lastmovearr[selectednames[i]][2]*(1-dampening))+(movearr[selectednames[i]][2])));
                }//end for i
            }else{//cluster in 3d
                for(int i=selectnamesnum;--i>=0;){
                    posarr[selectnames[i]][0]+=(currcool*((lastmovearr[selectnames[i]][0]*multdamp)+(movearr[selectnames[i]][0])));
                    posarr[selectnames[i]][1]+=(currcool*((lastmovearr[selectnames[i]][1]*multdamp)+(movearr[selectnames[i]][1])));
                    posarr[selectnames[i]][2]+=(currcool*((lastmovearr[selectnames[i]][2]*multdamp)+(movearr[selectnames[i]][2])));
                }//end for i
            }
        }else{
            if(cluster2dbutton.isSelected()){
                for(int i=elements;--i>=0;){
                    posarr[i][0]+=(currcool*((lastmovearr[i][0]*multdamp)+(movearr[i][0])));
                    posarr[i][1]+=(currcool*((lastmovearr[i][1]*multdamp)+(movearr[i][1])));
                    //posarr[i][2]+=(currcool*((lastmovearr[i][2]*(1-dampening))+(movearr[i][2])));
                }// end for i
            }else{//cluster in 3d
                for(int i=elements;--i>=0;){
                    posarr[i][0]+=(currcool*((lastmovearr[i][0]*multdamp)+(movearr[i][0])));
                    posarr[i][1]+=(currcool*((lastmovearr[i][1]*multdamp)+(movearr[i][1])));
                    posarr[i][2]+=(currcool*((lastmovearr[i][2]*multdamp)+(movearr[i][2])));
                }// end for i
            }
        }
        //return posarr;
    }// end domove
    */
    //--------------------------------------------------------------------------
    /*
    void getmovement(float[][] posarr, minattvals[] attvals, float[][] movement,int[] selectnames){
        //use the positions of all elements and their attraction/repulsion values to
        //calculate a movement vector for each (take into account the last movement).
        //repulsion doesn't have a specific value as all evalues below a certain point
        //are simply regarded as insignificant. therefore use a one formula for all
        //to compute the repulsive forces (see getrepulse)
        //inline all the getrepulse and getattract methods to reduce the number of method calls
        int i,j;
        int attnum=java.lang.reflect.Array.getLength(attvals);
        double[] currmoverep=new double[3];
        double[] currmoveatt=new double[3];
        float tmp;
        double totaldist=0;
        double totalmove=1;
        double distx,disty,distz;
        int selectnamesnum=java.lang.reflect.Array.getLength(selectnames);
        float repfac=repfactor;
        //float attfac=attfactor;
        //int attpow=attvalpow;
        int reppow=repvalpow;
        int hnum,qnum;
        float weight1=1,weight2=1;
        if((moveselectedonly.isSelected())&&(selectnamesnum>0)&&(selectnamesnum!=elements)){
            HashMap tmphash=new HashMap<Integer,Integer>((int)(selectnamesnum/0.8)+1,0.8f);
            Integer[] hashkeys=new Integer[elements];
            //hashkey[] hashkeys=new hashkey[elements];
            for(i=elements;--i>=0;){
                //hashkeys[i]=new hashkey(i);
                hashkeys[i]=new Integer(i);
            }
            if(cluster2dbutton.isSelected()){
                //no point in inlining this bit; I have very few selected sequences and therefore the gain should be marginal
                //cluster only the selected sequences in 2D
                for(i=selectnamesnum;--i>=0;){
                    //currmoveatt=getminattract(posarr[selectnames[i]],currmoveatt,minattract);
                    getminattract(posarr[selectnames[i]],currmoveatt,minattract);
                    weight1=1;
                    if(weights!=null){
                        weight1=weights[selectnames[i]];
                    }
                    tmphash.put(hashkeys[selectnames[i]],null);
                    movement[selectnames[i]][0]+=currmoveatt[0]*weight1;
                    movement[selectnames[i]][1]+=currmoveatt[1]*weight1;
                    for(j=elements;--j>=0;){
                        if(j==selectnames[i]){
                            continue;
                        }
                        //currmoverep=getrepulse2d(posarr[selectnames[i]],posarr[j],currmoverep,repvalpow, repfactor, rand);
                        getrepulse2d(posarr[selectnames[i]],posarr[j],currmoverep,repvalpow, repfactor, rand);
                        weight2=1;
                        if(weights!=null){
                            weight2=weights[selectnames[j]];
                        }
                        movement[selectnames[i]][0]+=currmoverep[0]*weight2;
                        movement[selectnames[i]][1]+=currmoverep[1]*weight2;
                    }//end for j
                }//end for i
                for(i=attnum;--i>=0;){
                    if(tmphash.containsKey(hashkeys[attvals[i].query])){
                        //currmoveatt=getattract2d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow,attfactor);
                        getattract2d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow,attfactor);
                        weight1=1;
                        weight2=1;
                        if(weights!=null){
                            weight1=weights[attvals[i].query];
                            weight2=weights[attvals[i].hit];
                        }
                        movement[attvals[i].query][0]+=currmoveatt[0]*weight2;
                        movement[attvals[i].query][1]+=currmoveatt[1]*weight2;
                        if(tmphash.containsKey(hashkeys[attvals[i].hit])){
                            movement[attvals[i].hit][0]-=currmoveatt[0]*weight1;
                            movement[attvals[i].hit][1]-=currmoveatt[1]*weight1;
                        }
                    }else if(tmphash.containsKey(hashkeys[attvals[i].hit])){
                        //currmoveatt=getattract2d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow,attfactor);
                        getattract2d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow,attfactor);
                        weight1=1;
                        if(weights!=null){
                            weight1=weights[attvals[i].query];
                        }
                        movement[attvals[i].hit][0]-=currmoveatt[0]*weight1;
                        movement[attvals[i].hit][1]-=currmoveatt[1]*weight1;
                    }
                }//end for i
                for(i=selectnamesnum;--i>=0;){
                    movement[selectnames[i]][0]/=elements;
                    movement[selectnames[i]][1]/=elements;
                    movement[selectnames[i]][2]=0;
                    totaldist=java.lang.Math.sqrt((movement[selectnames[i]][0]*movement[selectnames[i]][0])+(movement[selectnames[i]][1]*movement[selectnames[i]][1]));
                    if(totaldist>maxmove){
                        tmp=(float)(maxmove/totaldist);
                        movement[selectnames[i]][0]*=tmp;
                        movement[selectnames[i]][1]*=tmp;
                        //movement[selectednames[i]][2]*=maxmove/totaldist;
                    }
                }//end for i
            }else{
                //cluster only the selected sequences in 3D
                for(i=selectnamesnum;--i>=0;){
                    //currmoveatt=getminattract(posarr[selectnames[i]],currmoveatt,minattract);
                    getminattract(posarr[selectnames[i]],currmoveatt,minattract);
                    weight1=1;
                    if(weights!=null){
                        weight1=weights[selectnames[i]];
                    }
                    movement[selectnames[i]][0]+=currmoveatt[0]*weight1;
                    movement[selectnames[i]][1]+=currmoveatt[1]*weight1;
                    movement[selectnames[i]][2]+=currmoveatt[2]*weight1;
                    tmphash.put(hashkeys[selectnames[i]],null);
                    for(j=elements;--j>=0;){
                        if(j==selectnames[i]){
                            continue;
                        }
                        //currmoverep=getrepulse3d(posarr[selectnames[i]],posarr[j],currmoverep,repvalpow, repfactor,rand);
                        getrepulse3d(posarr[selectnames[i]],posarr[j],currmoverep,repvalpow, repfactor,rand);
                        weight2=1;
                        if(weights!=null){
                            weight2=weights[selectnames[j]];
                        }
                        movement[selectnames[i]][0]+=currmoverep[0]*weight2;
                        movement[selectnames[i]][1]+=currmoverep[1]*weight2;
                        movement[selectnames[i]][2]+=currmoverep[2]*weight2;
                    }//end for j
                }//end for i
                for(i=attnum;--i>=0;){
                    if(tmphash.containsKey(hashkeys[attvals[i].query])){
                        //currmoveatt=getattract3d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow,attfactor);
                        //no point in inlining this bit; I reuse the results a second time
                        getattract3d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow,attfactor);
                        weight1=1;
                        weight2=1;
                        if(weights!=null){
                            weight1=weights[attvals[i].query];
                            weight2=weights[attvals[i].hit];
                        }
                        movement[attvals[i].query][0]+=currmoveatt[0]*weight2;
                        movement[attvals[i].query][1]+=currmoveatt[1]*weight2;
                        movement[attvals[i].query][2]+=currmoveatt[2]*weight2;
                        if(tmphash.containsKey(hashkeys[attvals[i].hit])){
                            movement[attvals[i].hit][0]-=currmoveatt[0]*weight1;
                            movement[attvals[i].hit][1]-=currmoveatt[1]*weight1;
                            movement[attvals[i].hit][2]-=currmoveatt[2]*weight1;
                        }
                    }else if(tmphash.containsKey(hashkeys[attvals[i].hit])){
                        //currmoveatt=getattract3d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow,attfactor);
                        getattract3d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow,attfactor);
                        weight1=1;
                        if(weights!=null){
                            weight1=weights[attvals[i].query];
                        }
                        movement[attvals[i].hit][0]-=currmoveatt[0]*weight1;
                        movement[attvals[i].hit][1]-=currmoveatt[1]*weight1;
                        movement[attvals[i].hit][2]-=currmoveatt[2]*weight1;
                    }
                }//end for i
                for(i=selectnamesnum;--i>=0;){
                    movement[selectnames[i]][0]/=elements;
                    movement[selectnames[i]][1]/=elements;
                    movement[selectnames[i]][2]/=elements;
                    totaldist=java.lang.Math.sqrt((movement[selectnames[i]][0]*movement[selectnames[i]][0])+(movement[selectnames[i]][1]*movement[selectnames[i]][1])+(movement[selectnames[i]][2]*movement[selectnames[i]][2]));
                    if(totaldist>maxmove){
                        tmp=(float)(maxmove/totaldist);
                        movement[selectnames[i]][0]*=tmp;
                        movement[selectnames[i]][1]*=tmp;
                        movement[selectnames[i]][2]*=tmp;
                    }
                }//end for i
            }
        }else{//if no sequences were selected or all should be used
            if(cluster2dbutton.isSelected()){
                //cluster all in 2D
                if(weights!=null){
                    for(i=elements;--i>=0;){
                        weight1=weights[i];
                        movement[i][0]-=posarr[i][0]*minattract*weight1;
                        movement[i][1]-=posarr[i][1]*minattract*weight1;
                        for(j=i+1;j<elements;j++){
                            weight2=weights[j];
                            distx=posarr[j][0]-posarr[i][0];
                            disty=posarr[j][1]-posarr[i][1];
                            if(distx==0 && disty==0){
                                //if two points are at exactly the same position I need to add some random effect
                                movement[i][0]+=repfactor*(rand.nextDouble()-0.5)*0.001;
                                movement[i][1]+=repfactor*(rand.nextDouble()-0.5)*0.001;
                                movement[j][0]-=repfactor*(rand.nextDouble()-0.5)*0.001;
                                movement[j][1]-=repfactor*(rand.nextDouble()-0.5)*0.001;
                                //return movement;
                            }else{
                                totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty));
                                //here I scale the force**repvalpow
                                totalmove=1;
                                for(int r=repvalpow;--r>=0;){
                                    totalmove*=totaldist;
                                }
                                totalmove=(repfactor/totalmove);
                                movement[i][0]+=(-distx/totaldist)*totalmove*weight2;
                                movement[i][1]+=(-disty/totaldist)*totalmove*weight2;
                                movement[j][0]-=(-distx/totaldist)*totalmove*weight1;
                                movement[j][1]-=(-disty/totaldist)*totalmove*weight1;
                            }
                        }//end for j
                    }//end for i
                    for(i=attnum;--i>=0;){
                        getattract2d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow,attfactor);
                        weight1=weights[attvals[i].query];
                        weight2=weights[attvals[i].hit];
                        movement[attvals[i].query][0]+=currmoveatt[0]*weight2;
                        movement[attvals[i].query][1]+=currmoveatt[1]*weight2;
                        movement[attvals[i].hit][0]-=currmoveatt[0]*weight1;
                        movement[attvals[i].hit][1]-=currmoveatt[1]*weight1;
                    }//end for i
                }else{//if weights==null, the use a default weighting of 1
                    for(i=elements;--i>=0;){
                        movement[i][0]-=posarr[i][0]*minattract;
                        movement[i][1]-=posarr[i][1]*minattract;
                        for(j=i+1;j<elements;j++){
                            distx=posarr[j][0]-posarr[i][0];
                            disty=posarr[j][1]-posarr[i][1];
                            if(distx==0 && disty==0){
                                //if two points are at exactly the same position I need to add some random effect
                                movement[i][0]+=repfactor*(rand.nextDouble()-0.5)*0.001;
                                movement[i][1]+=repfactor*(rand.nextDouble()-0.5)*0.001;
                                movement[j][0]-=repfactor*(rand.nextDouble()-0.5)*0.001;
                                movement[j][1]-=repfactor*(rand.nextDouble()-0.5)*0.001;
                                //return movement;
                            }else{
                                totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty));
                                //here I scale the force**repvalpow
                                totalmove=1;
                                for(int r=repvalpow;--r>=0;){
                                    totalmove*=totaldist;
                                }
                                totalmove=(repfactor/totalmove);
                                movement[i][0]+=(-distx/totaldist)*totalmove;
                                movement[i][1]+=(-disty/totaldist)*totalmove;
                                movement[j][0]-=(-distx/totaldist)*totalmove;
                                movement[j][1]-=(-disty/totaldist)*totalmove;
                            }
                        }//end for j
                    }//end for i
                    for(i=attnum;--i>=0;){
                        getattract2d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow,attfactor);
                        movement[attvals[i].query][0]+=currmoveatt[0];
                        movement[attvals[i].query][1]+=currmoveatt[1];
                        movement[attvals[i].hit][0]-=currmoveatt[0];
                        movement[attvals[i].hit][1]-=currmoveatt[1];
                    }//end for i
                }//end else weights==1
                for(i=elements;--i>=0;){
                    movement[i][0]/=elements;
                    movement[i][1]/=elements;
                    movement[i][2]=0;
                    totaldist=java.lang.Math.sqrt((movement[i][0]*movement[i][0])+(movement[i][1]*movement[i][1]));
                    if(totaldist>maxmove){
                        tmp=(float)(maxmove/totaldist);
                        movement[i][0]*=tmp;
                        movement[i][1]*=tmp;
                    }
                }//end for i
            }else{
                //cluster all in 3D
                if(weights!=null){
                    for(i=elements;--i>=0;){
                        weight1=weights[i];
                        movement[i][0]-=posarr[i][0]*minattract*weight1;
                        movement[i][1]-=posarr[i][1]*minattract*weight1;
                        movement[i][2]-=posarr[i][2]*minattract*weight1;
                        for(j=i;++j<elements;){
                            weight2=weights[j];
                            distx=posarr[j][0]-posarr[i][0];
                            disty=posarr[j][1]-posarr[i][1];
                            distz=posarr[j][2]-posarr[i][2];
                            if(distx==0 && disty==0 && distz==0){
                                //if two points are at exactly the same position I need to add some random effect
                                movement[i][0]+=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[i][1]+=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[i][2]+=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[j][0]-=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[j][1]-=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[j][2]-=repfac*(rand.nextDouble()-0.5)*0.001;
                                //return movement;
                            }else{
                                totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty)+(distz*distz));
                                //here I scale the force**repvalpow
                                totalmove=1;
                                for(int r=reppow;--r>=0;){
                                    totalmove*=totaldist;
                                }
                                totalmove=(repfac/totalmove);
                                movement[i][0]+=(-distx/totaldist)*totalmove*weight2;
                                movement[i][1]+=(-disty/totaldist)*totalmove*weight2;
                                movement[i][2]+=(-distz/totaldist)*totalmove*weight2;
                                movement[j][0]-=(-distx/totaldist)*totalmove*weight1;
                                movement[j][1]-=(-disty/totaldist)*totalmove*weight1;
                                movement[j][2]-=(-distz/totaldist)*totalmove*weight1;
                            }
                        }//end for j
                    }//end for i
                    for(i=attnum;--i>=0;){
                        hnum=attvals[i].hit;
                        qnum=attvals[i].query;
                        weight1=weights[hnum];
                        weight2=weights[hnum];
                        distx=posarr[hnum][0]-posarr[qnum][0];
                        disty=posarr[hnum][1]-posarr[qnum][1];
                        distz=posarr[hnum][2]-posarr[qnum][2];
                        totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty)+(distz*distz));
                        //scale totalmove with distance**attvalpow
                        totalmove=1;
                        for(int r=attvalpow;--r>=0;){
                            totalmove*=totaldist;
                        }
                        totalmove=totalmove*attvals[i].att*attfactor;
                        if(attvals[i].att<0){
                            //in cae of repulsion I want to react inversely to the attractive forces
                            totalmove=-1/totalmove;
                        }
                        if(totaldist!=0){
                            movement[qnum][0]+=(distx/totaldist)*totalmove*weight2;
                            movement[qnum][1]+=(disty/totaldist)*totalmove*weight2;
                            movement[qnum][2]+=(distz/totaldist)*totalmove*weight2;
                            movement[hnum][0]-=(distx/totaldist)*totalmove*weight1;
                            movement[hnum][1]-=(disty/totaldist)*totalmove*weight1;
                            movement[hnum][2]-=(distz/totaldist)*totalmove*weight1;
                        }//else{
                    }//end for i
                }else{//use weights==1
                    for(i=elements;--i>=0;){
                        movement[i][0]-=posarr[i][0]*minattract;
                        movement[i][1]-=posarr[i][1]*minattract;
                        movement[i][2]-=posarr[i][2]*minattract;
                        for(j=i;++j<elements;){
                            distx=posarr[j][0]-posarr[i][0];
                            disty=posarr[j][1]-posarr[i][1];
                            distz=posarr[j][2]-posarr[i][2];
                            if(distx==0 && disty==0 && distz==0){
                                //if two points are at exactly the same position I need to add some random effect
                                movement[i][0]+=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[i][1]+=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[i][2]+=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[j][0]-=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[j][1]-=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[j][2]-=repfac*(rand.nextDouble()-0.5)*0.001;
                                //return movement;
                            }else{
                                totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty)+(distz*distz));
                                //here I scale the force**repvalpow
                                totalmove=1;
                                for(int r=reppow;--r>=0;){
                                    totalmove*=totaldist;
                                }
                                totalmove=(repfac/totalmove);
                                movement[i][0]+=(-distx/totaldist)*totalmove;
                                movement[i][1]+=(-disty/totaldist)*totalmove;
                                movement[i][2]+=(-distz/totaldist)*totalmove;
                                movement[j][0]-=(-distx/totaldist)*totalmove;
                                movement[j][1]-=(-disty/totaldist)*totalmove;
                                movement[j][2]-=(-distz/totaldist)*totalmove;
                            }
                        }//end for j
                    }//end for i
                    for(i=attnum;--i>=0;){
                        hnum=attvals[i].hit;
                        qnum=attvals[i].query;
                        distx=posarr[hnum][0]-posarr[qnum][0];
                        disty=posarr[hnum][1]-posarr[qnum][1];
                        distz=posarr[hnum][2]-posarr[qnum][2];
                        totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty)+(distz*distz));
                        //scale totalmove with distance**attvalpow
                        totalmove=1;
                        for(int r=attvalpow;--r>=0;){
                            totalmove*=totaldist;
                        }
                        totalmove=totalmove*attvals[i].att*attfactor;
                        if(attvals[i].att<0){
                            //in cae of repulsion I want to react inversely to the attractive forces
                            totalmove=-1/totalmove;
                        }
                        if(totaldist!=0){
                            movement[qnum][0]+=(distx/totaldist)*totalmove;
                            movement[qnum][1]+=(disty/totaldist)*totalmove;
                            movement[qnum][2]+=(distz/totaldist)*totalmove;
                            movement[hnum][0]-=(distx/totaldist)*totalmove;
                            movement[hnum][1]-=(disty/totaldist)*totalmove;
                            movement[hnum][2]-=(distz/totaldist)*totalmove;
                        }
                    }//end for i
                }//end else weights!=null
                for(i=elements;--i>=0;){
                    movement[i][0]/=elements;
                    movement[i][1]/=elements;
                    movement[i][2]/=elements;
                    totaldist=java.lang.Math.sqrt((movement[i][0]*movement[i][0])+(movement[i][1]*movement[i][1])+(movement[i][2]*movement[i][2]));
                    if(totaldist>maxmove){
                        tmp=(float)(maxmove/totaldist);
                        movement[i][0]*=tmp;
                        movement[i][1]*=tmp;
                        movement[i][2]*=tmp;
                    }
                }//end for i
            }//end clustering all in 3D
        }//end else cluster only the selected sequences
    }// end getmovement
    */
    //-------------------------------------------
    /*
    //static double[] getattract2d(float[] pos1, float[] pos2, float attval, double[] movement, int attvalpow, float attfactor){
    static void getattract2d(float[] pos1, float[] pos2, float attval, double[] movement, int attvalpow, float attfactor){
        //get the attractive forces for 2d only (forget Z-axis)
        //tmpattvals are between 0 and 1 (or ==2 for evalue==0) (o=no attraction, 1=max attraction)
        double distx=pos2[0]-pos1[0];
        double disty=pos2[1]-pos1[1];
        double totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty));
        //scale totalmove with distance**attvalpow
        double totalmove=1;
        for(int i=attvalpow;--i>=0;){
            totalmove*=totaldist;
        }
        totalmove=totalmove*attval*attfactor;
        if(attval<0){
            //in cae of repulsion I want to react inversely to the attractive forces
            totalmove=-1/totalmove;
        }
        if(totaldist!=0){
            movement[0]=(distx/totaldist)*totalmove;
            movement[1]=(disty/totaldist)*totalmove;
            movement[2]=0;
        }else{
            //repulsion values will differentially move them
            movement[0]=0;
            movement[1]=0;
            movement[2]=0;
        }
        //return movement;
    }//end getattract2d
    */
    //--------------------------------------------------------------------------
    /*
    //static double[] getattract3d(float[] pos1, float[] pos2, float attval,double[] movement,int attvalpow, float attfactor){
    static void getattract3d(float[] pos1, float[] pos2, float attval,double[] movement,int attvalpow, float attfactor){
        //similar to getrepulse but this is an attractive force that scales with distance**2 (the further away, the greater)
        //which way is pos1 going to move, given pos2
        //tmpattvals are between 0 and 1 (or ==2 for evalue==0) (o=no attraction, 1=max attraction)
        double distx=pos2[0]-pos1[0];
        double disty=pos2[1]-pos1[1];
        double distz=pos2[2]-pos1[2];
        double totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty)+(distz*distz));
        //scale totalmove with distance**attvalpow
        double totalmove=1;
        for(int i=attvalpow;--i>=0;){
            totalmove*=totaldist;
        }
        totalmove=totalmove*attval*attfactor;
        if(attval<0){
            //in cae of repulsion I want to react inversely to the attractive forces
            totalmove=-1/totalmove;
        }
        if(totaldist!=0){
            movement[0]=(distx/totaldist)*totalmove;
            movement[1]=(disty/totaldist)*totalmove;
            movement[2]=(distz/totaldist)*totalmove;
        }else{
            //the repulsion values will differentially move them
            movement[0]=0;
            movement[1]=0;
            movement[2]=0;
        }
        //return movement;
    }// end getattract
    
    //-------------------------------------------
    
    //static double[] getrepulse2d(float[] pos1,float[] pos2,double[] movement, int repvalpow, float repfactor, Random rand){
    static void getrepulse2d(float[] pos1,float[] pos2,double[] movement, int repvalpow, float repfactor, Random rand){
        //get the repulsion in 2d only
        double distx=pos2[0]-pos1[0];
        double disty=pos2[1]-pos1[1];
        double totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty));
        if(totaldist==0){
            //if two points are at exactly the same position I need to add some random effect
            movement[0]=repfactor*(rand.nextDouble()-0.5)*0.001;
            movement[1]=repfactor*(rand.nextDouble()-0.5)*0.001;
            movement[2]=0;//repfactor*(rand.nextDouble()-0.5)*0.001;
            return;
        }
        //here I scale the force**repvalpow
        double totalmove=1;
        for(int i=repvalpow;--i>=0;){
            totalmove*=totaldist;
        }
        totalmove=repfactor/totalmove;
        movement[0]=(-distx/totaldist)*totalmove;
        movement[1]=(-disty/totaldist)*totalmove;
        movement[2]=0;
        //return movement;
    }//end getrepulse2d
    
    //--------------------------------------------------------------------------
    
    //static double[] getrepulse3d(float[] pos1,float[] pos2,double[] movement, int repvalpow, float repfactor, Random rand){
    static void getrepulse3d(float[] pos1,float[] pos2,double[] movement, int repvalpow, float repfactor, Random rand){
        //given these two objects, which way will object 1 move? force scales with 1/distance**2
        double distx=pos2[0]-pos1[0];
        double disty=pos2[1]-pos1[1];
        double distz=pos2[2]-pos1[2];
        double totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty)+(distz*distz));
        if(totaldist==0){
            //if two points are at exactly the same position I need to add some random effect
            movement[0]=repfactor*(rand.nextDouble()-0.5)*0.001;
            movement[1]=repfactor*(rand.nextDouble()-0.5)*0.001;
            movement[2]=repfactor*(rand.nextDouble()-0.5)*0.001;
            return;
        }
        //here I scale the force**repvalpow
        double totalmove=1;
        for(int i=repvalpow;--i>=0;){
            totalmove*=totaldist;
        }
        totalmove=repfactor/totalmove;
        movement[0]=(-distx/totaldist)*totalmove;
        movement[1]=(-disty/totaldist)*totalmove;
        movement[2]=(-distz/totaldist)*totalmove;
        //return movement;
    }// end getrepulse
    
    //--------------------------------------------------------------------------
    
    //static double[] getminattract2d(float[] pos1,double[] movement,double minattract){
    static void getminattract2d(float[] pos1,double[] movement,double minattract){
        //get minimum attraction for 2d only
        movement[0]=-pos1[0]*minattract;
        movement[1]=-pos1[1]*minattract;
        movement[2]=0;
        //return movement;
    }//end getminattract2d
    
    //--------------------------------------------------------------------------
    
    //static double[] getminattract(float[] pos1, double[] movement,double minattract){
    static void getminattract(float[] pos1, double[] movement,double minattract){
        //which way is pos1 going to move if it is attracted by the origin
        movement[0]=-pos1[0]*minattract;
        movement[1]=-pos1[1]*minattract;
        movement[2]=-pos1[2]*minattract;
        //return movement;
    }// end getminattract
    */
    //--------------------------------------------------------------------------
    /*
    float[][] seedstart(float[][] inarr){
        //seed the positions array with random numbers([-1 to 1[)
        int arrsize=java.lang.reflect.Array.getLength(inarr);
        if(arrsize==0){
            return inarr;
        }
        int j;
        int arrdepth=java.lang.reflect.Array.getLength(inarr[0]);
        for(int i=0;i<arrdepth;i++){
            for(j=0;j<arrsize;j++){
                inarr[j][i]=rand.nextFloat()*2-1;
            }// end for j
        }// end for i
        return inarr;
    }// end seedstart
    */
    //--------------------------------------------------------------------------
    /*
    minattvals[] getattvals(minattvals[] inattvals,double minpval){
        //use for filtering attraction values [-1<0<1]; -1 and 1 are max repulse/attract
        int attnum=java.lang.reflect.Array.getLength(inattvals);
        ArrayList <minattvals>retvec=new ArrayList<minattvals>();
        for(int i=0;i<attnum;i++){
            if((inattvals[i].att>=minpval) || (inattvals[i].att<=-minpval)){
                retvec.add(inattvals[i]);
            }
        }//end for i
        //minattvals[] retarr=new minattvals[retvec.size()];
        minattvals[] retarr=(minattvals[])retvec.toArray(new minattvals[0]);
        //retvec.copyInto(retarr);
        return retarr;
    }//end getattvals
    
    //--------------------------------------------------------------------------
    
    minattvals[] getattvals(minhsp[] indata,double minpval){
        if(indata==null){//possible (if alternate data source was loaded)
            System.out.println("indata is null");
            return myattvals;
        }
        System.out.println("indata is size:"+java.lang.reflect.Array.getLength(indata));
        //compute the attraction values for all sequence pairs
        //values range from 0(no attraction) to 1(max); -1 denotes identity
        //indata is composed of one array of hsp objects
        //NOTE: possible multiple pvalues per hsp object (multiple hsp's for same sequence pair)
        int j;
        ArrayList<minattvals> tmpvec=new ArrayList<minattvals>();
        int datanum=java.lang.reflect.Array.getLength(indata);
        HashMap myhash=new HashMap(datanum);
        float newatt;
        String key;
        float maxattval=0;
        if(attvalcompcheckbox.isSelected()){
            attvalsimple=false;
        }else{
            attvalsimple=true;
        }
        minattvals curratt=null;
        maxvalfound=0;//init to zero, is assigned value in getattvalsimple or mult
        //NOTE: this is not necessarily a symmetrical array. compute all values
        //and then symmetrize computing the average values
        if(rescalepvaluescheckbox.isSelected()==false){
            //make the attraction values
            if(attvalsimple){
                for(int i=datanum;--i>=0;){
                    if(indata[i].query<indata[i].hit){
                        key=indata[i].query+"_"+indata[i].hit;
                    }else{
                        key=indata[i].hit+"_"+indata[i].query;
                    }
                    if(myhash.containsKey(key)){
                        curratt=(minattvals)myhash.get(key);
                        if(curratt.att==-1){
                            //in this case keep the -1
                        }else{
                            newatt=getattvalsimple(indata[i].val,elements,minpval);
                            if(newatt==-1){
                                curratt.att=-1;
                            }else{
                                newatt/=2;
                                curratt.att+=newatt;
                            }
                        }
                    }else{
                        //if I've never encountered this query-hit pair before
                        curratt=new minattvals();
                        if(indata[i].query<indata[i].hit){
                            curratt.query=indata[i].query;
                            curratt.hit=indata[i].hit;
                        }else{
                            curratt.hit=indata[i].query;
                            curratt.query=indata[i].hit;
                        }
                        curratt.att=getattvalsimple(indata[i].val,elements,minpval);
                        if(curratt.att!=-1){
                            curratt.att/=2;
                        }
                        if(curratt.att!=0){
                            myhash.put(key,curratt);
                            tmpvec.add(curratt);
                        }
                    }
                    if(curratt.att>maxattval){
                        maxattval=curratt.att;
                    }
                }//end for i
            }else{
                for(int i=0;i<datanum;i++){
                    if(indata[i].query<indata[i].hit){
                        key=indata[i].query+"_"+indata[i].hit;
                    }else{
                        key=indata[i].hit+"_"+indata[i].query;
                    }
                    if(myhash.containsKey(key)){
                        curratt=(minattvals)myhash.get(key);
                        if(curratt.att==-1){
                            //in this case keep the -1
                        }else{
                            newatt=getattvalmult(indata[i].val,elements,minpval);
                            if(newatt==-1){
                                curratt.att=-1;
                            }else{
                                newatt/=2;
                                curratt.att+=newatt;
                            }
                        }
                    }else{
                        //if I've never encountered this query-hit pair before
                        curratt=new minattvals();
                        if(indata[i].query<indata[i].hit){
                            curratt.query=indata[i].query;
                            curratt.hit=indata[i].hit;
                        }else{
                            curratt.hit=indata[i].query;
                            curratt.query=indata[i].hit;
                        }
                        curratt.att=getattvalmult(indata[i].val,elements,minpval);
                        if(curratt.att !=-1){
                            curratt.att/=2;
                        }
                        if(curratt.att!=0){
                            tmpvec.add(curratt);
                            myhash.put(key,curratt);
                        }
                    }
                    if(curratt.att>maxattval){
                        maxattval=curratt.att;
                    }
                }//end for i
            }
            //divide all vals by maxattval (-->range: 0-1)
            //standard, just divide all values by the maximum value
            //note, this does NOT symmetrize the attractions
            if(usescval==false){
                for(int i=tmpvec.size()-1;i>=0;i--){
                    if(((minattvals)tmpvec.get(i)).att==-1){
                        ((minattvals)tmpvec.get(i)).att=1;
                        //System.out.println(((minattvals)tmpvec.elementAt(i)).query+" "+((minattvals)tmpvec.elementAt(i)).hit+" :"+((minattvals)tmpvec.elementAt(i)).att);
                    }else{
                        ((minattvals)tmpvec.get(i)).att/=maxattval;
                        //System.out.println(((minattvals)tmpvec.elementAt(i)).query+" "+((minattvals)tmpvec.elementAt(i)).hit+" :"+((minattvals)tmpvec.elementAt(i)).att);
                    }
                }// end for i
                //System.out.println("maxattval"+maxattval+" offset="+0);
                p2attfactor=maxattval;
                p2attoffset=0;
            }else{//if using scval
                p2attfactor=1;
                p2attoffset=0;
            }
        }else{//if rescalepvaluecheckbox==true
            float minattval=java.lang.Float.MAX_VALUE;
            //rescale the attraction values to range from 0 to 1 (with the smallest positive non-zero value as zero.
            if(attvalsimple){
                for(int i=0;i<datanum;i++){
                    if(indata[i].query<indata[i].hit){
                        key=indata[i].query+"_"+indata[i].hit;
                    }else{
                        key=indata[i].hit+"_"+indata[i].query;
                    }
                    if(myhash.containsKey(key)){
                        curratt=(minattvals)myhash.get(key);
                        if(curratt.att==-1){
                            //in this case keep the -1
                        }else{
                            newatt=getattvalsimple(indata[i].val,elements,minpval);
                            if(newatt==-1){
                                curratt.att=-1;
                            }else{
                                newatt/=2;
                                curratt.att+=newatt;
                            }
                        }
                    }else{
                        //if I've never encountered this query-hit pair before
                        curratt=new minattvals();
                        if(indata[i].query<indata[i].hit){
                            curratt.query=indata[i].query;
                            curratt.hit=indata[i].hit;
                        }else{
                            curratt.hit=indata[i].query;
                            curratt.query=indata[i].hit;
                        }
                        curratt.att=getattvalsimple(indata[i].val,elements,minpval);
                        if(curratt.att!=-1){
                            curratt.att/=2;
                        }
                        if(curratt.att!=0){
                            myhash.put(key,curratt);
                            tmpvec.add(curratt);
                        }
                    }
                    if(curratt.att>maxattval){
                        maxattval=curratt.att;
                    }
                    if((curratt.att>0)&&(curratt.att<minattval)){
                        minattval=curratt.att;
                    }
                }//end for i
            }else{
                for(int i=0;i<datanum;i++){
                    if(indata[i].query<indata[i].hit){
                        key=indata[i].query+"_"+indata[i].hit;
                    }else{
                        key=indata[i].hit+"_"+indata[i].query;
                    }
                    if(myhash.containsKey(key)){
                        curratt=(minattvals)myhash.get(key);
                        if(curratt.att==-1){
                            //in this case keep the -1
                        }else{
                            newatt=getattvalmult(indata[i].val,elements,minpval);
                            if(newatt==-1){
                                curratt.att=-1;
                            }else{
                                newatt/=2;
                                curratt.att+=newatt;
                            }
                        }
                        
                    }else{
                        //if I've never encountered this query-hit pair before
                        curratt=new minattvals();
                        if(indata[i].query<indata[i].hit){
                            curratt.query=indata[i].query;
                            curratt.hit=indata[i].hit;
                        }else{
                            curratt.hit=indata[i].query;
                            curratt.query=indata[i].hit;
                        }
                        curratt.att=getattvalmult(indata[i].val,elements,minpval);
                        if(curratt.att!=-1){
                            curratt.att/=2;
                        }
                        if(curratt.att!=0){
                            myhash.put(key,curratt);
                            tmpvec.add(curratt);
                        }
                        
                    }
                    if(curratt.att>maxattval){
                        maxattval=curratt.att;
                    }
                    if((curratt.att>0)&&(curratt.att<minattval)){
                        minattval=curratt.att;
                    }
                }//end for i
            }
            //and divide all vals by maxattval and offset by minattval(-->range: 0-1)
            float divval=maxattval-minattval;
            for(int i=tmpvec.size()-1;i>=0;i--){
                if(((minattvals)tmpvec.get(i)).att==-1){
                    ((minattvals)tmpvec.get(i)).att=1;
                }else{
                    ((minattvals)tmpvec.get(i)).att=(((minattvals)tmpvec.get(i)).att-minattval)/divval;
                }
            }// end for i
            //System.out.println("maxattval"+maxattval+" offset="+minattval);
            p2attfactor=divval;
            p2attoffset=minattval;
        }
        if(usescval){
            minpvalbutton.setText("Use S/C values better than:");
        }
        //minattvals[] retarr=new minattvals[tmpvec.size()];
        //tmpvec.copyInto(retarr);
        minattvals[] retarr=(minattvals[])tmpvec.toArray(new minattvals[0]);
        System.out.println("attvals size="+java.lang.reflect.Array.getLength(retarr));
        //for(int i=java.lang.reflect.Array.getLength(retarr);--i>=0;){
        //    System.out.println("query="+retarr[i].query+" hit="+retarr[i].hit+" val="+retarr[i].att);
        //}//end for i
        return retarr;
    }// end getattvals
    */
    //------------------------------
    /*
    float getattvalsimple(double[] invec,int dbsize,double minpval){
        //System.out.println("simple");
        //this actually takes all data from vector and makes ONE number out of it
        //just use the BEST value!
        if(invec==null){
            return 0;
        }else if(java.lang.reflect.Array.getLength(invec)<1){//if I have no hits
            return 0;
        }
        double bestval=invec[0];//is a p-value (should be from 0 to 1)
        double currval;
        if(usescval){
            bestval=0;
            for(int i=java.lang.reflect.Array.getLength(invec)-1;i>=0;i--){
                currval=invec[i];
                if(currval>bestval){
                    bestval=currval;
                }
            }// end for i
            if(bestval<maxvalfound){//maxvalfound=worst accepted value (comes from P-values where larger=worse)
                maxvalfound=bestval;
            }
            currval=bestval;
            if(currval<minpval){//minpval also functions as minscoreval here
                //System.out.println(" currval="+currval+" is less than minpval="+minpval+" returning zero");
                return 0;
            }
        }else{
            for(int i=java.lang.reflect.Array.getLength(invec)-1;i>=1;i--){
                currval=invec[i];
                if(currval<bestval){
                    bestval=currval;
                }
            }// end for i
            if(bestval>maxvalfound){
                maxvalfound=bestval;
            }
            if(bestval==0){
                return -1;//this is identity
            }else if(bestval>1){//should never happen to p-values!
                return 0;
            }else if(bestval>minpval){//if this value is worse than permitted
                return 0;
            }
            //now all pvalues between 0 and 1
            currval=(-1*java.lang.Math.log(bestval));//ln10;//don't need it here as I convert all to relative attractions
        }
        return (float)(currval);
    }// end gettaval
    
    float getattvalmult(double[] invec,int dbsize,double minpval){
        //System.out.println("complex");
        //this actually takes all data from vector and makes ONE number out of it
        //new: multiply the pvalues of different hsp's
        if(invec==null){
            return 0;
        }else if(java.lang.reflect.Array.getLength(invec)<1){
            return 0;
        }
        double currval=invec[0];
        if(usescval){
            //then I am using score values (no logarithming at the end)
            currval=0;
            for(int i=java.lang.reflect.Array.getLength(invec);--i>=0;){//sum the values
                currval+=invec[i];
            }// end for i
            if(currval<maxvalfound){//maxvalfound=worst accepted value (comes from P-values where larger=worse)
                maxvalfound=currval;
            }
            if(currval<minpval){//minpval also functions as minscoreval here
                return 0;
            }
        }else{//then I am using P-values
            for(int i=java.lang.reflect.Array.getLength(invec);--i>=1;){
                currval*=invec[i];
            }// end for i
            if(currval>maxvalfound){
                maxvalfound=currval;
            }
            if(currval==0){
                return -1;//this is identity
            }else if(currval>1){//should never happen to p-values!
                return 0;
            }else if(currval>minpval){//if this value is worse than permitted
                return 0;
            }
            //now all pvalues between 0 and 1
            currval=(-1*java.lang.Math.log(currval));///ln10;
        }
        return (float)(currval);
    }// end gettaval
    */
    //--------------------------------------------------------------------------
    
    void updateselected(int[] tmpreg,boolean deselect){
        //if deselect==true I want to deselect all the sequences in the selectes region
        //as input this has 4 coordinates
        //what this does is to take all sequences that are within those coordinates
        //and select/deselect them. It also needs to update the selection in shownamedialog
        ArrayList <Integer> tmpselect=new ArrayList<Integer>();
        int minx=tmpreg[0];
        int maxx=tmpreg[2];
        int miny=tmpreg[1];
        int maxy=tmpreg[3];
        int tmpint;
        if(minx>maxx){
            tmpint=maxx;
            maxx=minx;
            minx=tmpint;
        }
        if(miny>maxy){
            tmpint=maxy;
            maxy=miny;
            miny=tmpint;
        }
        //now get all sequences that are within those 4 coords
        float[][] currpos=data.posarrtmp;
        for(int i=0;i<data.elements;i++){
            //see if this object is within the selected region
            if((currpos[i][0]>minx)&&(currpos[i][0]<maxx)&&(currpos[i][1]>miny)&&(currpos[i][1]<maxy)){
                tmpselect.add(new Integer(i));
            }
        }// end for i
        //now select those that are not present in selectednames and deselect those that are
        int tmparrsize=tmpselect.size();
        Integer[] tmparr=(Integer[]) tmpselect.toArray(new Integer[0]);
        tmpselect.clear();
        int j;
        int selectedelements=java.lang.reflect.Array.getLength(data.selectednames);
        boolean isnew;
        for(int i=0;i<tmparrsize;i++){
            tmpint=tmparr[i].intValue();
            isnew=true;
            for(j=0;j<selectedelements;j++){
                if(data.selectednames[j]==tmpint){//if this element was already selected
                    isnew=false;//do not add it to the new Vector
                    data.selectednames[j]=-1;//mark this as unselected
                }
            }// end for j
            if(isnew && (deselect==false)){//if this element should be selected
                tmpselect.add(tmparr[i]);
            }
        }// end for i
        for(int i=0;i<selectedelements;i++){
            //now add those that were selected before
            if(data.selectednames[i]>-1){
                tmpselect.add(new Integer(data.selectednames[i]));
            }
        }//end for i
        //now assign the new selecteds to selectelements and update shownamedialog
        data.selectednames=new int[tmpselect.size()];
        for(int i=tmpselect.size()-1;i>=0;i--){
            data.selectednames[i]=((Integer)tmpselect.get(i)).intValue();
        }// end for i
        if(shownames!=null){
            shownames.setselected(data.selectednames,(shownames.showall==false));
        }
        setclearbuttontext();
    }// end updateselected
    
    void updatetmpselected(int[] tmpreg){
        //as input this has 4 coordinates
        //what this does is to take all sequences that are within those coordinates
        //and select/deselect them. It also needs to update the selection in shownamedialog
        ArrayList <Integer>tmpselect=new ArrayList<Integer>();
        int minx=tmpreg[0];
        int maxx=tmpreg[2];
        int miny=tmpreg[1];
        int maxy=tmpreg[3];
        int tmpint;
        if(minx>maxx){
            tmpint=maxx;
            maxx=minx;
            minx=tmpint;
        }
        if(miny>maxy){
            tmpint=maxy;
            maxy=miny;
            miny=tmpint;
        }
        //now get all sequences that are within those 4 coords
        float[][] currpos=data.posarrtmp;
        for(int i=0;i<data.elements;i++){
            //see if this object is within the selected region
            if((currpos[i][0]>minx)&&(currpos[i][0]<maxx)&&(currpos[i][1]>miny)&&(currpos[i][1]<maxy)){
                tmpselect.add(new Integer(i));
            }
        }// end for i
        //now select those that are not present in selectednames and deselect those that are
        int j;
        int selectedelements=java.lang.reflect.Array.getLength(data.selectednames);
        int[] tmpdeselected=new int[selectedelements];
        boolean isnew;
        for(int i=0;i<tmpselect.size();i++){
            tmpint=((Integer)tmpselect.get(i)).intValue();
            isnew=true;
            for(j=0;j<selectedelements;j++){
                if(data.selectednames[j]==tmpint){
                    isnew=false;
                    tmpdeselected[j]=1;
                    break;
                }
            }//end for j
            if(isnew==false){
                tmpselect.remove(i);
                i--;
            }
        }//end for loop through vector
        //now add the already selected names
        for(int i=selectedelements-1;i>=0;i--){
            //now add those that were selected before
            if(tmpdeselected[i]!=1){
                tmpselect.add(0,new Integer(data.selectednames[i]));
            }
        }//end for i
        //now assign the new selecteds to tmpselectelements and update shownamedialog
        int[] tmpselectednames=new int[tmpselect.size()];
        for(int i=tmpselect.size()-1;i>=0;i--){
            tmpselectednames[i]=((Integer)tmpselect.get(i)).intValue();
        }// end for i
        if(shownames!=null){
            shownames.setselected(tmpselectednames,true);
        }
        //setclearbuttontext();
    }// end updatetmpselected
    
    //--------------------------------------------------------------------------
    
    void moveselected(int xint, int yint){
        //move selected sequences by movex and movy in 3d space
        //use the rotmtx to compute the x,y and z coord to move them by
        float x=(float)((float)xint/draw1.xscale);
        float y=(float)((float)yint/draw1.yscale);
        double[][]invrot=getinvrot(data.rotmtx);
        float[] move=new float[3];
        move[0]=(float)(invrot[0][0]*x+invrot[0][1]*y);
        move[1]=(float)(invrot[1][0]*x+invrot[1][1]*y);
        move[2]=(float)(invrot[2][0]*x+invrot[2][1]*y);
        //now move each of the selected seqs by that amount
        int selectedelements=java.lang.reflect.Array.getLength(data.selectednames);
        //for each of the selected sequences shift it's position in posarr by x,y,z coords.
        //selectednames is a int[] with the selected sequences indices
        for(int i=0;i<selectedelements;i++){
            data.myposarr[data.selectednames[i]][0]+=move[0];
            data.myposarr[data.selectednames[i]][1]+=move[1];
            data.myposarr[data.selectednames[i]][2]+=move[2];
        }// end for i
    }//end moveselected
    
    static double[][] getinvrot(double[][] m){
        //get the inverse of a 3x3 matrix inmtx
        double[][] retmtx={
            {((m[1][1]*m[2][2])-(m[2][1]*m[1][2])),-((m[0][1]*m[2][2])-(m[2][1]*m[0][2])),((m[0][1]*m[1][2])-(m[1][1]*m[0][2]))},
            {-((m[1][0]*m[2][2])-(m[2][0]*m[1][2])),((m[0][0]*m[2][2])-(m[2][0]*m[0][2])),-((m[0][0]*m[1][2])-(m[1][0]*m[0][2]))},
            {((m[1][0]*m[2][1])-(m[2][0]*m[1][1])),-((m[0][0]*m[2][1])-(m[2][0]*m[0][1])),((m[0][0]*m[1][1])-(m[1][0]*m[0][1]))}
        };
        return retmtx;
    }//end getinvrot
    
    //--------------------------------------------------------------------------
    
    void setclearbuttontext(){
        //if I have some sequences selected --. text=Clear selection
        if(java.lang.reflect.Array.getLength(data.selectednames)>0){
            button_select_all_or_clear.setText("Clear Selection");
        }else{
            button_select_all_or_clear.setText("Select All");
        }
    }//end setclearbuttontext
    
    //--------------------------------------------------------------------------
 /*
    Vector[][] gettmpblasthits(Vector[][] blasthits,int[] selectednames,double minpval){
        //only use hsp's coming from or pertaining to selected seqs
        int elements=java.lang.reflect.Array.getLength(blasthits);
        int selectedelements=java.lang.reflect.Array.getLength(selectednames);
        Vector[][] retarr=new Vector[elements][elements];
        Vector empty=new Vector();
        int j;
        int k;
        boolean isselected;
        for(int i=0;i<elements;i++){
            isselected=false;
            for(j=0;j<selectedelements;j++){
                if(i==selectednames[j]){
                    isselected=true;
                    break;
                }// end if
            }// end for j
            if(isselected){//if this is one of the selected sequences
                //add all the blast hits
                for(j=0;j<elements;j++){
                    if(blasthits[i][j]==null){
                        continue;
                    }
                    retarr[i][j]=new Vector();
                    for(k=0;k<blasthits[i][j].size();k++){
                        //check if this hsp is better than minpval
                        if(((proghsp)blasthits[i][j].elementAt(k)).pvalue<=minpval){
                            retarr[i][j].addElement(blasthits[i][j].elementAt(k));
                        }
                    }//end for k
                }// end for j
            }else{//if this is not one of the selected sequences
                //only add those blast hits going to selected sequences
                //initialize all to empty vector
                for(j=0;j<elements;j++){
                    retarr[i][j]=empty;
                }// end for j
                //now add the selected hits
                for(j=0;j<selectedelements;j++){
                    if(blasthits[i][j]==null){
                        continue;
                    }
                    retarr[i][selectednames[j]]=new Vector();
                    for(k=0;k<blasthits[i][selectednames[j]].size();k++){
                        //check if this hsp is better than minpval
                        if(((proghsp)blasthits[i][selectednames[j]].elementAt(k)).pvalue<=minpval){
                            retarr[i][selectednames[j]].addElement(blasthits[i][selectednames[j]].elementAt(k));
                        }
                    }//end for k
                }// end for j
            }// end else isselected
        }// end for i
        return retarr;
    }// end gettmpblasthits
  */
    
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    
    class drawpanel extends JPanel implements java.awt.print.Printable{
        
        public void drawpanel(){
            this.setOpaque(true);//should be set by default, but can't hurt to re-set
            this.setDoubleBuffered(true);//should be set by default, but can't hurt to re-set
        }
        
        void init(){
            setSize(graphpanel.getWidth(),graphpanel.getHeight());
            setBackground(new java.awt.Color(0.95f,0.95f,0.95f));
            data.posarr=new float[0][0];
            data.colorarr=makecolors(colornum);
            data.colorcutoffs=new float[colornum];
            if(data.usescval==true){
                for(int i=0;i<colornum;i++){
                    data.colorcutoffs[i]=(float)(((float)i/(float)(colornum))*data.p2attfactor);
                }
            }else{
                for(int i=0;i<colornum;i++){
                    data.colorcutoffs[i]=((float)i/(float)(colornum));
                }
            }
            refdbstring="";
            for(int i=0;i<java.lang.reflect.Array.getLength(data.referencedb);i++){
                refdbstring+=data.referencedb[i]+"; ";
            }//end for i
            //set the drawbox coordinates (a box of size 1 in x,y,z space)
            box=new double[8][3];
            box[0][0]=0;
            box[0][1]=0;
            box[0][2]=0;
            
            box[1][0]=1;
            box[1][1]=0;
            box[1][2]=0;
            
            box[2][0]=1;
            box[2][1]=1;
            box[2][2]=0;
            
            box[3][0]=0;
            box[3][1]=1;
            box[3][2]=0;
            
            box[4][0]=0;
            box[4][1]=0;
            box[4][2]=1;
            
            box[5][0]=1;
            box[5][1]=0;
            box[5][2]=1;
            
            box[6][0]=1;
            box[6][1]=1;
            box[6][2]=1;
            
            box[7][0]=0;
            box[7][1]=1;
            box[7][2]=1;
            
            boxdata=new double[8][3];
            originpos[0]=0;
            originpos[1]=0;
        }
        
        int[] originpos=new int[2];
        int colornum=10;
        //float[] colorcutoffs=new float[0];
        boolean drawbox=false;
        //float[][] posarr =new float[0][0];
        //float[][] posarrtmp={{1,1}};
        //int[][] drawarrtmp={{1,1}};
        HashMap frustration;//color code for the frustration of connections
        //java.awt.Color[] colorarr=new java.awt.Color[0];//makecolors(colornum);//make 10+1 colors in colorarr
        java.awt.Color[] frustcolorarr=new java.awt.Color[0];
        java.awt.Color bgcolor=new java.awt.Color(1f,1f,1f);
        java.awt.Color fgcolor=new java.awt.Color(0,0,0);
        java.awt.Color blastcirclecolor=new java.awt.Color(0,1f,0);
        //Vector[] draworder=new Vector[0];
        //ArrayList[] draworder=new ArrayList[0];//ArrayList is the same as a Vector but unsynchronized (faster)
        int xadd=25;
        int yadd=25;
        int drawwidth;
        int drawheight;
        int stereoangle=4;
        double[][] box;
        double[][] boxdata;
        double sinxz;
        double cosxz;
        double sinyz;
        double cosyz;
        double sinxy;
        double cosxy;
        double alphayz;
        double alphaxz;
        double alphaxy;
        double xscale=1;
        double yscale=1;
        boolean rotatez=false;
        double maxx=1;
        double maxy=1;
        double minx=1;
        double miny=1;
        double xmove=0;
        double ymove=0;
        double zmove=0;
        
        int[] draw=new int[4];//the coordinates of the field selected by the mouse
        java.awt.Color selectedcolor=new java.awt.Color(1f,0f,0f);
        java.awt.Font myfont=this.getFont();
        java.awt.Font smallfont=new java.awt.Font("Monospace", 0, 9);
        int fontsize=myfont.getSize();
        int fontwidth=(int)(fontsize/2);
        //int ovalsize=10;
        //int groupsize=4;
        //int dotsize=2;
        String refdbstring="";
        java.awt.Color blasthitcolor=new java.awt.Color(0.4f,0.4f,0.4f);
        selectclass tmpselectclass;
        //float zoomfactor=1;
        int xtranslate=0;
        int ytranslate=0;
        //ArrayList<int[][]> polygons=makepolygons.get(groupsize);
        int[][] tmpposarr;
        int[] xposarr;
        int[] yposarr;
        int posnum;
        
        //double[][] rotmtx={{1,0,0},{0,1,0},{0,0,1}};//the performed rotations
        double[][] tmprotmtx={{1,0,0},{0,1,0},{0,0,1}};//new double[3][3];//the current mouse rotation
        //double[][] myrotmtx={{1,0,0},{0,1,0},{0,0,1}};//new double[3][3];//both of the above together
        //----------------------------------------------------------------------
        
        java.awt.Color[] makecolors(int num){
            java.awt.Color[] retarr=new java.awt.Color[num];//was num+1
            int worstred=230;
            int worstgreen=230;
            int worstblue=230;
            int bestred=0;
            int bestgreen=0;
            int bestblue=0;
            //now compute the stepwise gradient
            float redstep=((float)(bestred-worstred))/((float)num);//were all three num+1
            float greenstep=((float)(bestgreen-worstgreen))/((float)num);
            float bluestep=((float)(bestblue-worstblue))/((float)num);
            for(int i=0;i<num;i++){
                retarr[i]=new java.awt.Color((int)(worstred+(i*redstep)),(int)(worstgreen+(i*greenstep)),(int)(worstblue+(i*bluestep)));
            }
            //retarr[num-1]=new java.awt.Color(bestred,bestgreen,bestblue);
            return retarr;
        }// end makecolors
        
        //----------------------------------------------------------------------
        
        public int print(java.awt.Graphics g, java.awt.print.PageFormat pf, int pi) throws java.awt.print.PrinterException {
            double xscalel=pf.getImageableWidth()/((double)graphpanel.getWidth());
            double yscalel=pf.getImageableHeight()/((double)graphpanel.getHeight());
            if (pi >= 1) {
                return java.awt.print.Printable.NO_SUCH_PAGE;
            }
            java.awt.Graphics2D g2d=(java.awt.Graphics2D) g;
            java.awt.geom.AffineTransform trans=java.awt.geom.AffineTransform.getScaleInstance(xscalel,yscalel);
            g2d.translate(pf.getImageableX(), pf.getImageableY());
            g2d.scale(xscalel,yscalel);
            graphpanel.setDoubleBuffered(false);
            paintComponent(g2d);
            graphpanel.setDoubleBuffered(true);
            return java.awt.print.Printable.PAGE_EXISTS;
        }
        
        //----------------------------------------------------------------------
        
        public void paintComponent(java.awt.Graphics g){
            //super.paintComponent(g);
            //long time=System.currentTimeMillis();
            g.setFont(myfont);
            if(stereocheckboxmenuitem.isSelected()){
                //draw the graph in stereo vision
                //draw the first picture as normal
                drawwidth=(int)((graphpanel.getWidth()-xadd-xadd)/2);
                drawheight=(int)(graphpanel.getHeight()-yadd-yadd);
                g.translate(xtranslate,ytranslate);
                g.setColor(bgcolor);
                g.fillRect(xadd-xtranslate,yadd-ytranslate,drawwidth+drawwidth,drawheight);
                drawwidth*=data.zoomfactor;
                drawheight*=data.zoomfactor;
                g.setColor(fgcolor);
                if(repaint!=null){//error loadign or similar
                    g.drawString(repaint,(drawwidth/2)-xtranslate, (drawheight/2)-ytranslate);
                }else{
                    //now rotate the graph by the stereoangle along the x-axis
                    mousemove[0]+=((float)stereoangle/360f)*drawwidth;
                    getangles();
                    if(drawbox==false){
                        drawdata(g);
                    }else{
                        drawdata(g);
                        drawbox(g);
                    }
                    //now rotate back to what it was before
                    mousemove[0]-=((float)stereoangle/360f)*(drawwidth);
                }
                if(showorigcheckbox.isSelected()){
                    g.setColor(java.awt.Color.red);
                    g.drawString("X(0,0,0)", originpos[0], originpos[1]);
                    //g.fillRect(originpos[0], originpos[1],5,5);
                }
                g.setColor(java.awt.Color.black);
                if(showinfocheckbox.isSelected()){
                    g.drawString(" INFO: blast="+data.blastpath+" refdb="+refdbstring,(colornum+9)*fontwidth-xtranslate,fontsize-ytranslate);
                    g.setFont(smallfont);
                    g.drawString("x:"+minx+" to "+maxx,fontwidth-xtranslate,graphpanel.getHeight()-fontsize-3-ytranslate);
                    g.drawString("y:"+miny+" to "+maxy,fontwidth-xtranslate,graphpanel.getHeight()-ytranslate-3);
                    g.drawString((float)((int)(data.myrotmtx[0][0]*100))/100+","+(float)((int)(data.myrotmtx[0][1]*100))/100+","+(float)((int)(data.myrotmtx[0][2]*100))/100,xadd-xtranslate,yadd+fontsize-ytranslate);
                    g.drawString((float)((int)(data.myrotmtx[1][0]*100))/100+","+(float)((int)(data.myrotmtx[1][1]*100))/100+","+(float)((int)(data.myrotmtx[1][2]*100))/100,xadd-xtranslate,yadd+2*fontsize-ytranslate);
                    g.drawString((float)((int)(data.myrotmtx[2][0]*100))/100+","+(float)((int)(data.myrotmtx[2][1]*100))/100+","+(float)((int)(data.myrotmtx[2][2]*100))/100,xadd-xtranslate,yadd+3*fontsize-ytranslate);
                }
                g.drawString("Round: "+data.rounds,(int)(graphpanel.getWidth()/2)-xtranslate,graphpanel.getHeight()-fontsize-ytranslate);
                //docalc="true";
                recalc=true;
                //and draw the second image
                drawwidth=(int)((graphpanel.getWidth()-xadd-xadd)/2);
                g.translate(xtranslate+drawwidth,ytranslate);
                drawwidth*=data.zoomfactor;
                g.setColor(fgcolor);
                if(repaint!=null){
                    //g.drawString(repaint,(drawwidth/2)-xtranslate+drawwidth, (drawheight/2)-ytranslate);
                }else{
                    getangles();
                    if(drawbox==false){
                        drawdata(g);
                    }else{
                        drawdata(g);
                        drawbox(g);
                    }
                }
            }else{//don't draw in stereo
                drawwidth=(int)(graphpanel.getWidth()-xadd-xadd);
                drawheight=(int)(graphpanel.getHeight()-yadd-yadd);
                g.translate(xtranslate,ytranslate);
                g.setColor(bgcolor);
                g.fillRect(xadd-xtranslate,yadd-ytranslate,drawwidth,drawheight);
                drawwidth*=data.zoomfactor;
                drawheight*=data.zoomfactor;
                g.setColor(fgcolor);
                if(repaint!=null){
                    g.drawString(repaint,(drawwidth/2)-xtranslate, (drawheight/2)-ytranslate);
                }else{
                    getangles();
                    if(drawbox==false){
                        drawdata(g);
                    }else{
                        drawdata(g);
                        drawbox(g);
                    }
                }
                if(showorigcheckbox.isSelected()){
                    g.setColor(java.awt.Color.red);
                    g.drawString("X(0,0,0)", originpos[0], originpos[1]);
                    //g.fillRect(originpos[0], originpos[1],5,5);
                }
                g.setColor(java.awt.Color.black);
                if(showinfocheckbox.isSelected()){
                    g.drawString(" INFO: blast="+data.blastpath+" refdb="+refdbstring,(colornum+9)*fontwidth-xtranslate,fontsize-ytranslate);
                    g.setFont(smallfont);
                    g.drawString("x:"+minx+" to "+maxx,fontwidth-xtranslate,graphpanel.getHeight()-fontsize-3-ytranslate);
                    g.drawString("y:"+miny+" to "+maxy,fontwidth-xtranslate,graphpanel.getHeight()-ytranslate-3);
                    g.drawString((float)((int)(data.myrotmtx[0][0]*100))/100+","+(float)((int)(data.myrotmtx[0][1]*100))/100+","+(float)((int)(data.myrotmtx[0][2]*100))/100,xadd-xtranslate,yadd+fontsize-ytranslate);
                    g.drawString((float)((int)(data.myrotmtx[1][0]*100))/100+","+(float)((int)(data.myrotmtx[1][1]*100))/100+","+(float)((int)(data.myrotmtx[1][2]*100))/100,xadd-xtranslate,yadd+2*fontsize-ytranslate);
                    g.drawString((float)((int)(data.myrotmtx[2][0]*100))/100+","+(float)((int)(data.myrotmtx[2][1]*100))/100+","+(float)((int)(data.myrotmtx[2][2]*100))/100,xadd-xtranslate,yadd+3*fontsize-ytranslate);
                }
                g.drawString("Round: "+data.rounds,(int)(graphpanel.getWidth()/2)-xtranslate,graphpanel.getHeight()-fontsize-ytranslate);
                //docalc="true";
                recalc=true;
            }
            //System.out.println("drawtime="+(System.currentTimeMillis()-time));
        }// end paintcomponent
        
        //----------------------------------------------------------------------
        
        void getangles(){
            //set alphaxy alphaxz and alphayz and their sin and cos (the rotation angles)
            //first compute the angle xy by taking the last rotaion angles and the last mouse movement
            //tmpvec contains movement of x,y,and z by the mouse)
            double[] tmpvec={java.lang.Math.PI/2*((double)mousemove[0]/(drawwidth/2)),java.lang.Math.PI/2*((double)mousemove[1]/(drawwidth/2)),0};
            xmove=tmpvec[0]+0.001;//avoid div0 errors
            ymove=tmpvec[1]+0.001;
            //if rotatez the calculate the rotation angle for z axis
            if((rotatez==true)||(cluster2dbutton.isSelected())){
                zmove=ymove+xmove;
                tmprotmtx[0][0]=java.lang.Math.cos(zmove);
                tmprotmtx[0][1]=-java.lang.Math.sin(zmove);
                tmprotmtx[0][2]=0;
                tmprotmtx[1][0]=java.lang.Math.sin(zmove);
                tmprotmtx[1][1]=java.lang.Math.cos(zmove);
                tmprotmtx[1][2]=0;
                tmprotmtx[2][0]=0;
                tmprotmtx[2][1]=0;
                tmprotmtx[2][2]=1;
            }else{//get the parameters for the xy rotation
                if(xmove==0 && ymove==0){
                    return;
                }
                double rotz=java.lang.Math.atan(ymove/xmove)-java.lang.Math.PI/2;//the amount I have to rotate in xy
                double rotx=java.lang.Math.sqrt(xmove*xmove+ymove*ymove);//the amount to rotate by
                if((xmove<0)){
                    rotx=-rotx;
                }
                double cosrotz=java.lang.Math.cos(rotz);
                double sinrotz=java.lang.Math.sin(rotz);
                double cosrotx=java.lang.Math.cos(rotx);
                double sinrotx=java.lang.Math.sin(rotx);
                tmprotmtx[0][0]=((cosrotz*cosrotz)+(-sinrotz*cosrotx*-sinrotz));
                tmprotmtx[0][1]=((cosrotz*sinrotz)+(-sinrotz*cosrotx*cosrotz));
                tmprotmtx[0][2]=(-sinrotz*-sinrotx);
                tmprotmtx[1][0]=((sinrotz*cosrotz)+(cosrotz*cosrotx*-sinrotz));
                tmprotmtx[1][1]=((sinrotz*sinrotz)+(cosrotz*cosrotx*cosrotz));
                tmprotmtx[1][2]=(cosrotz*-sinrotx);
                tmprotmtx[2][0]=(sinrotx*-sinrotz);
                tmprotmtx[2][1]=(sinrotx*cosrotz);
                tmprotmtx[2][2]=cosrotx;
            }//end xy rotation
            xmove=0;
            ymove=0;
            zmove=0;
            //now multiply tmprotmtx to rotmtx and save in myrotmtx
            data.myrotmtx[0][0]=(tmprotmtx[0][0]*data.rotmtx[0][0])+(tmprotmtx[0][1]*data.rotmtx[1][0])+(tmprotmtx[0][2]*data.rotmtx[2][0]);
            data.myrotmtx[0][1]=(tmprotmtx[0][0]*data.rotmtx[0][1])+(tmprotmtx[0][1]*data.rotmtx[1][1])+(tmprotmtx[0][2]*data.rotmtx[2][1]);
            data.myrotmtx[0][2]=(tmprotmtx[0][0]*data.rotmtx[0][2])+(tmprotmtx[0][1]*data.rotmtx[1][2])+(tmprotmtx[0][2]*data.rotmtx[2][2]);
            data.myrotmtx[1][0]=(tmprotmtx[1][0]*data.rotmtx[0][0])+(tmprotmtx[1][1]*data.rotmtx[1][0])+(tmprotmtx[1][2]*data.rotmtx[2][0]);
            data.myrotmtx[1][1]=(tmprotmtx[1][0]*data.rotmtx[0][1])+(tmprotmtx[1][1]*data.rotmtx[1][1])+(tmprotmtx[1][2]*data.rotmtx[2][1]);
            data.myrotmtx[1][2]=(tmprotmtx[1][0]*data.rotmtx[0][2])+(tmprotmtx[1][1]*data.rotmtx[1][2])+(tmprotmtx[1][2]*data.rotmtx[2][2]);
            data.myrotmtx[2][0]=(tmprotmtx[2][0]*data.rotmtx[0][0])+(tmprotmtx[2][1]*data.rotmtx[1][0])+(tmprotmtx[2][2]*data.rotmtx[2][0]);
            data.myrotmtx[2][1]=(tmprotmtx[2][0]*data.rotmtx[0][1])+(tmprotmtx[2][1]*data.rotmtx[1][1])+(tmprotmtx[2][2]*data.rotmtx[2][1]);
            data.myrotmtx[2][2]=(tmprotmtx[2][0]*data.rotmtx[0][2])+(tmprotmtx[2][1]*data.rotmtx[1][2])+(tmprotmtx[2][2]*data.rotmtx[2][2]);
        }// end getangles
        
        //-----------------------------------
        
        void drawbox(java.awt.Graphics g){
            //draw a box onto the screen showing the current orientation of worldspace
            //do the xz rotation
            double minxl=data.myrotmtx[0][0]*box[0][0]+data.myrotmtx[0][1]*box[0][1]+data.myrotmtx[0][2]*box[0][2];
            double maxxl=minxl;
            double minyl=data.myrotmtx[1][0]*box[0][0]+data.myrotmtx[1][1]*box[0][1]+data.myrotmtx[1][2]*box[0][2];
            double maxyl=minyl;
            boxdata[3][0]=data.myrotmtx[0][0]*box[3][0]+data.myrotmtx[0][1]*box[3][1]+data.myrotmtx[0][2]*box[3][2];
            boxdata[3][1]=data.myrotmtx[1][0]*box[3][0]+data.myrotmtx[1][1]*box[3][1]+data.myrotmtx[1][2]*box[3][2];
            boxdata[3][2]=data.myrotmtx[2][0]*box[3][0]+data.myrotmtx[2][1]*box[3][1]+data.myrotmtx[2][2]*box[3][2];
            if(boxdata[3][0]<minxl){
                minxl=boxdata[3][0];
            }else if(boxdata[3][0]>maxxl){
                maxxl=boxdata[3][0];
            }
            if(boxdata[3][1]<minyl){
                minyl=boxdata[3][1];
            }else if(boxdata[3][1]>maxyl){
                maxyl=boxdata[3][1];
            }
            boxdata[2][0]=data.myrotmtx[0][0]*box[2][0]+data.myrotmtx[0][1]*box[2][1]+data.myrotmtx[0][2]*box[2][2];
            boxdata[2][1]=data.myrotmtx[1][0]*box[2][0]+data.myrotmtx[1][1]*box[2][1]+data.myrotmtx[1][2]*box[2][2];
            boxdata[2][2]=data.myrotmtx[2][0]*box[2][0]+data.myrotmtx[2][1]*box[2][1]+data.myrotmtx[2][2]*box[2][2];
            if(boxdata[2][0]<minxl){
                minxl=boxdata[2][0];
            }else if(boxdata[2][0]>maxxl){
                maxxl=boxdata[2][0];
            }
            if(boxdata[2][1]<minyl){
                minyl=boxdata[2][1];
            }else if(boxdata[2][1]>maxyl){
                maxyl=boxdata[2][1];
            }
            boxdata[7][0]=data.myrotmtx[0][0]*box[7][0]+data.myrotmtx[0][1]*box[7][1]+data.myrotmtx[0][2]*box[7][2];
            boxdata[7][1]=data.myrotmtx[1][0]*box[7][0]+data.myrotmtx[1][1]*box[7][1]+data.myrotmtx[1][2]*box[7][2];
            boxdata[7][2]=data.myrotmtx[2][0]*box[7][0]+data.myrotmtx[2][1]*box[7][1]+data.myrotmtx[2][2]*box[7][2];
            if(boxdata[7][0]<minxl){
                minxl=boxdata[7][0];
            }else if(boxdata[7][0]>maxxl){
                maxxl=boxdata[7][0];
            }
            if(boxdata[7][1]<minyl){
                minyl=boxdata[7][1];
            }else if(boxdata[7][1]>maxyl){
                maxyl=boxdata[7][1];
            }
            double xoffset=-minxl;
            double yoffset=-minyl;
            double xfac=50;
            double yfac=50;
            g.setColor(bgcolor);
            g.fillRect(xadd,yadd,(int)(xfac*1.5),(int)(yfac*1.5));
            g.setColor(fgcolor);
            //now I have the coordinates; now draw the box
            g.drawLine((int)((boxdata[0][0]+xoffset)*xfac+xadd),(int)((boxdata[0][1]+yoffset)*yfac+yadd),(int)((boxdata[3][0]+xoffset)*xfac+xadd),(int)((boxdata[3][1]+yoffset)*yfac+yadd));
            g.drawLine((int)((boxdata[2][0]+xoffset)*xfac+xadd),(int)((boxdata[2][1]+yoffset)*yfac+yadd),(int)((boxdata[3][0]+xoffset)*xfac+xadd),(int)((boxdata[3][1]+yoffset)*yfac+yadd));
            g.drawLine((int)((boxdata[3][0]+xoffset)*xfac+xadd),(int)((boxdata[3][1]+yoffset)*yfac+yadd),(int)((boxdata[7][0]+xoffset)*xfac+xadd),(int)((boxdata[7][1]+yoffset)*yfac+yadd));
            g.setColor(java.awt.Color.red);
            g.drawString("Y",(int)((boxdata[0][0]+xoffset)*xfac+xadd),(int)((boxdata[0][1]+yoffset)*yfac+yadd));
            g.drawString("X",(int)((boxdata[2][0]+xoffset)*xfac+xadd),(int)((boxdata[2][1]+yoffset)*yfac+yadd));
            g.drawString("0,0,0",(int)((boxdata[3][0]+xoffset)*xfac+xadd),(int)((boxdata[3][1]+yoffset)*yfac+yadd));
            g.drawString("Z",(int)((boxdata[7][0]+xoffset)*xfac+xadd),(int)((boxdata[7][1]+yoffset)*yfac+yadd));
        }// end drawbox
        
        //---------------------------------------
        
        void makedrawdata(){
            //use the data from posarr and the rotation angles xz and yz to
            double xfac=1;
            double yfac=1;
            double xoffset;
            double yoffset;
            float[][] tposarrtmp=data.posarrtmp;
            float[][] posarr=data.posarr;
            double[][] myrotmtx=data.myrotmtx;
            //compute the new positions and save them to posarrtmp
            for(int i=data.elements;--i>=0;){
                //same as for boxdata
                tposarrtmp[i][0]=(float)(myrotmtx[0][0]*posarr[i][0]+myrotmtx[0][1]*posarr[i][1]+myrotmtx[0][2]*posarr[i][2]);
                tposarrtmp[i][1]=(float)(myrotmtx[1][0]*posarr[i][0]+myrotmtx[1][1]*posarr[i][1]+myrotmtx[1][2]*posarr[i][2]);
                data.drawarrtmp[i][0]=(int)tposarrtmp[i][0];
                data.drawarrtmp[i][1]=(int)tposarrtmp[i][1];
                //I don't need to calculate the 3rd dimention
                //tposarrtmp[i][2]=(float)(myrotmtx[2][0]*posarr[i][0]+myrotmtx[2][1]*posarr[i][1]+myrotmtx[2][2]*posarr[i][2]);
                //now add perspective using posarrtmp[i][2]
                //maxdiff is the maximum value for any coordinate in the system
                //posarrtmp[i][0]=posarrtmp[i][0]*maxdiff/(maxdiff+posarrtmp[i][2]);
                //posarrtmp[i][1]=posarrtmp[i][1]*maxdiff/(maxdiff+posarrtmp[i][2]);
            }// end for i
            int selectedelements=java.lang.reflect.Array.getLength(data.selectednames);
            if((zoom)&&(selectedelements>0)){//if zoom=true
                //only get maxx and maxy and minx and miny from selected sequences
                //I know I have at least one selected element (checked at zoombuttonactionperformed)
                maxx=tposarrtmp[data.selectednames[0]][0];
                maxy=tposarrtmp[data.selectednames[0]][1];
                minx=maxx;
                miny=maxy;
                for(int i=selectedelements;--i>=0;){//for (int i=1;i<selectedelements;i++){
                    if(maxx<tposarrtmp[data.selectednames[i]][0]){
                        maxx=tposarrtmp[data.selectednames[i]][0];
                    }else if(minx>tposarrtmp[data.selectednames[i]][0]){
                        minx=tposarrtmp[data.selectednames[i]][0];
                    }
                    if(maxy<tposarrtmp[data.selectednames[i]][1]){
                        maxy=tposarrtmp[data.selectednames[i]][1];
                    }else if(miny>tposarrtmp[data.selectednames[i]][1]){
                        miny=tposarrtmp[data.selectednames[i]][1];
                    }
                }// end for i
                xfac=drawwidth/(maxx-minx);
                yfac=drawheight/(maxy-miny);
                xoffset=(-minx);
                yoffset=(-miny);
                if(maxx-minx==0){
                    System.out.println("isZero (in Zoom)!");
                }
                for(int i=0;i<data.elements;i++){
                    tposarrtmp[i][0]=(float)(((tposarrtmp[i][0]+xoffset)*xfac)+xadd);
                    tposarrtmp[i][1]=(float)(((tposarrtmp[i][1]+yoffset)*yfac)+yadd);
                    data.drawarrtmp[i][0]=(int)tposarrtmp[i][0];
                    data.drawarrtmp[i][1]=(int)tposarrtmp[i][1];
                    //I don't need to calculate the 3rd dimention
                }// end for i
                originpos[0]=(int)((xoffset*xfac)+xadd);
                originpos[1]=(int)((yoffset*yfac)+yadd);
            }else{//if(zoom==false){//if I don't want to zoom in on selected sequences
                maxx=tposarrtmp[0][0];
                maxy=tposarrtmp[0][1];
                minx=maxx;
                miny=maxy;
                for(int i=data.elements;--i>=0;){//for (int i=1;i<edelements;i++){
                    if(maxx<tposarrtmp[i][0]){
                        maxx=tposarrtmp[i][0];
                    }else if(minx>tposarrtmp[i][0]){
                        minx=tposarrtmp[i][0];
                    }
                    if(maxy<tposarrtmp[i][1]){
                        maxy=tposarrtmp[i][1];
                    }else if(miny>tposarrtmp[i][1]){
                        miny=tposarrtmp[i][1];
                    }
                }// end for i
                if(maxx-minx==0){
                    //System.out.println("isZero!");
                }
                xfac=drawwidth/(maxx-minx);
                yfac=drawheight/(maxy-miny);
                xoffset=(-minx);
                yoffset=(-miny);
                for(int i=data.elements;--i>=0;){
                    tposarrtmp[i][0]=(float)(((tposarrtmp[i][0]+xoffset)*xfac)+xadd);
                    tposarrtmp[i][1]=(float)(((tposarrtmp[i][1]+yoffset)*yfac)+yadd);
                    data.drawarrtmp[i][0]=(int)tposarrtmp[i][0];
                    data.drawarrtmp[i][1]=(int)tposarrtmp[i][1];
                    //I don't need to calculate the 3rd dimention
                }// end for i
                originpos[0]=(int)((xoffset*xfac)+xadd);
                originpos[1]=(int)((yoffset*yfac)+yadd);
            }//end if zoom==false
            xscale=xfac;
            yscale=yfac;
        }// end makedrawdata
        
        //---------------------------------------
        
        void drawdata(java.awt.Graphics g1d){
            java.awt.Graphics2D g=(java.awt.Graphics2D) g1d;
            if(antialiasingcheckboxmenuitem.isSelected()){
                g.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING,java.awt.RenderingHints.VALUE_ANTIALIAS_ON);
            }else{
                g.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING,java.awt.RenderingHints.VALUE_ANTIALIAS_OFF);
            }
            float halfdot=(float)(data.dotsize/2);
            float halfoval=(float)(data.ovalsize/2);
            float halfgroup=(float)(data.groupsize/2);
            //old, slow for many elements
            int[][] tposarrtmp=data.drawarrtmp;
            int elements=data.elements;//java.lang.reflect.Array.getLength(tposarrtmp);//was posarr 19.1.04
            int tmpint;
            int dotsize=data.dotsize;
            int ovalsize=data.ovalsize;

            if(elements>0){
                makedrawdata();
                g.setColor(fgcolor);
                g.drawRect(xadd-xtranslate,yadd-ytranslate,graphpanel.getWidth()-(2*xadd),graphpanel.getHeight()-(2*yadd));
                if(checkbox_show_connections.isSelected()){
                    //here i am looking for a speedup by eliminating the multiple "for" loops
                    //to this effect the drawing order is computed in an earlier step. here I just have to loop through
                    //the 2d array and draw the elements in the according color
                    //int colornum=java.lang.reflect.Array.getLength(colorarr);
                    if(colorfrustrationcheckbox.isSelected()==false){
                        //color the lines by their blast P-value
                        int j;
                        int[] currdraw;
                        int vecsize;
                        if(data.draworder==null || java.lang.reflect.Array.getLength(data.draworder)<1){
                            data.draworder=getdraworder(data.myattvals,colornum);
                        }
                        for(int i=0;i<colornum;i++){
                            vecsize=data.draworder[i].size();
                            g.setColor(data.colorarr[i]);
                            g.drawString("-",(5+i)*fontwidth-xtranslate,fontsize-ytranslate);
                            //draw all the vector elements
                            for(j=vecsize;--j>=0;){
                                currdraw=(int[])data.draworder[i].get(j);
                                //g.drawLine((int)tposarrtmp[currdraw[0]][0],(int)tposarrtmp[currdraw[0]][1],(int)tposarrtmp[currdraw[1]][0],(int)tposarrtmp[currdraw[1]][1]);
                                g.drawLine(tposarrtmp[currdraw[0]][0],tposarrtmp[currdraw[0]][1],tposarrtmp[currdraw[1]][0],tposarrtmp[currdraw[1]][1]);
                            }// end for j
                        }// end for i
                        g.setColor(fgcolor);
                        g.drawString("Worst",-xtranslate,fontsize-ytranslate);
                        g.drawString("Best",((colornum+5)*fontwidth)-xtranslate,fontsize-ytranslate);
                    }else{//if colorfrustrationcheckbox.isSelected()
                        //color the lines by their "frustration" (i.e. too long or too short for P-value)
                        frustration=getfrustration(data.myattvals,data.posarr);//get the frustration for each line
                        int j;
                        //int[] currdraw;
                        int vecsize;
                        String key;
                        int[] tmparr;
                        minattvals currfrust;
                        if(java.lang.reflect.Array.getLength(data.draworder)<1){
                            data.draworder=getdraworder(data.myattvals,colornum);//draw the lines in the same order as before
                        }
                        for(int i=0;i<colornum;i++){
                            vecsize=data.draworder[i].size();
                            //draw all the vector elements
                            for(j=0;j<vecsize;j++){
                                tmparr=(int[])data.draworder[i].get(j);
                                //if(tmparr[0]<tmparr[1]){
                                key=tmparr[0]+"_"+tmparr[1];
                                //}else{
                                //    key=tmparr[1]+"_"+tmparr[0];
                                //}
                                if(frustration.containsKey(key)){
                                    currfrust=(minattvals)frustration.get(key);
                                    if(currfrust.att>0){
                                        //then the coloring should be in blue(i.e. too short for attval)
                                        g.setColor(new java.awt.Color((1-currfrust.att),(1-currfrust.att),1));
                                    }else{
                                        //then the coloring should be in red(i.e. too long for attval)
                                        g.setColor(new java.awt.Color(1,(1+currfrust.att),(1+currfrust.att)));
                                    }
                                    //g.drawLine((int)tposarrtmp[currfrust.query][0],(int)tposarrtmp[currfrust.query][1],(int)tposarrtmp[currfrust.hit][0],(int)tposarrtmp[currfrust.hit][1]);
                                    g.drawLine(tposarrtmp[currfrust.query][0],tposarrtmp[currfrust.query][1],tposarrtmp[currfrust.hit][0],tposarrtmp[currfrust.hit][1]);
                                }else{
                                    System.err.println("no value for frustration key '"+key+"'");
                                }
                            }// end for j
                        }// end for i
                        g.setColor(fgcolor);
                        g.drawString("Worst",-xtranslate,fontsize-ytranslate);
                        g.drawString("Best",((colornum+5)*fontwidth)-xtranslate,fontsize-ytranslate);
                    }
                }
                if(checkbox_show_names.isSelected()){
                    int selectednamesnum=java.lang.reflect.Array.getLength(data.selectednames);
                    String[] namearr=data.namearr;
                    if(selectednamesnum==0){
                        for(int i=0;i<elements;i++){
                            g.drawString(String.valueOf(i)+"-"+namearr[i],(int)tposarrtmp[i][0],(int)tposarrtmp[i][1]);
                        }// end for i
                    }else{
                        for(int i=0;i<selectednamesnum;i++){
                            g.drawString(String.valueOf(i)+"-"+namearr[data.selectednames[i]],(int)tposarrtmp[data.selectednames[i]][0],(int)tposarrtmp[data.selectednames[i]][1]);
                        }// end for i
                    }
                }else if(checkbox_show_numbers.isSelected()){
                    for(int i=0;i<elements;i++){
                        g.drawString(String.valueOf(i),(int)tposarrtmp[i][0],(int)tposarrtmp[i][1]);
                    }// end for i
                }
                if(dotsfirst==true){
                    //now draw the sequence dots
                    if(lengthcolormenuitem.isSelected()){
                        float[] seqlengths=data.seqlengths;
                        for(int i=0;i<elements;i++){
                            g.setColor(new java.awt.Color(1-seqlengths[i],1-seqlengths[i],seqlengths[i]));
                            g.fillOval((int)(tposarrtmp[i][0]-halfdot),(int)(tposarrtmp[i][1]-halfdot),dotsize,dotsize);
                        }
                    }else{
                        g.setColor(fgcolor);
                        for(int i=0;i<elements;i++){
                            g.fillOval((int)(tposarrtmp[i][0]-halfdot),(int)(tposarrtmp[i][1]-halfdot),dotsize,dotsize);
                        }
                    }
                }
                //draw the sequences from older selection steps
                //now draw the elements in the selectvec
                for(int i=selectvec.size()-1;i>=0;i--){
                    tmpselectclass=(selectclass)selectvec.elementAt(i);
                    g.setColor(tmpselectclass.color);
                    int tmpelements=java.lang.reflect.Array.getLength(tmpselectclass.selectednames);
                    for(int j=0;j<tmpelements;j++){
                        g.fillOval((int)(tposarrtmp[tmpselectclass.selectednames[j]][0]-halfoval),(int)(tposarrtmp[tmpselectclass.selectednames[j]][1]-halfoval),ovalsize,ovalsize);
                    }
                }//end for i
                //draw the sequence groups
                if(data.showseqgroups){
                    seqgroup mygroup;
                    for(int i=data.seqgroupsvec.size()-1;i>=0;i--){
                        mygroup=(seqgroup)data.seqgroupsvec.elementAt(i);
                        if(mygroup.hide==true){
                            continue;
                        }
                        g.setColor(mygroup.color);
                        if(mygroup.type==0){
                            tmpint=((seqgroup)data.seqgroupsvec.elementAt(i)).size;
                            halfgroup=(float)tmpint/2;
                            for(int j=java.lang.reflect.Array.getLength(mygroup.sequences)-1;j>=0;j--){
                                if(mygroup.sequences[j]<elements){
                                    g.fillOval((int)(tposarrtmp[mygroup.sequences[j]][0]-halfgroup),(int)(tposarrtmp[mygroup.sequences[j]][1]-halfgroup),tmpint,tmpint);
                                }else{
                                    System.err.println("sequence number "+mygroup.sequences[j]+" is not present in file; removing entry from group");
                                    mygroup.remove(j);
                                }
                            }//end for j
                        }else{
                            tmpposarr=mygroup.polygon;//((int[][])polygons.elementAt(mygroup.polygon));
                            posnum=tmpposarr[2][0];//the number of points
                            xposarr=new int[posnum];
                            yposarr=new int[posnum];
                            for(int j=java.lang.reflect.Array.getLength(mygroup.sequences)-1;j>=0;j--){
                                if(mygroup.sequences[j]<elements){
                                    for(int k=0;k<posnum;k++){
                                        xposarr[k]=(int)(tmpposarr[0][k]+tposarrtmp[mygroup.sequences[j]][0]);
                                        yposarr[k]=(int)(tmpposarr[1][k]+tposarrtmp[mygroup.sequences[j]][1]);
                                    }//end for k
                                    g.fillPolygon(xposarr,yposarr,posnum);
                                }else{
                                    System.err.println("sequence number "+mygroup.sequences[j]+" is not present in file; removing entry from group");
                                    mygroup.remove(j);
                                }
                            }//end for j
                        }
                    }//end for i
                }//end if showseqgroups
                if(groupseqs!=null){//draw the sequences from the currently selected group
                    g.setColor(groupseqscolor);
                    for(int i=java.lang.reflect.Array.getLength(groupseqs);--i>=0;){
                        g.fillOval((int)(tposarrtmp[groupseqs[i]][0]-halfoval),(int)(tposarrtmp[groupseqs[i]][1]-halfoval),ovalsize,ovalsize);
                    }//end for i
                }
                //draw the current selection
                if(clusterconf==null){
                    g.setColor(selectedcolor);
                    for(int i=java.lang.reflect.Array.getLength(data.selectednames)-1;i>=0;i--){
                        g.fillOval((int)(tposarrtmp[data.selectednames[i]][0]-halfoval),(int)(tposarrtmp[data.selectednames[i]][1]-halfoval),ovalsize,ovalsize);
                    }//end for selected
                }else{
                    for(int i=java.lang.reflect.Array.getLength(data.selectednames)-1;i>=0;i--){
                        if (clusterconf[i]<0){
                            clusterconf[i]=0;
                        }
                        g.setColor(new java.awt.Color((int)(selectedcolor.getRed()*clusterconf[i]),(int)(selectedcolor.getGreen()*clusterconf[i]),(int)(selectedcolor.getBlue()*clusterconf[i])));
                        g.fillOval((int)(tposarrtmp[data.selectednames[i]][0]-halfoval),(int)(tposarrtmp[data.selectednames[i]][1]-halfoval),ovalsize,ovalsize);
                    }//end for selected
                }
                g.setColor(fgcolor);
                for(int i=java.lang.reflect.Array.getLength(data.selectednames)-1;i>=0;i--){
                    g.drawOval((int)(tposarrtmp[data.selectednames[i]][0]-halfoval),(int)(tposarrtmp[data.selectednames[i]][1]-halfoval),ovalsize,ovalsize);
                }//end for selected
                if(java.lang.reflect.Array.getLength(blastselectseqs)>0){
                    //draw these sequences as green dots
                    g.setColor(blasthitcolor);
                    if(showblasthitnamescheckbox.isSelected()){
                        //write the names to screen
                        for(int i=java.lang.reflect.Array.getLength(blastselectseqs)-1;i>=0;i--){
                            g.drawString(String.valueOf(blastselectseqs[i]),(int)tposarrtmp[blastselectseqs[i]][0],(int)tposarrtmp[blastselectseqs[i]][1]);
                        }//end for i
                    }//end if show names for blasthits
                    g.setColor(blastcirclecolor);
                    for(int i=java.lang.reflect.Array.getLength(blastselectseqs)-1;i>=0;i--){
                        g.fillOval((int)(tposarrtmp[blastselectseqs[i]][0]-halfoval),(int)(tposarrtmp[blastselectseqs[i]][1]-halfoval),ovalsize,ovalsize);
                    }//end for selected
                    //now draw the blast query
                    g.setColor(java.awt.Color.magenta);
                    g.fillOval((int)(tposarrtmp[blastselectseqs[0]][0]-halfoval),(int)(tposarrtmp[blastselectseqs[0]][1]-halfoval),ovalsize,ovalsize);
                    //now draw the black outline
                    g.setColor(fgcolor);
                    for(int i=java.lang.reflect.Array.getLength(blastselectseqs)-1;i>=0;i--){
                        g.drawOval((int)(tposarrtmp[blastselectseqs[i]][0]-halfoval),(int)(tposarrtmp[blastselectseqs[i]][1]-halfoval),ovalsize,ovalsize);
                    }//end for selected
                }
                if(dotsfirst==false){
                    //now draw the sequence dots
                    if(lengthcolormenuitem.isSelected()){
                        float[] seqlengths=data.seqlengths;
                        for(int i=0;i<elements;i++){
                            g.setColor(new java.awt.Color(1-seqlengths[i],1-seqlengths[i],seqlengths[i]));
                            g.fillOval((int)(tposarrtmp[i][0]-halfdot),(int)(tposarrtmp[i][1]-halfdot),dotsize,dotsize);
                        }
                    }else{
                        g.setColor(fgcolor);
                        for(int i=0;i<elements;i++){
                            g.fillOval((int)(tposarrtmp[i][0]-halfdot),(int)(tposarrtmp[i][1]-halfdot),dotsize,dotsize);
                        }
                    }
                }
            }// end if elements>0
            if(button_select_move.isSelected()&&mouse_is_pressed){
                g.setColor(java.awt.Color.orange);
                if(currmousepos[0]<selectstart[0]){
                    draw[0]=currmousepos[0];
                    draw[2]=selectstart[0];
                }else{
                    draw[0]=selectstart[0];
                    draw[2]=currmousepos[0];
                }
                if(currmousepos[1]<selectstart[1]){
                    draw[1]=currmousepos[1];
                    draw[3]=selectstart[1];
                }else{
                    draw[1]=selectstart[1];
                    draw[3]=currmousepos[1];
                }
                g.drawRect(draw[0],draw[1],(draw[2]-draw[0]),(draw[3]-draw[1]));
            }
        }// end drawdata
        
        //----------------------------------------------------------------------
        
        ArrayList[] getdraworder(minattvals[] myattvals,int colornum){
            //bin the lines connecting datapoints by the color they are assigned (makes subsequent drawing quicker)
            int elements=java.lang.reflect.Array.getLength(myattvals);
            //System.out.println("attvalelements="+elements);
            ArrayList[] retarr=new ArrayList[colornum];
            for(int i=0;i<colornum;i++){
                retarr[i]=new ArrayList<int[]>();
            }
            int[] mydraw=new int[2];//connection between sequences i & j
            //for(int i=0;i<elements;i++){
            float[] colorcutoffs=data.colorcutoffs;
            for(int i=elements;--i>=0;){
                //System.out.println("attvals "+myattvals[i].query+";"+myattvals[i].hit+"="+myattvals[i].att);
                if((myattvals[i].att==0)||(myattvals[i].att<=colorcutoffs[0])){//if I am below the lowest to draw
                    continue;
                }else if(myattvals[i].att<=colorcutoffs[1]){//from 0 to 1
                    mydraw=new int[2];
                    mydraw[0]=myattvals[i].query;
                    mydraw[1]=myattvals[i].hit;
                    //retarr[0].addElement(mydraw);
                    retarr[0].add(mydraw);
                }else if(myattvals[i].att<=colorcutoffs[2]){//from 1 to 2
                    mydraw=new int[2];
                    mydraw[0]=myattvals[i].query;
                    mydraw[1]=myattvals[i].hit;
                    retarr[1].add(mydraw);
                }else if(myattvals[i].att<=colorcutoffs[3]){//from 2 to 3
                    mydraw=new int[2];
                    mydraw[0]=myattvals[i].query;
                    mydraw[1]=myattvals[i].hit;
                    retarr[2].add(mydraw);
                }else if(myattvals[i].att<=colorcutoffs[4]){//from 3 to 4
                    mydraw=new int[2];
                    mydraw[0]=myattvals[i].query;
                    mydraw[1]=myattvals[i].hit;
                    retarr[3].add(mydraw);
                }else if(myattvals[i].att<=colorcutoffs[5]){//from 4 to 5
                    mydraw=new int[2];
                    mydraw[0]=myattvals[i].query;
                    mydraw[1]=myattvals[i].hit;
                    retarr[4].add(mydraw);
                }else if(myattvals[i].att<=colorcutoffs[6]){//from 5 to 6
                    mydraw=new int[2];
                    mydraw[0]=myattvals[i].query;
                    mydraw[1]=myattvals[i].hit;
                    retarr[5].add(mydraw);
                }else if(myattvals[i].att<=colorcutoffs[7]){//from 6 to 7
                    mydraw=new int[2];
                    mydraw[0]=myattvals[i].query;
                    mydraw[1]=myattvals[i].hit;
                    retarr[6].add(mydraw);
                }else if(myattvals[i].att<=colorcutoffs[8]){//from 7 to 8
                    mydraw=new int[2];
                    mydraw[0]=myattvals[i].query;
                    mydraw[1]=myattvals[i].hit;
                    retarr[7].add(mydraw);
                }else if(myattvals[i].att<=colorcutoffs[9]){//from 8 to 9
                    mydraw=new int[2];
                    mydraw[0]=myattvals[i].query;
                    mydraw[1]=myattvals[i].hit;
                    retarr[8].add(mydraw);
                }else if(myattvals[i].att>colorcutoffs[9]){//better than 9
                    mydraw=new int[2];
                    mydraw[0]=myattvals[i].query;
                    mydraw[1]=myattvals[i].hit;
                    retarr[9].add(mydraw);
                }
            }// end for i
            return retarr;
        }// end getdraworder
        
        //----------------------------------------------------------------------
        
        HashMap getfrustration(minattvals[] attvals,float[][] posarr){
            //get wether each connection is longer or shorter than expected based on the attraction value
            int seqnum=java.lang.reflect.Array.getLength(attvals);
            HashMap retarr=new HashMap((int)(seqnum/0.8)+1,0.8f);
            minattvals[] minattarr=new minattvals[seqnum];
            int q,h;
            float tmpfloat;
            float att;
            float tmpx,tmpy,tmpz;
            float avgatt=0;
            float avgdist=0;
            float maxval=0;
            int counter=0;
            float tmpval;
            String key;
            for(int i=0;i<seqnum;i++){
                att=attvals[i].att;
                q=attvals[i].query;
                h=attvals[i].hit;
                key=q+"_"+h;
                if(att>0){//else don't do the calculations as I won't be drawing a line
                    counter++;
                    tmpx=posarr[q][0]-posarr[h][0];
                    tmpy=posarr[q][1]-posarr[h][1];
                    tmpz=posarr[q][2]-posarr[h][2];
                    tmpval=(float)java.lang.Math.sqrt(tmpx*tmpx+tmpy*tmpy+tmpz*tmpz);
                    minattarr[i]=new minattvals(q,h,tmpval);
                    retarr.put(key,minattarr[i]);
                    avgdist+=tmpval;
                    avgatt+=att;
                }
            }//end for i
            avgdist/=counter;
            avgatt/=counter;
            tmpfloat=avgatt*avgdist;
            //note: the higher the attval, the smaller the dist!
            for(int i=0;i<seqnum;i++){
                att=attvals[i].att;
                q=attvals[i].query;
                h=attvals[i].hit;
                if(att>0){
                    minattarr[i].att=(tmpfloat-(att*minattarr[i].att));
                    if(maxval<minattarr[i].att){
                        maxval=minattarr[i].att;
                    }
                    if(maxval<-minattarr[i].att){
                        maxval=-minattarr[i].att;
                    }
                }
            }//end for i
            //System.out.println("avgatt="+avgatt+" avgdist="+avgdist);
            for(int i=0;i<seqnum;i++){
                if(attvals[i].att>0){
                    minattarr[i].att/=maxval;
                }//else I have a null entry
            }//end fo ri
            return retarr;
        }//end getfrustration
        
        //-----------------------------------
        
    }// end class drawpanel
    
    class computethread extends java.lang.Thread{
        
        public computethread(clustermain_graphics parent){
            this.parent=parent;
            this.didrun=false;
            this.stop=false;
        }
        
        boolean stop=true;
        boolean didrun=false;
        String tmpstr="";
        float tmpcool=1;
        clustermain_graphics parent;
        
@Override
        public void run(){
            this.didrun=true;
            data.roundsdone=0;
            while (stop==false){
                data.rounds++;
                if(data.roundslimit!=-1){
                    data.roundsdone++;
                    button_start_stop_resume.setText("STOP ("+data.roundsdone+"/"+data.roundslimit+")");
                    if(data.roundsdone>=data.roundslimit){
                        stop=true;
                        synchronized(parent){
                            parent.notify();
                        }
                    }
                }
                //first see whether the main window has the focus
                if(graphpanel.isFocusOwner()){
                    //mainwindow has the focus, then see whether I am done drawing
                    try{
                        //while(docalc.equals("false")){
                        while(recalc==false){
                            this.sleep(100);
                        }
                    }catch (InterruptedException ie){
                        System.err.println("Interrupted sleep in computethread");
                    }
                }
                //myposarr=recluster3d();
                if(data.changedvals){
                    updatevals();
                    if(data.hidebelowold!=data.hidebelow){
                        //draw1.draworder=new Vector[0];
                        data.draworder=new ArrayList[0];
                        data.hidebelowold=data.hidebelow;
                    }
                    data.changedvals=false;
                }//end if changedvals
                ClusterMethods.recluster3d(parent.data);
                if((level==0)&&(savepos)){
                    try{
                        PrintWriter outwriter=new PrintWriter(new BufferedWriter(new FileWriter("positionfile.dat")));
                        //now save the current variables to file
                        outwriter.println("sequences: "+java.lang.reflect.Array.getLength(data.myposarr));
                        outwriter.println("values: ");
                        outwriter.println("minattract "+data.minattract);
                        outwriter.println("maxmove "+data.maxmove);
                        outwriter.println("cooling "+data.cooling);
                        outwriter.println("currcool "+data.currcool);
                        outwriter.println("mineval "+data.mineval);
                        outwriter.println("minpval "+data.minpval);
                        outwriter.println("repfactor "+data.repfactor);
                        outwriter.println("attfactor "+data.attfactor);
                        outwriter.println("hidebelow "+data.hidebelow);
                        outwriter.println("dampening "+data.dampening);
                        outwriter.println("rounds "+data.rounds);
                        printout.saveseqpos(outwriter,data.myposarr);
                        outwriter.flush();
                        outwriter.close();
                    }catch (IOException e){
                        System.err.println("unable to save positions to positionfile.dat");
                        e.printStackTrace();
                    }
                }//end if savepos
                data.posarr=data.myposarr;
                tmpcool=(((float)((int)(data.currcool*100000)))/100000);
                if(myoptionswindow!=null){
                    myoptionswindow.currcoolfield.setText(String.valueOf(tmpcool));
                }
                if(tmpcool<=1e-5){
                    stop=true;
                    button_start_stop_resume.setText("DONE (absolute zero)");
                }
                if(data.rounds%skiprounds==0 && data.nographics==false){
                    //docalc="false";//don't do any further calculations until the drawing is done!
                    recalc=false;
                    repaint();
                }
            }// end while
            button_start_stop_resume.setText("Resume");
            button_start_stop_resume.setEnabled(true);
            parent.repaint();
        }// end run
        
    }// end class computethread

    /*
    
    class getmovethread extends java.lang.Thread{
        
        public getmovethread(float[][] myposarr,minattvals[] myattvals,float[][] mymovearr,int myi,int cpu,HashMap selectnamehash,int[] selectnames,String syncon, clustertest parent){
            this.done=false;
            this.posarr=myposarr;
            this.attvals=myattvals;
            this.movearr=mymovearr;
            this.myi=myi;
            this.cpu=cpu;
            this.doselected=true;
            this.tmphash=selectnamehash;
            this.parent=parent;
            this.syncon=syncon;
            this.selectnames=selectnames;
        }
        
        public getmovethread(float[][] myposarr,minattvals[] myattvals,float[][] mymovearr,int myi,int cpu,String syncon,clustertest parent){
            this.done=false;
            this.posarr=myposarr;
            this.attvals=myattvals;
            this.movearr=mymovearr;
            this.myi=myi;
            this.cpu=cpu;
            this.doselected=false;
            this.tmphash=null;
            this.parent=parent;
            this.syncon=syncon;
        }
        
        boolean done;
        boolean doselected;
        float[][] posarr;
        minattvals[] attvals;
        float[][] movearr;
        int myi;
        int cpu;
        HashMap tmphash=null;
        clustertest parent;
        String syncon;
        int[] selectnames;//a local copy of parent.selectednames, as that may change during a calculation
        
@Override
        public void run(){
            //System.out.println("start="+start);
            //use the positions of all elements and their attraction/repulsion values to
            //calculate a movement vector for each (take into account the last movement).
            //repulsion doesn't have a specific value as all evalues below a certain point
            //are simply regarded as insignificant. therefore use a one formula for all
            //to compute the repulsive forces (see getrepulse)
            int j;
            int k;
            int elements=java.lang.reflect.Array.getLength(posarr);
            double[] currmoverep=new double[3];
            double[] currmoveatt=new double[3];
            double totaldist=0;
            //double totalmove=0;
            //double avgmove=0;
            //double minattelems=minattract*elements;
            int attnum=java.lang.reflect.Array.getLength(attvals);
            if(doselected){
                //System.out.println("getmoveselected");
                int selectnamesnum=java.lang.reflect.Array.getLength(selectnames);
                //now get from where to where I should do my calculations
                int start,end;//,attstart,attend;
                if(myi==(cpu-1)){
                    start=myi*(int)(selectnamesnum/cpu);
                    end=selectnamesnum;
                    //attstart=myi*(int)(attnum/cpu);
                    //attend=attnum;
                }else{
                    start=myi*(int)(selectnamesnum/cpu);
                    end=(myi+1)*(int)(selectnamesnum/cpu);
                    //attstart=myi*(int)(attnum/cpu);
                    //attend=(myi+1)*(int)(attnum/cpu);
                }
                if(parent.cluster2dbutton.isSelected()){
                    //System.out.println("clusterselected2D");
                    //cluster only the selected sequences in 2D
                    for(int i=start;i<end;i++){
                        //currmoveatt=getminattract(posarr[selectnames[i]],currmoveatt,minattract);
                        getminattract(posarr[selectnames[i]],currmoveatt,minattract);
                        movearr[selectnames[i]][0]+=currmoveatt[0];
                        movearr[selectnames[i]][1]+=currmoveatt[1];
                        for(j=elements;--j>=0;){
                            //currmoverep=getrepulse2d(posarr[selectnames[i]],posarr[j],currmoverep,repvalpow,repfactor,rand);
                            getrepulse2d(posarr[selectnames[i]],posarr[j],currmoverep,repvalpow,repfactor,rand);
                            movearr[selectnames[i]][0]+=currmoverep[0];
                            movearr[selectnames[i]][1]+=currmoverep[1];
                        }//end for j
                    }//end for i
                    //now add the attraction values, but only for the query or hit sequences in my part of the selectnames array (assigned in recluster3d)
                    for(int i=attnum;--i>=0;){
                        if((tmphash.containsKey(String.valueOf(attvals[i].query)))&&(((Integer)tmphash.get(String.valueOf(attvals[i].query))).intValue()==myi)){
                            //currmoveatt=getattract2d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow, attfactor);
                            getattract2d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow, attfactor);
                            movearr[attvals[i].query][0]+=currmoveatt[0];
                            movearr[attvals[i].query][1]+=currmoveatt[1];
                            //movement[attvals[i].query][2]+=currmoveatt[2];
                            if((tmphash.containsKey(String.valueOf(attvals[i].hit)))&&(((Integer)tmphash.get(String.valueOf(attvals[i].hit))).intValue()==myi)){
                                movearr[attvals[i].hit][0]-=currmoveatt[0];
                                movearr[attvals[i].hit][1]-=currmoveatt[1];
                                //movement[attvals[i].hit][2]-=currmoveatt[2];
                            }
                        }else if((tmphash.containsKey(String.valueOf(attvals[i].hit)))&&(((Integer)tmphash.get(String.valueOf(attvals[i].hit))).intValue()==myi)){
                            //currmoveatt=getattract2d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow, attfactor);
                            getattract2d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow, attfactor);
                            movearr[attvals[i].hit][0]-=currmoveatt[0];
                            movearr[attvals[i].hit][1]-=currmoveatt[1];
                            //movement[attvals[i].hit][2]-=currmoveatt[2];
                        }
                    }//end for i
                    //double totaldistsq;
                    for(int i=start;i<end;i++){
                        movearr[selectnames[i]][0]/=elements;
                        movearr[selectnames[i]][1]/=elements;
                        movearr[selectnames[i]][2]=0;
                        totaldist=java.lang.Math.sqrt((movearr[selectnames[i]][0]*movearr[selectnames[i]][0])+(movearr[selectnames[i]][1]*movearr[selectnames[i]][1]));
                        if(totaldist>maxmove){
                            movearr[selectnames[i]][0]*=maxmove/totaldist;
                            movearr[selectnames[i]][1]*=maxmove/totaldist;
                        }
                    }//end for i
                }else{
                    //cluster only the selected sequences in 3D
                    //System.out.println("clusterselected3D");
                    for(int i=start;i<end;i++){
                        //currmoveatt=getminattract(posarr[selectnames[i]],currmoveatt,minattract);
                        getminattract(posarr[selectnames[i]],currmoveatt,minattract);
                        movearr[selectnames[i]][0]+=currmoveatt[0];
                        movearr[selectnames[i]][1]+=currmoveatt[1];
                        movearr[selectnames[i]][2]+=currmoveatt[2];
                        for(j=elements;--j>=0;){
                            //currmoverep=getrepulse3d(posarr[selectnames[i]],posarr[j],currmoverep,repvalpow,repfactor,rand);
                            getrepulse3d(posarr[selectnames[i]],posarr[j],currmoverep,repvalpow,repfactor,rand);
                            movearr[selectnames[i]][0]+=currmoverep[0];
                            movearr[selectnames[i]][1]+=currmoverep[1];
                            movearr[selectnames[i]][2]+=currmoverep[2];
                        }//end for j
                    }//end for i
                    //now add the attraction values, but only for the query or hit sequences in my part of the selectnames array (assigned in recluster3d)
                    for(int i=attnum;--i>=0;){
                        if((tmphash.containsKey(String.valueOf(attvals[i].query)))&&(((Integer)tmphash.get(String.valueOf(attvals[i].query))).intValue()==myi)){
                            //currmoveatt=getattract3d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow, attfactor);
                            getattract3d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow, attfactor);
                            movearr[attvals[i].query][0]+=currmoveatt[0];
                            movearr[attvals[i].query][1]+=currmoveatt[1];
                            movearr[attvals[i].query][2]+=currmoveatt[2];
                            if((tmphash.containsKey(String.valueOf(attvals[i].hit)))&&(((Integer)tmphash.get(String.valueOf(attvals[i].hit))).intValue()==myi)){
                                movearr[attvals[i].hit][0]-=currmoveatt[0];
                                movearr[attvals[i].hit][1]-=currmoveatt[1];
                                movearr[attvals[i].hit][2]-=currmoveatt[2];
                            }
                        }else if((tmphash.containsKey(String.valueOf(attvals[i].hit)))&&(((Integer)tmphash.get(String.valueOf(attvals[i].hit))).intValue()==myi)){
                            //currmoveatt=getattract3d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow, attfactor);
                            getattract3d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow, attfactor);
                            movearr[attvals[i].hit][0]-=currmoveatt[0];
                            movearr[attvals[i].hit][1]-=currmoveatt[1];
                            movearr[attvals[i].hit][2]-=currmoveatt[2];
                        }
                    }//end for i
                    //double totaldistsq;
                    for(int i=start;i<end;i++){
                        movearr[selectnames[i]][0]/=elements;
                        movearr[selectnames[i]][1]/=elements;
                        movearr[selectnames[i]][2]/=elements;
                        totaldist=java.lang.Math.sqrt((movearr[selectnames[i]][0]*movearr[selectnames[i]][0])+(movearr[selectnames[i]][1]*movearr[selectnames[i]][1])+(movearr[selectnames[i]][2]*movearr[selectnames[i]][2]));
                        if(totaldist>maxmove){
                            movearr[selectnames[i]][0]*=maxmove/totaldist;
                            movearr[selectnames[i]][1]*=maxmove/totaldist;
                            movearr[selectnames[i]][2]*=maxmove/totaldist;
                        }
                    }//end for i
                }
            }else{//if no sequences were selected or all should be used
                //now get from where to where I should do my calculations
                int start,end;//,attstart,attend;
                if(myi==(cpu-1)){
                    //special case, do everything from here to end to avoid rounding errors
                    start=myi*(int)(elements/cpu);
                    end=elements;
                    //attstart=myi*(int)(attnum/cpu);
                    //attend=attnum;
                }else{
                    start=myi*(int)(elements/cpu);
                    end=(myi+1)*(int)(elements/cpu);
                    //attstart=myi*(int)(attnum/cpu);
                    //attend=(myi+1)*(int)(attnum/cpu);
                }
                if(cluster2dbutton.isSelected()){
                    //cluster all in 2D
                    //System.out.println("cluster2D");
                    for(int i=start;i<end;i++){
                        //currmoveatt=getminattract(posarr[i],currmoveatt,minattract);
                        getminattract(posarr[i],currmoveatt,minattract);
                        movearr[i][0]+=currmoveatt[0];
                        movearr[i][1]+=currmoveatt[1];
                        for(j=elements;--j>=0;){
                            //currmoverep=getrepulse2d(posarr[i],posarr[j],currmoverep,repvalpow,repfactor,rand);
                            getrepulse2d(posarr[i],posarr[j],currmoverep,repvalpow,repfactor,rand);
                            movearr[i][0]+=currmoverep[0];
                            movearr[i][1]+=currmoverep[1];
                        }//end for j
                    }//end for i
                    for(int i=attnum;--i>=0;){
                        if(attvals[i].query>=start&&attvals[i].query<end){
                            //currmoveatt=getattract2d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow, attfactor);
                            getattract2d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow, attfactor);
                            movearr[attvals[i].query][0]+=currmoveatt[0];
                            movearr[attvals[i].query][1]+=currmoveatt[1];
                        }
                        if(attvals[i].hit>=start && attvals[i].hit<end){
                            //currmoveatt=getattract2d(posarr[attvals[i].hit],posarr[attvals[i].query],attvals[i].att,currmoveatt,attvalpow, attfactor);
                            getattract2d(posarr[attvals[i].hit],posarr[attvals[i].query],attvals[i].att,currmoveatt,attvalpow, attfactor);
                            movearr[attvals[i].hit][0]+=currmoveatt[0];
                            movearr[attvals[i].hit][1]+=currmoveatt[1];
                        }
                    }//end for i
                    //double totaldistsq;
                    for(int i=start;i<end;i++){
                        movearr[i][0]/=elements;
                        movearr[i][1]/=elements;
                        movearr[i][2]=0;
                        totaldist=java.lang.Math.sqrt((movearr[i][0]*movearr[i][0])+(movearr[i][1]*movearr[i][1]));
                        if(totaldist>maxmove){
                            movearr[i][0]*=maxmove/totaldist;
                            movearr[i][1]*=maxmove/totaldist;
                        }
                    }//end for i
                }else{
                    //cluster all in 3D
                    //System.out.println("cluster3D");
                    for(int i=start;i<end;i++){
                        //currmoveatt=getminattract(posarr[i],currmoveatt,minattract);
                        getminattract(posarr[i],currmoveatt,minattract);
                        movearr[i][0]+=currmoveatt[0];
                        movearr[i][1]+=currmoveatt[1];
                        movearr[i][2]+=currmoveatt[2];
                        for(j=0;j<elements;j++){
                            //currmoverep=getrepulse3d(posarr[i],posarr[j],currmoverep,repvalpow,repfactor,rand);
                            getrepulse3d(posarr[i],posarr[j],currmoverep,repvalpow,repfactor,rand);
                            movearr[i][0]+=currmoverep[0];
                            movearr[i][1]+=currmoverep[1];
                            movearr[i][2]+=currmoverep[2];
                        }//end for j
                    }//end for i
                    for(int i=attnum;--i>=0;){
                        if(attvals[i].query>=start&&attvals[i].query<end){
                            //currmoveatt=getattract3d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow, attfactor);
                            getattract3d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow, attfactor);
                            movearr[attvals[i].query][0]+=currmoveatt[0];
                            movearr[attvals[i].query][1]+=currmoveatt[1];
                            movearr[attvals[i].query][2]+=currmoveatt[2];
                        }
                        if(attvals[i].hit>=start && attvals[i].hit<end){
                            //currmoveatt=getattract3d(posarr[attvals[i].hit],posarr[attvals[i].query],attvals[i].att,currmoveatt,attvalpow, attfactor);
                            getattract3d(posarr[attvals[i].hit],posarr[attvals[i].query],attvals[i].att,currmoveatt,attvalpow, attfactor);
                            movearr[attvals[i].hit][0]+=currmoveatt[0];
                            movearr[attvals[i].hit][1]+=currmoveatt[1];
                            movearr[attvals[i].hit][2]+=currmoveatt[2];
                        }
                    }//end for i
                    //double totaldistsq;
                    for(int i=start;i<end;i++){
                        movearr[i][0]/=elements;
                        movearr[i][1]/=elements;
                        movearr[i][2]/=elements;
                        totaldist=java.lang.Math.sqrt((movearr[i][0]*movearr[i][0])+(movearr[i][1]*movearr[i][1])+(movearr[i][2]*movearr[i][2]));
                        if(totaldist>maxmove){
                            movearr[i][0]*=maxmove/totaldist;
                            movearr[i][1]*=maxmove/totaldist;
                            movearr[i][2]*=maxmove/totaldist;
                        }
                    }//end for i
                }//end clustering all in 3D
            }//end if not doselected
            this.done=true;
            //synchronized(parent){
            //    parent.notify();
            //}// end synchronized parent
            synchronized(syncon){
                syncon.notify();
            }// end synchronized parent
        }// end run
        
    }//end class getmovethread
    */
}//end class
