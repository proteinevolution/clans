package clans.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Vector;
import java.util.concurrent.CancellationException;

import javax.swing.SwingWorker;

import clans.gui.Shapes;
import clans.model.ClusterData;
import clans.model.SequenceGroup;
import clans.model.microarray.Replicates;
import clans.model.proteins.AminoAcidSequence;
import clans.model.proteins.MinimalAttractionValue;
import clans.model.proteins.MinimalHsp;

public class ClusterDataLoadHelper {

	public int autosaveIntervalMinutes = -1;
	public java.io.File file=null;
	public MinimalHsp[] blasthits=null;
    public AminoAcidSequence[] inaln=null;
    public MinimalAttractionValue[] attvals=null;
    public float[][] posarr=null;
    public float[] weights=null;
    public int ovalsize=4;
    public int dotsize=2;
    public int groupsize=4;
    public double minpval=1;
    public double cooling=1;
    public double currcool=1;
    public float attfactor=10;
    public int attvalpow=1;
    public float repfactor=10;
    public int repvalpow=1;
    public float dampening=0.2f;
    public float zoom=1;
    public double minattract=1;
    public float maxmove=0.1f;
    public double pval=-1;
    public double[][] rotmtx={{1,0,0},{0,1,0},{0,0,1}};//default
    public boolean cluster2d=false;
    public boolean showinfo=false;
    public boolean usescval=false;
    public boolean complexatt=true;
    public java.util.Vector<SequenceGroup> seqgroupsvec=new Vector<SequenceGroup>();
    public String blastpath="blastall -p blastp";
    public String formatdbpath="formatdb";
    public java.util.ArrayList<File> mapfiles=new java.util.ArrayList<File>();
    public java.util.ArrayList<File> lookupfiles=new java.util.ArrayList<File>();
    public java.util.Vector<Replicates> affyfiles=null;
    public java.awt.Color[] colorarr=null;
    public float[] colorcutoffs=null;
    public boolean usefoldchange=false;
    public boolean avgfoldchange=false;
    public String namesdmp_file=null;
    public String nodesdmp_file=null;
    public int rounds=0;
 
	public boolean parse_params_block(BufferedReader inread) throws IOException{
		String[] tmparr;
		String inline;

        while (((inline = inread.readLine()) != null) && (inline.equalsIgnoreCase("</param>") == false)) {
            //skip empty lines
            if (inline.length() == 0) {
                continue;
            }
            tmparr = inline.split("=");
            if (tmparr.length != 2) {
                System.err.println("ERROR reading savefile for line '" + inline + "'");
                return false;
            }
            if (tmparr[0].equalsIgnoreCase("autosaveinterval")) {
                try {
                    autosaveIntervalMinutes = Integer.parseInt(tmparr[1]);
                } catch (NumberFormatException ne) {
                    System.err.println("ERROR parsing integer from " + tmparr[1]);
                    return false;
                }
            } else if (tmparr[0].equalsIgnoreCase("pval")) {
                try {
                    pval = Double.parseDouble(tmparr[1]);
                } catch (NumberFormatException ne) {
                    System.err.println("ERROR parsing double from " + tmparr[1]);
                    return false;
                }
            } else if (tmparr[0].equalsIgnoreCase("usescval")) {
                if (tmparr[1].startsWith("t") || tmparr[1].startsWith("T")) {
                    usescval = true;
                } else {
                    usescval = false;
                }
            } else if (tmparr[0].equalsIgnoreCase("complexatt")) {
                if (tmparr[1].startsWith("t") || tmparr[1].startsWith("T")) {
                    complexatt = true;
                } else {
                    complexatt = false;
                }
            } else if (tmparr[0].equalsIgnoreCase("cooling")) {
                try {
                    cooling = Double.parseDouble(tmparr[1]);
                } catch (NumberFormatException ne) {
                    System.err.println("ERROR parsing double from " + tmparr[1]);
                    return false;
                }
            } else if (tmparr[0].equalsIgnoreCase("currcool")) {
                try {
                    currcool = Double.parseDouble(tmparr[1]);
                } catch (NumberFormatException ne) {
                    System.err.println("ERROR parsing double from " + tmparr[1]);
                    return false;
                }
            } else if (tmparr[0].equalsIgnoreCase("maxmove")) {
                try {
                    maxmove = Float.parseFloat(tmparr[1]);
                } catch (NumberFormatException ne) {
                    System.err.println("ERROR parsing float from " + tmparr[1]);
                    return false;
                }
            } else if (tmparr[0].equalsIgnoreCase("cluster2d")) {
                if (tmparr[1].equalsIgnoreCase("true") || tmparr[1].startsWith("t") || tmparr[1].startsWith("T")) {
                    cluster2d = true;
                } else {
                    cluster2d = false;
                }
            } else if (tmparr[0].equalsIgnoreCase("usefoldchange")) {
                if (tmparr[1].equalsIgnoreCase("true") || tmparr[1].startsWith("t") || tmparr[1].startsWith("T")) {
                    usefoldchange = true;
                } else {
                    usefoldchange = false;
                }
            } else if (tmparr[0].equalsIgnoreCase("avgfoldchange")) {
                if (tmparr[1].equalsIgnoreCase("true") || tmparr[1].startsWith("t") || tmparr[1].startsWith("T")) {
                    avgfoldchange = true;
                } else {
                    avgfoldchange = false;
                }
            } else if (tmparr[0].equalsIgnoreCase("showinfo")) {
                if (tmparr[1].equalsIgnoreCase("true") || tmparr[1].startsWith("t") || tmparr[1].startsWith("T")) {
                    showinfo = true;
                } else {
                    showinfo = false;
                }
            } else if (tmparr[0].equalsIgnoreCase("attfactor")) {
                try {
                    attfactor = Float.parseFloat(tmparr[1]);
                } catch (NumberFormatException ne) {
                    System.err.println("ERROR parsing float from " + tmparr[1]);
                    return false;
                }
            } else if (tmparr[0].equalsIgnoreCase("repfactor")) {
                try {
                    repfactor = Float.parseFloat(tmparr[1]);
                } catch (NumberFormatException ne) {
                    System.err.println("ERROR parsing float from " + tmparr[1]);
                    return false;
                }
            } else if (tmparr[0].equalsIgnoreCase("attvalpow")) {
                try {
                    attvalpow = (int) Float.parseFloat(tmparr[1]);
                } catch (NumberFormatException ne) {
                    System.err.println("ERROR parsing float from " + tmparr[1]);
                    return false;
                }
            } else if (tmparr[0].equalsIgnoreCase("repvalpow")) {
                try {
                    repvalpow = (int) Float.parseFloat(tmparr[1]);
                } catch (NumberFormatException ne) {
                    System.err.println("ERROR parsing float from " + tmparr[1]);
                    return false;
                }
            } else if (tmparr[0].equalsIgnoreCase("dampening")) {
                try {
                    dampening = Float.parseFloat(tmparr[1]);
                } catch (NumberFormatException ne) {
                    System.err.println("ERROR parsing float from " + tmparr[1]);
                    return false;
                }
            } else if (tmparr[0].equalsIgnoreCase("zoom")) {
                try {
                    zoom = Float.parseFloat(tmparr[1]);
                } catch (NumberFormatException ne) {
                    System.err.println("ERROR parsing float from " + tmparr[1]);
                    return false;
                }
            } else if (tmparr[0].equalsIgnoreCase("minattract")) {
                try {
                    minattract = Double.parseDouble(tmparr[1]);
                } catch (NumberFormatException ne) {
                    System.err.println("ERROR parsing double from " + tmparr[1]);
                    return false;
                }
            } else if (tmparr[0].equalsIgnoreCase("dotsize")) {
                try {
                    dotsize = (int) Float.parseFloat(tmparr[1]);
                } catch (NumberFormatException ne) {
                    System.err.println("ERROR parsing float from " + tmparr[1]);
                    return false;
                }
            } else if (tmparr[0].equalsIgnoreCase("ovalsize")) {
                try {
                    ovalsize = (int) Float.parseFloat(tmparr[1]);
                } catch (NumberFormatException ne) {
                    System.err.println("ERROR parsing float from " + tmparr[1]);
                    return false;
                }
            } else if (tmparr[0].equalsIgnoreCase("groupsize")) {
                try {
                    groupsize = (int) Float.parseFloat(tmparr[1]);
                } catch (NumberFormatException ne) {
                    System.err.println("ERROR parsing float from " + tmparr[1]);
                    return false;
                }
            } else if (tmparr[0].equalsIgnoreCase("blastpath")) {
                blastpath = tmparr[1];
            } else if (tmparr[0].equalsIgnoreCase("formatdbpath")) {
                formatdbpath = tmparr[1];
            } else if (tmparr[0].equalsIgnoreCase("functionmapfile")) {
                //old version: only one file should be present in these save-files
                java.io.File tmpfile = new java.io.File(tmparr[1]);
                if (tmpfile.canRead()) {
                    mapfiles.add(tmpfile);
                    lookupfiles.add(null);
                } else {
                    System.err.println("Warning: cannot read from '" + tmparr[1] + "' setting mapfile to null");
                }
            } else if (tmparr[0].equalsIgnoreCase("functionlookupfile")) {
                java.io.File tmpfile = new java.io.File(tmparr[1]);
                if (tmpfile.canRead()) {
                    lookupfiles.add(tmpfile);
                } else {
                    System.err.println("Warning: cannot read from '" + tmparr[1] + "' setting lookupfile to null");
                }
            } else if (tmparr[0].equalsIgnoreCase("colorarr")) {
                tmparr = tmparr[1].split("\\):");
                java.awt.Color[] colorarr = new java.awt.Color[tmparr.length];
                try {
                    String[] tmptmp;
                    int red, green, blue;
                    for (int i = 0; i < tmparr.length; i++) {
                        tmptmp = tmparr[i].substring(1).split(";");//remove the "(" and split the numbers
                        red = Integer.parseInt(tmptmp[0]);
                        green = Integer.parseInt(tmptmp[1]);
                        blue = Integer.parseInt(tmptmp[2]);
                        colorarr[i] = new java.awt.Color(red, green, blue);
                    }//end for i
                    this.colorarr = colorarr;
                } catch (NumberFormatException ne) {
                    System.err.println("Warning: unable to parse color array from '" + inline + "'");
                }
            } else if (tmparr[0].equalsIgnoreCase("colorcutoffs")) {
                tmparr = tmparr[1].split(";");
                float[] cutoffs = new float[tmparr.length];
                try {
                    for (int i = 0; i < tmparr.length; i++) {
                        cutoffs[i] = Float.parseFloat(tmparr[i]);
                    }//end for i
                    colorcutoffs = cutoffs;
                } catch (NumberFormatException ne) {
                    System.err.println("Warning: unable to parse float[] from '" + inline + "'");
                }
            }else if (tmparr[0].equalsIgnoreCase("namesdmp_file")) {
                namesdmp_file = tmparr[1];
            } else if (tmparr[0].equalsIgnoreCase("nodesdmp_file")) {
                nodesdmp_file = tmparr[1];
            } else if (tmparr[0].equalsIgnoreCase("rounds_done")) {
                try {
                    rounds = Integer.parseInt(tmparr[1]);
                } catch (NumberFormatException ne) {
                    System.err.println("ERROR parsing int from " + tmparr[1]);
                    return false;
                }
            }
            
        }//end while </param>
	return true;
	}
	
	public boolean parse_function_block(BufferedReader buffered_file_handle) throws IOException {
    	String[] tmp;
        String inline;
        while (((inline = buffered_file_handle.readLine()) != null) && (inline.equalsIgnoreCase("</function>") == false)) {
            //skip empty lines
            if (inline.length() == 0) {
                continue;
            }
            tmp = inline.split("';'");
            if (tmp.length != 2) {
                System.err.println("ERROR readung mapping data on line '" + inline + "'");
                return false;
                }
            
            java.io.File tmpfile = new java.io.File(tmp[0]);
            if (tmpfile.canRead()) {
                mapfiles.add(tmpfile);
            } else {
                System.err.println("WARNING: cannot read from '" + tmp[0] + "' skipping entry");
            }
            if (tmp[1].equalsIgnoreCase("NONE")) {
                lookupfiles.add(null);
            } else {
                tmpfile = new java.io.File(tmp[1]);
                if (tmpfile.canRead()) {
                    lookupfiles.add(tmpfile);
                } else {
                    System.err.println("WARNING: cannot read from '" + tmp[1] + "' skipping entry");
                }
            }

        }
        return true;

	}
	
	public boolean parse_affyfiles_block(BufferedReader buffered_file_handle) throws IOException {
		affyfiles = new Vector<Replicates>();
        Replicates rep = null;
        String[] tmparr;
        String inline;
        
        while (((inline = buffered_file_handle.readLine()) != null) && (inline.equalsIgnoreCase("</affyfiles>") == false)) {
            //skip empty lines
            if (inline.length() == 0) {
                continue;
            }
            inline = inline.trim();
            if (inline.equalsIgnoreCase("<")) {
                rep = new Replicates();
            } else if (inline.equalsIgnoreCase(">")) {
                affyfiles.add(rep);
            } else {
                tmparr = inline.split("=");
                if (tmparr[0].equalsIgnoreCase("replicates")) {
                    try {
                        rep.replicates = Integer.parseInt(tmparr[1]);
                    } catch (NumberFormatException ne) {
                        System.err.println("ERROR: unable to parse int from '" + tmparr[1] + "'");
                        return false;
                    }
                } else if (tmparr[0].equalsIgnoreCase("wtreplicates")) {
                    try {
                        rep.wtreplicates = Integer.parseInt(tmparr[1]);
                    } catch (NumberFormatException ne) {
                        System.err.println("ERROR: unable to parse int from '" + tmparr[1] + "'");
                        return false;
                    }
                } else if (tmparr[0].equalsIgnoreCase("name")) {
                    rep.name = tmparr[1];
                } else if (tmparr[0].equalsIgnoreCase("wtname")) {
                    rep.wtname = tmparr[1];
                } else if (tmparr[0].equalsIgnoreCase("abbreviation")) {
                    rep.abbreviation = tmparr[1];
                } else if (tmparr[0].equalsIgnoreCase("replicate")) {
                    tmparr = tmparr[1].split("';'");
                    int repnum = tmparr.length;
                    if (repnum != rep.replicates) {
                        System.err.println("unequal number of elements found in replicate; found:" + repnum + " expected:" + rep.replicates);
                    }
                    rep.replicate = new File[repnum];
                    for (int i = 0; i < repnum; i++) {
                        rep.replicate[i] = new File(tmparr[i]);
                        if (rep.replicate[i].canRead() == false) {
                            System.err.println("WARNING: cannot read file '" + tmparr[i] + "'");
                        }
                    }//end for i
                } else if (tmparr[0].equalsIgnoreCase("wtreplicate")) {
                    tmparr = tmparr[1].split("';'");
                    int repnum = tmparr.length;
                    if (repnum != rep.wtreplicates) {
                        System.err.println("unequal number of elements found in wtreplicate; found:" + repnum + " expected:" + rep.wtreplicates);
                    }
                    rep.wtreplicate = new File[repnum];
                    for (int i = 0; i < repnum; i++) {
                        rep.wtreplicate[i] = new File(tmparr[i]);
                        if (rep.wtreplicate[i].canRead() == false) {
                            System.err.println("WARNING: cannot read file '" + tmparr[i] + "'");
                        }
                    }//end for i
                }
            }
        }
        
        return true;
	}
	
	public boolean parse_rotmtx_block(BufferedReader buffered_file_handle) throws IOException {
		String tmpstr = "";
		String inline;
        while (((inline = buffered_file_handle.readLine()) != null) && (inline.equalsIgnoreCase("</rotmtx>") == false)) {
            //skip empty lines
            if (inline.length() == 0) {
                continue;
            }
            tmpstr += inline;
        }
        String[] tmparr = tmpstr.split(";", 0);
        if (tmparr.length != 9) {
            System.err.println("ERROR reading rotation matrix values from '" + tmpstr + "'");
            return false;
            }
        
        double[][] temp_rotmtx = new double[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                try {
                    temp_rotmtx[i][j] = Double.parseDouble(tmparr[(i * 3) + j]);
                } catch (NumberFormatException ne) {
                    System.err.println("unable to parse double from '" + tmparr[(i * 3) + j] + "'");
                    return false;
                }
            }
        }

        rotmtx = temp_rotmtx;
        return true;
	}
	
	public boolean parse_seq_block(BufferedReader buffered_file_handle, int expected_sequences) throws IOException {
		int counter = 0;
        String currname = "";
        String currseq = "";
        String inline;
        
        while (((inline = buffered_file_handle.readLine()) != null) && (inline.equalsIgnoreCase("</seq>") == false)) {
            //skip empty lines
            if (inline.length() == 0) {
                continue;
            }
            //read the sequence data
            if (inline.startsWith(">")) {
                if (currname.length() > 0) {
                    inaln[counter] = new AminoAcidSequence();
                    inaln[counter].name = currname;
                    inaln[counter].seq = currseq;
                    counter++;
                }
                currname = inline.substring(1).trim();
                currseq = "";
            } else {
                currseq += inline;
            }
        }
        
        inaln[counter] = new AminoAcidSequence();
        inaln[counter].name = currname;
        inaln[counter].seq = currseq;
        counter++;
        if (counter != expected_sequences) {
            System.err.println("ERROR, not found the number of specified sequences, expected:" + expected_sequences + " found:" + counter);
            return false;
        }
        return true;
	}

	public boolean parse_weight_block(BufferedReader buffered_file_handle, int expected_sequences) throws IOException {
		String currname = "";
        String currseq = "";
        String inline;
        
        if (expected_sequences < 1) {
            System.err.println("ERROR; <weight>...</weight> statements have to come AFTER <seq>...</seq> statement");
            javax.swing.JOptionPane.showMessageDialog(new javax.swing.JFrame(), "ERROR; <weight>...</weight> statements have to come AFTER <seq>...</seq> statement");
            return false;
        }
        weights = new float[expected_sequences];
        //now make a Hash that links the weights values to the sequence names
        HashMap<String, Float> tmphash = new HashMap<String, Float>();
        //and now read the data
        while (((inline = buffered_file_handle.readLine()) != null) && (inline.equalsIgnoreCase("</weight>") == false)) {
            //skip empty lines
            if (inline.length() == 0) {
                continue;
            }
            //read the sequence data
            if (inline.startsWith(">")) {
                if (currname.length() > 0) {
                    if (tmphash.containsKey(currname)) {
                        try {
                            tmphash.put(currname, new Float(Float.parseFloat(currseq.trim())));
                        } catch (NumberFormatException ne) {
                            System.err.println("unable to parse number from '" + currseq + "'");
                            return false;
                        }
                    } else {//if this name is unknown
                        System.err.println("WARNING: name '" + currname + "' is not a sequence name.");
                    }
                }
                currname = inline.substring(1).trim();
                currseq = "";
            } else {
                currseq += inline;
            }
        }//end seq part
        if (currname.length() > 0) {
            if (tmphash.containsKey(currname)) {
                try {
                    tmphash.put(currname, new Float(Float.parseFloat(currseq.trim())));
                } catch (NumberFormatException ne) {
                    System.err.println("unable to parse number from '" + currseq + "'");
                    return false;
                }
            } else {//if this name is unknown
                System.err.println("WARNING: name '" + currname + "' is not a sequence name.");
            }
        }
        //and now convert the tmphash data to a weights array
        //for that loop through the sequence names and assign the corresponding weights
        for (int i = 0; i < expected_sequences; i++) {
            if (tmphash.containsKey(inaln[i].name)) {
                weights[i] = ((Float) tmphash.get(inaln[i].name)).floatValue();
            } else {//default weight of 1
                weights[i] = 1;
            }
        }//end for i
        
        return true;
	}

	public boolean parse_pos_block(BufferedReader buffered_file_handle, int expected_sequences) throws IOException {
		int counter = 0;
        int mypos;
        String[] tmparr;
        String inline;
        
        while (((inline = buffered_file_handle.readLine()) != null) && (inline.equalsIgnoreCase("</pos>") == false)) {
            //skip empty lines
            if (inline.length() == 0) {
                continue;
            }
            tmparr = inline.trim().split("\\s+");
            if (tmparr.length != 4) {
                System.err.println("ERROR reading positions from " + inline);
                return false;
            }
            //else
            try {
                mypos = Integer.parseInt(tmparr[0]);
                posarr[mypos][0] = Float.parseFloat(tmparr[1]);
                posarr[mypos][1] = Float.parseFloat(tmparr[2]);
                posarr[mypos][2] = Float.parseFloat(tmparr[3]);
                counter++;
            } catch (NumberFormatException ne) {
                System.err.println("ERROR, unable to parse int,double,double,double from " + inline);
                return false;
            }
        }//end while pos
        if (counter != expected_sequences) {
            System.err.println("ERROR, not found the necessary number of positions");
            return false;
        }		
        return true;
	}
	
	/**
	 * Convenience method that calls {@code parse_hsp_block(BufferedReader, SwingWorker)} with SwingWorker {@code null}. This
	 * method should be used in no-GUI mode.
	 */
	public boolean parse_hsp_block(BufferedReader buffered_file_handle) throws IOException{
		return parse_hsp_block(buffered_file_handle, null);
	}

	/**
	 * Parses a CLANS file HSP block.
	 * 
	 * @param buffered_file_handle
	 *            The file handle with the first HSP block line ready.
	 * @param worker
	 *            If not null, this worker is used on a regular basis to check whether loading should be canceled. This
	 *            is used by the GUI to act on user aborts.
	 * @return true if the parsing was completed without errors, false else.
	 * @throws IOException
	 * @throws CancellationException
	 *             If worker is cancelled and loading should therefore be cancelled. In headless (no-GUI) mode this
	 *             Exception will not occur.
	 */
	public boolean parse_hsp_block(BufferedReader buffered_file_handle, SwingWorker<Void, Integer> worker)
            throws IOException, CancellationException {
        String[] split_result;
        String tmpstr;
        HashMap<MinimalHsp, MinimalHsp> hspHash = new HashMap<MinimalHsp, MinimalHsp>();
        String inline;
        int characters_for_round_number = 9;
        String format_for_round_number = "%" + characters_for_round_number + "s";
        
        if (pval < 0) {
            pval = 1;
        }
        
        int count = 0;
        String lastline = "";
        while (((inline = buffered_file_handle.readLine()) != null) && (inline.equalsIgnoreCase("</hsp>") == false)) {

            // skip empty lines
            if (inline.length() == 0) {
                continue;
            }

            count++;

            if (count % 10000 == 0) {
                System.out.print(".");

                // every 10000 HSPs as tradeoff between time wasted on checking and responsiveness of the cancelation
                ClusterData.checkWorkerStatus(worker); // stop saving if told so by user (only in GUI mode)

                if (count % 500000 == 0) {
                    System.out.println(String.format(format_for_round_number, count));
                }
            }

            split_result = inline.split(":");
            if (inline.equalsIgnoreCase("DONE")) {
                System.err.println("manually truncated file, returning results so far");
                return false;
            }
            
            if (split_result.length != 2) {
                System.err.println("ERROR parsing HSP data from " + inline + " right after line '" + lastline + "'");
                return false;
            }

            tmpstr = split_result[1]; // String with one or multiple p-values separated by any whitespace characters
            split_result = split_result[0].split("\\s+"); // the query and subject identification numbers as Strings
            
            if (split_result.length != 2) {
                System.err.println("ERROR, wrong hsp line " + inline);
                return false;
            }
            
            try {
                int myi = Integer.parseInt(split_result[0]);
                int myj = Integer.parseInt(split_result[1]);

                if (myi == myj) { // Connections to itself are useless
                    continue;
                }

                split_result = tmpstr.split("\\s+");
                if (split_result.length == 0) {
                    continue;
                }

                MinimalHsp newHSP = new MinimalHsp(myi, myj);
                if (hspHash.containsKey(newHSP)) {
                    MinimalHsp oldHSP = hspHash.get(newHSP);
                    for (int i = 0; i < split_result.length; i++) {
                        oldHSP.addpval(Double.parseDouble(split_result[i]));
                    }
                } else {
                    newHSP.val = new double[split_result.length];
                    for (int i = 0; i < split_result.length; i++) {
                        newHSP.val[i] = Double.parseDouble(split_result[i]);
                    }
                    hspHash.put(newHSP, newHSP);
                }
            } catch (NumberFormatException Ne) {
                System.err.println("ERROR, unable to parse int int: double ... from " + inline);
                return false;
            }
            
            lastline = inline;
        }
        
        System.out.println(String.format(format_for_round_number, count));
        System.out.println("total HSPs: " + count);
        
        if (inline.equalsIgnoreCase("</hsp>")) {
            blasthits = hspHash.values().toArray(new MinimalHsp[0]);
        } else {
            System.err.println("ERROR reading truncated file; <hsp> tag not closed by </hsp>");
        }
        
        return true;
	}

	public boolean parse_att_block(BufferedReader buffered_file_handle) throws IOException {
		blasthits = null;
        Vector<MinimalAttractionValue> tmpvec = new Vector<MinimalAttractionValue>();
        String[] tmparr;
        String tmpstr = "";
        String inline;

        System.out.println("reading attraction value data");
        
        if (pval < 0) {
            pval = 0;
        }
        
        int elements, counter = 0;
        while (((inline = buffered_file_handle.readLine()) != null) && (inline.equalsIgnoreCase("</att>") == false)) {
            //skip empty lines
            if (inline.length() == 0) {
                continue;
            }
            inline = inline.trim();
            tmparr = inline.split("[\\s+]");
            if (counter % 1000 == 0) {
                System.out.print(counter);
            } else if (counter % 10 == 0) {
                System.out.print(".");
            }
            elements = tmparr.length;
            if (elements != 3) {
                System.err.println("unable to parse attraction value information from '" + inline + "'");
            }
            try {
                int query = Integer.parseInt(tmparr[0]);
                int hit = Integer.parseInt(tmparr[1]);
                float att = Float.parseFloat(tmparr[2]);
                tmpvec.addElement(new MinimalAttractionValue(query, hit, att));
            } catch (NumberFormatException ne) {
                System.err.println("unable to parse float for " + tmpstr + " in " + inline);
            }
            counter++;
        }
        attvals = new MinimalAttractionValue[tmpvec.size()];
        tmpvec.copyInto(attvals);
        System.out.println("done reading matrix data");
        return true;
	}

	public boolean parse_seqgroups_block(BufferedReader buffered_file_handle) throws IOException {
		 //while I am reading the sequence groups
        String[] tmparr;
        String tmpstr = "";
        SequenceGroup currgroup = null;
        int tmpval;
        String inline;
        
        while (((inline = buffered_file_handle.readLine()) != null) && (inline.equalsIgnoreCase("</seqgroups>") == false)) {
            //skip empty lines
            if (inline.length() == 0) {
                continue;
            }
            tmparr = inline.split("=", 2);
            if (tmparr.length != 2) {
                System.err.println("ERROR reading from savefile on line '" + inline + "'");
                return false;
            }
            if (tmparr[0].equalsIgnoreCase("name")) {
                if (currgroup != null) {
                    seqgroupsvec.addElement(currgroup);
                }
                currgroup = new SequenceGroup();
                currgroup.name = tmparr[1];
                currgroup.type = 0;
            } else if (inline.startsWith("type=")) {
                try {
                    tmpstr = inline.substring(5).trim();
                    currgroup.type = Integer.parseInt(tmpstr);
                } catch (NumberFormatException ne) {
                    System.err.println("ERROR, unable to parse int from '" + tmpstr + "'");
                }
            } else if (inline.startsWith("size=")) {
                try {
                    tmpstr = inline.substring(5).trim();
                    currgroup.size = Integer.parseInt(tmpstr);
                } catch (NumberFormatException ne) {
                    System.err.println("ERROR, unable to parse int from '" + tmpstr + "'");
                }
            } else if (inline.startsWith("hide=")) {
                try {
                    tmpstr = inline.substring(5).trim();
                    tmpval = Integer.parseInt(tmpstr);
                    if (tmpval == 1) {
                        currgroup.hide = true;
                    } else {
                        currgroup.hide = false;
                    }
                } catch (NumberFormatException ne) {
                    System.err.println("ERROR, unable to parse int from '" + tmpstr + "'");
                }
            } else if (tmparr[0].equalsIgnoreCase("color")) {
                tmparr = tmparr[1].split(";");
                try {
                    int red = Integer.parseInt(tmparr[0]);
                    int green = Integer.parseInt(tmparr[1]);
                    int blue = Integer.parseInt(tmparr[2]);
                    if(tmparr.length<4){
                        currgroup.color = new java.awt.Color(red, green, blue);
                    }else{
                        int alpha=Integer.parseInt(tmparr[3]);
                        currgroup.color = new java.awt.Color(red, green, blue,alpha);
                    }
                } catch (NumberFormatException ne) {
                    System.err.println("ERROR parsing numbers from '" + inline + "'");
                    return false;
                }

            } else if (tmparr[0].equalsIgnoreCase("numbers")) {
                
                // allow empty groups
                if (tmparr[1].length() == 0) {
                    continue;
                }
                
                tmparr = tmparr[1].split(";");
                int num = tmparr.length;
                int[] retarr = new int[num];

                try {
                    for (int i = 0; i < num; i++) {
                        retarr[i] = Integer.parseInt(tmparr[i]);
                    }
                } catch (NumberFormatException ne) {
                    System.err.println("ERROR parsing numbers from '" + inline + "'");
                    return false;
                }
                
                currgroup.sequences = retarr;

                // as "numbers" is the last property of seqgroups, we now also create the groups icon
                currgroup.polygon = Shapes.get(currgroup.type, currgroup.size);
                
            } else {
                System.err.println("Error reading savefile in line" + inline);
                return false;
            }
        }

        if (currgroup != null) {
            seqgroupsvec.addElement(currgroup);
        }
        
        return true;
	}

	public boolean parse_mtx_block(BufferedReader buffered_file_handle, int expected_sequences) throws IOException {
		//load a matrix of attraction values
        blasthits = null;
        System.out.println("reading matrix");
        int counter = 0;
        String[] tmparr;
        HashMap<MinimalAttractionValue, MinimalAttractionValue> tmphash = new HashMap<MinimalAttractionValue, MinimalAttractionValue>();
        String inline;
        
        while (((inline = buffered_file_handle.readLine()) != null) && (inline.equalsIgnoreCase("</mtx>") == false)) {
            //skip empty lines
            if (inline.length() == 0) {
                continue;
            }
            tmparr = inline.trim().split(";");
            if (tmparr.length != expected_sequences) {
                System.err.println("ERROR reading positions from " + inline + "; expecting " + expected_sequences + " values");
                return false;
            }
            try {
                for (int i = 0; i < expected_sequences; i++) {
                    float tmpval = Float.parseFloat(tmparr[i]);
                    if (tmpval != 0) {
                        MinimalAttractionValue curratt = new MinimalAttractionValue(counter, i);
                        if (tmphash.containsKey(curratt)) {
                            tmphash.get(curratt).att += tmpval / 2;
                        } else {
                            curratt.att = tmpval / 2;
                            tmphash.put(curratt, curratt);
                        }
                    }
                }//end for i
            } catch (NumberFormatException ne) {
                System.err.println("ERROR, unable to parse float array from " + inline);
                return false;
            }
            counter++;
        }//end while pos
        attvals = (MinimalAttractionValue[]) (tmphash.values().toArray(new MinimalAttractionValue[0]));
        System.out.println("done reading matrix:" + counter);
        if (counter != expected_sequences) {
            System.err.println("ERROR, not found the necessary number of positions");
            return false;
        }
        
        return true;		
	}	
}
