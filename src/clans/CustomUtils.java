/*
 * customutils.java
 *
 * Created on September 19, 2003, 3:43 PM
 */
package clans;

import java.io.*;
import java.util.*;

/**
 *
 * @author  tancred
 */
public class CustomUtils {

    /** Creates a new instance of customutils */
    public CustomUtils() {
    }
    //--------------------------------------------------------------------------

    static Vector loadgroups(Vector groupsvec, String[] namearr, File infile) {
        //read the seqgroup information from a separate file.
        //the file has to either contain sequence names or numbers. if numbers,
        //the names of the sequences have to be specified as well.
        //read the groups, assign sequence names where necessary and find those names in namearr
        //then assign groups based on the position of the sequence in namearr
        try {
            BufferedReader inread = new BufferedReader(new FileReader(infile));
            String firstline = inread.readLine();
            inread.close();
            if (firstline.startsWith("<ids>")) {
                //I want to read cluster data not group data
                System.out.println("reading cluster data");
                return loadclusters(groupsvec, namearr, infile);
            }//else do the usual read
            inread = new BufferedReader(new FileReader(infile));
            Vector tmpvec = new Vector();
            HashMap mynamehash = new HashMap();
            String[] tmparr;
            String tmpstr = "";
            int elements = java.lang.reflect.Array.getLength(namearr);
            //System.out.println("parentnames="+elements);
            for (int i = 0; i < elements; i++) {
                tmparr = namearr[i].split("\\s+");
                mynamehash.put(tmparr[0], new Integer(i));
                //System.out.println("adding '"+tmparr[0]+" as "+i);
            }//end for i
            String[] mynamearr = null;
            int maxnum = 0;
            String inline;
            while ((inline = inread.readLine()) != null) {
                //skip empty lines
                if (inline.length() == 0) {
                    continue;
                }
                if (inline.equalsIgnoreCase("<seq>")) {
                    //read the sequence names used in this file
                    while ((inline = inread.readLine()) != null) {
                        //skip empty lines
                        if (inline.length() == 0) {
                            continue;
                        }
                        if (inline.equalsIgnoreCase("</seq>")) {
                            //exit this while loop
                            break;
                        } else if (inline.startsWith(">")) {
                            tmpstr = inline.substring(1).trim();
                            tmparr = tmpstr.split("\\s+");
                            tmpvec.addElement(tmparr[0]);
                        }
                    }
                    maxnum = tmpvec.size();
                    mynamearr = new String[maxnum];
                    tmpvec.copyInto(mynamearr);
                } else if (inline.equalsIgnoreCase("<seqgroups>")) {
                    //read the groups defined in this file; format is CLANS
                    if (maxnum == 0) {
                        System.err.println("ERROR; missing sequence names");
                        return groupsvec;
                    }
                    int count = 0;
                    java.awt.Color color = java.awt.Color.red;
                    String name = null;
                    int type = 0;
                    int size = 0, tmpval;
                    boolean hide = false;
                    while ((inline = inread.readLine()) != null) {
                        //skip empty lines
                        if (inline.length() == 0) {
                            continue;
                        }
                        if (inline.equalsIgnoreCase("</seqgroups>")) {
                            //exit this while loop
                            break;
                        }
                        if (inline.startsWith("name=")) {
                            name = inline.substring(5).trim();
                            type = 0;
                        } else if (inline.startsWith("type=")) {
                            try {
                                tmpstr = inline.substring(5).trim();
                                type = Integer.parseInt(tmpstr);
                            } catch (NumberFormatException ne) {
                                System.err.println("ERROR, unable to parse int from '" + tmpstr + "'");
                            }
                        } else if (inline.startsWith("size=")) {
                            try {
                                tmpstr = inline.substring(5).trim();
                                size = Integer.parseInt(tmpstr);
                            } catch (NumberFormatException ne) {
                                System.err.println("ERROR, unable to parse int from '" + tmpstr + "'");
                            }
                        } else if (inline.startsWith("hide=")) {
                            try {
                                tmpstr = inline.substring(5).trim();
                                tmpval = Integer.parseInt(tmpstr);
                                if (tmpval == 1) {
                                    hide = true;
                                } else {
                                    hide = false;
                                }
                            } catch (NumberFormatException ne) {
                                System.err.println("ERROR, unable to parse int from '" + tmpstr + "'");
                            }
                        } else if (inline.startsWith("numbers=")) {
                            count++;
                            seqgroup mygroup = new seqgroup();
                            if (name != null) {
                                mygroup.name = name;
                                name = null;
                            } else {
                                mygroup.name = String.valueOf(count);
                            }
                            mygroup.type = type;
                            mygroup.color = color;
                            mygroup.hide = hide;
                            mygroup.size = size;
                            mygroup.polygon = makepolygons.get(mygroup.type, mygroup.size);
                            color = java.awt.Color.red;
                            //now convert the numbers you read to the numbers used in namearr
                            tmparr = inline.substring(8).trim().split(";");
                            elements = java.lang.reflect.Array.getLength(tmparr);
                            int[] tmp = new int[elements];
                            //System.out.println("in:"+inline);
                            try {
                                for (int i = 0; i < elements; i++) {
                                    tmpstr = tmparr[i];
                                    tmp[i] = Integer.parseInt(tmpstr);
                                }//end for i
                            } catch (NumberFormatException ne) {
                                System.err.println("Unable to parse int from '" + tmpstr + "' in line '" + inline + "'");
                                continue;
                            }
                            tmpvec.clear();
                            for (int i = 0; i < elements; i++) {
                                tmpstr = mynamearr[tmp[i]];
                                //System.out.print(tmpstr+"/");
                                if (mynamehash.containsKey(tmpstr)) {
                                    tmpvec.addElement(mynamehash.get(tmpstr));
                                    //System.out.println("\t"+tmpstr+"-->"+((Integer)mynamehash.get(tmpstr)).intValue()+":"+namearr[((Integer)mynamehash.get(tmpstr)).intValue()]+";");
                                } else {
                                    System.err.println("no correspondence found for '" + tmpstr + "' in current graph; skipping entry " + java.lang.reflect.Array.getLength(namearr));
                                }
                            }//end for i
                            //System.out.println();
                            elements = tmpvec.size();
                            mygroup.sequences = new int[elements];
                            //System.out.print("out:");
                            for (int i = 0; i < elements; i++) {
                                mygroup.sequences[i] = ((Integer) tmpvec.elementAt(i)).intValue();
                                //System.out.print(mygroup.sequences[i]+";");
                            }//end for i
                            //System.out.println();
                            if (elements > 0) {
                                groupsvec.addElement(mygroup);
                            }
                        } else if (inline.startsWith("color=")) {
                            tmpstr = inline.substring(6).trim();
                            tmparr = tmpstr.split(";");
                            if ((java.lang.reflect.Array.getLength(tmparr) != 3) && (java.lang.reflect.Array.getLength(tmparr) != 4)) {
                                System.err.println("ERROR reading color info from line '" + inline + "'");
                            }
                            int red, green, blue;
                            try {
                                tmpstr = tmparr[0];
                                red = Integer.parseInt(tmpstr);
                                tmpstr = tmparr[1];
                                green = Integer.parseInt(tmpstr);
                                tmpstr = tmparr[2];
                                blue = Integer.parseInt(tmpstr);
                                if(java.lang.reflect.Array.getLength(tmparr)<4){
                                    color = new java.awt.Color(red, green, blue);
                                }else{
                                    tmpstr=tmparr[3];
                                    int alpha=Integer.parseInt(tmpstr);
                                    color=new java.awt.Color(red, green, blue, alpha);
                                }
                            } catch (NumberFormatException ne) {
                                System.err.println("unable to convert '" + tmpstr + "' to int in line '" + inline + "'");
                                color = java.awt.Color.red;
                            }
                        }
                    }
                }
            }//end while
            inread.close();
        } catch (IOException ioe) {
            System.err.println("IOERROR; unable to read from file '" + infile.getAbsolutePath() + "'");
        }
        return groupsvec;
    }//end loadgroups

    static Vector loadclusters(Vector groupsvec, String[] namesarr, File infile) {
        //read the cluster data generated by the iterative clustering approach
        //(append groups to groupsvec)
        //System.out.println("in loadcluster groupsvec.size=" + groupsvec.size());
        HashMap<String, Integer> mynamehash = new HashMap<String, Integer>();
        String[] tmparr;
        String tmpstr = "";
        int elements = java.lang.reflect.Array.getLength(namesarr);
        //System.out.println("parentnames="+elements);
        for (int i = 0; i < elements; i++) {
            tmparr = namesarr[i].split("\\s+");
            if (mynamehash.containsKey(tmparr[0])) {
                System.err.println("ERROR: non unique identifier in parent: name '" + tmparr[0] + "' already defined as -->'" + mynamehash.get(tmparr[0]) + "'");
                return groupsvec;
            } else {
                mynamehash.put(tmparr[0], new Integer(i));
            }
        }//end for i
        //System.out.println("parent short names:" + mynamehash.size());
        try {
            BufferedReader inread = new BufferedReader(new FileReader(infile));
            String inline;
            HashMap<String, Integer> clusternamehash = new HashMap<String, Integer>();
            String numberstring = "";
            String namestring = "";
            while ((inline = inread.readLine()) != null) {
                if (inline.equalsIgnoreCase("<ids>")) {
                    System.out.println("reading cluster IDS");
                    boolean stopids = false;
                    while (stopids == false) {
                        inline = inread.readLine();
                        if (inline == null) {
                            System.err.println("ERROR, reached EOF before end of <ids>");
                            return groupsvec;
                        } else if (inline.equalsIgnoreCase("</ids>")) {
                            stopids = true;
                        } else {
                            tmparr = inline.split("\\s+", 2);
                            numberstring = tmparr[0];
                            namestring = tmparr[1];
                            int spaceidx = namestring.indexOf(" ", 1);
                            int quoteidx = namestring.indexOf("\"", 1);
                            if (spaceidx < quoteidx) {
                                namestring = namestring.substring(1, spaceidx);
                            } else {
                                namestring = namestring.substring(1, quoteidx);
                            }
                            if (mynamehash.containsKey(namestring)) {
                                clusternamehash.put(numberstring, mynamehash.get(namestring));
                            } else {
                                System.err.println("WARNING: name '" + namestring + "' is not known in parent");
                                clusternamehash.put(numberstring, null);
                            }
                        }
                    }//end while stopids==false
                    System.out.println("DONE reading cluster IDS; found:" + clusternamehash.size());
                } else if (inline.equalsIgnoreCase("<clusters>")) {
                    tmpstr = inread.readLine();//next line contains offset=number (substring position 7)
                    float offset = 0;
                    try {
                        tmpstr = tmpstr.substring(7);
                        offset = Float.parseFloat(tmpstr);
                    } catch (NumberFormatException ne) {
                        System.err.println("ERROR trying to parse float from '" + tmpstr + "'");
                        return groupsvec;
                    }
                    System.out.println("reading clusters; offset=" + offset);
                    boolean stopclusters = false;
                    while (stopclusters ==false) {
                        inline = inread.readLine();
                        //System.out.println("reading clusters line '"+inline+"'");
                        if (inline == null) {
                            System.err.println("ERROR, reached EOF before end of <clusters>");
                            return groupsvec;
                        } else if (inline.equalsIgnoreCase("</clusters>")) {
                            //System.out.println("STOPCLUSTERS");
                            stopclusters = true;
                        } else if (inline.equalsIgnoreCase("<cluster>")) {
                            //System.out.print(".");
                            seqgroup mycluster = new seqgroup();
                            int clustersize = -1;
                            String nameline = inread.readLine();
                            String sizeline = inread.readLine();
                            String elementline = inread.readLine();
                            tmpstr = inread.readLine();//this should be </cluster>
                            if ((tmpstr.equalsIgnoreCase("</cluster>") == false) || (nameline.startsWith("name=") == false) || (sizeline.startsWith("size=") == false) || (elementline.startsWith("elements=") == false)) {
                                System.err.println("ERROR parsing data from cluster:");
                                System.err.println("<cluster>");
                                System.err.println(nameline);
                                System.err.println(sizeline);
                                System.err.println(elementline);
                                System.err.println(tmpstr);
                            } else {//all is well
                                //System.out.println("read cluster:");
                                //System.out.println("\t"+nameline);
                                //System.out.println("\t"+sizeline);
                                //System.out.println("\t"+elementline);
                                mycluster.name = nameline.substring(5) + "_of:" + offset;
                                try {
                                    clustersize = Integer.parseInt(sizeline.substring(5));
                                } catch (NumberFormatException ne) {
                                    System.err.println("ERROR trying to parse int from line '" + sizeline + "'");
                                    return groupsvec;
                                }
                                tmparr = elementline.substring(9).split(";");
                                if (java.lang.reflect.Array.getLength(tmparr) != clustersize) {
                                    System.err.println("WARNING: unequal cluster size: should=" + clustersize + " is=" + java.lang.reflect.Array.getLength(tmparr));
                                    clustersize = java.lang.reflect.Array.getLength(tmparr);
                                }
                                mycluster.sequences = new int[clustersize];
                                Integer myint;
                                for (int i = 0; i < clustersize; i++) {
                                    if ((myint = clusternamehash.get(tmparr[i])) != null) {
                                        mycluster.sequences[i] = myint.intValue();
                                    } else {
                                        System.err.println("WARNING: undefined name for '" + tmparr[i] + "'");
                                    }
                                }//end for i
                                //System.out.println("adding cluster '"+mycluster.name+"' with "+java.lang.reflect.Array.getLength(mycluster.sequences)+" elements");
                                groupsvec.add(mycluster);
                            }
                        }else{
                            System.out.println("read unknown line '"+inline+"'");
                        }
                    }//end while stopclusters==false
                }//else don't do anything
            }//end while reading file
            inread.close();
        } catch (IOException ioe) {
            System.err.println("IOERROR; unable to read from file '" + infile.getAbsolutePath() + "'");
        }
        System.out.println("read " + groupsvec.size() + " cluster entries");
        return groupsvec;
    }//end loadclusters

    //--------------------------------------------------------------------------
    static void saverun(saverunobject in, String[] namesarr,boolean nographics) {
        //save tha data to the selected file
        int i, j;
        try {
            PrintWriter outwrite = new PrintWriter(new BufferedWriter(new FileWriter(in.file)));
            int seqs = java.lang.reflect.Array.getLength(in.inaln);
            outwrite.println("sequences=" + seqs);
            outwrite.println("<param>");
            outwrite.println("maxmove=" + in.maxmove);
            outwrite.println("pval=" + in.pval);
            outwrite.println("usescval=" + in.usescval);
            outwrite.println("complexatt=" + in.complexatt);
            outwrite.println("cooling=" + in.cooling);
            outwrite.println("currcool=" + in.currcool);
            outwrite.println("attfactor=" + in.attfactor);
            outwrite.println("attvalpow=" + in.attvalpow);
            outwrite.println("repfactor=" + in.repfactor);
            outwrite.println("repvalpow=" + in.repvalpow);
            outwrite.println("dampening=" + in.dampening);
            outwrite.println("minattract=" + in.minattract);
            outwrite.println("cluster2d=" + in.cluster2d);
            outwrite.println("blastpath=" + in.blastpath);
            outwrite.println("formatdbpath=" + in.formatdbpath);
            outwrite.println("showinfo=" + in.showinfo);
            outwrite.println("zoom=" + in.zoom);
            outwrite.println("dotsize=" + in.dotsize);
            outwrite.println("ovalsize=" + in.ovalsize);
            outwrite.println("groupsize=" + in.groupsize);
            outwrite.println("usefoldchange=" + in.usefoldchange);
            outwrite.println("avgfoldchange=" + in.avgfoldchange);
            
            if(in.namesdmp_file!=null && in.nodesdmp_file!=null){
                outwrite.println("namesdmp_file="+in.namesdmp_file);
                outwrite.println("nodesdmp_file="+in.nodesdmp_file);
            }
            
            outwrite.print("colorcutoffs=");
            for (i = 0; i < java.lang.reflect.Array.getLength(in.colorcutoffs); i++) {
                outwrite.print(in.colorcutoffs[i] + ";");
            }//end for i
            outwrite.println();
            
            outwrite.print("colorarr=");
            for (i = 0; i < java.lang.reflect.Array.getLength(in.colorarr); i++) {
                java.awt.Color tmp = in.colorarr[i];
                outwrite.print("(" + tmp.getRed() + ";" + tmp.getGreen() + ";" + tmp.getBlue() + "):");
            }//end for i
            outwrite.println();
            
            outwrite.println("rounds_done=" + in.rounds);
            outwrite.println("</param>");
            
            if ((in.mapfiles!=null) && (in.mapfiles.size() > 0)) {
                outwrite.println("<function>");
                int num = in.mapfiles.size();
                for (i = 0; i < num; i++) {
                    if ((in.lookupfiles!=null) && (in.lookupfiles.get(i) != null)) {
                        outwrite.println(((File) in.mapfiles.get(i)).getAbsolutePath() + "';'" + ((File) in.lookupfiles.get(i)).getAbsolutePath());
                    } else {
                        outwrite.println(((File) in.mapfiles.get(i)).getAbsolutePath() + "';'NONE");
                    }
                }//end for i
                outwrite.println("</function>");
            }
            if (in.affyfiles != null) {
                outwrite.println("<affyfiles>");
                int repnum = in.affyfiles.size();
                replicates rep;
                for (i = 0; i < repnum; i++) {
                    rep = (replicates) (in.affyfiles.get(i));
                    outwrite.println("<");
                    outwrite.println("abbreviation=" + rep.abbreviation);
                    outwrite.println("replicates=" + rep.replicates);
                    outwrite.println("wtreplicates=" + rep.wtreplicates);
                    outwrite.println("name=" + rep.name);
                    outwrite.println("wtname=" + rep.wtname);
                    outwrite.print("replicate=");
                    for (j = java.lang.reflect.Array.getLength(rep.replicate); --j >= 0;) {
                        outwrite.print(rep.replicate[j].getAbsolutePath() + "';'");
                    }//end for j
                    outwrite.println();
                    outwrite.print("wtreplicate=");
                    for (j = java.lang.reflect.Array.getLength(rep.wtreplicate); --j >= 0;) {
                        outwrite.print(rep.wtreplicate[j].getAbsolutePath() + "';'");
                    }//end for j
                    outwrite.println();
                    outwrite.println(">");
                }//end for i
                outwrite.println("</affyfiles>");
            }
            outwrite.println("<rotmtx>");
            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    outwrite.print(in.rotmtx[i][j] + ";");
                }//end for j
                outwrite.println();
            }//end for i
            outwrite.println("</rotmtx>");
            //first write the sequences to file
            outwrite.println("<seq>");
            for (i = 0; i < seqs; i++) {
                outwrite.println(">" + namesarr[i]);
                outwrite.println(in.inaln[i].seq);
            }
            outwrite.println("</seq>");
            //write the sequence weights
            if (in.weights != null) {
                outwrite.println("<weight>");
                int weightsnum = java.lang.reflect.Array.getLength(in.weights);
                for (i = 0; i < weightsnum; i++) {
                    outwrite.println(">" + namesarr[i]);
                    outwrite.println(in.weights[i]);
                }//end for i
                outwrite.println("</weight>");
            }
            //write the sequence groups
            if ((in.seqgroupsvec!=null) && (in.seqgroupsvec.size() > 0)) {
                outwrite.println("<seqgroups>");
                seqgroup mygroup;
                for (i = 0; i < in.seqgroupsvec.size(); i++) {
                    mygroup = (seqgroup) in.seqgroupsvec.elementAt(i);
                    if (java.lang.reflect.Array.getLength(mygroup.sequences) == 0) {
                        System.err.println("WARNING: seqgroup " + mygroup.name + " has zero elements; skipping save");
                        continue;
                    }
                    outwrite.println("name=" + mygroup.name);
                    outwrite.println("type=" + mygroup.type);
                    outwrite.println("size=" + mygroup.size);
                    if (mygroup.hide == true) {
                        outwrite.println("hide=1");
                    } else {
                        outwrite.println("hide=0");
                    }
                    outwrite.println("color=" + mygroup.color.getRed() + ";" + mygroup.color.getGreen() + ";" + mygroup.color.getBlue()+ ";" + mygroup.color.getAlpha());
                    outwrite.print("numbers=");
                    for (j = java.lang.reflect.Array.getLength(mygroup.sequences) - 1; j >= 0; j--) {
                        outwrite.print(mygroup.sequences[j] + ";");
                    }//end for j
                    outwrite.println();
                }//end for i
                outwrite.println("</seqgroups>");
            }
            //next write the sequence positions
            outwrite.println("<pos>");
            for (i = 0; i < seqs; i++) {
                outwrite.println(i + " " + in.posarr[i][0] + " " + in.posarr[i][1] + " " + in.posarr[i][2]);
            }//end for i
            outwrite.println("</pos>");
            if (in.blasthits != null) {
                //next write the blast hsp results
                outwrite.println("<hsp>");
                int hspnum = java.lang.reflect.Array.getLength(in.blasthits);
                int tmpsize = 0;
                for (i = 0; i < hspnum; i++) {
                    outwrite.print(in.blasthits[i].query + " " + in.blasthits[i].hit + ":");
                    tmpsize = java.lang.reflect.Array.getLength(in.blasthits[i].val);
                    for (j = 0; j < tmpsize; j++) {
                        outwrite.print(in.blasthits[i].val[j] + " ");
                    }//end for j
                    outwrite.println();
                }//end for i
                outwrite.println("</hsp>");
            } else if (in.attvals != null) {
                //write the attvals instead of the hsp's
                outwrite.println("<att>");
                int elements = java.lang.reflect.Array.getLength(in.attvals);
                for (i = 0; i < elements; i++) {
                    outwrite.println(in.attvals[i].query + " " + in.attvals[i].hit + " " + in.attvals[i].att);
                }
                outwrite.println("</att>");
            } else {
                System.err.println("unable to print attraction values or list of HSP's to file");
                javax.swing.JOptionPane.showMessageDialog(new javax.swing.JFrame(), "unable to write attraction values or list of HSP's to file; no data");
            }
            outwrite.close();
        } catch (IOException e) {
            System.err.println("IOError writing to " + in.file.getName());
            if(nographics==false){
                javax.swing.JOptionPane.showMessageDialog(new javax.swing.JFrame(), "IOERROR writing to '" + in.file.getName() + "'");
            }
        }
    }//end saverun

    //--------------------------------------------------------------------------
    static void saveattvalstofile (File savefile, minattvals[] attvals, boolean nographics){
        try {
            PrintWriter outwrite = new PrintWriter(new BufferedWriter(new FileWriter(savefile)));
            minattvals myatt;
            for(int i=java.lang.reflect.Array.getLength(attvals);--i>=0;){
                myatt=attvals[i];
                outwrite.println(myatt.query+" "+myatt.hit+" "+myatt.att);
            }//end for i
            outwrite.close();
        } catch (IOException e) {
            System.err.println("IOError writing to " + savefile.getName());
            if(nographics==false){
                javax.swing.JOptionPane.showMessageDialog(new javax.swing.JFrame(), "IOERROR writing to '" + savefile.getName() + "'");
            }
        }
    }

    //--------------------------------------------------------------------------
    static saverunobject loadrun(File infile) {
        //load stuff from a savefile
        //!!!NOTE!!!: this was edited so as to combine hsp's with the same query-hit combination irrespective of which sequence is the query and which the hit
        //this is a valid approach in this case, as I later on anyways symmetrize the sequence interactions
        saverunobject myrun = new saverunobject();
        myrun.file = null;//if myrun has a filename all was read ok
        Vector nullvector = new Vector();
        try {
            BufferedReader inread = new BufferedReader(new FileReader(infile));
            System.out.println("loading data from '" + infile.getAbsolutePath() + "'");
            String inline;
            int seqs = -1;
            while ((inline = inread.readLine()) != null) {
                inline = inline.trim();
                if (inline.length() == 0) {
                    continue;
                } else if (inline.startsWith("#")) {//skip comment lines
                    continue;
                }
                if (inline.startsWith("sequences=")) {
                    try {
                        seqs = Integer.parseInt(inline.substring(10));
                        System.out.println("sequences=" + seqs);
                    } catch (NumberFormatException ne) {
                        System.err.println("Error parsing int from '" + inline.substring(10) + "'");
                    }
                    myrun.inaln = new AminoAcidSequence[seqs];
                    myrun.posarr = new float[seqs][3];
                    myrun.blasthits = new minhsp[0];
                    continue;
                } else if (seqs != -1) {
                    if (inline.equalsIgnoreCase("<param>")) {
                        String[] tmparr;
                        String tmpstr;
                        while (((inline = inread.readLine()) != null) && (inline.equalsIgnoreCase("</param>") == false)) {
                            //skip empty lines
                            if (inline.length() == 0) {
                                continue;
                            }
                            tmparr = inline.split("=");
                            if (java.lang.reflect.Array.getLength(tmparr) != 2) {
                                System.err.println("ERROR reading savefile for line '" + inline + "'");
                                return myrun;
                            }
                            if (tmparr[0].equalsIgnoreCase("pval")) {
                                try {
                                    myrun.pval = Double.parseDouble(tmparr[1]);
                                } catch (NumberFormatException ne) {
                                    System.err.println("ERROR parsing double from " + tmparr[1]);
                                    return myrun;
                                }
                            } else if (tmparr[0].equalsIgnoreCase("usescval")) {
                                if (tmparr[1].startsWith("t") || tmparr[1].startsWith("T")) {
                                    myrun.usescval = true;
                                } else {
                                    myrun.usescval = false;
                                }
                            } else if (tmparr[0].equalsIgnoreCase("complexatt")) {
                                if (tmparr[1].startsWith("t") || tmparr[1].startsWith("T")) {
                                    myrun.complexatt = true;
                                } else {
                                    myrun.complexatt = false;
                                }
                            } else if (tmparr[0].equalsIgnoreCase("cooling")) {
                                try {
                                    myrun.cooling = Double.parseDouble(tmparr[1]);
                                } catch (NumberFormatException ne) {
                                    System.err.println("ERROR parsing double from " + tmparr[1]);
                                    return myrun;
                                }
                            } else if (tmparr[0].equalsIgnoreCase("currcool")) {
                                try {
                                    myrun.currcool = Double.parseDouble(tmparr[1]);
                                } catch (NumberFormatException ne) {
                                    System.err.println("ERROR parsing double from " + tmparr[1]);
                                    return myrun;
                                }
                            } else if (tmparr[0].equalsIgnoreCase("maxmove")) {
                                try {
                                    myrun.maxmove = Float.parseFloat(tmparr[1]);
                                } catch (NumberFormatException ne) {
                                    System.err.println("ERROR parsing float from " + tmparr[1]);
                                    return myrun;
                                }
                            } else if (tmparr[0].equalsIgnoreCase("cluster2d")) {
                                if (tmparr[1].equalsIgnoreCase("true") || tmparr[1].startsWith("t") || tmparr[1].startsWith("T")) {
                                    myrun.cluster2d = true;
                                } else {
                                    myrun.cluster2d = false;
                                }
                            } else if (tmparr[0].equalsIgnoreCase("usefoldchange")) {
                                if (tmparr[1].equalsIgnoreCase("true") || tmparr[1].startsWith("t") || tmparr[1].startsWith("T")) {
                                    myrun.usefoldchange = true;
                                } else {
                                    myrun.usefoldchange = false;
                                }
                            } else if (tmparr[0].equalsIgnoreCase("avgfoldchange")) {
                                if (tmparr[1].equalsIgnoreCase("true") || tmparr[1].startsWith("t") || tmparr[1].startsWith("T")) {
                                    myrun.avgfoldchange = true;
                                } else {
                                    myrun.avgfoldchange = false;
                                }
                            } else if (tmparr[0].equalsIgnoreCase("showinfo")) {
                                if (tmparr[1].equalsIgnoreCase("true") || tmparr[1].startsWith("t") || tmparr[1].startsWith("T")) {
                                    myrun.showinfo = true;
                                } else {
                                    myrun.showinfo = false;
                                }
                            } else if (tmparr[0].equalsIgnoreCase("attfactor")) {
                                try {
                                    myrun.attfactor = Float.parseFloat(tmparr[1]);
                                } catch (NumberFormatException ne) {
                                    System.err.println("ERROR parsing float from " + tmparr[1]);
                                    return myrun;
                                }
                            } else if (tmparr[0].equalsIgnoreCase("repfactor")) {
                                try {
                                    myrun.repfactor = Float.parseFloat(tmparr[1]);
                                } catch (NumberFormatException ne) {
                                    System.err.println("ERROR parsing float from " + tmparr[1]);
                                    return myrun;
                                }
                            } else if (tmparr[0].equalsIgnoreCase("attvalpow")) {
                                try {
                                    myrun.attvalpow = (int) Float.parseFloat(tmparr[1]);
                                } catch (NumberFormatException ne) {
                                    System.err.println("ERROR parsing float from " + tmparr[1]);
                                    return myrun;
                                }
                            } else if (tmparr[0].equalsIgnoreCase("repvalpow")) {
                                try {
                                    myrun.repvalpow = (int) Float.parseFloat(tmparr[1]);
                                } catch (NumberFormatException ne) {
                                    System.err.println("ERROR parsing float from " + tmparr[1]);
                                    return myrun;
                                }
                            } else if (tmparr[0].equalsIgnoreCase("dampening")) {
                                try {
                                    myrun.dampening = Float.parseFloat(tmparr[1]);
                                } catch (NumberFormatException ne) {
                                    System.err.println("ERROR parsing float from " + tmparr[1]);
                                    return myrun;
                                }
                            } else if (tmparr[0].equalsIgnoreCase("zoom")) {
                                try {
                                    myrun.zoom = Float.parseFloat(tmparr[1]);
                                } catch (NumberFormatException ne) {
                                    System.err.println("ERROR parsing float from " + tmparr[1]);
                                    return myrun;
                                }
                            } else if (tmparr[0].equalsIgnoreCase("minattract")) {
                                try {
                                    myrun.minattract = Double.parseDouble(tmparr[1]);
                                } catch (NumberFormatException ne) {
                                    System.err.println("ERROR parsing double from " + tmparr[1]);
                                    return myrun;
                                }
                            } else if (tmparr[0].equalsIgnoreCase("dotsize")) {
                                try {
                                    myrun.dotsize = (int) Float.parseFloat(tmparr[1]);
                                } catch (NumberFormatException ne) {
                                    System.err.println("ERROR parsing float from " + tmparr[1]);
                                    return myrun;
                                }
                            } else if (tmparr[0].equalsIgnoreCase("ovalsize")) {
                                try {
                                    myrun.ovalsize = (int) Float.parseFloat(tmparr[1]);
                                } catch (NumberFormatException ne) {
                                    System.err.println("ERROR parsing float from " + tmparr[1]);
                                    return myrun;
                                }
                            } else if (tmparr[0].equalsIgnoreCase("groupsize")) {
                                try {
                                    myrun.groupsize = (int) Float.parseFloat(tmparr[1]);
                                } catch (NumberFormatException ne) {
                                    System.err.println("ERROR parsing float from " + tmparr[1]);
                                    return myrun;
                                }
                            } else if (tmparr[0].equalsIgnoreCase("blastpath")) {
                                myrun.blastpath = tmparr[1];
                            } else if (tmparr[0].equalsIgnoreCase("formatdbpath")) {
                                myrun.formatdbpath = tmparr[1];
                            } else if (tmparr[0].equalsIgnoreCase("functionmapfile")) {
                                //old version: only one file should be present in these save-files
                                java.io.File tmpfile = new java.io.File(tmparr[1]);
                                if (tmpfile.canRead()) {
                                    myrun.mapfiles.add(tmpfile);
                                    myrun.lookupfiles.add(null);
                                } else {
                                    System.err.println("Warning: cannot read from '" + tmparr[1] + "' setting mapfile to null");
                                }
                            } else if (tmparr[0].equalsIgnoreCase("functionlookupfile")) {
                                java.io.File tmpfile = new java.io.File(tmparr[1]);
                                if (tmpfile.canRead()) {
                                    myrun.lookupfiles.add(tmpfile);
                                } else {
                                    System.err.println("Warning: cannot read from '" + tmparr[1] + "' setting lookupfile to null");
                                }
                            } else if (tmparr[0].equalsIgnoreCase("colorarr")) {
                                tmparr = tmparr[1].split("\\):");
                                java.awt.Color[] colorarr = new java.awt.Color[java.lang.reflect.Array.getLength(tmparr)];
                                try {
                                    String[] tmptmp;
                                    int red, green, blue;
                                    for (int i = 0; i < java.lang.reflect.Array.getLength(tmparr); i++) {
                                        tmptmp = tmparr[i].substring(1).split(";");//remove the "(" and split the numbers
                                        red = Integer.parseInt(tmptmp[0]);
                                        green = Integer.parseInt(tmptmp[1]);
                                        blue = Integer.parseInt(tmptmp[2]);
                                        colorarr[i] = new java.awt.Color(red, green, blue);
                                    }//end for i
                                    myrun.colorarr = colorarr;
                                } catch (NumberFormatException ne) {
                                    System.err.println("Warning: unable to parse color array from '" + inline + "'");
                                }
                            } else if (tmparr[0].equalsIgnoreCase("colorcutoffs")) {
                                tmparr = tmparr[1].split(";");
                                float[] cutoffs = new float[java.lang.reflect.Array.getLength(tmparr)];
                                try {
                                    for (int i = 0; i < java.lang.reflect.Array.getLength(tmparr); i++) {
                                        cutoffs[i] = Float.parseFloat(tmparr[i]);
                                    }//end for i
                                    myrun.colorcutoffs = cutoffs;
                                } catch (NumberFormatException ne) {
                                    System.err.println("Warning: unable to parse float[] from '" + inline + "'");
                                }
                            }else if (tmparr[0].equalsIgnoreCase("namesdmp_file")) {
                                myrun.namesdmp_file = tmparr[1];
                            } else if (tmparr[0].equalsIgnoreCase("nodesdmp_file")) {
                                myrun.nodesdmp_file = tmparr[1];
                            } else if (tmparr[0].equalsIgnoreCase("rounds_done")) {
                                try {
                                    myrun.rounds = Integer.parseInt(tmparr[1]);
                                } catch (NumberFormatException ne) {
                                    System.err.println("ERROR parsing int from " + tmparr[1]);
                                    return myrun;
                                }
                            }
                        }//end while </param>
                    } else if (inline.equalsIgnoreCase("<function>")) {
                        File mapfile = null;
                        File lookupfile = null;
                        String[] tmp;
                        while (((inline = inread.readLine()) != null) && (inline.equalsIgnoreCase("</function>") == false)) {
                            //skip empty lines
                            if (inline.length() == 0) {
                                continue;
                            }
                            tmp = inline.split("';'");
                            if (java.lang.reflect.Array.getLength(tmp) != 2) {
                                System.err.println("ERROR readung mapping data on line '" + inline + "'");
                            } else {
                                java.io.File tmpfile = new java.io.File(tmp[0]);
                                if (tmpfile.canRead()) {
                                    myrun.mapfiles.add(tmpfile);
                                } else {
                                    System.err.println("WARNING: cannot read from '" + tmp[0] + "' skipping entry");
                                }
                                if (tmp[1].equalsIgnoreCase("NONE")) {
                                    myrun.lookupfiles.add(null);
                                } else {
                                    tmpfile = new java.io.File(tmp[1]);
                                    if (tmpfile.canRead()) {
                                        myrun.lookupfiles.add(tmpfile);
                                    } else {
                                        System.err.println("WARNING: cannot read from '" + tmp[1] + "' skipping entry");
                                    }
                                }
                            }
                        }//end while reading
                    } else if (inline.equalsIgnoreCase("<affyfiles>")) {
                        myrun.affyfiles = new Vector();
                        replicates rep = null;
                        String[] tmparr;
                        while (((inline = inread.readLine()) != null) && (inline.equalsIgnoreCase("</affyfiles>") == false)) {
                            //skip empty lines
                            if (inline.length() == 0) {
                                continue;
                            }
                            inline = inline.trim();
                            if (inline.equalsIgnoreCase("<")) {
                                rep = new replicates();
                            } else if (inline.equalsIgnoreCase(">")) {
                                myrun.affyfiles.add(rep);
                            } else {
                                tmparr = inline.split("=");
                                if (tmparr[0].equalsIgnoreCase("replicates")) {
                                    try {
                                        rep.replicates = Integer.parseInt(tmparr[1]);
                                    } catch (NumberFormatException ne) {
                                        System.err.println("ERROR: unable to parse int from '" + tmparr[1] + "'");
                                        return myrun;
                                    }
                                } else if (tmparr[0].equalsIgnoreCase("wtreplicates")) {
                                    try {
                                        rep.wtreplicates = Integer.parseInt(tmparr[1]);
                                    } catch (NumberFormatException ne) {
                                        System.err.println("ERROR: unable to parse int from '" + tmparr[1] + "'");
                                        return myrun;
                                    }
                                } else if (tmparr[0].equalsIgnoreCase("name")) {
                                    rep.name = tmparr[1];
                                } else if (tmparr[0].equalsIgnoreCase("wtname")) {
                                    rep.wtname = tmparr[1];
                                } else if (tmparr[0].equalsIgnoreCase("abbreviation")) {
                                    rep.abbreviation = tmparr[1];
                                } else if (tmparr[0].equalsIgnoreCase("replicate")) {
                                    tmparr = tmparr[1].split("';'");
                                    int repnum = java.lang.reflect.Array.getLength(tmparr);
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
                                    int repnum = java.lang.reflect.Array.getLength(tmparr);
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
                    } else if (inline.equalsIgnoreCase("<rotmtx>")) {
                        String tmpstr = "";
                        while (((inline = inread.readLine()) != null) && (inline.equalsIgnoreCase("</rotmtx>") == false)) {
                            //skip empty lines
                            if (inline.length() == 0) {
                                continue;
                            }
                            tmpstr += inline;
                        }
                        String[] tmparr = tmpstr.split(";", 0);
                        if (java.lang.reflect.Array.getLength(tmparr) != 9) {
                            System.err.println("ERROR reading rotation matrix values from '" + tmpstr + "'");
                        } else {
                            double[][] retmtx = new double[3][3];
                            boolean error = false;
                            for (int i = 0; i < 3; i++) {
                                for (int j = 0; j < 3; j++) {
                                    try {
                                        retmtx[i][j] = Double.parseDouble(tmparr[(i * 3) + j]);
                                    } catch (NumberFormatException ne) {
                                        System.err.println("unable to parse double from '" + tmparr[(i * 3) + j] + "'");
                                        i = 4;
                                        error = true;
                                    }
                                }
                            }
                            if (error == false) {
                                myrun.rotmtx = retmtx;
                            }
                        }
                    } else if (inline.equalsIgnoreCase("<seq>")) {
                        int counter = 0;
                        String currname = "";
                        String currseq = "";
                        while (((inline = inread.readLine()) != null) && (inline.equalsIgnoreCase("</seq>") == false)) {
                            //skip empty lines
                            if (inline.length() == 0) {
                                continue;
                            }
                            //read the sequence data
                            if (inline.startsWith(">")) {
                                if (currname.length() > 0) {
                                    myrun.inaln[counter] = new AminoAcidSequence();
                                    myrun.inaln[counter].name = currname;
                                    myrun.inaln[counter].seq = currseq;
                                    counter++;
                                }
                                currname = inline.substring(1).trim();
                                currseq = "";
                            } else {
                                currseq += inline;
                            }
                        }//end seq part
                        myrun.inaln[counter] = new AminoAcidSequence();
                        myrun.inaln[counter].name = currname;
                        myrun.inaln[counter].seq = currseq;
                        counter++;
                        if (counter != seqs) {
                            System.err.println("ERROR, not found the number of specified sequences, expected:" + seqs + " found:" + counter);
                            return myrun;
                        }
                    } else if (inline.equalsIgnoreCase("<weight>")) {
                        int counter = 0;
                        String currname = "";
                        String currseq = "";
                        if (seqs < 1) {
                            System.err.println("ERROR; <weight>...</weight> statements have to come AFTER <seq>...</seq> statement");
                            javax.swing.JOptionPane.showMessageDialog(new javax.swing.JFrame(), "ERROR; <weight>...</weight> statements have to come AFTER <seq>...</seq> statement");
                            return myrun;
                        }
                        myrun.weights = new float[seqs];
                        //now make a Hash that links the weights values to the sequence names
                        HashMap tmphash = new HashMap();
                        //and now read the data
                        while (((inline = inread.readLine()) != null) && (inline.equalsIgnoreCase("</weight>") == false)) {
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
                                            return myrun;
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
                                    return myrun;
                                }
                            } else {//if this name is unknown
                                System.err.println("WARNING: name '" + currname + "' is not a sequence name.");
                            }
                        }
                        //and now convert the tmphash data to a weights array
                        //for that loop through the sequence names and assign the corresponding weights
                        for (int i = 0; i < seqs; i++) {
                            if (tmphash.containsKey(myrun.inaln[i].name)) {
                                myrun.weights[i] = ((Float) tmphash.get(myrun.inaln[i].name)).floatValue();
                            } else {//default weight of 1
                                myrun.weights[i] = 1;
                            }
                        }//end for i
                    } else if (inline.equalsIgnoreCase("<pos>")) {
                        int counter = 0;
                        int mypos;
                        String[] tmparr;
                        while (((inline = inread.readLine()) != null) && (inline.equalsIgnoreCase("</pos>") == false)) {
                            //skip empty lines
                            if (inline.length() == 0) {
                                continue;
                            }
                            tmparr = inline.trim().split("\\s+");
                            if (java.lang.reflect.Array.getLength(tmparr) != 4) {
                                System.err.println("ERROR reading positions from " + inline);
                                return myrun;
                            }
                            //else
                            try {
                                mypos = Integer.parseInt(tmparr[0]);
                                myrun.posarr[mypos][0] = Float.parseFloat(tmparr[1]);
                                myrun.posarr[mypos][1] = Float.parseFloat(tmparr[2]);
                                myrun.posarr[mypos][2] = Float.parseFloat(tmparr[3]);
                                counter++;
                            } catch (NumberFormatException ne) {
                                System.err.println("ERROR, unable to parse int,double,double,double from " + inline);
                                return myrun;
                            }
                        }//end while pos
                        if (counter != seqs) {
                            System.err.println("ERROR, not found the necessary number of positions");
                            return myrun;
                        }
                    } else if (inline.equalsIgnoreCase("<hsp>")) {
                        //Vector tmpvec=new Vector();
                        String[] tmparr;
                        String tmpstr;
                        minhsp currhsp;
                        HashMap hsphash = new HashMap();
                        String hspkey = "";
                        boolean addvals = false;
                        if (myrun.pval < 0) {
                            myrun.pval = 1;
                        }
                        int count = 0;
                        String lastline = "";
                        while (((inline = inread.readLine()) != null) && (inline.equalsIgnoreCase("</hsp>") == false)) {
                            //skip empty lines
                            if (inline.length() == 0) {
                                continue;
                            }
                            if (count % 1000 == 0) {
                                System.out.print(".");
                                if (count % 100000 == 0) {
                                    System.out.print(count);
                                }
                            }
                            count++;
                            tmparr = inline.split(":");
                            int myi, myj;
                            if (inline.equalsIgnoreCase("DONE")) {
                                System.err.println("manually truncated file, returning results so far");
                                return myrun;
                            }
                            if (java.lang.reflect.Array.getLength(tmparr) != 2) {
                                System.err.println("ERROR parsing HSP data from " + inline + " right after line '" + lastline + "'");
                                return myrun;
                            }
                            //else
                            try {
                                tmpstr = tmparr[1];
                                tmparr = tmparr[0].split("\\s+");
                                if (java.lang.reflect.Array.getLength(tmparr) != 2) {
                                    System.err.println("ERROR, wrong hsp line " + inline);
                                    return myrun;
                                }
                                myi = Integer.parseInt(tmparr[0]);
                                myj = Integer.parseInt(tmparr[1]);
                                if (myi == myj) {
                                    continue;
                                }
                                //if(myi>myj){
                                //    int tmp=myi;
                                //    myi=myj;
                                //    myj=tmp;
                                //}
                                tmparr = tmpstr.split("\\s+");
                                if (java.lang.reflect.Array.getLength(tmparr) > 0) {
                                    hspkey = myi + "_" + myj;
                                    if (hsphash.containsKey(hspkey)) {
                                        currhsp = (minhsp) hsphash.get(hspkey);
                                        for (int i = 0; i < java.lang.reflect.Array.getLength(tmparr); i++) {
                                            currhsp.addpval(Double.parseDouble(tmparr[i]));
                                        }//end for i
                                    } else {
                                        currhsp = new minhsp();
                                        currhsp.query = myi;
                                        currhsp.hit = myj;
                                        currhsp.val = new double[java.lang.reflect.Array.getLength(tmparr)];
                                        for (int i = 0; i < java.lang.reflect.Array.getLength(tmparr); i++) {
                                            currhsp.val[i] = Double.parseDouble(tmparr[i]);
                                        }//end for i
                                    }
                                } else {
                                    continue;
                                }
                                hsphash.put(hspkey, currhsp);
                            } catch (NumberFormatException Ne) {
                                System.err.println("ERROR, unable to parse int int: double ... from " + inline);
                                return myrun;
                            }
                            lastline = inline;
                        }//end while hsp
                        if (inline.equalsIgnoreCase("</hsp>")) {
                            myrun.blasthits = (minhsp[]) hsphash.values().toArray(new minhsp[0]);
                        } else {
                            System.err.println("ERROR reading truncated file; <hsp> tag not closed by </hsp>");
                        }
                    } else if (inline.equalsIgnoreCase("<att>")) {
                        myrun.blasthits = null;
                        Vector tmpvec = new Vector();
                        //read in the attvals
                        System.out.println("reading attraction value data");
                        String[] tmparr;
                        String tmpstr = "";
                        if (myrun.pval < 0) {
                            myrun.pval = 0;
                        }
                        int elements, firstelements = -1, counter = 0;
                        while (((inline = inread.readLine()) != null) && (inline.equalsIgnoreCase("</att>") == false)) {
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
                            elements = java.lang.reflect.Array.getLength(tmparr);
                            if (elements != 3) {
                                System.err.println("unable to parse attraction value information from '" + inline + "'");
                            }
                            minattvals curratt = new minattvals();
                            try {
                                curratt.query = Integer.parseInt(tmparr[0]);
                                curratt.hit = Integer.parseInt(tmparr[1]);
                                curratt.att = Float.parseFloat(tmparr[2]);
                            } catch (NumberFormatException ne) {
                                System.err.println("unable to parse float for " + tmpstr + " in " + inline);
                            }
                            tmpvec.addElement(curratt);
                            counter++;
                        }
                        myrun.attvals = new minattvals[tmpvec.size()];
                        tmpvec.copyInto(myrun.attvals);
                        System.out.println("done reading matrix data");
                    } else if (inline.equalsIgnoreCase("<seqgroups>")) {
                        //while I am reading the sequence groups
                        String[] tmparr;
                        String tmpstr = "";
                        seqgroup currgroup = null;
                        int tmpval;
                        while (((inline = inread.readLine()) != null) && (inline.equalsIgnoreCase("</seqgroups>") == false)) {
                            //skip empty lines
                            if (inline.length() == 0) {
                                continue;
                            }
                            tmparr = inline.split("=");
                            if (java.lang.reflect.Array.getLength(tmparr) != 2) {
                                System.err.println("ERROR reading from savefile on line '" + inline + "'");
                                return myrun;
                            }
                            if (tmparr[0].equalsIgnoreCase("name")) {
                                if (currgroup != null) {
                                    myrun.seqgroupsvec.addElement(currgroup);
                                }
                                currgroup = new seqgroup();
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
                                    if(java.lang.reflect.Array.getLength(tmparr)<4){
                                        currgroup.color = new java.awt.Color(red, green, blue);
                                    }else{
                                        int alpha=Integer.parseInt(tmparr[3]);
                                        currgroup.color = new java.awt.Color(red, green, blue,alpha);
                                    }
                                } catch (NumberFormatException ne) {
                                    System.err.println("ERROR parsing numbers from '" + inline + "'");
                                    return myrun;
                                }
                            } else if (tmparr[0].equalsIgnoreCase("numbers")) {
                                tmparr = tmparr[1].split(";");
                                int num = java.lang.reflect.Array.getLength(tmparr);
                                int[] retarr = new int[num];
                                currgroup.polygon = makepolygons.get(currgroup.type, currgroup.size);
                                try {
                                    for (int i = 0; i < num; i++) {
                                        retarr[i] = Integer.parseInt(tmparr[i]);
                                    }//end for i
                                } catch (NumberFormatException ne) {
                                    System.err.println("ERROR parsing numbers from '" + inline + "'");
                                    return myrun;
                                }
                                currgroup.sequences = retarr;
                            } else {
                                System.err.println("Error reading savefile in line" + inline);
                                return myrun;
                            }
                        }//end while !=/seqgroups
                        if (currgroup != null) {
                            myrun.seqgroupsvec.addElement(currgroup);
                        }
                    } else if (inline.equalsIgnoreCase("<mtx>")) {
                        //load a matrix of attraction values
                        myrun.blasthits = null;
                        System.out.println("reading matrix");
                        int counter = 0;
                        int mypos;
                        String[] tmparr;
                        HashMap tmphash = new HashMap();
                        String key;
                        minattvals curratt;
                        float tmpval;
                        while (((inline = inread.readLine()) != null) && (inline.equalsIgnoreCase("</mtx>") == false)) {
                            //skip empty lines
                            if (inline.length() == 0) {
                                continue;
                            }
                            tmparr = inline.trim().split(";");
                            if (java.lang.reflect.Array.getLength(tmparr) != seqs) {
                                System.err.println("ERROR reading positions from " + inline + "; expecting " + seqs + " values");
                                return myrun;
                            }
                            try {
                                for (int i = 0; i < seqs; i++) {
                                    tmpval = Float.parseFloat(tmparr[i]);
                                    if (tmpval != 0) {
                                        if (counter < i) {
                                            key = counter + "_" + i;
                                        } else {
                                            key = i + "_" + counter;
                                        }
                                        if (tmphash.containsKey(key)) {
                                            curratt = (minattvals) tmphash.get(key);
                                            curratt.att += tmpval / 2;
                                        } else {
                                            curratt = new minattvals();
                                            curratt.query = counter;
                                            curratt.hit = i;
                                            curratt.att = tmpval / 2;
                                            tmphash.put(key, curratt);
                                        }
                                    }
                                }//end for i
                            } catch (NumberFormatException ne) {
                                System.err.println("ERROR, unable to parse float array from " + inline);
                                return myrun;
                            }
                            counter++;
                        }//end while pos
                        myrun.attvals = (minattvals[]) (tmphash.values().toArray(new minattvals[0]));
                        System.out.println("done reading matrix:" + counter);
                        if (counter != seqs) {
                            System.err.println("ERROR, not found the necessary number of positions");
                            return myrun;
                        }
                    } else {
                        System.err.println("ERROR, wrong format! unknown specs on line " + inline);
                        return myrun;
                    }
                } else {
                    System.out.println("assuming BioLayout format");
                    inread.close();
                    myrun = loadrunbiolayout(infile);
                    return myrun;
                    //System.err.println("currently not implemented");
                    //return null;
                }
            }//end while reading from file
            inread.close();
        } catch (IOException e) {
            System.err.println("IOError unable to read from " + infile.getAbsolutePath());
            return myrun;
        }
        if (myrun.posarr == null) {
            //give it random values
            int num = java.lang.reflect.Array.getLength(myrun.inaln);
            Random rand = new Random(System.currentTimeMillis());
            myrun.posarr = new float[num][3];
            for (int i = 0; i < num; i++) {
                myrun.posarr[i][0] = rand.nextFloat();
                myrun.posarr[i][1] = rand.nextFloat();
                myrun.posarr[i][2] = rand.nextFloat();
            }//end for i
        }
        myrun.file = infile;//marker for successful read
        return myrun;
    }//end loadrun

    //--------------------------------------------------------------------------
    static saverunobject loadrunbiolayout(File infile) {
        //this is supposed to read data from a biolayout input file
        saverunobject myrun = new saverunobject();
        float attval;
        String name1;
        String name2;
        ArrayList namelist = new ArrayList();
        ArrayList datlist = new ArrayList();
        try {
            BufferedReader inread = new BufferedReader(new FileReader(infile));
            String inline;
            String[] tmparr;
            nameval mynameval;//see bottom of this file for class (saves name & attval bot biolayout info)
            int vals = 0;
            int seq1num, seq2num, seqnum = 0;
            HashMap nameshash = new HashMap();
            while ((inline = inread.readLine()) != null) {
                if (inline.length() == 0) {
                    continue;
                } else if (inline.startsWith("//")) {
                    continue;
                } else {//if this is contact info
                    //split the line on spaces either in two or three
                    tmparr = inline.split("\\s+");
                    if (java.lang.reflect.Array.getLength(tmparr) == 2) {
                        if ((vals != 2) && (vals != 0)) {
                            System.out.println("ERROR found 2 elements (expected 3) on line '" + inline + "'");
                            System.err.println("Mixed 2 and 3 element entries in file " + infile.getName() + " the graph layout WILL misrepresent distances");
                            System.err.println("NEVER mix 2 and 3 element entries!");
                        }
                        vals = 2;
                        //then assign the connection a strength of 1
                        name1 = tmparr[0];
                        name2 = tmparr[1];
                        attval = 1;
                        // now see if the names are already known, else, add as new names
                        if (nameshash.containsKey(name1) == false) {
                            //System.out.println("new name "+name1);
                            namelist.add(name1);
                            seq1num = seqnum;
                            nameshash.put(name1, new Integer(seq1num));
                            seqnum++;
                        } else {
                            seq1num = ((Integer) nameshash.get(name1)).intValue();
                        }
                        if (nameshash.containsKey(name2) == false) {
                            //System.out.println("new name "+name2);
                            namelist.add(name2);
                            seq2num = seqnum;
                            nameshash.put(name2, new Integer(seq2num));
                            seqnum++;
                        } else {
                            seq2num = ((Integer) nameshash.get(name2)).intValue();
                        }
                        datlist.add(new minattvals(seq1num, seq2num, attval));
                    } else if (java.lang.reflect.Array.getLength(tmparr) == 3) {
                        if ((vals != 3) && (vals != 0)) {
                            System.out.println("ERROR found 3 elements (expected 2) on line '" + inline + "'");
                            System.err.println("Mixed 2 and 3 element entries in file " + infile.getName() + " the graph layout WILL misrepresent distances");
                            System.err.println("NEVER mix 2 and 3 element entries!");
                        }
                        vals = 3;
                        name1 = tmparr[0];
                        name2 = tmparr[1];
                        try {
                            attval = Float.parseFloat(tmparr[2]);
                        } catch (NumberFormatException ne) {
                            System.err.println("unable to parse float from '" + tmparr[2] + "' in '" + inline + "'");
                            return myrun;
                        }
                        // now see if the names are already known, else, add as new names
                        if (nameshash.containsKey(name1) == false) {
                            //System.out.println("new name "+name1);
                            namelist.add(name1);
                            seq1num = seqnum;
                            nameshash.put(name1, new Integer(seq1num));
                            seqnum++;
                        } else {
                            seq1num = ((Integer) nameshash.get(name1)).intValue();
                        }
                        if (nameshash.containsKey(name2) == false) {
                            //System.out.println("new name "+name2);
                            namelist.add(name2);
                            seq2num = seqnum;
                            nameshash.put(name2, new Integer(seq2num));
                            seqnum++;
                        } else {
                            seq2num = ((Integer) nameshash.get(name2)).intValue();
                        }
                        datlist.add(new minattvals(seq1num, seq2num, attval));
                    } else {
                        System.err.println("ERROR: line '" + inline + "' does not correspond to known format");
                        return myrun;
                    }
                }//end if not empty or // line
            }//end while reading from file
            inread.close();
        } catch (IOException ioe) {
            System.err.println("IOERROR reading from " + infile.getName());
            return myrun;
        }
        //now make an array out of the arrayList
        myrun.attvals = (minattvals[]) datlist.toArray(new minattvals[0]);
        java.util.Random rand = new java.util.Random(System.currentTimeMillis());
        myrun.inaln = new AminoAcidSequence[namelist.size()];
        myrun.posarr = new float[namelist.size()][3];
        for (int i = namelist.size(); --i >= 0;) {
            myrun.inaln[i] = new AminoAcidSequence(((String) namelist.get(i)), "");
            myrun.posarr[i][0] = rand.nextFloat();
            myrun.posarr[i][1] = rand.nextFloat();
            myrun.posarr[i][2] = rand.nextFloat();
        }//end for i
        //now normalize the attraction values and symmetrize the array of attvals
        myrun.file = infile;//marker for successful read
        return myrun;
    }//end loadrunbiolayout

    //--------------------------------------------------------------------------
    static saverunobject loadrunalternate(File infile) {
        System.out.println("reding matrixdata");
        //load data from a file with fasta format sequence input (sequence length can be null)
        //and a matrix with pairwise similarity values
        //format:
        //sequences=number_of_sequences
        //<seqs>
        //sequences (fasta format, sequences optional, names mandatory)
        //</seqs>
        //<mtx>
        //matrix_with_attracion-values
        //</mtx>
        saverunobject myrun = new saverunobject();
        myrun.file = null;//if myrun has a filename all was read ok
        try {
            BufferedReader inread = new BufferedReader(new FileReader(infile));
            String inline;
            int seqs = -1;
            while ((inline = inread.readLine()) != null) {
                inline = inline.trim();
                if (inline.length() == 0) {
                    continue;
                } else if (inline.startsWith("#")) {//skip comment lines
                    continue;
                }
                if (inline.startsWith("sequences=")) {
                    try {
                        seqs = Integer.parseInt(inline.substring(10));
                    } catch (NumberFormatException ne) {
                        System.err.println("Error parsing int from '" + inline.substring(10) + "'");
                    }
                    myrun.inaln = new AminoAcidSequence[seqs];
                    myrun.posarr = new float[seqs][3];
                    myrun.blasthits = null;
                    continue;
                } else if (seqs != -1) {
                    inline = inline.trim();
                    if (inline.equalsIgnoreCase("<seqs>")) {
                        System.out.println("reading sequences");
                        int counter = 0;
                        String currname = "";
                        String currseq = "";
                        while (((inline = inread.readLine().trim()) != null) && (inline.equalsIgnoreCase("</seqs>") == false)) {
                            //System.out.println("inline='"+inline+"'");
                            //skip empty lines
                            if (inline.length() == 0) {
                                continue;
                            }
                            //read the sequence data
                            if (inline.startsWith(">")) {
                                if (currname.length() > 0) {
                                    myrun.inaln[counter] = new AminoAcidSequence();
                                    myrun.inaln[counter].name = currname;
                                    myrun.inaln[counter].seq = currseq;
                                    counter++;
                                }
                                currname = inline.substring(1);
                                currseq = "";
                            } else {
                                currseq += inline;
                            }
                        }//end seq part
                        myrun.inaln[counter] = new AminoAcidSequence();
                        myrun.inaln[counter].name = currname;
                        myrun.inaln[counter].seq = currseq;
                        counter++;
                        if (counter != seqs) {
                            System.err.println("ERROR, not found the number of specified sequences");
                            return myrun;
                        }
                        System.out.println("done reading sequences:" + seqs);
                    } else if (inline.equalsIgnoreCase("<mtx>")) {
                        System.out.println("reading matrix");
                        int counter = 0;
                        int mypos;
                        String[] tmparr;
                        HashMap tmphash = new HashMap();
                        String key;
                        minattvals curratt;
                        float tmpval;
                        while (((inline = inread.readLine()) != null) && (inline.equalsIgnoreCase("</mtx>") == false)) {
                            //skip empty lines
                            if (inline.length() == 0) {
                                continue;
                            }
                            tmparr = inline.trim().split("\\s+");
                            if (java.lang.reflect.Array.getLength(tmparr) != seqs) {
                                System.err.println("ERROR reading positions from " + inline + "; expecting " + seqs + " values");
                                return myrun;
                            }
                            try {
                                for (int i = 0; i < seqs; i++) {
                                    tmpval = Float.parseFloat(tmparr[i]);
                                    if (tmpval != 0) {
                                        if (counter < i) {
                                            key = counter + "_" + i;
                                        } else {
                                            key = i + "_" + counter;
                                        }
                                        if (tmphash.containsKey(key)) {
                                            curratt = (minattvals) tmphash.get(key);
                                            curratt.att += tmpval / 2;
                                        } else {
                                            curratt = new minattvals();
                                            curratt.query = counter;
                                            curratt.hit = i;
                                            curratt.att = tmpval / 2;
                                            tmphash.put(key, curratt);
                                        }
                                    }
                                }//end for i
                            } catch (NumberFormatException ne) {
                                System.err.println("ERROR, unable to parse float array from " + inline);
                                return myrun;
                            }
                            counter++;
                        }//end while pos
                        myrun.attvals = (minattvals[]) (tmphash.values().toArray(new minattvals[0]));
                        System.out.println("done reading matrix:" + counter);
                        if (counter != seqs) {
                            System.err.println("ERROR, not found the necessary number of positions");
                            return myrun;
                        }
                    } else if (inline.equalsIgnoreCase("<seqgroups>")) {
                        //while I am reading the sequence groups
                        String[] tmparr;
                        String tmpstr;
                        seqgroup currgroup = null;
                        while (((inline = inread.readLine()) != null) && (inline.equalsIgnoreCase("</seqgroups>") == false)) {
                            //skip empty lines
                            if (inline.length() == 0) {
                                continue;
                            }
                            tmparr = inline.split("=");
                            if (java.lang.reflect.Array.getLength(tmparr) != 2) {
                                System.err.println("ERROR reading from savefile on line '" + inline + "'");
                                return myrun;
                            }
                            if (tmparr[0].equalsIgnoreCase("name")) {
                                if (currgroup != null) {
                                    myrun.seqgroupsvec.addElement(currgroup);
                                }
                                currgroup = new seqgroup();
                                currgroup.name = tmparr[1];
                            } else if (tmparr[0].equalsIgnoreCase("color")) {
                                tmparr = tmparr[1].split(";");
                                try {
                                    int red = Integer.parseInt(tmparr[0]);
                                    int green = Integer.parseInt(tmparr[1]);
                                    int blue = Integer.parseInt(tmparr[2]);
                                    currgroup.color = new java.awt.Color(red, green, blue);
                                } catch (NumberFormatException ne) {
                                    System.err.println("ERROR parsing numbers from '" + inline + "'");
                                    return myrun;
                                }
                            } else if (tmparr[0].equalsIgnoreCase("numbers")) {
                                tmparr = tmparr[1].split(";");
                                int num = java.lang.reflect.Array.getLength(tmparr);
                                int[] retarr = new int[num];
                                try {
                                    for (int i = 0; i < num; i++) {
                                        retarr[i] = Integer.parseInt(tmparr[i]);
                                    }//end for i
                                } catch (NumberFormatException ne) {
                                    System.err.println("ERROR parsing numbers from '" + inline + "'");
                                    return myrun;
                                }
                                currgroup.sequences = retarr;
                            } else {
                                System.err.println("Error reading savefile in line" + inline);
                                return myrun;
                            }
                        }//end while !=/seqgroups
                        if (currgroup != null) {
                            myrun.seqgroupsvec.addElement(currgroup);
                        }
                    } else {
                        System.err.println("ERROR, wrong format! unknown keyword on line " + inline);
                        return myrun;
                    }
                } else {
                    System.err.println("ERROR, wrong format! Fist line has to start with: sequences=number_of_sequences");
                    return myrun;
                }
            }//end while reading from file
            inread.close();
        } catch (IOException e) {
            System.err.println("IOError unable to read from " + infile.getName());
            return myrun;
        }
        myrun.file = infile;
        return myrun;
    }//end loadrunalternate
}//end class

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
class nameval {
    //used to save biolayout info

    public nameval() {
    }

    public nameval(String name, float val) {
        this.name = name;
        this.val = val;
    }
    String name = "";
    float val = 0;
}//end class

