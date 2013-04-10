package clans;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import java.util.Vector;

/**
 *
 * @author tancred
 */
public class ClusterData {

    //variables initialized on creation
    MinimalHsp[] blasthits=null;
    AminoAcidSequence[] sequences=null;
    String[] sequence_names=null;
    HashMap<String, Integer> nameshash=null;
    double eval=-1;
    double pval=-1;
    float scval=-1;
    int verbose=-1;
    int cpu=-1;
    boolean savepos=false;
    String cmd=null;
    String blastpath=null;
    boolean addblastvbparam=false;
    String formatdbpath=null;
    String[] referencedb=null;
    StringBuffer errbuff=null;
    
    private String input_filename=null;

    //variables initialized later on
    boolean nographics=false;
    boolean complexatt=true;
    int seqnum=0;
    float[] seqlengths=null;
    getmovethread[] movethreads=null;
    boolean usescval=false;
    minattvals[] myattvals=null;
    
    //variables I use as part of the clustering
    int rounds=0;
    float[][] myposarr=null;
    float[][] posarr=null;
    boolean cluster2d=false;
    float maxmove=0.1f;
    double pvalue_threshold=1;
    double mineval=1;
    float hidebelow=0;
    float hidebelowold=0;
    double cooling=1;
    double currcool=1;
    saverunobject saveddata=null;
    float attfactor=10f;
    float repfactor=0.5f;
    int attvalpow=1;
    int repvalpow=1;
    float dampening=0.2f;
    double minattract=1;
    float[] weights=null;
    ArrayList<File> mapfiles=null;
    ArrayList<File> lookupfiles=null;
    java.util.Vector<replicates> affyfiles=null;
    boolean usefoldchange=false;
    boolean avgfoldchange=false;
    String namesdmp_file="not_spcified";
    String nodesdmp_file="not_specified";
    float zoomfactor=1;
    boolean showinfo=false;
    int[] selectednames=new int[0];
    int[] selectnames=new int[0];
    int selnamenum=0;
    float[][] lastmovearr=null;
    float[][] mymovearr=null;
    float[][] posarrtmp=null;
    int[][] drawarrtmp=null;
    ArrayList<ArrayList<int[]>> draworder=null;
    static int dimensions=3;
    int elements=-1;
    double[][] rotmtx={{1,0,0},{0,1,0},{0,0,1}};//the performed rotations
    double[][] myrotmtx={{1,0,0},{0,1,0},{0,0,1}};//new double[3][3];//both of the above together
    minattvals[] orgattvals=null;
    boolean attvalsimple=false;
    boolean rescalepvalues=false;
    double maxvalfound=0;
    float p2attfactor=1;
    float p2attoffset=0;
    int ovalsize=4;
    int dotsize=2;
    int groupsize=4;
    java.util.Vector<seqgroup> seqgroupsvec=new java.util.Vector<seqgroup>();
    ArrayList <int[][]>polygons=null;
    boolean showseqgroups=false;
    boolean changedvals=false;
    java.awt.Color[] colorarr=null;
    float[] colorcutoffs=null;
    int roundsdone=0;
    int roundslimit=-1;
    boolean moveselectedonly=false;

    public ClusterData(MinimalHsp[] blasthits, AminoAcidSequence[] sequences, String[] namearr, HashMap<String, Integer> nameshash, 
    		double eval, double pval, float scval, int verbose, int cpu, boolean savepos, String cmd, String blastpath,
    		boolean addblastvbparam, String formatdbpath, String[] referencedb, StringBuffer errbuff, 
    		String input_filename) {
        
    	this.sequences=ClusterMethods.remove_gaps_from_sequences(sequences);
    	
        this.movethreads=new getmovethread[cpu];
        
    	this.blasthits=blasthits;
    	this.sequence_names=namearr;
        this.nameshash=nameshash;
        this.eval=eval;
        
        if (pval != -1) {
        	this.pval = pval;	
        }
        
        this.scval=scval;
        this.verbose=verbose;
        this.cpu=cpu;
        this.savepos=savepos;
        this.cmd=cmd;
        this.blastpath=blastpath;
        this.addblastvbparam=addblastvbparam;
        this.formatdbpath=formatdbpath;
        this.referencedb=referencedb;
        this.errbuff=errbuff;
        
        this.seqnum=namearr.length;
        
        this.setInputFilename(input_filename);
    }
    
    /**
     * Sets the filename as an absolute one independent of whether the input filename is relative or absolute.
     * 
     * @param new_filename a relative or absolute input filename
     */
    public void setInputFilename(String new_filename) {
    	    	
    	this.input_filename = new_filename;
        if (!new File(new_filename).isAbsolute()) {
        	// the directory java was executed from is used as basis for all relative filenames
        	this.input_filename = System.getProperty("user.dir") + System.getProperty("file.separator") + new_filename;
        }
    }
    
    /**
     * 
     * @return the absolute filename of the file associated with this clustering
     */
    public String getAbsoluteInputfileName() {
    	return this.input_filename;
    }
    
    /**
     * 
     * @return the base name of the file associated with this clustering
     */
    public String getBaseInputfileName() {
    	return new File(this.input_filename).getName();
    }
    
    
    public void resetDrawOrder() {
    	draworder = new ArrayList<ArrayList<int[]>>();
    }
    
    //--------------------------------------------------------------------------
	static saverunobject load_run_from_file(File infile) {
	    //load stuff from a savefile
	    //!!!NOTE!!!: this was edited so as to combine hsp's with the same query-hit combination irrespective of which sequence is the query and which the hit
	    //this is a valid approach in this case, as I later on anyways symmetrize the sequence interactions
	    saverunobject myrun = new saverunobject();
	    myrun.file = null;//if myrun has a filename all was read ok
	    try {
	        BufferedReader buffered_file_handle = new BufferedReader(new FileReader(infile));
	        
	        System.out.println("LOADING data from '" + infile.getAbsolutePath() + "'");
	        
	        String inline;
	        int expected_sequences = -1;
	        while ((inline = buffered_file_handle.readLine()) != null) {
	            inline = inline.trim();
	            if (inline.length() == 0) {
	                continue;
	            } else if (inline.startsWith("#")) {//skip comment lines
	                continue;
	            }
	            if (inline.startsWith("sequences=")) {
	                try {
	                    expected_sequences = Integer.parseInt(inline.substring(10));
	                    System.out.println("sequences=" + expected_sequences);
	                } catch (NumberFormatException ne) {
	                    System.err.println("Error parsing int from '" + inline.substring(10) + "'");
	                }
	                myrun.inaln = new AminoAcidSequence[expected_sequences];
	                myrun.posarr = new float[expected_sequences][3];
	                myrun.blasthits = new MinimalHsp[0];
	                continue;
	            } else if (expected_sequences != -1) {

	            	if (inline.equalsIgnoreCase("<param>")) {
	                	if (!myrun.parse_params_block(buffered_file_handle)) {
	                		System.err.println("WARNING: could not parse <param> block");
	                		return myrun;
	                	}
	                	
	                } else if (inline.equalsIgnoreCase("<function>")) {
	                	if (!myrun.parse_function_block(buffered_file_handle)) {
	                		System.err.println("WARNING: could not parse <function> block");
	                		return myrun;
	                	}
                
	                } else if (inline.equalsIgnoreCase("<affyfiles>")) {
	                    if (!myrun.parse_affyfiles_block(buffered_file_handle)) {
	                    	System.err.println("WARNING: could not parse <affyfiles> block");
	                    	return myrun;
	                    }
	                	
	                } else if (inline.equalsIgnoreCase("<rotmtx>")) {
	                    if (!myrun.parse_rotmtx_block(buffered_file_handle)) {
	                    	System.err.println("WARNING: could not parse <rotmtx> block");
	                    	return myrun;
	                    }
	                	
	                } else if (inline.equalsIgnoreCase("<seq>")) {
	                	if (!myrun.parse_seq_block(buffered_file_handle, expected_sequences)) {
	                		System.err.println("WARNING: could not parse <seq> block");
	                    	return myrun;
	                	}
	                    
	                } else if (inline.equalsIgnoreCase("<weight>")) {
	                	if (!myrun.parse_weight_block(buffered_file_handle, expected_sequences)) {
	                		System.err.println("WARNING: could not parse <weight> block");
	                    	return myrun;
	                	}
	                	
	                } else if (inline.equalsIgnoreCase("<pos>")) {
	                	if (!myrun.parse_pos_block(buffered_file_handle, expected_sequences)) {
	                		System.err.println("WARNING: could not parse <pos> block");
	                    	return myrun;
	                	}
	                    
	                } else if (inline.equalsIgnoreCase("<hsp>")) {
	                	if (!myrun.parse_hsp_block(buffered_file_handle)) {
	                		System.err.println("WARNING: could not parse <hsp> block");
	                    	return myrun;
	                	}
	                	
	                } else if (inline.equalsIgnoreCase("<att>")) {
	                  	if (!myrun.parse_att_block(buffered_file_handle)) {
	                		System.err.println("WARNING: could not parse <att> block");
	                    	return myrun;
	                	}
	                  	
	                } else if (inline.equalsIgnoreCase("<seqgroups>")) {
	                  	if (!myrun.parse_seqgroups_block(buffered_file_handle)) {
	                		System.err.println("WARNING: could not parse <seqgroups> block");
	                    	return myrun;
	                	}
	                  	
	                } else if (inline.equalsIgnoreCase("<mtx>")) {
	                	if (!myrun.parse_mtx_block(buffered_file_handle, expected_sequences)) {
	                		System.err.println("WARNING: could not parse <mtx> block");
	                    	return myrun;
	                	}
	                	
	                } else {
	                    System.err.println("ERROR, wrong format! unknown specs on line " + inline);
	                    return myrun;
	                }
	            } else {
	                System.out.println("assuming BioLayout format");
	                buffered_file_handle.close();
	                myrun = load_biolayout_file(infile);
	                return myrun;
	            }
	        }

	        buffered_file_handle.close();

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
	        }
	    }
	    
	    myrun.file = infile; //marker for successful read
	    return myrun;
	}
	
	public void load_clans_file(String input_filename){
    	
        saverunobject loaded_data = ClusterData.load_run_from_file(new java.io.File(input_filename));
        
        this.input_filename = input_filename;
        
        if(loaded_data.file == null){
            System.err.println("ERROR reading saved data from '" +input_filename + "'; aborting read");
            return;
        }
      
        System.out.println("File loaded:" + input_filename);

       sequences=ClusterMethods.remove_gaps_from_sequences(loaded_data.inaln);

       myposarr=loaded_data.posarr;
       blasthits=loaded_data.blasthits;
       usescval=loaded_data.usescval;

        if(blasthits==null){
            //first time I load myattvals; cannot be anything else; don't need to sync
           myattvals=loaded_data.attvals;
        }
        
       complexatt=loaded_data.complexatt;
       maxmove=loaded_data.maxmove;
       pvalue_threshold=loaded_data.pval;
       cooling=loaded_data.cooling;
       currcool=loaded_data.currcool;
       attfactor=loaded_data.attfactor;
       repfactor=loaded_data.repfactor;
       attvalpow=loaded_data.attvalpow;
       repvalpow=loaded_data.repvalpow;
       dampening=loaded_data.dampening;
       minattract=loaded_data.minattract;
       weights=loaded_data.weights;
       mapfiles=loaded_data.mapfiles;
       lookupfiles=loaded_data.lookupfiles;
       affyfiles=loaded_data.affyfiles;
       usefoldchange=loaded_data.usefoldchange;
       avgfoldchange=loaded_data.avgfoldchange;
       namesdmp_file=loaded_data.namesdmp_file;
       nodesdmp_file=loaded_data.nodesdmp_file;

        //be careful not to overwrite any blastpath and formatdbpath setting passed via command line
        if(blastpath.equals("blastall -p blastp")){//if it was not changed via command line
           blastpath=loaded_data.blastpath;
        }
        if(formatdbpath.equals("formatdb -p T")){//if it was not changed via command line
           formatdbpath=loaded_data.formatdbpath;
        }
        
       zoomfactor=loaded_data.zoom;
       cluster2d=loaded_data.cluster2d;
       showinfo=loaded_data.showinfo;
       int number_of_sequences =sequences.length;
       nameshash=new HashMap<String, Integer>((int)(number_of_sequences/0.75)+1,(float)0.75);//holds info about which name is which array number
       sequence_names=new String[number_of_sequences];
        for(int i=0;i<number_of_sequences;i++){
           sequence_names[i]=sequences[i].name.trim();
           sequences[i].name = "sequence"+i;
           nameshash.put(sequences[i].name,new Integer(i));
        }
       elements=sequence_names.length;
       selectednames=new int[0];
       posarr=myposarr;
       lastmovearr=new float[elements][ClusterData.dimensions];
       mymovearr=new float[elements][ClusterData.dimensions];
       posarrtmp=new float[elements][ClusterData.dimensions];
       drawarrtmp=new int[elements][ClusterData.dimensions];
       resetDrawOrder();
        
       myrotmtx=loaded_data.rotmtx;
        for (int i = 0; i < 3; i ++){
        	for (int j = 0; j < 3; j ++){
        		rotmtx[i][j]=myrotmtx[i][j];
            }
        }
        
       orgattvals=null;
        //first time I load myattvals; don't need to sync as nothing else can be using this yet
       
       compute_attraction_values();
       
       dotsize=loaded_data.dotsize;
       ovalsize=loaded_data.ovalsize;
       groupsize=loaded_data.groupsize;
       polygons=makepolygons.get(groupsize);
       seqgroupsvec=loaded_data.seqgroupsvec;
        if(seqgroupsvec.size()>0){
           showseqgroups=true;
        }
       changedvals=true;

        if(loaded_data.colorarr!=null){
            System.out.println("setting colorarr");
           colorarr=loaded_data.colorarr;
        }
        if(loaded_data.colorcutoffs!=null){
            System.out.println("setting colorcutoffs");
           colorcutoffs=loaded_data.colorcutoffs;
        }
        
       rounds = loaded_data.rounds;
        
       System.out.println("seqnum="+number_of_sequences);
       seqlengths=new float[number_of_sequences];
        float maxlength=0;
        for(int i=0;i<number_of_sequences;i++){
           seqlengths[i]=sequences[i].seq.length();
            if(seqlengths[i]>maxlength){
                maxlength=seqlengths[i];
            }
        }
        for(int i=0;i<number_of_sequences;i++){
           seqlengths[i]/=maxlength;
        }
    }
    
    //--------------------------------------------------------------------------
	static saverunobject load_biolayout_file(File infile) {
	    //this is supposed to read data from a biolayout input file
	    saverunobject myrun = new saverunobject();
	    float attval;
	    String name1;
	    String name2;
	    ArrayList<String> namelist = new ArrayList<String>();
	    ArrayList<minattvals> datlist = new ArrayList<minattvals>();
	    try {
	        BufferedReader inread = new BufferedReader(new FileReader(infile));
	        String inline;
	        String[] tmparr;
	        int vals = 0;
	        int seq1num, seq2num, seqnum = 0;
	        HashMap<String, Integer> nameshash = new HashMap<String, Integer>();
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
	                        inread.close();
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
	                    inread.close();
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
	static saverunobject load_matrix_file(File infile) {
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
	                        inread.close();
	                        return myrun;
	                    }
	                    System.out.println("done reading sequences:" + seqs);
	                } else if (inline.equalsIgnoreCase("<mtx>")) {
	                    System.out.println("reading matrix");
	                    int counter = 0;
	                    String[] tmparr;
	                    HashMap<String, minattvals> tmphash = new HashMap<String, minattvals>();
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
	                            inread.close();
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
	                            inread.close();
	                            return myrun;
	                        }
	                        counter++;
	                    }//end while pos
	                    myrun.attvals = (minattvals[]) (tmphash.values().toArray(new minattvals[0]));
	                    System.out.println("done reading matrix:" + counter);
	                    if (counter != seqs) {
	                        System.err.println("ERROR, not found the necessary number of positions");
	                        inread.close();
	                        return myrun;
	                    }
	                } else if (inline.equalsIgnoreCase("<seqgroups>")) {
	                    //while I am reading the sequence groups
	                    String[] tmparr;
	                    seqgroup currgroup = null;
	                    while (((inline = inread.readLine()) != null) && (inline.equalsIgnoreCase("</seqgroups>") == false)) {
	                        //skip empty lines
	                        if (inline.length() == 0) {
	                            continue;
	                        }
	                        tmparr = inline.split("=");
	                        if (java.lang.reflect.Array.getLength(tmparr) != 2) {
	                            System.err.println("ERROR reading from savefile on line '" + inline + "'");
	                            inread.close();
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
	                                inread.close();
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
	                                inread.close();
	                                return myrun;
	                            }
	                            currgroup.sequences = retarr;
	                        } else {
	                            System.err.println("Error reading savefile in line" + inline);
	                            inread.close();
	                            return myrun;
	                        }
	                    }//end while !=/seqgroups
	                    if (currgroup != null) {
	                        myrun.seqgroupsvec.addElement(currgroup);
	                    }
	                } else {
	                    System.err.println("ERROR, wrong format! unknown keyword on line " + inline);
	                    inread.close();
	                    return myrun;
	                }
	            } else {
	                System.err.println("ERROR, wrong format! Fist line has to start with: sequences=number_of_sequences");
	                inread.close();
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

	public void save_to_file(java.io.File output_file){
		saverunobject myrun=new saverunobject();
	    myrun.file=output_file;
	    myrun.inaln=sequences;
	    myrun.blasthits=blasthits;
	    myrun.attvals=myattvals;
	    myrun.posarr=myposarr;
	    myrun.maxmove=maxmove;
	    myrun.pval=pvalue_threshold;
	    myrun.usescval=usescval;
	    if(attvalsimple){
	        myrun.complexatt=false;
	    }else{
	        myrun.complexatt=true;
	    }
	    myrun.rotmtx=rotmtx;
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
	    myrun.dotsize=dotsize;
	    myrun.ovalsize=ovalsize;
	    myrun.groupsize=groupsize;
	    myrun.mapfiles=mapfiles;
	    myrun.lookupfiles=lookupfiles;
	    myrun.usefoldchange=usefoldchange;
	    myrun.avgfoldchange=avgfoldchange;
	    myrun.affyfiles=affyfiles;
	    myrun.namesdmp_file=namesdmp_file;
	    myrun.nodesdmp_file=nodesdmp_file;
	    if(cluster2d){
	        myrun.cluster2d=true;
	    }else{
	        myrun.cluster2d=false;
	    }
	    myrun.colorarr=colorarr;
	    myrun.colorcutoffs=colorcutoffs;
	    
	    myrun.rounds = rounds;
	    
	    ClusterData.saverun(myrun, sequence_names, nographics);
	    myrun=null;
	}

	public void append_groups_or_clusters_from_file(File infile) {
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
	            append_clusters_from_file(infile);
	        }//else do the usual read
	        inread = new BufferedReader(new FileReader(infile));
	        Vector<String> tmpvec_strings = new Vector<String>();
	        Vector<Integer> tmpvec_integers = new Vector<Integer>();
	        HashMap<String, Integer> mynamehash = new HashMap<String, Integer>();
	        String[] tmparr;
	        String tmpstr = "";
	        int elements = java.lang.reflect.Array.getLength(sequence_names);
	        //System.out.println("parentnames="+elements);
	        for (int i = 0; i < elements; i++) {
	            tmparr = sequence_names[i].split("\\s+");
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
	                        tmpvec_strings.addElement(tmparr[0]);
	                    }
	                }
	                maxnum = tmpvec_strings.size();
	                mynamearr = new String[maxnum];
	                tmpvec_strings.copyInto(mynamearr);
	            } else if (inline.equalsIgnoreCase("<seqgroups>")) {
	                //read the groups defined in this file; format is CLANS
	                if (maxnum == 0) {
	                    System.err.println("ERROR; missing sequence names");
	                    inread.close();
	                    return;
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
	                        tmpvec_strings.clear();
	                        for (int i = 0; i < elements; i++) {
	                            tmpstr = mynamearr[tmp[i]];
	                            //System.out.print(tmpstr+"/");
	                            if (mynamehash.containsKey(tmpstr)) {
	                                tmpvec_integers.addElement(mynamehash.get(tmpstr));
	                                //System.out.println("\t"+tmpstr+"-->"+((Integer)mynamehash.get(tmpstr)).intValue()+":"+namearr[((Integer)mynamehash.get(tmpstr)).intValue()]+";");
	                            } else {
	                                System.err.println("no correspondence found for '" + tmpstr + "' in current graph; skipping entry " + java.lang.reflect.Array.getLength(sequence_names));
	                            }
	                        }//end for i
	                        //System.out.println();
	                        elements = tmpvec_strings.size();
	                        mygroup.sequences = new int[elements];
	                        //System.out.print("out:");
	                        for (int i = 0; i < elements; i++) {
	                            mygroup.sequences[i] = ((Integer) tmpvec_integers.elementAt(i)).intValue();
	                            //System.out.print(mygroup.sequences[i]+";");
	                        }//end for i
	                        //System.out.println();
	                        if (elements > 0) {
	                            this.seqgroupsvec.addElement(mygroup);
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
	}//end loadgroups

	private void append_clusters_from_file(File infile) {
	    //read the cluster data generated by the iterative clustering approach
	    //(append groups to groupsvec)
	    //System.out.println("in loadcluster groupsvec.size=" + groupsvec.size());
	    HashMap<String, Integer> mynamehash = new HashMap<String, Integer>();
	    String[] tmparr;
	    String tmpstr = "";
	    int elements = java.lang.reflect.Array.getLength(sequence_names);
	    //System.out.println("parentnames="+elements);
	    for (int i = 0; i < elements; i++) {
	        tmparr = sequence_names[i].split("\\s+");
	        if (mynamehash.containsKey(tmparr[0])) {
	            System.err.println("ERROR: non unique identifier in parent: name '" + tmparr[0] + "' already defined as -->'" + mynamehash.get(tmparr[0]) + "'");
	            return;
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
	                        inread.close();
	                        return;
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
	                    inread.close();
	                    return;
	                }
	                System.out.println("reading clusters; offset=" + offset);
	                boolean stopclusters = false;
	                while (stopclusters ==false) {
	                    inline = inread.readLine();
	                    //System.out.println("reading clusters line '"+inline+"'");
	                    if (inline == null) {
	                        System.err.println("ERROR, reached EOF before end of <clusters>");
	                        inread.close();
	                        return;
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
	                                inread.close();
	                                return;
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
	                            this.seqgroupsvec.add(mycluster);
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
	    System.out.println("read " + this.seqgroupsvec.size() + " cluster entries");
	}//end loadclusters

	//--------------------------------------------------------------------------
	private static void saverun(saverunobject in, String[] namesarr,boolean nographics) {
	    //save tha data to the selected file
	    int i, j;
	    try {
	        PrintWriter outwrite = new PrintWriter(new BufferedWriter(new FileWriter(in.file)));

	        outwrite.println("sequences=" + in.inaln.length);
	        
	        outwrite.println("<param>");
	        
	        outwrite.println("attfactor=" + in.attfactor);
	        outwrite.println("attvalpow=" + in.attvalpow);
	        outwrite.println("avgfoldchange=" + in.avgfoldchange);
	        
	        outwrite.println("blastpath=" + in.blastpath);
	        
	        outwrite.println("cluster2d=" + in.cluster2d);
	        
	        outwrite.print("colorarr=");
	        for (i = 0; i < java.lang.reflect.Array.getLength(in.colorarr); i++) {
	            java.awt.Color tmp = in.colorarr[i];
	            outwrite.print("(" + tmp.getRed() + ";" + tmp.getGreen() + ";" + tmp.getBlue() + "):");
	        }//end for i
	        outwrite.println();
	        
	        outwrite.print("colorcutoffs=");
	        for (i = 0; i < java.lang.reflect.Array.getLength(in.colorcutoffs); i++) {
	            outwrite.print(in.colorcutoffs[i] + ";");
	        }//end for i
	        outwrite.println();

	        outwrite.println("complexatt=" + in.complexatt);
	        outwrite.println("cooling=" + in.cooling);
	        outwrite.println("currcool=" + in.currcool);
	        
	        outwrite.println("dampening=" + in.dampening);
	        outwrite.println("dotsize=" + in.dotsize);
	        
	        outwrite.println("formatdbpath=" + in.formatdbpath);
	        
	        outwrite.println("groupsize=" + in.groupsize);
	        
	        outwrite.println("maxmove=" + in.maxmove);
	        outwrite.println("minattract=" + in.minattract);
	        
	        if(in.namesdmp_file!=null && in.nodesdmp_file!=null){
	            outwrite.println("namesdmp_file="+in.namesdmp_file);
	            outwrite.println("nodesdmp_file="+in.nodesdmp_file);
	        }

	        outwrite.println("ovalsize=" + in.ovalsize);
	        
	        outwrite.println("pval=" + in.pval);
	        	        
	        outwrite.println("repfactor=" + in.repfactor);
	        outwrite.println("repvalpow=" + in.repvalpow);
	        outwrite.println("rounds_done=" + in.rounds);

	        outwrite.println("showinfo=" + in.showinfo);

	        outwrite.println("usefoldchange=" + in.usefoldchange);
	        outwrite.println("usescval=" + in.usescval);

	        outwrite.println("zoom=" + in.zoom);
	        
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
	        for (i = 0; i < in.inaln.length; i++) {
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
	                
	                Arrays.sort(mygroup.sequences);
	                for (j=0; j < mygroup.sequences.length; j++) {
	                    outwrite.print(mygroup.sequences[j] + ";");
	                }//end for j
	                outwrite.println();
	            }//end for i
	            outwrite.println("</seqgroups>");
	        }
	        //next write the sequence positions
	        outwrite.println("<pos>");
	        for (i = 0; i < in.inaln.length; i++) {
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
	public void save_attraction_values_to_file(File savefile){
	    try {
	        PrintWriter outwrite = new PrintWriter(new BufferedWriter(new FileWriter(savefile)));
	        minattvals myatt;
	        for(int i=java.lang.reflect.Array.getLength(myattvals);--i>=0;){
	            myatt=myattvals[i];
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
	
	public void compute_attraction_values(){
	
		if(blasthits==null){//possible (if alternate data source was loaded)
	        System.out.println("blasthits is null");
	        return;
	    }
		
	    System.out.println("blasthits is size:"+java.lang.reflect.Array.getLength(blasthits));
	    ArrayList<minattvals> tmpvec=new ArrayList<minattvals>();
	    int datanum=java.lang.reflect.Array.getLength(blasthits);
	    HashMap<String, minattvals> myhash=new HashMap<String, minattvals>(datanum);
	    float newatt;
	    String key;
	    float maxattval=0;
	    minattvals curratt=null;
	    maxvalfound=0;//init to zero, is assigned value in getattvalsimple or mult
	    //NOTE: this is not necessarily a symmetrical array. compute all values
	    //and then symmetrize computing the average values

	    if(rescalepvalues==false){
	        //make the attraction values
	        if(attvalsimple){
	            for(int i=datanum;--i>=0;){
	                if(blasthits[i].query<blasthits[i].hit){
	                    key=blasthits[i].query+"_"+blasthits[i].hit;
	                }else{
	                    key=blasthits[i].hit+"_"+blasthits[i].query;
	                }
	                if(myhash.containsKey(key)){
	                    curratt=(minattvals)myhash.get(key);
	                    if(curratt.att==-1){
	                        //in this case keep the -1
	                    }else{
	                        newatt=ClusterMethods.getattvalsimple(blasthits[i].val,elements,pvalue_threshold,this);
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
	                    if(blasthits[i].query<blasthits[i].hit){
	                        curratt.query=blasthits[i].query;
	                        curratt.hit=blasthits[i].hit;
	                    }else{
	                        curratt.hit=blasthits[i].query;
	                        curratt.query=blasthits[i].hit;
	                    }
	                    curratt.att=ClusterMethods.getattvalsimple(blasthits[i].val,elements,pvalue_threshold,this);
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
	                if(blasthits[i].query<blasthits[i].hit){
	                    key=blasthits[i].query+"_"+blasthits[i].hit;
	                }else{
	                    key=blasthits[i].hit+"_"+blasthits[i].query;
	                }
	                if(myhash.containsKey(key)){
	                    curratt=(minattvals)myhash.get(key);
	                    if(curratt.att==-1){
	                        //in this case keep the -1
	                    }else{
	                        newatt=ClusterMethods.getattvalmult(blasthits[i].val,elements,pvalue_threshold,this);
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
	                    if(blasthits[i].query<blasthits[i].hit){
	                        curratt.query=blasthits[i].query;
	                        curratt.hit=blasthits[i].hit;
	                    }else{
	                        curratt.hit=blasthits[i].query;
	                        curratt.query=blasthits[i].hit;
	                    }
	                    curratt.att=ClusterMethods.getattvalmult(blasthits[i].val,elements,pvalue_threshold,this);
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
	                if(blasthits[i].query<blasthits[i].hit){
	                    key=blasthits[i].query+"_"+blasthits[i].hit;
	                }else{
	                    key=blasthits[i].hit+"_"+blasthits[i].query;
	                }
	                if(myhash.containsKey(key)){
	                    curratt=(minattvals)myhash.get(key);
	                    if(curratt.att==-1){
	                        //in this case keep the -1
	                    }else{
	                        newatt=ClusterMethods.getattvalsimple(blasthits[i].val, elements, pvalue_threshold, this);
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
	                    if(blasthits[i].query<blasthits[i].hit){
	                        curratt.query=blasthits[i].query;
	                        curratt.hit=blasthits[i].hit;
	                    }else{
	                        curratt.hit=blasthits[i].query;
	                        curratt.query=blasthits[i].hit;
	                    }
	                    curratt.att=ClusterMethods.getattvalsimple(blasthits[i].val, elements, pvalue_threshold, this);
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
	                if(blasthits[i].query<blasthits[i].hit){
	                    key=blasthits[i].query+"_"+blasthits[i].hit;
	                }else{
	                    key=blasthits[i].hit+"_"+blasthits[i].query;
	                }
	                if(myhash.containsKey(key)){
	                    curratt=(minattvals)myhash.get(key);
	                    if(curratt.att==-1){
	                        //in this case keep the -1
	                    }else{
	                        newatt=ClusterMethods.getattvalmult(blasthits[i].val, elements, pvalue_threshold, this);
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
	                    if(blasthits[i].query<blasthits[i].hit){
	                        curratt.query=blasthits[i].query;
	                        curratt.hit=blasthits[i].hit;
	                    }else{
	                        curratt.hit=blasthits[i].query;
	                        curratt.query=blasthits[i].hit;
	                    }
	                    curratt.att = ClusterMethods.getattvalmult(blasthits[i].val, elements, pvalue_threshold, this);
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
	    myattvals = (minattvals[])tmpvec.toArray(new minattvals[0]);
	    System.out.println("attvals size=" + myattvals.length);
	}// end getattvals

	public void initialize(){
	    // compute "attraction" values for all sequence pairs from the hsp objects
		// and initialize the positions randomly
		
	    rounds=0;
	    currcool=1;
	    
	    zoomfactor = 1;
	    this.reset_rotmtx();
	    
		mymovearr = null;
	    lastmovearr = null;       
	
	    elements = sequence_names.length;
	    
	    if(elements == 0){ // no elements means nothing to do
	    	myposarr = new float[0][0];
	    	posarr=myposarr;
	    	return;
	    }
	    
	    myposarr=new float[elements][dimensions];
	    posarrtmp=new float[elements][dimensions];
	    drawarrtmp=new int[elements][dimensions];
	    mymovearr=new float[elements][dimensions];
	    lastmovearr=new float[elements][dimensions];
	    for(int i = 0; i < elements; i++){
	        mymovearr[i][0]=0;
	        mymovearr[i][1]=0;
	        mymovearr[i][2]=0;
	    	
	        lastmovearr[i][0]=0;
	        lastmovearr[i][1]=0;
	        lastmovearr[i][2]=0;
	    }
	    
	    // compute the "attraction values"
	    if(blasthits != null){
	        if(myattvals == null){
	            //synchronized(myattvals){//myattvals is null here; cannot sync on it
	            compute_attraction_values();
	            //}
	        }
	    }
	    
	    for(int i = 0; i < myposarr.length; i++){
	    	for(int j = 0; j < myposarr[j].length; j++){
	    		myposarr[i][j] = ClusterMethods.rand.nextFloat() * 2 - 1;
	        }
	    }

	    posarr=myposarr;
	}

	private void reset_rotmtx() {
		for (int i=0; i<3; i++) {
			for (int j=0; j<3; j++) {
				if (i == j) {
					rotmtx[i][j] = 1; 
					myrotmtx[i][j] = 1;
				} else {
					rotmtx[i][j] = 0;
					myrotmtx[i][j] = 0;
				}
			}
		}
	}
}