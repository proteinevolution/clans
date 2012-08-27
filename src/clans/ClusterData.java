package clans;

/**
 *
 * @author tancred
 */
public class ClusterData {

    //variables initialized on creation
    minhsp[] blasthits=null;
    AminoAcidSequence[] sequences=null;
    String[] namearr=null;
    java.util.HashMap nameshash=null;
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
    String input_filename=null;

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
    double minpval=1;
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
    java.util.ArrayList<String> mapfiles=null;
    java.util.ArrayList<String> lookupfiles=null;
    java.util.Vector<String> affyfiles=null;
    boolean usefoldchange=false;
    boolean avgfoldchange=false;
    String namesdmp_file="not_spcified";
    String nodesdmp_file="not_specified";
    float zoomfactor=1;
    boolean showinfo=true;
    int[] selectednames=new int[0];
    int[] selectnames=new int[0];
    int selnamenum=0;
    float[][] lastmovearr=null;
    float[][] mymovearr=null;
    float[][] posarrtmp=null;
    int[][] drawarrtmp=null;
    java.util.ArrayList[] draworder=null;
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
    java.util.Vector seqgroupsvec=new java.util.Vector();
    java.util.ArrayList <int[][]>polygons=null;
    boolean showseqgroups=false;
    boolean changedvals=false;
    java.awt.Color[] colorarr=null;
    float[] colorcutoffs=null;
    int roundsdone=0;
    int roundslimit=-1;
    boolean moveselectedonly=false;

    public ClusterData(minhsp[] blasthits, AminoAcidSequence[] sequences, String[] namearr, java.util.HashMap nameshash, 
    		double eval, double pval, float scval, int verbose, int cpu, boolean savepos, String cmd, String blastpath,
    		boolean addblastvbparam, String formatdbpath, String[] referencedb, StringBuffer errbuff, 
    		String input_filename) {
        
    	this.sequences=ClusterMethods.remove_gaps_from_sequences(sequences);
    	
        this.movethreads=new getmovethread[cpu];
        
    	this.blasthits=blasthits;
    	this.namearr=namearr;
        this.nameshash=nameshash;
        this.eval=eval;
        this.pval=pval;
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
        this.input_filename=input_filename;
        this.seqnum=namearr.length;
    }
    
    public void load_from_file(String input_filename){
    	
        saverunobject loaded_data = CustomUtils.load_run_from_file(new java.io.File(input_filename));
        
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
       minpval=loaded_data.pval;
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
       nameshash=new java.util.HashMap((int)(number_of_sequences/0.75)+1,(float)0.75);//holds info about which name is which array number
       namearr=new String[number_of_sequences];
        for(int i=0;i<number_of_sequences;i++){
           namearr[i]=sequences[i].name.trim();
           sequences[i].name = "sequence"+i;
           nameshash.put(sequences[i].name,new Integer(i));
        }
       elements=namearr.length;
       selectednames=new int[0];
       posarr=myposarr;
       lastmovearr=new float[elements][ClusterData.dimensions];
       mymovearr=new float[elements][ClusterData.dimensions];
       posarrtmp=new float[elements][ClusterData.dimensions];
       drawarrtmp=new int[elements][ClusterData.dimensions];
       draworder=new java.util.ArrayList[0];
        
       myrotmtx=loaded_data.rotmtx;
        for (int i = 0; i < 3; i ++){
        	for (int j = 0; j < 3; j ++){
        		rotmtx[i][j]=myrotmtx[i][j];
            }
        }
        
       orgattvals=null;
        //first time I load myattvals; don't need to sync as nothing else can be using this yet
       myattvals=ClusterMethods.compute_attraction_values(blasthits,minpval,this);
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
    
    public void save_to_file(java.io.File output_file){
    	saverunobject myrun=new saverunobject();
        myrun.file=output_file;
        myrun.inaln=sequences;
        myrun.blasthits=blasthits;
        myrun.attvals=myattvals;
        myrun.posarr=myposarr;
        myrun.maxmove=maxmove;
        myrun.pval=minpval;
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
        
        CustomUtils.saverun(myrun, namearr, nographics);
        myrun=null;
    }
    
}
