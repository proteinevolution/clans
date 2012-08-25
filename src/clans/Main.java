/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package clans;
/*
 * main.java
 *
 **Copyright (C) 2004 Tancred Frickey
 *Distributed under the GNU General Public Licence
 *This program is free software; you can redistribute it and/or modify
 *it under the terms of the GNU General Public License as published by
 *the Free Software Foundation.
 *
 *This program is distributed in the hope that it will be useful,
 *but WITHOUT ANY WARRANTY; without even the implied warranty of
 *MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *GNU General Public License for more details.
 *
 */
import java.io.*;
import java.util.*;
/**
 *
 * @author tancred
 */
public class Main {

    public static void main(String[] args) {
        System.out.println("new CLANS");
        if(java.lang.reflect.Array.getLength(args)==0){
            docalc=false;
            //printinfo();
            //System.exit(1);
        }else if(java.lang.reflect.Array.getLength(args)==1){
            loadsaved=args[0];
            docalc=false;
            args=new String[0];
        }
        String[] inargs=args;
        boolean allok=true;
        //System.out.println("checking args");
        allok=checkargs(args);//look for -conf and other pre-conf settings
        if(allok==false){
            System.err.println("unable to check args.");
            return;
        }
        //System.out.println("reading conf");
        File testfile=new File(conffilename);
        if(testfile.exists()){
            allok=readconf(conffilename);
            if(allok==false){
                System.err.println("unable to read conffile "+conffilename+"; using defaults");
                //errbuff.append("unable to read conffile '"+conffilename+"'; using defaults");
            }
        }else{
            System.err.println("Warning: "+conffilename+" does not exists, will be using default and command-line options only");
        }
        //System.out.println("reading args");
        allok=readargs(inargs);
        if(allok==false){
            System.err.println("unable to read args; exiting program.");
            printinfo();
            return;
        }else{
            //System.out.println("doing check");
            allok=docheck();
        }
        if(allok){
            makeend();
            return;
        }else{
            System.err.println("Error in docheck().");
            return;
        }
    }//end main
    
    static StringBuffer errbuff=new StringBuffer();
    static String conffilename="clans.conf";//name and location of the config file
    static String infilename="stdin";//default input
    static String cmd="";//command to prepend to anything the system executes (i.e. nice -19)
    //static String blastpath="blastall -p blastp -I T -T T ";//command necessary to start blast (if path location is needed don't forget it)
    static String blastpath="blastp ";//command necessary to start blast (if path location is needed don't forget it)
    static int blastblocks=50;//the number of sequences to do at the same time in one blast
    //static String blastpath="blastp ";//command necessary to start blast (if path location is needed don't forget it)
    static boolean addblastvbparam=true;//check to see whether I have more sequences than blast would normally return hits for
    //static String formatdbpath="formatdb -p T";//command needed to execute formatdb
    static String formatdbpath="makeblastdb -dbtype prot";//command needed to execute formatdb (since blast+ 2.2.26 -dbtype is no longer an optional entry)
    //static String formatdbpath="makeblastdb -dbtype prot";//command needed to execute formatdb (or makeblastdb for blast+)
    static String[] referencedb=new String[0];//holds the databases to blast against to generate psiblast profiles
    static boolean skipcheckdone=true; //check for a DONE in tmpblasthsp and then skip all further checks (if false)
    static double eval=10;//default maximum evalue to accept for hsp
    static double pval=0.1;//default maximum pvalue to accept for hsp
    static float coverage=(float)0;//necessary minimal coverage for blast hsp
    static float scval=(float)-1;//necessary minimal blast score/collumn for blast hsp
    static float ident=(float)0;//necessary minimal identity to query for blast hsp
    static float minconf=6;//what is the minimal confidence value to take as significant
    static int verbose=1;//verbose mode
    static int cpu=1; //how many cpu's should I use
    static boolean useallrounds=false;//do you want to use all or oly the last psiblast round for computing confidences?
    static String saveblastname="tmpblasthsp.txt";//name of blast savefile
    static boolean readblast=true;//do I want to read blasts from a savefile?
    static boolean lowmem=false;//do I want to temporarily save results to hdd?
    static boolean savepos=false;//do I want to save my sequence positions in 3d?
    static boolean docalc=true;//do I want to do calculations or just open a clustertest window?
    static boolean nographics=false;//do I want to NOT start the graphical interface?
    static String savetoname=null;//define a savefile to save results to (only if used in conjunction with -dorounds and -load)
    static int dorounds=-1;//how many rounds to cluster by (only if used in conjunction with -load)
    //--------------------------------------------------------------------------
    //variables used for adding new sequences to an already present dataset
    static String olddata="";
    static String newseqs="";
    static String loadsaved=null;
    static boolean enrichseqs=false;
    static double gatherseqseval=1e-10;
    static double rmseqseval=1e-25;
    static int maxenrichseqsnum=-1;
    static int exhaustive=1;//if I add new sequences do the fast version and only look for one way blast hits!
    //--------------------------------------------------------------------------
    
    static void printinfo(){
        //print the program arguments to stdout
        System.out.println("USAGE: java -jar clans.jar [options]");
        System.out.println("If a outOfMemoryError occurs, try running it via java -Xmx###m -jar programname options");
        System.out.println("where ### is the number of megabytes of ram you are willing to allocate to this process");
        System.out.println("-------------------OPTIONS--------------------");
        System.out.println("-? print this information");
        System.out.println("-conf name of configuration file (def: clans.conf)");
        System.out.println("-infile name of input file");
        System.out.println("-load name-of-savefile");
        System.out.println("-rounds (int) (def:"+dorounds+") how many rounds to cluster for (only used in conjunction with -load)");
        System.out.println("-saveto String where to save the results to (only used in conjunction with -load and -rounds)");
        System.out.println("-loadalt name-of-alternate-format-savefile");
        System.out.println("-lowmem t/F (doesn't do much at the moment)");
        System.out.println("-cmd String to prepend to all commands (i.e nice (unix) cmd (windows)) (def: \"\")");
        System.out.println("-blastpath \"path to blast [options]\" (def: "+blastpath+")");
        System.out.println("-blastblocks \"number of sequences to search in one blast\" (def: "+blastblocks+")");
        System.out.println("-addblastvb \"check to see whether you have more than the default 250/500 sequences blast returns hits for\" (def:"+addblastvbparam+")");
        System.out.println("-formatdbpath \"path to formatdb executable\" (def: "+formatdbpath+")");
        System.out.println("-referencedb \"databases to psiblast against to generate the profile\"");
        System.out.println("-skipcheckdone \"check for a DONE in the tmpblasthsp file and skip further checks if found\" (def:"+skipcheckdone+")");
        System.out.println("-eval maximum e-value to collect hits for (def: 10)");
        System.out.println("-pval maximum p-value to collect hits for (def: 0.1)");
        System.out.println("-scval minimum score/query_length to collect hits for (def: -1);(if >=0 it disables the P-value and E-value filters and returns the scores normalized from 0-1 as attraction values!)");
        System.out.println("-verbose verbosity of the program (def:1)");
        System.out.println("-cpu number of CPU to use (def: 1)");
        System.out.println("-readblast T/f read former blast results");
        System.out.println("-savepos t/F do you want to save the positions in 3D of the sequences during clustering");
        System.out.println("-docalc T/f do calculation or just load interface");
        System.out.println("-nographics t/F only do the blast runs; no graphical interface, results are saved to tmpblasthsp.txt");
        System.out.println("-olddata name of savefile (def: \"\")");
        System.out.println("-newseqs filename with new sequences to add to cluster (need to specify olddata) (def: \"\")");
        System.out.println("-enrichseqs t/F take newseqs as such or enrich with close relatives (blast/psiblast)");
        System.out.println("-gatherseqseval gather sequences for blast hits up to this evalue (in enrichment) (def: 1e-10)");
        System.out.println("-rmseqseval remove sequences more similar than rmseqseval (in enrichment) (def: 1e-20)");
        System.out.println("-maxenrichseqsnum get at most this number of sequences for each query (in enrichment) (def: unlimited)");
        System.out.println("-exhaustive (number) 0=one way search; 1=backvalidation; 2=redo all blast runs (def: 1)");
        System.out.println("            when adding sequences, calculate the pairwise blast values by(see above)");
        System.out.println("-------------------OPTIONS--------------------");
    }//end printinfo
    
    //--------------------------------------------------------------------------
    
    static void printargs(){
        //print the program arguments to stdout
        System.out.println("------------------SETTINGS-------------------");
        System.out.println("conffile="+conffilename);
        System.out.println("infilename="+infilename);
        System.out.println("loadname="+loadsaved);
        if(loadsaved!=null && dorounds>=0){
            System.out.println("rounds="+dorounds);
            System.out.println("saveto="+savetoname);
        }
        System.out.println("cmd="+cmd);
        System.out.println("blastpath="+blastpath);
        System.out.println("blastblocks="+blastblocks);
        System.out.println("addblastvb="+addblastvbparam);
        System.out.println("formatdbpath="+formatdbpath);
        System.out.print("referencedb: ");
        for(int i=0;i<java.lang.reflect.Array.getLength(referencedb);i++){
            System.out.print(referencedb[i]+"; ");
        }//end for i
        System.out.println();
        System.out.println("skipcheckdone="+skipcheckdone);
        System.out.println("eval="+String.valueOf(eval));
        System.out.println("pval="+String.valueOf(pval));
        System.out.println("scval="+String.valueOf(scval));
        System.out.println("verbose="+String.valueOf(verbose));
        System.out.println("cpu="+String.valueOf(cpu));
        System.out.println("lowmem="+lowmem);
        System.out.println("savepos="+savepos);
        System.out.println("docalc="+docalc);
        System.out.println("nographics="+nographics);
        System.out.println("readblast="+readblast);
        System.out.println("olddata="+olddata);
        System.out.println("newseqs="+newseqs);
        System.out.println("enrichseqs="+enrichseqs);
        if(enrichseqs){
            System.out.println("gatherseqseval="+gatherseqseval);
            System.out.println("rmseqseval="+rmseqseval);
            System.out.println("maxenrichseqsnum="+maxenrichseqsnum);
        }
        System.out.println("exhaustive="+exhaustive+"; 0=one way search; 1=backvalidation; 2=redo all blast runs");
        System.out.println("------------------SETTINGS-------------------");
    }//end printargs
    
    //--------------------------------------------------------------------------
    
    static void makeend(){
        //finalize the program; print the output to file, clean up, wait for all to finish, etc.
        
    }// end makeend
    
    //--------------------------------------------------------------------------
    
    static boolean docheck(){
        //does the actual computational parts of the program
        if(verbose>0){
            printargs();
        }
        if(docalc){
            //System.out.println("in docalc");
            if(olddata.length()==0){//if no olddata was defined
                aaseq[] inaln=readaln.fastaread(infilename);
                if(verbose>3){
                    System.out.println("sequences read:");
                    for(int i=0;i<java.lang.reflect.Array.getLength(inaln);i++){
                        System.out.println(i+" "+inaln[i].name);
                        System.out.println(i+" "+inaln[i].seq);
                    }
                }//end verbose 3
                int seqnum=java.lang.reflect.Array.getLength(inaln);
                if(seqnum<=1){
                    System.err.println("One or less sequences read, nothing to do.");
                    return true;
                }else{//if I have at least two seqs
                    //now set up a vector array that will hold the blast hsp's
                    HashMap nameshash=new HashMap((int)(seqnum/0.84),(float)0.85);//holds info about which name is which array number
                    String[] namearr=new String[seqnum];
                    for(int i=0;i<seqnum;i++){
                        namearr[i]=inaln[i].name;
                        inaln[i].name=new String("sequence"+i);//assign internal names to the sequences
                        nameshash.put(inaln[i].name,new Integer(i));
                        inaln[i].seq=inaln[i].seq.toUpperCase();
                    }
                    //now see whether I am using blast or blast+
                    boolean isblastplus=true;
                    if(blastpath.contains("blastall") || blastpath.contains("blastpgp")){
                        isblastplus=false;
                    }
                    searchblastv2 mysearchblast=new searchblastv2(errbuff,addblastvbparam,isblastplus,cpu,blastblocks,cmd,blastpath,formatdbpath,eval,pval,coverage,scval,ident,verbose,saveblastname,readblast,nameshash);
                    minhsp[] blasthits=mysearchblast.gethits(inaln);
                    mysearchblast=null;
                    if(nographics==false){
                        //System.out.println("starting clustertest");
                        System.out.println("...reading data");
                        clusterdata myclusterdata=new clusterdata(blasthits,inaln,namearr,nameshash,eval,pval,scval,verbose,cpu,savepos,cmd,blastpath,addblastvbparam,formatdbpath,referencedb,errbuff,loadsaved);
                        myclusterdata.nographics=false;
                        myclusterdata.roundslimit=dorounds;//set the limit of how often to run this
                        clustermain_graphics myclusterer=new clustermain_graphics(myclusterdata);
                        myclusterer.setVisible(true);
                        //clustertest myclustertest=new clustertest(blasthits,inaln,namearr,nameshash,eval,pval,scval,verbose,cpu,savepos,cmd,blastpath,addblastvbparam,formatdbpath,referencedb,errbuff,loadsaved);
                        //myclustertest.setVisible(true);
                    }else{
                        System.out.println("DONE. To visualize results restart program with the -nographics F option.");
                    }
                }//end else one or less seqs
            }else if(newseqs.length()>0){
                java.util.Random rand=new java.util.Random(System.currentTimeMillis());
                System.out.println("reading old data");
                saverunobject readdata=customutils.loadrun(new java.io.File(olddata));
                System.out.println("reading new sequences");
                aaseq[] newaln=readaln.read(newseqs);
                if(enrichseqs){
                    if(java.lang.reflect.Array.getLength(referencedb)==0){
                        System.err.println("ERROR, no referencedb specified, unable to enrich dataset, skipping.");
                    }else{
                        System.out.println("starting sequences="+java.lang.reflect.Array.getLength(newaln));
                        enrichutils myenrich=new enrichutils();
                        newaln=myenrich.enrich(newaln,cmd,blastpath,formatdbpath,referencedb,cpu,gatherseqseval,rmseqseval,maxenrichseqsnum);
                        System.out.println("enriched sequences="+java.lang.reflect.Array.getLength(newaln));
                    }
                }
                //now copy the position of the points and assign random positions to the new sequences
                int newelements=java.lang.reflect.Array.getLength(newaln);
                int readelements=java.lang.reflect.Array.getLength(readdata.inaln);
                int allelements=newelements+readelements;
                String[] allnamearr=new String[allelements];
                aaseq[] allaln=new aaseq[allelements];
                HashMap allnameshash=new HashMap();
                float[][] allposarr=new float[allelements][3];
                for(int i=0;i<readelements;i++){
                    allnamearr[i]=readdata.inaln[i].name;
                    readdata.inaln[i].name=new String("sequence"+i);
                    allnameshash.put(readdata.inaln[i].name,new Integer(i));
                    allaln[i]=readdata.inaln[i];
                    allposarr[i][0]=readdata.posarr[i][0];
                    allposarr[i][1]=readdata.posarr[i][1];
                    allposarr[i][2]=readdata.posarr[i][2];
                }//end for i
                int[] newnumarr=new int[newelements];//will hold the numbers of the new sequences
                for(int i=0;i<newelements;i++){
                    allnamearr[readelements+i]=newaln[i].name;
                    newaln[i].name=new String("sequence"+(readelements+i));
                    allnameshash.put(newaln[i].name, new Integer(readelements+i));
                    newnumarr[i]=readelements+i;
                    allaln[readelements+i]=newaln[i];
                    allposarr[readelements+i][0]=rand.nextFloat();
                    allposarr[readelements+i][1]=rand.nextFloat();
                    allposarr[readelements+i][2]=rand.nextFloat();
                }//end for i
                searchblast mysearchblast=new searchblast(errbuff,addblastvbparam);
                double mypval=readdata.pval;
                float mymaxmove=readdata.maxmove;
                minhsp[] newblasthits=mysearchblast.gethits(readdata.inaln,readdata.blasthits,newaln,cmd,formatdbpath,blastpath,cpu,eval,pval,coverage,scval,ident,verbose,allnameshash,useallrounds,lowmem,referencedb,exhaustive,readblast,true);
                //now keep all of the old data and add the new data
                ArrayList addblasthits=new ArrayList();
                for(int i=java.lang.reflect.Array.getLength(newblasthits);--i>=0;){
                    if(newblasthits[i].query>=readelements || newblasthits[i].hit>=readelements){
                        //then I want to add this one
                        addblasthits.add(newblasthits[i]);
                    }
                }//end for i
                //now I know which of the "new" blast hits to add
                int oldnum=java.lang.reflect.Array.getLength(readdata.blasthits);
                minhsp[] blasthits=new minhsp[oldnum+addblasthits.size()];
                System.arraycopy(readdata.blasthits,0,blasthits,0,oldnum);
                for(int i=addblasthits.size();--i>=0;){
                    blasthits[oldnum+i]=(minhsp)addblasthits.get(i);
                }//end for i
                //doen adding the blast hit data
                newaln=null;
                readdata=null;
                if(nographics==false){
                    System.out.println("...reading data");
                    clusterdata myclusterdata=new clusterdata(new minhsp[0],new aaseq[0],new String[0],new HashMap(),eval,pval,scval,verbose,cpu,savepos,cmd,blastpath,addblastvbparam,formatdbpath,referencedb,errbuff,loadsaved);
                    myclusterdata.nographics=false;
                    myclusterdata.roundslimit=dorounds;//set the limit of how often to run this
                    clustermain_graphics myclusterer=new clustermain_graphics(myclusterdata);
                    myclusterer.initaddedseqs(blasthits,allaln,allnamearr,allnameshash,newnumarr,allposarr,mymaxmove,mypval,true);
                    readdata=null;
                    myclusterer.setVisible(true);
                    //clustertest myclustertest=new clustertest(new minhsp[0],new aaseq[0],new String[0],new HashMap(),eval,pval,scval,verbose,cpu,savepos,cmd,blastpath,addblastvbparam,formatdbpath,referencedb,errbuff,loadsaved);
                    //myclustertest.initaddedseqs(blasthits,allaln,allnamearr,allnameshash,newnumarr,allposarr,mymaxmove,mypval,true);
                    //readdata=null;
                    //myclustertest.setVisible(true);
                }else{
                    System.out.println("DONE. To visualize results restart program with the -nographics F option");
                }
            }else{
                System.out.println("Reading old data from "+olddata);
                saverunobject readdata=customutils.loadrun(new java.io.File(olddata));
                int seqnum=java.lang.reflect.Array.getLength(readdata.inaln);
                HashMap nameshash=new HashMap((int)(seqnum/0.74),(float)0.75);//holds info about which name is which array number
                String[] namearr=new String[seqnum];
                for(int i=0;i<seqnum;i++){
                    namearr[i]=readdata.inaln[i].name;
                    readdata.inaln[i].name=new String("sequence"+i);//assign internal names to the sequences
                    nameshash.put(readdata.inaln[i].name,new Integer(i));
                    readdata.inaln[i].seq=readdata.inaln[i].seq.toUpperCase();
                }
                if(nographics==false){
                    System.out.println("...reading data");
                    clusterdata myclusterdata=new clusterdata(new minhsp[0],new aaseq[0],new String[0],new HashMap(),eval,pval,scval,verbose,cpu,savepos,cmd,blastpath,addblastvbparam,formatdbpath,referencedb,errbuff,loadsaved);
                    myclusterdata.nographics=false;
                    myclusterdata.roundslimit=dorounds;//set the limit of how often to run this
                    clustermain_graphics myclusterer=new clustermain_graphics(myclusterdata);
                    myclusterer.initaddedseqs(readdata.blasthits,readdata.inaln,namearr,nameshash,new int[0],readdata.posarr,readdata.maxmove,readdata.pval,false);
                    readdata=null;
                    myclusterer.setVisible(true);
                    //clustertest myclustertest=new clustertest(new minhsp[0],new aaseq[0],new String[0],new HashMap(),eval,pval,scval,verbose,cpu,savepos,cmd,blastpath,addblastvbparam,formatdbpath,referencedb,errbuff,loadsaved);
                    //myclustertest.initaddedseqs(readdata.blasthits,readdata.inaln,namearr,nameshash,new int[0],readdata.posarr,readdata.maxmove,readdata.pval,false);
                    //readdata=null;
                    //myclustertest.setVisible(true);
                }else{
                    System.out.println("DONE. To visualize results restart program with the -nographics F option");
                }
            }
        }else{//if docalc=false; i.e. the -load option was set
            //now see whether I want to recluster this dataset
            if(dorounds>=0 && savetoname!=null){
                // run CLANS in command line mode. No gui will be started. Results will be saved to a file.
                System.out.println("LOADING data from '"+loadsaved+"' and running in non-graphical mode");
                clusterdata myclusterdata=new clusterdata(new minhsp[0],new aaseq[0],new String[0],new HashMap(),eval,pval,scval,verbose,cpu,savepos,cmd,blastpath,addblastvbparam,formatdbpath,referencedb,errbuff,loadsaved);
                myclusterdata.roundslimit=dorounds;//set the limit of how often to run this
                clustermain_nographics myclusterer=new clustermain_nographics(myclusterdata);
                if(myclusterer.data.loadsaved!=null){
                    System.out.println("loading data from "+myclusterer.data.loadsaved);
                    clustermethods.loaddata(myclusterer.data);
                }
//                myclusterer.initgraph();//initialize the clustering
                myclusterer.startstopthread();//start the thread
                int waittime=15000;//15 seconds
                synchronized(myclusterer){
                    while(myclusterer.mythread.stop==false){
                        System.out.println("done clustering round "+myclusterer.data.roundsdone);
                        try{
                            myclusterer.wait(waittime);
                        }catch(InterruptedException ie){
                            System.err.println("ERROR, interrupted wait in Main\n");
                            System.exit(-5);
                        }
                    }
                }
                File savefile=new File(savetoname);
                clustermethods.savetofile(savefile,myclusterer.data);
                System.out.println("done clustering, saving results to file '"+savefile.getAbsolutePath()+"'");
                //myclusterer.dispose();
            }else{//if the reclustering is NOT the case
                if(nographics==false){//just load the file as usual and display the results
                    //clustertest myclustertest=new clustertest(new minhsp[0],new aaseq[0],new String[0],new HashMap(),eval,pval,scval,verbose,cpu,savepos,cmd,blastpath,addblastvbparam,formatdbpath,referencedb,errbuff,loadsaved);
                    //myclustertest.setVisible(true);
                    System.out.println("LOADING data from '"+loadsaved+"'");
                    clusterdata myclusterdata=new clusterdata(new minhsp[0],new aaseq[0],new String[0],new HashMap(),eval,pval,scval,verbose,cpu,savepos,cmd,blastpath,addblastvbparam,formatdbpath,referencedb,errbuff,loadsaved);
                    myclusterdata.nographics=false;
                    myclusterdata.roundslimit=dorounds;//set the limit of how often to run this
                    clustermain_graphics myclusterer=new clustermain_graphics(myclusterdata);
                    myclusterer.setVisible(true);

                }else{
                    System.out.println("Nothing to do!, try starting the program with the -nographics option set to false");
                }
            }
        }
        return true;
    }//end docheck
    
    //--------------------------------------------------------------------------
    
    //--------------------------------------------------------------------------
    //-------------------------setup stuff--------------------------------------
    //--------------------------------------------------------------------------
    
    static boolean readargs(String[] args){
        //read the arguments from the command line and those that are passed from the readconf method
        int argsize=java.lang.reflect.Array.getLength(args);
        int i=0;
        while(i<argsize){
            if(args[i].equals("?")||args[i].equals("-?")){
                printinfo();
                System.exit(0);
            }
            if(args[i].equalsIgnoreCase("-conf")||args[i].equalsIgnoreCase("-c")){
                //this shopuld have been read in checkargs, so here just skip it
                i++;
                if((i)<argsize){
                    //do nothing and increase i
                }else{
                    System.err.println("Error reading -conf, missing argument.");
                    return false;
                }
                i++;
                continue;
            }
            if((args[i].equalsIgnoreCase("-infile"))||(args[i].equalsIgnoreCase("-i"))||(args[i].equalsIgnoreCase("-in"))){
                i++;
                if((i)<argsize){
                    infilename=args[i];
                }else{
                    System.err.println("Error reading -infile, missing argument.");
                    return false;
                }
                i++;
                continue;
            }//end infile
            if((args[i].equalsIgnoreCase("-load"))||(args[i].equalsIgnoreCase("-l"))){
                i++;
                if((i)<argsize){
                    loadsaved=args[i];
                    docalc=false;
                }else{
                    System.err.println("Error reading -load, missing argument.");
                    return false;
                }
                i++;
                continue;
            }//end load
            if((args[i].equalsIgnoreCase("-saveto"))){
                i++;
                if((i)<argsize){
                    savetoname=args[i];
                }else{
                    System.err.println("Error reading -saveto, missing argument.");
                    return false;
                }
                i++;
                continue;
            }//end saveto
            if((args[i].equalsIgnoreCase("-referencedb"))||(args[i].equalsIgnoreCase("-refdb"))){
                i++;
                if((i)<argsize){
                    referencedb=args[i].split("\\s+",0);
                }else{
                    System.err.println("Error reading -referencedb, missing argument.");
                    return false;
                }
                i++;
                continue;
            }//end referencedb
            if(args[i].equalsIgnoreCase("-cmd")){
                int quotesfound=0;
                cmd="";
                if((i+1)<argsize){
                    if(args[i+1].indexOf("\"")==-1){// if the next elem is not in quotes
                        i++;
                        cmd=args[i];
                        i++;
                        continue;
                    }else{
                        while (quotesfound<2){
                            i++;
                            if(i>=argsize){
                                System.err.println("Error reading -cmd, missing argument.");
                                return false;
                            }
                            String curr=args[i];
                            int qindex;
                            if((qindex=curr.indexOf("\""))>-1){
                                quotesfound+=1;
                                if(quotesfound==1){
                                    cmd+=curr.substring(qindex+1);
                                }
                                if(quotesfound==2){
                                    cmd+=" "+curr.substring(0,qindex);
                                }
                                //now check if the second quote is somewhere in command
                                if((qindex=cmd.indexOf("\""))>-1){
                                    quotesfound+=1;
                                    cmd=cmd.substring(0,cmd.indexOf("\""));
                                    continue;
                                }
                                continue;
                            }// end if index1
                            if(quotesfound==1){
                                cmd=cmd+" "+curr;
                            }
                        }// end while quotesfound
                    }// end if firstelem had quotes
                }else{
                    System.err.println("Error reading -cmd, missing argument.");
                    return false;
                }
                i++;
                continue;
            }// end if -cmd
            if((args[i].equalsIgnoreCase("-blastpath"))||(args[i].equalsIgnoreCase("-blast"))){
                int quotesfound=0;
                blastpath="";
                if(args[i+1].indexOf("\"")==-1){// if the next elem is not in quotes
                    i++;
                    blastpath=args[i];
                    i++;
                    continue;
                }else{
                    while (quotesfound<2){
                        i++;
                        if(i>=argsize){
                            System.err.println("Error reading -blastpath, missing argument.");
                            return false;
                        }
                        String curr=args[i];
                        int qindex;
                        if((qindex=curr.indexOf("\""))>-1){//if this element has quotes
                            quotesfound++;
                            if(quotesfound==1){//if this is the first quote I found
                                blastpath=curr.substring(qindex+1);//everything from the quote to the end
                            }
                            if(quotesfound==2){//if I found the seconf quote
                                blastpath+=" "+curr.substring(0,qindex);//add the rest
                            }
                            //now check if the second quote is somewhere in command
                            if((qindex=blastpath.indexOf("\""))>-1){
                                quotesfound+=1;
                                blastpath=blastpath.substring(0,blastpath.indexOf("\""));
                                continue;
                            }
                        }else{//if this element does not contain a quote
                            if(quotesfound==1){//but I already did find one before
                                blastpath+=" "+curr;
                            }
                        }//end else quote found
                    }// end while quotesfound
                }// end if firstelem had quotes
                i++;
                continue;
            }// end if -blastpath
            if((args[i].equalsIgnoreCase("-addblastvb"))){
                i++;
                if(i<argsize){
                    if(args[i].equalsIgnoreCase("FALSE")||args[i].equalsIgnoreCase("F")){
                        addblastvbparam=false;
                    }else{
                        addblastvbparam=true;//default
                    }
                }else{
                    System.err.println("Error reading -addblastvb, missing argument.");
                    return false;
                }
                i++;
                continue;
            }// end readblast
            if((args[i].equalsIgnoreCase("-formatdbpath"))||(args[i].equalsIgnoreCase("-fdb"))){
                int quotesfound=0;
                formatdbpath="";
                if(args[i+1].indexOf("\"")==-1){// if the next elem is not in quotes
                    i++;
                    formatdbpath=args[i];
                    i++;
                    continue;
                }else{
                    while (quotesfound<2){
                        i++;
                        if(i>=argsize){
                            System.err.println("Error reading -formatdbpath, missing argument.");
                            return false;
                        }
                        String curr=args[i];
                        int qindex;
                        if((qindex=curr.indexOf("\""))>-1){
                            quotesfound+=1;
                            if(quotesfound==1){
                                formatdbpath+=curr.substring(qindex+1);
                            }
                            if(quotesfound==2){
                                formatdbpath+=" "+curr.substring(0,qindex);
                            }
                            //now check if the second quote is somewhere in formatdbpath
                            if((qindex=formatdbpath.indexOf("\""))>-1){
                                quotesfound+=1;
                                formatdbpath=formatdbpath.substring(0,formatdbpath.indexOf("\""));
                                continue;
                            }
                            continue;
                        }// end if index1
                        if(quotesfound==1){
                            formatdbpath=formatdbpath+" "+curr;
                        }
                    }// end while quotesfound
                    i++;
                    continue;
                }// end if firstelem had quotes
            }// end if -formatdbpath
            if((args[i].equalsIgnoreCase("-blastblocks"))||(args[i].equalsIgnoreCase("-bb"))){
                i++;
                if((i)<argsize){
                    try{
                        blastblocks=Integer.parseInt(args[i]);
                    }catch (NumberFormatException e){
                        System.err.println("unable to parse int from '"+args[i]+"' in -blastblocks");
                        return false;
                    }
                }else{
                    System.err.println("Error reading -blastblocks, missing argument.");
                    return false;
                }
                i++;
                continue;
            }//end eval
            if((args[i].equalsIgnoreCase("-eval"))||(args[i].equalsIgnoreCase("-e"))){
                i++;
                if((i)<argsize){
                    try{
                        eval=Double.parseDouble(args[i]);
                        if(eval<pval){
                            pval=eval;
                        }
                    }catch (NumberFormatException e){
                        System.err.println("unable to parse double from '"+args[i]+"' in -eval");
                        return false;
                    }
                }else{
                    System.err.println("Error reading -eval, missing argument.");
                    return false;
                }
                i++;
                continue;
            }//end eval
            if((args[i].equalsIgnoreCase("-pval"))||(args[i].equalsIgnoreCase("-p"))){
                i++;
                if((i)<argsize){
                    try{
                        pval=Double.parseDouble(args[i]);
                    }catch (NumberFormatException e){
                        System.err.println("unable to parse double from '"+args[i]+"' in -pval");
                        return false;
                    }
                }else{
                    System.err.println("Error reading -pval, missing argument.");
                    return false;
                }
                i++;
                continue;
            }//end pval
            if((args[i].equalsIgnoreCase("-scval"))||(args[i].equalsIgnoreCase("-sc"))){
                i++;
                if((i)<argsize){
                    try{
                        scval=Float.parseFloat(args[i]);
                    }catch (NumberFormatException e){
                        System.err.println("unable to parse float from '"+args[i]+"' in -pval");
                        return false;
                    }
                }else{
                    System.err.println("Error reading -pval, missing argument.");
                    return false;
                }
                i++;
                continue;
            }//end scval
            if((args[i].equalsIgnoreCase("-verbose"))||(args[i].equalsIgnoreCase("-v"))){
                i++;
                if((i)<argsize){
                    try{
                        verbose=Integer.parseInt(args[i]);
                    }catch (NumberFormatException e){
                        System.err.println("unable to parse int from "+args[i]+" in -verbose.");
                        return false;
                    }
                }else{
                    System.err.println("Error reading -verbose, missing argument.");
                    return false;
                }
                i++;
                continue;
            }// end in -verbose||-v
            if((args[i].equalsIgnoreCase("-cpu"))){
                i++;
                if((i)<argsize){
                    try{
                        cpu=Integer.parseInt(args[i]);
                    }catch (NumberFormatException e){
                        System.err.println("unable to parse int from "+args[i]+" in -cpu.");
                        return false;
                    }
                }else{
                    System.err.println("Error reading -cpu, missing argument.");
                    return false;
                }
                i++;
                continue;
            }// end in -cpu
            if((args[i].equalsIgnoreCase("-dorounds"))){
                i++;
                if((i)<argsize){
                    try{
                        dorounds=Integer.parseInt(args[i]);
                    }catch (NumberFormatException e){
                        System.err.println("unable to parse int from "+args[i]+" in -dorounds.");
                        return false;
                    }
                }else{
                    System.err.println("Error reading -dorounds, missing argument.");
                    return false;
                }
                i++;
                continue;
            }// end in -dorounds
            if((args[i].equalsIgnoreCase("-readblast"))){
                i++;
                if(i<argsize){
                    if(args[i].equalsIgnoreCase("FALSE")||args[i].equalsIgnoreCase("F")){
                        readblast=false;
                    }else{
                        readblast=true;//default
                    }
                }else{
                    System.err.println("Error reading -readblast, missing argument.");
                    return false;
                }
                i++;
                continue;
            }// end readblast
            if((args[i].equalsIgnoreCase("-lowmem"))){
                i++;
                if(i<argsize){
                    if(args[i].equalsIgnoreCase("TRUE")||args[i].equalsIgnoreCase("T")){
                        lowmem=true;
                    }else{
                        lowmem=false;//default
                    }
                }else{
                    System.err.println("Error reading -lowmem, missing argument.");
                    return false;
                }
                i++;
                continue;
            }// end lowmem
            if((args[i].equalsIgnoreCase("-savepos"))){
                i++;
                if(i<argsize){
                    if(args[i].equalsIgnoreCase("TRUE")||args[i].equalsIgnoreCase("T")){
                        savepos=true;
                    }else{
                        savepos=false;//default
                    }
                }else{
                    System.err.println("Error reading -savepos, missing argument.");
                    return false;
                }
                i++;
                continue;
            }// end savepos
            if((args[i].equalsIgnoreCase("-docalc"))){
                i++;
                if(i<argsize){
                    if(args[i].equalsIgnoreCase("TRUE")||args[i].equalsIgnoreCase("T")){
                        docalc=true;
                    }else{
                        docalc=false;
                    }
                }else{
                    System.err.println("Error reading -docalc, missing argument.");
                    return false;
                }
                i++;
                continue;
            }// end docalc
            if((args[i].equalsIgnoreCase("-nographics"))){
                i++;
                if(i<argsize){
                    if(args[i].equalsIgnoreCase("TRUE")||args[i].equalsIgnoreCase("T")){
                        nographics=true;
                    }else{
                        nographics=false;
                    }
                }else{
                    System.err.println("Error reading -nographics, missing argument.");
                    return false;
                }
                i++;
                continue;
            }// end nographics
            if((args[i].equalsIgnoreCase("-exhaustive"))){
                i++;
                if((i)<argsize){
                    try{
                        exhaustive=Integer.parseInt(args[i]);
                    }catch (NumberFormatException e){
                        System.err.println("unable to parse int from "+args[i]+" in -exhaustive.");
                        return false;
                    }
                    if (exhaustive>2){
                        exhaustive=2;
                    }else if(exhaustive<0){
                        exhaustive=0;
                    }
                }else{
                    System.err.println("Error reading -exhaustive, missing argument.");
                    return false;
                }
                i++;
                continue;
            }// end in -cpu
            if((args[i].equalsIgnoreCase("-olddata"))){
                i++;
                if((i)<argsize){
                    olddata=args[i];
                }else{
                    System.err.println("Error reading -olddata, missing argument.");
                    return false;
                }
                i++;
                continue;
            }//end olddata
            if((args[i].equalsIgnoreCase("-newseqs"))){
                i++;
                if((i)<argsize){
                    newseqs=args[i];
                }else{
                    System.err.println("Error reading -newseqs, missing argument.");
                    return false;
                }
                i++;
                continue;
            }//end newseqs
            if((args[i].equalsIgnoreCase("-enrichseqs"))){
                i++;
                if(i<argsize){
                    if(args[i].equalsIgnoreCase("true")||args[i].equalsIgnoreCase("t")){
                        enrichseqs=true;
                    }else{
                        enrichseqs=false;
                    }
                }else{
                    System.err.println("Error reading -enrichseqs, missing argument.");
                    return false;
                }
                i++;
                continue;
            }//end enrichseqs
            if(args[i].equalsIgnoreCase("-gatherseqseval")){
                i++;
                if(i<argsize){
                    try{
                        gatherseqseval=Double.parseDouble(args[i]);
                    }catch (NumberFormatException e){
                        System.err.println("ERROR: unable to parse double from "+args[i]+" in -gatherseqseval");
                        return false;
                    }
                }else{
                    System.err.println("Error reading -gatherseqseval, missing argument.");
                    return false;
                }
                i++;
                continue;
            }//end gatherseqseval
            if(args[i].equalsIgnoreCase("-rmseqseval")){
                i++;
                if(i<argsize){
                    try{
                        rmseqseval=Double.parseDouble(args[i]);
                    }catch (NumberFormatException e){
                        System.err.println("ERROR: unable to parse double from "+args[i]+" in -rmseqseval");
                        return false;
                    }
                }else{
                    System.err.println("Error reading -rmseqseval, missing argument.");
                    return false;
                }
                i++;
                continue;
            }//end rmseqseval
            if(args[i].equalsIgnoreCase("-maxenrichseqsnum")){
                i++;
                if(i<argsize){
                    try{
                        gatherseqseval=Integer.parseInt(args[i]);
                    }catch (NumberFormatException e){
                        System.err.println("ERROR: unable to parse int from "+args[i]+" in -maxenrichseqsnum");
                        return false;
                    }
                }else{
                    System.err.println("Error reading -maxenrichseqsnum, missing argument.");
                    return false;
                }
                i++;
                continue;
            }//end maxenrichseqsnum
            //if I get here I have an unknown parameter
            System.err.println("unknown option "+args[i]);
            return false;
        }// end while i<argsize
        return true;
    }// end readargs
    
    //--------------------------------------------------------------------------
    
    static boolean readconf(String filename){
        //read the configuration file
        try{
            BufferedReader infile=new BufferedReader(new FileReader(filename));
            String inline;
            int enddata=0;
            while((inline=infile.readLine())!=null){//while I am not at EOF
                inline=inline.trim();
                if((enddata=inline.indexOf("#"))>-1){//if this is a line with a comment on it
                    if(enddata>1){//if I have some data on this line
                        inline=inline.substring(0,enddata);
                        if((readargs(inline.split("\\s",0)))==false){
                            System.err.println("Error reading on line "+inline);
                            return false;
                        }
                    }else{
                        continue;
                    }
                }else{//if this line has no comment on it
                    if(inline.length()>0){
                        if((readargs(inline.split("\\s",0)))==false){
                            System.err.println("Error reading on line "+inline);
                            return false;
                        }
                    }
                }
            }// end while
        }catch (IOException e){
            System.err.println("IOError reading from "+filename);
            return false;
        }
        return true;
    }// end readconf
    
    //--------------------------------------------------------------------------
    
    static boolean checkargs(String [] args){
        //look for any pre-conffile settings
        int argsize=java.lang.reflect.Array.getLength(args);
        for(int i=0;i<argsize;i++){
            if((args[i].equalsIgnoreCase("-conf"))||(args[i].equalsIgnoreCase("-c"))){
                if((i+1)<argsize){
                    conffilename=args[i+1];
                }else{
                    return false;
                }
            }// end in -conf || -c
            if((args[i].equalsIgnoreCase("-verbose"))||(args[i].equalsIgnoreCase("-v"))){
                if((i+1)<argsize){
                    try{
                        verbose=Integer.parseInt(args[i+1]);
                    }catch (NumberFormatException e){
                        System.err.println("unable to parse int from "+args[i+1]+" in -verbose.");
                        return false;
                    }
                }else{
                    return false;
                }
            }// end in -verbose||-v
        }// end for i
        return true;
    }//end checkargs
    
    
}
