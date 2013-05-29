/*
 * viewblasthitsutils.java
 *
 * Created on November 3, 2003, 12:25 PM
 */
package clans;
import java.util.*;
import java.io.*;
/**
 *
 * @author  tancred
 */
public class viewblasthitsutils {
    
    /** Creates a new instance of viewblasthitsutils */
    public viewblasthitsutils() {
    }
    
    hsp[] gethsps(int referenceseqnum,AminoAcidSequence[] inaln,String cmd,String formatdbpath,String blastpath,boolean addblastvbparam,String[] referencedb,double mineval,double minpval){
        //get all the blast hits to this sequence.
        Vector<hsp[]> retvec = new Vector<hsp[]>();
        String basename=String.valueOf(System.currentTimeMillis());
        String blastcommand="";
        String tmpfilestring=basename+".query";
        String tmpcheckfile=basename+".chkpnt";
        String dbstring=basename;
        System.out.println("passed formatdbpath="+formatdbpath);
        Runtime rt=Runtime.getRuntime();
        int seqnum=inaln.length;
        if(addblastvbparam==true){
            //see whether I have more sequences than would normally be returned by blast
            if(seqnum>250){//250 is the default used by blast
                if(blastpath.contains("blastall")||blastpath.contains("blastpgp")){//non-blast+ versions
                    if(blastpath.indexOf("-b ")>-1 || blastpath.indexOf("-v ")>-1){
                        //if either of these is set, then don't re-set them
                    }else{
                        //if neither -v nor -b was set in blastpath, set them to allow (on average) one hit per sequence
                        blastpath+=" -v "+seqnum+" -b "+seqnum;
                    }
                }else if(blastpath.contains("blastp") || blastpath.contains("psiblast")){//potential blast+ versions
                    if(blastpath.indexOf("-num_descriptions ")>-1 || blastpath.indexOf("-num_alignments ")>-1){
                        //if either of these is set, then don't re-set them
                    }else{
                        //if neither -v nor -b was set in blastpath, set them to allow (on average) one hit per sequence
                        blastpath+=" -num_descriptions "+seqnum+" -num_alignments "+seqnum;
                    }
                }
            }
        }
        HashMap<String, Integer> nameshash = new HashMap<String, Integer>();
        for(int i=0;i<seqnum;i++){
            nameshash.put(inaln[i].name, new Integer(i));
        }//end for i
        //format the database
        try{
            PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter(basename)));
            for(int i=0;i<seqnum;i++){
                outwrite.println(">"+inaln[i].name);
                outwrite.println(inaln[i].seq);
            }//end for i
            outwrite.close();
        }catch (IOException ie){
            System.err.println("IOERROR writing to "+basename);
            return new hsp[0];
        }
        try{
            String formatcmd="";
            if(cmd.length()>0){
                formatcmd=cmd+" ";
            }
            if(formatdbpath.contains("makeblastdb")){//blast+ version
                formatcmd+=formatdbpath+" -in "+basename;
            }else{//former blast version
                formatcmd+=formatdbpath+" -i "+basename;
            }
            StringBuffer inout=new StringBuffer();
            StringBuffer errout=new StringBuffer();
            System.out.println("trying "+formatcmd);
            BufferedReader perr;
            BufferedReader pin;
            threadstreamreader perrread;
            threadstreamreader pinread;
            Process p=rt.exec(formatcmd);
            perr=new BufferedReader(new InputStreamReader(p.getErrorStream()));//process error stream
            pin=new BufferedReader(new InputStreamReader(p.getInputStream()));//output of process
            //note: threadstreamreader closes the input bufferedreader
            //threadstreamreaders are used because otherwise deadlock can occur if the output of the process is larger than the
            //buffersize java allocates to the bufferedreaer.
            perrread=new threadstreamreader(perr,errout);
            pinread=new threadstreamreader(pin,inout);
            try{
                perrread.start();
                pinread.start();
                p.waitFor();
                if(p.exitValue()!=0){
                    System.err.println("non-perfect exit for "+formatcmd);
                    System.err.print(errout.toString());
                }
                while (pinread.done==false){
                    synchronized(inout){
                        try{
                            inout.wait(10l);
                        }catch(InterruptedException e){
                            System.err.println("interrupted sleep in formatting database");
                            e.printStackTrace();
                        }
                    }
                }
            }catch (InterruptedException e){
                System.err.println("Interrupted process "+formatcmd);
                return new hsp[0];
            }
            perrread=null;//.clear();
            pinread=null;//.clear();
            System.out.println("Done formatdb for "+basename);
        }catch (IOException ie){
            System.err.println("ERROR formatting database "+basename);
            return new hsp[0];
        }
        //run against either the referencedb's or directly against the present sequences
        if((referencedb.length>0)&&(blastpath.indexOf("blastpgp")>-1)){
            System.out.println("profile psiblast run for "+basename);
            Vector<String> tmpvec=new Vector<String>();//supposed to hold the command parameters
            //first do the psiblast runs agains all specified databases
            //first put all database names into one string
            dbstring=basename;
            for(int i=0;i<referencedb.length;i++){
                dbstring+=" "+referencedb[i];
            }//end for i
            blastcommand=blastpath;
            String[] cmdarr=blastpath.split("\\s+");
            tmpvec.addElement("-T");
            tmpvec.addElement("T");
            tmpvec.addElement("-e");
            tmpvec.addElement(String.valueOf(mineval));
            tmpvec.addElement("-i");
            tmpvec.addElement(tmpfilestring);
            tmpvec.addElement("-C");
            tmpvec.addElement(tmpcheckfile);
            tmpvec.addElement("-d");
            tmpvec.addElement(dbstring);
            for(int i=cmdarr.length;i>0;i--){
                tmpvec.add(0,cmdarr[i-1]);
            }
            if(cmd.length()>0){
                blastcommand=cmd+" "+blastcommand;
                tmpvec.add(0, cmd);
            }
            cmdarr=new String[tmpvec.size()];
            tmpvec.copyInto(cmdarr);
            StringBuffer inout=new StringBuffer();
            StringBuffer errout=new StringBuffer();
            System.out.println("for "+inaln[referenceseqnum].name+" trying command="+blastcommand);
            try{
                //write the query to a file instead of using stdin; !stdin gives problems!
                PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter(tmpfilestring)));
                outwrite.println(">"+inaln[referenceseqnum].name);
                outwrite.println(inaln[referenceseqnum].seq.replaceAll("-",""));
                outwrite.close();
                BufferedReader perr;
                BufferedReader pin;
                threadstreamreader perrread;
                threadstreamreader pinread;
                Process p=rt.exec(cmdarr);
                perr=new BufferedReader(new InputStreamReader(p.getErrorStream()));//process error stream
                pin=new BufferedReader(new InputStreamReader(p.getInputStream()));//output of process
                //note: threadstreamreader closes the input bufferedreader
                //threadstreamreaders are used because otherwise deadlock can occur if the output of the process is larger than the
                //buffersize java allocates to the bufferedreaer.
                perrread=new threadstreamreader(perr,errout);
                pinread=new threadstreamreader(pin,inout);
                try{
                    perrread.start();
                    pinread.start();
                    p.waitFor();
                    if(p.exitValue()!=0){
                        System.err.println("non-perfect exit from blast for "+inaln[referenceseqnum].name);
                    }
                    while (pinread.done==false){
                        //System.out.println("waiting for process reader to finish");
                        synchronized(inout){
                            try{
                                inout.wait(10l);
                            }catch(InterruptedException e){
                                System.err.println("interrupted sleep in blastthread");
                                e.printStackTrace();
                            }
                        }
                    }
                }catch (InterruptedException e){
                    System.err.println("Interrupted process "+blastcommand+" for "+inaln[referenceseqnum].name);
                    return new hsp[0];
                }
                perrread=null;//.clear();
                pinread=null;//.clear();
            }catch (IOException ioe){
                System.err.println("IOError in "+blastcommand+" for "+inaln[referenceseqnum].name);
                for(int i=0;i<cmdarr.length;i++){
                    System.err.println(cmdarr[i]);
                }
                System.out.println("inoutsize="+inout.length());
                return new hsp[0];
            }
            System.err.print(errout.toString());
            //now the psiblast run should have generated a checkpoint file
            //then do the last round using the psiblast checkpoint file against blastdb
            System.out.println("actual psiblast run on "+inaln[referenceseqnum].name);
            blastcommand=blastpath+" -d "+basename+" -T T -e "+mineval+" -i "+tmpfilestring+" -R "+tmpcheckfile;
        }else if((referencedb.length>0)&&(blastpath.indexOf("psiblast")>-1)){
            System.out.println("profile psiblast run for "+basename);
            Vector<String> tmpvec = new Vector<String>();// supposed to hold the command parameters
            //first do the psiblast runs agains all specified databases
            //first put all database names into one string
            dbstring=basename;
            for(int i=0;i<referencedb.length;i++){
                dbstring+=" "+referencedb[i];
            }//end for i
            blastcommand=blastpath;
            String[] cmdarr=blastpath.split("\\s+");
            tmpvec.addElement("-html");
            tmpvec.addElement("-show_gis");
            tmpvec.addElement("-evalue");
            tmpvec.addElement(String.valueOf(mineval));
            tmpvec.addElement("-query");
            tmpvec.addElement(tmpfilestring);
            tmpvec.addElement("-out_pssm");
            tmpvec.addElement(tmpcheckfile);
            tmpvec.addElement("-db");
            tmpvec.addElement(dbstring);
            for(int i=cmdarr.length;i>0;i--){
                tmpvec.add(0,cmdarr[i-1]);
            }
            if(cmd.length()>0){
                blastcommand=cmd+" "+blastcommand;
                tmpvec.add(0, cmd);
            }
            cmdarr=new String[tmpvec.size()];
            tmpvec.copyInto(cmdarr);
            StringBuffer inout=new StringBuffer();
            StringBuffer errout=new StringBuffer();
            System.out.println("for "+inaln[referenceseqnum].name+" trying command="+blastcommand);
            try{
                //write the query to a file instead of using stdin; !stdin gives problems!
                PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter(tmpfilestring)));
                outwrite.println(">"+inaln[referenceseqnum].name);
                outwrite.println(inaln[referenceseqnum].seq.replaceAll("-",""));
                outwrite.close();
                BufferedReader perr;
                BufferedReader pin;
                threadstreamreader perrread;
                threadstreamreader pinread;
                Process p=rt.exec(cmdarr);
                perr=new BufferedReader(new InputStreamReader(p.getErrorStream()));//process error stream
                pin=new BufferedReader(new InputStreamReader(p.getInputStream()));//output of process
                //note: threadstreamreader closes the input bufferedreader
                //threadstreamreaders are used because otherwise deadlock can occur if the output of the process is larger than the
                //buffersize java allocates to the bufferedreaer.
                perrread=new threadstreamreader(perr,errout);
                pinread=new threadstreamreader(pin,inout);
                try{
                    perrread.start();
                    pinread.start();
                    p.waitFor();
                    if(p.exitValue()!=0){
                        System.err.println("non-perfect exit from blast for "+inaln[referenceseqnum].name);
                    }
                    while (pinread.done==false){
                        //System.out.println("waiting for process reader to finish");
                        synchronized(inout){
                            try{
                                inout.wait(10l);
                            }catch(InterruptedException e){
                                System.err.println("interrupted sleep in blastthread");
                                e.printStackTrace();
                            }
                        }
                    }
                }catch (InterruptedException e){
                    System.err.println("Interrupted process "+blastcommand+" for "+inaln[referenceseqnum].name);
                    return new hsp[0];
                }
                perrread=null;//.clear();
                pinread=null;//.clear();
            }catch (IOException ioe){
                System.err.println("IOError in "+blastcommand+" for "+inaln[referenceseqnum].name);
                for(int i=0;i<cmdarr.length;i++){
                    System.err.println(cmdarr[i]);
                }
                System.out.println("inoutsize="+inout.length());
                return new hsp[0];
            }
            System.err.print(errout.toString());
            //now the psiblast run should have generated a checkpoint file
            //then do the last round using the psiblast checkpoint file against blastdb
            System.out.println("actual psiblast run on "+inaln[referenceseqnum].name);
            blastcommand=blastpath+" -db "+basename+" -html -show_gis -evalue "+mineval+" -query "+tmpfilestring+" -in_pssm "+tmpcheckfile;
        }else{//if I only want to use the sequences database it is a lot easier
            if(blastpath.contains("blastall")){//non blast+ version
                blastcommand=blastpath+" -d "+basename+" -T T -e "+mineval+" -i "+tmpfilestring;//input will be stdin//or alternately a infile
            }else{//blast+ version of the programs
                blastcommand=blastpath+" -db "+basename+" -html -show_gis -evalue "+mineval+" -query "+tmpfilestring;//input will be stdin//or alternately a infile
            }
        }
        //now I need to run the blast (if referencedb I have a checkpoint file, else not.
        if(cmd.length()>0){
            blastcommand=cmd+" "+blastcommand;
        }
        StringBuffer inout=new StringBuffer();
        StringBuffer errout=new StringBuffer();
        System.out.println("for "+inaln[referenceseqnum].name+" trying command="+blastcommand);
        String myblast="";
        try{
            //write the query to a file instead of using stdin; !stdin gives problems!
            PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter(tmpfilestring)));
            outwrite.println(">"+inaln[referenceseqnum].name);
            outwrite.println(inaln[referenceseqnum].seq.replaceAll("-",""));
            outwrite.close();
            BufferedReader perr;
            BufferedReader pin;
            threadstreamreader perrread;
            threadstreamreader pinread;
            Process p=rt.exec(blastcommand);
            perr=new BufferedReader(new InputStreamReader(p.getErrorStream()));//process error stream
            pin=new BufferedReader(new InputStreamReader(p.getInputStream()));//output of process
            //note: threadstreamreader closes the input bufferedreader
            //threadstreamreaders are used because otherwise deadlock can occur if the output of the process is larger than the
            //buffersize java allocates to the bufferedreaer.
            perrread=new threadstreamreader(perr,errout);
            pinread=new threadstreamreader(pin,inout);
            try{
                perrread.start();
                pinread.start();
                p.waitFor();
                if(p.exitValue()!=0){
                    System.err.println("non-perfect exit from blast for "+inaln[referenceseqnum].name);
                }
                System.out.println("blast done for "+inaln[referenceseqnum].name);
                while (pinread.done==false){
                    //System.out.println("waiting for process reader to finish");
                    synchronized(inout){
                        try{
                            inout.wait(10l);
                        }catch(InterruptedException e){
                            System.err.println("interrupted sleep in blastthread");
                            e.printStackTrace();
                        }
                    }
                }
            }catch (InterruptedException e){
                System.err.println("Interrupted process "+blastcommand+" for "+inaln[referenceseqnum].name);
                return new hsp[0];
            }
            myblast=inout.toString();
            perrread=null;//.clear();
            pinread=null;//.clear();
        }catch (IOException ioe){
            System.err.println("IOError in "+blastcommand+" for "+inaln[referenceseqnum].name);
            return new hsp[0];
        }
        //parse the hsp's from myblast
        retvec = hspget.get(myblast, retvec, mineval, minpval, 0f, 0f, 0f, 0, nameshash, false);
        System.out.println("hsp's for "+inaln[referenceseqnum].name+" = "+retvec.size());
        //now output the hsp elements
        hsp[] retarr=new hsp[retvec.size()];
        retvec.copyInto(retarr);
        //for(int i=0;i<retarr.length;i++){
        //    System.out.println("hit:"+retarr[i].hname+" pval="+retarr[i].pvalue);
        //}
        //now remove the temporary files
        File rmfile=new File(basename);
        if(rmfile.exists()){
            rmfile.delete();
        }
        rmfile=new File(basename+".query");
        if(rmfile.exists()){
            rmfile.delete();
        }
        rmfile=new File(basename+".chkpnt");
        if(rmfile.exists()){
            rmfile.delete();
        }
        rmfile=new File(basename+".psd");
        if(rmfile.exists()){
            rmfile.delete();
        }
        rmfile=new File(basename+".phr");
        if(rmfile.exists()){
            rmfile.delete();
        }
        rmfile=new File(basename+".pin");
        if(rmfile.exists()){
            rmfile.delete();
        }
        rmfile=new File(basename+".pnd");
        if(rmfile.exists()){
            rmfile.delete();
        }
        rmfile=new File(basename+".pni");
        if(rmfile.exists()){
            rmfile.delete();
        }
        rmfile=new File(basename+".psi");
        if(rmfile.exists()){
            rmfile.delete();
        }
        rmfile=new File(basename+".psq");
        if(rmfile.exists()){
            rmfile.delete();
        }
        rmfile=new File("formatdb.log");
        if(rmfile.exists()){
            rmfile.delete();
        }
        return retarr;
    }//end gethsps
    
}//end class
