/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package clans.io;
import java.util.*;
import java.io.*;

import clans.model.proteins.HighScoringSegmentPair;
/**
 *
 * @author tancred
 */
public class HspReaderVersion2 {

    public static ArrayList<HighScoringSegmentPair> get(String tmpoutfile,double eval,double pval,float coverage,float scval,float ident,int verbose){
        ArrayList<HighScoringSegmentPair> retlist=new ArrayList<HighScoringSegmentPair>();
        //System.out.println("reading hsp data for file '"+tmpoutfile+"'");
        //this should open a BLAST xml output file and read the results
        try{
            BufferedReader inread=new BufferedReader(new FileReader(tmpoutfile));
            //first see whether the file completed all right
            String inline;
            double effspace=-1;
            double dbsize=-1;
            while((inline=inread.readLine())!=null){
                if(inline.contains("<Statistics_db-len>")){
                    inline=inline.substring(inline.indexOf(">")+1,inline.indexOf("</")).trim();
                    try{
                       dbsize=Double.parseDouble(inline);
                    }catch(NumberFormatException ne){
                        System.err.println("ERROR parsing number for effective search space from String '"+inline+"'");
                    }
                }
                if(inline.contains("<Statistics_eff-space>")){
                    //System.out.println("found effspace line '"+inline+"'");
                    inline=inline.substring(inline.indexOf(">")+1,inline.indexOf("</")).trim();
                    try{
                        effspace=Double.parseDouble(inline);
                    }catch(NumberFormatException ne){
                        System.err.println("ERROR parsing number for effective search space from String '"+inline+"'");
                    }
                    break;
                }
            }//end while reading first time
            //System.out.println("effspace="+effspace);
            if(effspace<0){
                System.err.println("WARNING: could not retrieve effective sequence space from blast file: effspace="+effspace);
                inread.close();
                return null;
            }else if(effspace==0){
                //sometimes happens with the legacy blast
                effspace=dbsize;//bad approximation, but better than nothing at all!
                //the new blast+ however, sometimes assigns negative values to the database size (roll-over error?)
                if(effspace<0){
                    effspace=-effspace;
                }
            }
            double checkpval=pval;
            if(pval<0){
            	//if pval was not set
            	checkpval=eval/effspace;
            }
            inread.close();
            inread=new BufferedReader(new FileReader(tmpoutfile));
            //now read and parse the hsp's from the file
            //<Iteration_query-def> encompasses the new query name
            //<Iteration_query-len> gives the query length
            //<Hit> to </Hit> covers one HSP
            //<Hit_def> to </Hit_def> gives the hit name
            //<Hsp_evalue> gives the evalue
            //<Hsp_query-from> and
            //<Hsp_query-to> give the coverage of the HSP to the query
            //<Hsp_identity> gives the identities
            //<Hsp_bit-score> gives the score
            String qname="";
            String hname="";
            int qlength=-1;
            int hlength=-1;
            int qstart=-1,qend=-1;
            double evalue=-1;
            //int identities;
            float score=-1;
            while((inline=inread.readLine())!=null){
                if(inline.contains("<Iteration_query-def>")){
                    qname=inline.substring(inline.indexOf(">")+1,inline.indexOf("</"));
                }else if(inline.contains("<Hit_def>")){
                    hname=inline.substring(inline.indexOf(">")+1,inline.indexOf("</"));
                }else if(inline.contains("<Hit_len>")){
                    try{
                        hlength=Integer.parseInt(inline.substring(inline.indexOf(">")+1,inline.indexOf("</")));
                    }catch (NumberFormatException ne){
                        System.err.println("ERROR, unable to parse int from '"+inline.substring(inline.indexOf(">")+1,inline.indexOf("</"))+"'");
                        hlength=-1;
                    }
                }else if(inline.contains("<Hit_def>")){
                    hname=inline.substring(inline.indexOf(">")+1,inline.indexOf("</"));
                }else if(inline.contains("<Iteration_query-len>")){
                    try{
                        qlength=Integer.parseInt(inline.substring(inline.indexOf(">")+1,inline.indexOf("</")));
                    }catch (NumberFormatException ne){
                        System.err.println("ERROR, unable to parse int from '"+inline.substring(inline.indexOf(">")+1,inline.indexOf("</"))+"'");
                        qlength=-1;
                    }
                }else if(inline.contains("<Hsp_evalue>")){
                    try{
                        evalue=Double.parseDouble(inline.substring(inline.indexOf(">")+1,inline.indexOf("</")));
                    }catch (NumberFormatException ne){
                        System.err.println("ERROR, unable to parse double from '"+inline.substring(inline.indexOf(">")+1,inline.indexOf("</"))+"'");
                        evalue=-1;
                    }
                }else if(inline.contains("<Hsp_bit-score>")){
                    try{
                        score=Float.parseFloat(inline.substring(inline.indexOf(">")+1,inline.indexOf("</")));
                    }catch (NumberFormatException ne){
                        System.err.println("ERROR, unable to parse float from '"+inline.substring(inline.indexOf(">")+1,inline.indexOf("</"))+"'");
                        score=-1;
                    }
                }else if(inline.contains("<Hsp_query-from>")){
                    try{
                        qstart=Integer.parseInt(inline.substring(inline.indexOf(">")+1,inline.indexOf("</")));
                    }catch (NumberFormatException ne){
                        System.err.println("ERROR, unable to parse int from '"+inline.substring(inline.indexOf(">")+1,inline.indexOf("</"))+"'");
                        qstart=-1;
                    }
                }else if(inline.contains("<Hsp_query-to>")){
                    try{
                        qend=Integer.parseInt(inline.substring(inline.indexOf(">")+1,inline.indexOf("</")));
                    }catch (NumberFormatException ne){
                        System.err.println("ERROR, unable to parse int from '"+inline.substring(inline.indexOf(">")+1,inline.indexOf("</"))+"'");
                        qend=-1;
                    }
                }else if(inline.contains("</Hsp>")){
                    //System.out.println("reached end of HSP: qlength="+qlength+" evalue="+evalue+" score="+score+" qstart="+qstart+" qend="+qend);
                    if(qlength>0 && evalue>=0 && score>0 && qstart>=0 && qend>0){
                        //System.out.println("\tpassing filter1; evalue="+evalue+" effspace="+effspace+" scval="+scval);
                        HighScoringSegmentPair newhsp=new HighScoringSegmentPair();
                        int uselength=qlength;
                        if(hlength<qlength){
                            uselength=hlength;
                        }
                        //float mycover=(float)(qend-qstart)/(float)uselength;
                        float myscval=(float)score/(float)uselength;
                        double mypval=evalue/effspace;
                        if(scval>0){//if I want to filter by bit-score per column
                            if(myscval>scval){//if I want to check for score/column
                                //System.out.println("\tpassing filter 2 myscval="+myscval);
                                newhsp.value=myscval;
                                newhsp.qname=qname;
                                newhsp.hname=hname;
                                retlist.add(newhsp);
                            }//else{//else don't add it
                                //System.out.println(myscval+" is smaller than "+scval+" not adding!");
                            //}
                        }else if(checkpval>=0 && mypval<=checkpval){//only look at this if the scval was NOT selected
                            //System.out.println("\tpassing filter 2 mypval="+mypval);
                            newhsp.qname=qname;
                            newhsp.hname=hname;
                            newhsp.value=mypval;
                            retlist.add(newhsp);
                        }
                    }
                    //qlength=-1;//don't set qlength back to -1 as I only get this value once per sequence not once per hsp
                    evalue=-1;
                    score=-1;
                    qstart=-1;
                    qend=-1;
                }
            }//end while reading second round
            inread.close();
        }catch (IOException ioe){
            System.err.println("IOERROR reading from '"+tmpoutfile+"'");
            return null;
        }
        return retlist;
    }//end get

}
