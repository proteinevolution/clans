package clans;
import java.io.*;
import java.util.*;
/**
 *
 * @author  tancred
 */
public class readsave {
    
    public static int[] checkblast(String filename,int seqnum,boolean skipcheck){
        //see if the data in tmpblasthsp is complete
        //return a two element array containing
        int[] retarr={0,0};
        HashMap<String, Object> nameshash=new HashMap<String, Object>();
        try{
            BufferedReader inread=new BufferedReader(new FileReader(filename));
            String inline;
            String tmpstr="";
            while ((inline=inread.readLine())!=null){
                if(inline.equalsIgnoreCase("DONE")){
                    if(skipcheck==false){
                        retarr[0]=seqnum;
                        retarr[1]=seqnum;
                        inread.close();
                        return retarr;
                    }
                }else if(inline.startsWith("hsp:")){
                    //parse the first number following this
                    if(inline.indexOf(";")>0){
                        tmpstr=inline.substring(5,inline.indexOf(";"));
                        try{
                            retarr[1]=Integer.parseInt(tmpstr);
                            nameshash.put(tmpstr,null);
                        }catch (NumberFormatException ne){
                            System.err.println("unable to parse int from "+tmpstr);
                            System.exit(1);
                        }
                    }
                }
            }//end while
            inread.close();
            retarr[0]=nameshash.size();
        }catch (IOException e){
            System.err.println("IOERROR reading from "+filename);
            System.exit(1);
        }//end catch
        return retarr;
    }//end checkblast
    
    //--------------------------------------------------------------------------
    
    public static MinimalHsp[] blast(String filename,double cutoff){
        //read the hsp data from a blast savefile
        MinimalHsp[] retarr=new MinimalHsp[0];
        try{
            BufferedReader inread=new BufferedReader(new FileReader(filename));
            String inline=inread.readLine();
            if(inline==null){
                System.err.println("Error reading former blast results from file: "+filename+"; NO DATA!; restart the program with the option \"-readblast F\" or remove that file.");
                System.exit(1);
            }
            //first line should be the number of array elements to create
            int seqs=0;
            try{
                seqs=Integer.parseInt(inline.substring(0,inline.indexOf(" sequences")));
            }catch (NumberFormatException num){
                System.err.println("NumberFormatError, Wrong Format!");
                System.exit(0);
            }
            MinimalHsp currhsp=new MinimalHsp();
            int ival=-1;
            int jval=-1;
            int kval=-1;
            int hspcount=-1;
            double pval=-1;
            String[] tmparr;
            HashMap hsphash=new HashMap();
            String hspkey;
            while ((inline=inread.readLine())!=null){
                if(inline.startsWith("#")){
                    continue;
                }
                if(inline.startsWith("hsp: ")){
                    hspcount++;
                    if(hspcount%10000==0){
                        System.out.print(hspcount+";");
                    }
                    if(ival!=jval && pval<=cutoff){
                        if((ival>-1)&&(jval>-1)&&(kval>-1)){
                            hspkey=ival+"_"+jval;
                            if(hsphash.containsKey(hspkey)==false){
                                currhsp.val=new double[1];
                                currhsp.val[0]=pval;
                                hsphash.put(hspkey,currhsp);
                            }else{
                                ((MinimalHsp)hsphash.get(hspkey)).addpval(pval);
                            }
                        }
                    }
                    tmparr=(inline.substring(5)).split(";",0);
                    if(java.lang.reflect.Array.getLength(tmparr)==3){
                        try{
                            ival=Integer.parseInt(tmparr[0]);
                            jval=Integer.parseInt(tmparr[1]);
                            kval=Integer.parseInt(tmparr[2]);
                        }catch (NumberFormatException nume){
                            System.err.println("unable to parse correct numbers from "+inline);
                            System.exit(0);
                        }
                        currhsp=new MinimalHsp();
                        currhsp.query=ival;
                        currhsp.hit=jval;
                    }else{
                        currhsp=new MinimalHsp();
                        pval=1;
                    }
                }else if(inline.startsWith("\tqname: ")){
                    //currhsp.qname=inline.substring(8);
                }else if(inline.startsWith("\tqseq: ")){
                    //currhsp.qseq=inline.substring(7);
                }else if(inline.startsWith("\thname: ")){
                    //currhsp.hname=inline.substring(8);
                }else if(inline.startsWith("\thseq: ")){
                    //currhsp.hseq=inline.substring(7);
                }else if(inline.startsWith("\tqstart: ")){
                }else if(inline.startsWith("\tqend: ")){
                }else if(inline.startsWith("\thstart: ")){
                }else if(inline.startsWith("\thend: ")){
                }else if(inline.startsWith("\tevalue: ")){
                }else if(inline.startsWith("\tvalue: ")||inline.startsWith("\tpvalue: ")){
                    //changed from "pvalue" to "value" as I want to also remember "scores per column" (saved as negative values)
                    //a frequent error on re-starting the blast run is that a "sequences" is present on the value line at some point
                    //if present, simply ignore it
                    if(inline.indexOf("sequences")>0){
                        inline=inline.substring(0,inline.indexOf("sequences"));
                    }
                    try{
                        pval=Double.parseDouble(inline.substring(8));
                    }catch (NumberFormatException num){
                        System.err.println("unable to parse number from "+inline);
                        System.exit(0);
                    }
                }else if(inline.startsWith("\tscval: ")){
                }else if(inline.startsWith("\tcoverage: ")){
                }else if(inline.startsWith("\tident: ")){
                }else if(inline.startsWith("\tdblength: ")){
                }else if(inline.startsWith("\tdblengtheff: ")){
                }else if(inline.startsWith("\tqlength: ")){
                }else if(inline.startsWith("\tqlengtheff: ")){
                }else if(inline.equalsIgnoreCase("DONE")){
                    System.out.println("done reading");
                    break;
                }else{
                    System.err.println("unknown line:'"+inline+"'");
                }
            }//end while reading
            if(pval<=cutoff){
                if((ival>-1)&&(jval>-1)&&(kval>-1)){
                    hspkey=ival+"_"+jval;
                    if(hsphash.containsKey(hspkey)==false){
                        currhsp.val=new double[1];
                        currhsp.val[0]=pval;
                        hsphash.put(hspkey,currhsp);
                    }else{
                        ((MinimalHsp)hsphash.get(hspkey)).addpval(pval);
                    }
                }
            }
            int elements=hsphash.size();
            retarr=new MinimalHsp[elements];
            System.out.println("hsp's="+hsphash.size());
            retarr=(MinimalHsp[])(hsphash.values().toArray(retarr));
        }catch (IOException e){
            System.out.println("IOerror reading blastfile "+filename);
            System.exit(0);
        }
        System.out.println("done reading");
        //now get all the minhsp elements form the hash
        return retarr;
    }//end blast
    
    //--------------------------------------------------------------------------
    
    static void parse_intermediate_results(ClusterData data){
        try{
            BufferedReader inread=new BufferedReader(new FileReader(data.getIntermediateResultfileName()));
            String inline=inread.readLine();
            int sequences=0;
            try{
                sequences=Integer.parseInt(inline.substring(11));
            }catch (NumberFormatException nume){
                System.out.println("unable to parse number of sequences from first line:"+inline);
            }
            float[][] posarr=new float[sequences][3];
            boolean readpos=false;
            String[] tmpstr;
            int ipos;
            while ((inline=inread.readLine())!=null){
                
                if (inline.startsWith("rounds")) {
                    try {
                        String[] x = inline.split(" ");
                        data.rounds = Integer.parseInt(x[1]);
                    } catch (NumberFormatException num) {
                        System.err.println("couldn't parse rounds line (" + inline + ") from "
                                + data.getIntermediateResultfileName());
                    }
                    continue;
                }
                
                if(inline.equals("positions:")){
                    readpos=true;
                    continue;
                }
                if(readpos){
                    tmpstr=inline.split(";");
                    try{
                        ipos=Integer.parseInt(tmpstr[0]);
                        posarr[ipos][0]=Float.parseFloat(tmpstr[1]);
                        posarr[ipos][1]=Float.parseFloat(tmpstr[2]);
                        posarr[ipos][2]=Float.parseFloat(tmpstr[3]);
                    }catch (NumberFormatException num){
                        System.out.println("unable to parse positions data from "+inline);
                        System.out.println(tmpstr[0]);
                        System.out.println(tmpstr[1]);
                        System.out.println(tmpstr[2]);
                        System.out.println(tmpstr[3]);
                        inread.close();
                        data.myposarr = posarr;
                        return;
                    }
                }//end if readpos
            }//end while reading
            inread.close();
            data.myposarr = posarr;
            return;
        }catch (IOException e){
            System.err.println("unable to read from " + data.getIntermediateResultfileName());
            e.printStackTrace();
            data.myposarr = new float[0][0];
            return;
        }
    }//end readpos
    
}
