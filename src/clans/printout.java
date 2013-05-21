/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package clans;
import java.io.*;
import java.util.*;
/**
 *
 * @author  tancred
 */
public class printout {
    
    /** Creates a new instance of printaln */
    public printout() {
    }
    
    static void printaln(AminoAcidSequence[] inaln,String outfilename, String oformat, int blocksize){
        //print the working alignment to file outfilename in the specified format
        if(oformat.equalsIgnoreCase("clustal")){
            printclustal(inaln,outfilename,blocksize);
        }else{
            printfasta(inaln,outfilename);
        }
    }//end printaln
    
    //--------------------------------------------------------------------------
    
    static void printconf(AminoAcidSequence[] inaln, double[][] confarr, String outfilename, String oformat, int blocksize){
        //print the position specific confidences for each sequence
        if(oformat.equalsIgnoreCase("clustal")){
            printclustal(inaln,confarr,outfilename,blocksize);
        }else{
            printfasta(inaln,confarr,outfilename);
        }
    }//end printconf
    
    //--------------------------------------------------------------------------
    
    static void printfiltered(AminoAcidSequence[] inaln, double[][] confarr, String outfilename, String oformat, int blocksize, float filtercutoff){
        //print the input alignment filtered for confidence.
        //replace low confidence residues by dashes
        int elements=java.lang.reflect.Array.getLength(inaln);
        int seqlength=inaln[0].seq.length();
        AminoAcidSequence[] tmpaln=new AminoAcidSequence[elements];
        for(int i=0;i<elements;i++){
            tmpaln[i]=new AminoAcidSequence();
            tmpaln[i].name=inaln[i].name;
            tmpaln[i].seq=inaln[i].seq;
        }//end for i
        int maxnamelength=0;
        String[] namearr=new String[elements];
        int j;
        StringBuffer tmpstr=new StringBuffer(blocksize+3);
        int position=0;
        for(int i=0;i<elements;i++){
            tmpstr.setLength(0);
            for(j=0;j<seqlength;j++){
                if(confarr[i][j]>=filtercutoff){
                    tmpstr.append(tmpaln[i].seq.charAt(j));
                }else{
                    tmpstr.append("-");
                }
            }//end for j
            tmpaln[i].seq=tmpstr.toString();
        }//end for i
        if(oformat.equalsIgnoreCase("clustal")){
            printclustal(tmpaln,outfilename,blocksize);
        }else{
            printfasta(tmpaln,outfilename);
        }
    }//end printfiltered
    
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    
    static void printfasta(AminoAcidSequence[] inaln, String outfilename){
        //print the alignment in fasta format
        try{
            PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter(outfilename)));
            int elements=java.lang.reflect.Array.getLength(inaln);
            for(int i=0;i<elements;i++){
                outwrite.println(">"+inaln[i].name);
                outwrite.println(inaln[i].seq);
            }//end for i
            outwrite.close();
        }catch (IOException e){
            System.err.println("IOERROR in printfiltered");
            e.printStackTrace();
        }
    }//end printfasta
    
    static void printfasta(AminoAcidSequence[] inaln, File outfile){
        //print the alignment in fasta format
        try{
            PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter(outfile)));
            int elements=java.lang.reflect.Array.getLength(inaln);
            for(int i=0;i<elements;i++){
                outwrite.println(">"+inaln[i].name);
                outwrite.println(inaln[i].seq);
            }//end for i
            outwrite.close();
        }catch (IOException e){
            System.err.println("IOERROR in printfiltered");
            e.printStackTrace();
        }
    }//end printfasta
    
    static void printfasta(AminoAcidSequence[] inaln, double[][] confarr, String outfilename){
        //print the confarr in fasta format
        try{
            PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter(outfilename)));
            int elements=java.lang.reflect.Array.getLength(inaln);
            if(elements<1){
                System.err.println("ERROR: empty alignment");
                return;
            }
            int seqlength=inaln[0].seq.length();
            int j;
            StringBuffer tmp=new StringBuffer();
            for(int i=0;i<elements;i++){
                outwrite.println(">"+inaln[i].name);
                for(j=0;j<seqlength;j++){
                    outwrite.print(((float)(int)(confarr[i][j]*100))/100+";");
                }//end for j
                outwrite.println();
            }//end for i
            outwrite.close();
        }catch (IOException e){
            System.err.println("IOERROR in printfiltered");
            e.printStackTrace();
        }
    }//end printfasta
    
    //--------------------------------------------------------------------------
    
    static void printclustal(AminoAcidSequence[] inaln, String outfilename, int blocksize){
        //print the alignment in clustal format
        try{
            PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter(outfilename)));
            int elements=java.lang.reflect.Array.getLength(inaln);
            if(elements<1){
                System.err.println("ERROR: empty alignment");
                return;
            }
            int position=0;
            int seqlength=inaln[0].seq.length();
            int maxnamelength=10;
            int i;
            for(i=0;i<elements;i++){
                if(inaln[i].name.length()>maxnamelength){
                    maxnamelength=inaln[i].name.length();
                }
            }//end for i
            StringBuffer tmpstr=new StringBuffer(maxnamelength+3);
            String[] namearr=new String[elements];
            for(i=0;i<elements;i++){
                tmpstr.setLength(0);
                tmpstr.append(inaln[i].name);
                while (tmpstr.length()<maxnamelength){
                    tmpstr.append(" ");
                }
                namearr[i]=tmpstr.toString();
            }//end for i
            outwrite.println("CLUSTAL");
            outwrite.println();
            outwrite.println();
            while ((position+blocksize)<=seqlength){
                for(i=0;i<elements;i++){
                    outwrite.println(namearr[i]+"  "+inaln[i].seq.substring(position,position+blocksize));
                }//end for i
                outwrite.println();
                position+=blocksize;
            }//end while
            for(i=0;i<elements;i++){
                tmpstr.setLength(0);
                tmpstr.append(inaln[i].name);
                while (tmpstr.length()<maxnamelength){
                    tmpstr.append(" ");
                }
                outwrite.println(tmpstr.toString()+"  "+inaln[i].seq.substring(position));
            }//end for i
            outwrite.println();
            outwrite.close();
        }catch (IOException e){
            System.err.println("IOERROR in printfiltered");
            e.printStackTrace();
        }
    }//end printclustal
    
    static void printclustal(AminoAcidSequence[] inaln, double[][] confarr, String outfilename, int blocksize){
        //print the confarr in clustal format
        try{
            PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter(outfilename)));
            int elements=java.lang.reflect.Array.getLength(inaln);
            if(elements<1){
                System.err.println("ERROR: empty alignment");
                return;
            }
            int position=0;
            int seqlength=inaln[0].seq.length();
            int maxnamelength=10;
            int i,j;
            for(i=0;i<elements;i++){
                if(inaln[i].name.length()>maxnamelength){
                    maxnamelength=inaln[i].name.length();
                }
            }//end for i
            StringBuffer tmpstr=new StringBuffer(maxnamelength+3);
            String[] namearr=new String[elements];
            for(i=0;i<elements;i++){
                tmpstr.setLength(0);
                tmpstr.append(inaln[i].name);
                while (tmpstr.length()<maxnamelength){
                    tmpstr.append(" ");
                }
                namearr[i]=tmpstr.toString();
            }//end for i
            outwrite.println("CLUSTAL");
            outwrite.println();
            outwrite.println();
            while ((position+blocksize)<=seqlength){
                for(i=0;i<elements;i++){
                    outwrite.print(namearr[i]+"  ");
                    for(j=position;j<position+blocksize;j++){
                        outwrite.print(((float)(int)(confarr[i][j]*100))/100+";");
                    }//end for j
                    outwrite.println();
                }//end for i
                outwrite.println();
                position+=blocksize;
            }//end while
            for(i=0;i<elements;i++){
                outwrite.print(namearr[i]+"  ");
                for(j=position;j<seqlength;j++){
                    outwrite.print(confarr[i][j]+";");
                }//end for j
                outwrite.println();
            }//end for i
            outwrite.println();
            outwrite.close();
        }catch (IOException e){
            System.err.println("IOERROR in printfiltered");
            e.printStackTrace();
        }
    }//end printclustal
    
    //--------------------------------------------------------------------------
    
    static void saveblast(String outfilename, Vector[][] blastvals){
        try{
            PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter(outfilename)));
            int vecsize=java.lang.reflect.Array.getLength(blastvals);
            int i,j,k;
            int tmpvecsize;
            Vector tmpvec;
            hsp tmphsp;
            outwrite.println(vecsize+" sequences");
            for(i=0;i<vecsize;i++){
                for(j=0;j<vecsize;j++){
                    //print the hsp data for this sequence comparison
                    tmpvec=blastvals[i][j];
                    tmpvecsize=tmpvec.size();
                    for(k=0;k<tmpvecsize;k++){
                        tmphsp=(hsp)tmpvec.elementAt(k);
                        //now print the hsp data
                        outwrite.println("hsp: "+i+";"+j+";"+k);
                        //outwrite.println("\tqname: "+tmphsp.qname);
                        //outwrite.println("\tqseq: "+tmphsp.qseq);
                        //outwrite.println("\thname: "+tmphsp.hname);
                        //outwrite.println("\thseq: "+tmphsp.hseq);
                        //outwrite.println("\tqstart: "+tmphsp.qstart);
                        //outwrite.println("\tqend: "+tmphsp.qend);
                        //outwrite.println("\thstart: "+tmphsp.hstart);
                        //outwrite.println("\thend: "+tmphsp.hend);
                        outwrite.println("\tvalue: "+tmphsp.value);
                    }//end for k
                }//end for j
            }//end for i
            outwrite.close();
        }catch (IOException e){
            System.out.println("Error printing blast results to "+outfilename);
        }
    }//end saveblast
    
    static void saveblastappend(PrintWriter outwrite, Vector[] blastvals, int i){
        int vecsize=java.lang.reflect.Array.getLength(blastvals);
        int j,k;
        int tmpvecsize;
        Vector tmpvec;
        hsp tmphsp;
        for(j=0;j<vecsize;j++){
            //print the hsp data for this sequence comparison
            tmpvec=blastvals[j];
            if(tmpvec!=null){
                tmpvecsize=tmpvec.size();
                for(k=0;k<tmpvecsize;k++){
                    tmphsp=(hsp)tmpvec.elementAt(k);
                    //now print the hsp data
                    outwrite.println("hsp: "+i+";"+j+";"+k);
                    //outwrite.println("\tqname: "+tmphsp.qname);
                    //outwrite.println("\tqseq: "+tmphsp.qseq);
                    //outwrite.println("\thname: "+tmphsp.hname);
                    //outwrite.println("\thseq: "+tmphsp.hseq);
                    //outwrite.println("\tqstart: "+tmphsp.qstart);
                    //outwrite.println("\tqend: "+tmphsp.qend);
                    //outwrite.println("\thstart: "+tmphsp.hstart);
                    //outwrite.println("\thend: "+tmphsp.hend);
                    outwrite.println("\tvalue: "+tmphsp.value);
                }//end for k
            }
        }//end for j
    }//end saveblastappend
    
    //--------------------------------------------------------------------------
    
    static void save_semicolon_delimited_positions(PrintWriter outwrite, float[][] positions) {
        outwrite.println("positions:");
        for (int i = 0; i < positions.length; i++) {
            outwrite.println(i + ";" + positions[i][0] + ";" + positions[i][1] + ";" + positions[i][2]);
        }
        outwrite.flush();
    }
    
}
