/*
 * utils.java
 *
 * Created on January 11, 2006, 12:00 PM
 */
package clans;
import java.io.*;
import java.util.*;

/**
 *
 * @author  tancred
 */
public class affyutils {
    
    /** Creates a new instance of utils */
    public affyutils() {
    }
    
    public static int load(File infile, affydialog parent){
        //save all the possible parameters of the parent object to a file
        try{
            BufferedReader inread=new BufferedReader(new FileReader(infile));
            String inline;
            String[] tmparr;
            String tmpstr;
            int tmpint;
            while((inline=inread.readLine())!=null){
                if(inline.startsWith("usefoldchange=")){
                    tmpstr=inline.substring(14);
                    if(tmpstr.startsWith("true")){
                        parent.parent.usefoldchange=true;
                    }else{
                        parent.parent.usefoldchange=false;
                    }
                }else if(inline.startsWith("avgfoldchange=")){
                    tmpstr=inline.substring(14);
                    if(tmpstr.startsWith("true")){
                        parent.parent.avgfoldchange=true;
                    }else{
                        parent.parent.avgfoldchange=false;
                    }
                }else if(inline.startsWith("wtfiles=")){
                    tmparr=inline.substring(8).split(";");
                    tmpint=java.lang.reflect.Array.getLength(tmparr);
                    parent.wtfiles=new File[tmpint];
                    for(int i=0;i<tmpint;i++){
                        parent.wtfiles[i]=new File(tmparr[i]);
                    }//end for i
                }else if(inline.startsWith("<data>")){
                    replicates myrep=null;
                    while((inline=inread.readLine()).equals("</data>")==false){
                        if(inline==null){
                            break;
                        }
                        if(inline.startsWith("name=")){
                            if(myrep!=null){
                                parent.datavec.addElement(myrep);
                            }
                            myrep=new replicates();
                            myrep.name=inline.substring(5);
                            //System.out.println("name="+myrep.name);
                        }else if(inline.startsWith("wtname=")){
                            myrep.wtname=inline.substring(7);
                            //System.out.println("wtname="+myrep.wtname);
                        }else if(inline.startsWith("replicate=")){
                            tmparr=inline.substring(10).split(";");
                            tmpint=java.lang.reflect.Array.getLength(tmparr);
                            myrep.replicate=new File[tmpint];
                            for(int i=0;i<tmpint;i++){
                                myrep.replicate[i]=new File(tmparr[i]);
                            }//end for i
                            myrep.replicates=tmpint;
                            
                        }else if(inline.startsWith("wtreplicate=")){
                            tmparr=inline.substring(12).split(";");
                            tmpint=java.lang.reflect.Array.getLength(tmparr);
                            myrep.wtreplicate=new File[tmpint];
                            for(int i=0;i<tmpint;i++){
                                myrep.wtreplicate[i]=new File(tmparr[i]);
                            }//end for i
                            myrep.wtreplicates=tmpint;
                            
                        }else{
                            System.err.println("unknown line:'"+inline+"'");
                        }
                    }//end while reading
                    if(myrep!=null){
                        parent.datavec.addElement(myrep);
                    }
                    //System.out.println("done reading data: "+parent.datavec.size()+" elements");
                }else{
                    System.err.println("unknown line:'"+inline+"'");
                }
            }//end while reading
            inread.close();
        }catch (IOException ioe){
            System.err.println("IOError writing to "+infile.getName());
            return -1;
        }
        return 0;
    }//end save
    
    //--------------------------------------------------------------------------
    
    static float getglobalrelstdev(Vector datvec){
        //calculate a value giving me the average noise ratio over this dataset
        //then return a float that, when multiplied with the average value for each sequence
        //tells me wether my datapoint is outside the "noise" range
        int datnum=datvec.size();
        float retval=0;
        float mystdev;
        float myavg,tmp;
        datapoint[] datarr;
        int conditions;
        int subtract;
        for(int i=datnum;--i>=0;){
            datarr=(datapoint[])datvec.get(i);
            myavg=0;
            mystdev=0;
            //use the slower version (2-pass) as that limits roundoff errors
            conditions=java.lang.reflect.Array.getLength(datarr);
            subtract=0;
            for(int j=conditions;--j>=0;){
                if(datarr[j]==null){//possible if one data set is missing some entries
                    subtract++;
                }else{
                    myavg+=datarr[j].value;
                }
            }//end for j
            myavg/=(conditions-subtract);
            for(int j=conditions;--j>=0;){
                if(datarr[j]!=null){
                    tmp=datarr[j].value-myavg;
                    mystdev+=(tmp*tmp);
                }
            }
            mystdev/=(conditions-1-subtract);
            mystdev=(float)java.lang.Math.sqrt(mystdev);
            retval+=mystdev/myavg;
        }//end for i
        retval/=datnum;
        return retval;
    }//end getglobalrelstdev
    
    //--------------------------------------------------------------------------
    
    static Vector getfoldchange(Vector datvec, boolean useavg){
        int datnum=datvec.size();
        datapoint[] tmparr;
        int arrsize;
        if(useavg){
            int valint;
            //the amount of data doesn't change
            for(int i=0;i<datnum;i++){
                tmparr=(datapoint[])datvec.elementAt(i);
                arrsize=java.lang.reflect.Array.getLength(tmparr);
                for(int j=0;j<arrsize;j++){
                    if(tmparr[j]!=null){
                        valint=java.lang.reflect.Array.getLength(tmparr[j].values);
                        for(int k=0;k<valint;k++){
                            tmparr[j].values[k]/=tmparr[j].wtval;
                        }//end for k
                        tmparr[j].value/=tmparr[j].wtval;
                        tmparr[j].stdev/=tmparr[j].wtval;
                    }
                }//end for j
            }//end for i
        }else{
            int valint,wtint;
            float[] newarr;
            float avgval=0;
            //the datanum increases
            for(int i=0;i<datnum;i++){
                tmparr=(datapoint[])datvec.elementAt(i);
                arrsize=java.lang.reflect.Array.getLength(tmparr);
                avgval=0;
                for(int j=0;j<arrsize;j++){
                    if(tmparr[j]!=null){
                        valint=java.lang.reflect.Array.getLength(tmparr[j].values);
                        wtint=java.lang.reflect.Array.getLength(tmparr[j].wtvalues);
                        newarr=new float[valint*wtint];
                        for(int k=0;k<valint;k++){
                            for(int l=0;l<wtint;l++){
                                newarr[k*wtint+l]=tmparr[j].values[k]/tmparr[j].wtvalues[l];
                                avgval+=newarr[k*wtint+l];
                            }//end for l
                        }//end for k
                        tmparr[j].value=avgval/(valint*wtint);
                        tmparr[j].values=newarr;
                        tmparr[j].stdev=-1;
                    }
                }//end for j
            }//end for i
        }
        return datvec;
    }//end getfoldchange
        
    //--------------------------------------------------------------------------
    
    static Vector readdata(Vector datavec){
        //load the data into a vector of datapoint[] and try to use less memory
        //at the end, each vector element should correspond to all the conditions for one spot-set(i.e. gene)
        int datnum=datavec.size();
        File tmpfile;
        replicates myreplicate;
        HashMap nameshash=new HashMap();
        Vector tmpvec=new Vector();
        indata tmpdat;
        datapoint mypoint;
        double sumval,sumsqval;
        datapoint[] tmparr;
        Vector retvec=new Vector();
        int unknownval=0;
        int presentval=1;
        int absentval=-1;
        float marginalval=0.5f;
        for(int i=0;i<datnum;i++){
            myreplicate=(replicates)datavec.elementAt(i);
            System.out.println("reading "+myreplicate.name+"/"+myreplicate.wtname);
            //read the data
            nameshash.clear();
            for(int j=0;j<myreplicate.replicates;j++){
                tmpfile=myreplicate.replicate[j];
                if((tmpvec=readfile(tmpfile,absentval,marginalval,presentval,unknownval))==null){
                    tmpvec=new Vector();
                }
                for(int k=tmpvec.size();--k>=0;){
                    tmpdat=(indata)tmpvec.elementAt(k);
                    if(nameshash.containsKey(tmpdat.name)){
                        //add this value to an existing datapoint
                        mypoint=((datapoint)nameshash.get(tmpdat.name));
                        mypoint.values[j]=(float)tmpdat.val;
                        mypoint.datpresences[j]=tmpdat.presentval;
                    }else{
                        //create a new datapoint and pad with the right number of empty elements
                        mypoint=new datapoint();
                        mypoint.name=tmpdat.name;
                        mypoint.info=myreplicate.name+";"+myreplicate.wtname;
                        mypoint.values=new float[myreplicate.replicates];
                        mypoint.datpresences=new float[myreplicate.replicates];
                        mypoint.wtvalues=new float[myreplicate.wtreplicates];
                        mypoint.wtpresences=new float[myreplicate.wtreplicates];
                        for(int l=myreplicate.replicates;--l>=0;){
                            mypoint.values[l]=0;
                            mypoint.datpresences[l]=unknownval;
                        }
                        for(int l=myreplicate.wtreplicates;--l>=0;){
                            mypoint.wtvalues[l]=0;
                            mypoint.wtpresences[l]=unknownval;
                        }//end for l
                        mypoint.values[j]=(float)tmpdat.val;
                        mypoint.datpresences[j]=tmpdat.presentval;
                        nameshash.put(tmpdat.name,mypoint);
                    }
                }//end for k
            }//end for j
            //now read the wt-data
            for(int j=0;j<myreplicate.wtreplicates;j++){
                tmpfile=myreplicate.wtreplicate[j];
                if((tmpvec=readfile(tmpfile,absentval,marginalval,presentval,unknownval))==null){
                    tmpvec=new Vector();
                }
                for(int k=tmpvec.size();--k>=0;){
                    tmpdat=(indata)tmpvec.elementAt(k);
                    if(nameshash.containsKey(tmpdat.name)){
                        //add this value to an existing datapoint
                        mypoint=((datapoint)nameshash.get(tmpdat.name));
                        mypoint.wtvalues[j]=(float)tmpdat.val;
                        mypoint.wtpresences[j]=tmpdat.presentval;
                    }else{
                        //forget about this point; should never happen!
                        System.err.println("ERROR found point with no data, but reference values:"+tmpdat.name+"("+myreplicate.name+"/"+myreplicate.wtname+")");
                    }
                }//end for k
            }//end for j
            System.out.println();//needed to get a new line after reading the files
            //now calculate the avg and stdev for each datapoint
            String[] names=(String[])(nameshash.keySet().toArray(new String[0]));
            java.util.Arrays.sort(names);
            int namesnum=java.lang.reflect.Array.getLength(names);
            tmparr=new datapoint[namesnum];
            for(int j=namesnum;--j>=0;){
                mypoint=(datapoint)nameshash.get(names[j]);
                sumval=0;
                sumsqval=0;
                for(int k=0;k<myreplicate.replicates;k++){
                    mypoint.presence+=mypoint.datpresences[k];
                    sumval+=mypoint.values[k];
                    sumsqval+=(mypoint.values[k]*mypoint.values[k]);
                }//end for k
                mypoint.presence/=myreplicate.replicates;
                mypoint.value=(float)(sumval/myreplicate.replicates);
                mypoint.stdev=(float)java.lang.Math.sqrt((sumsqval/myreplicate.replicates)-((mypoint.value)*(mypoint.value)));
                sumval=0;
                sumsqval=0;
                for(int k=0;k<myreplicate.wtreplicates;k++){
                    mypoint.wtpresence+=mypoint.wtpresences[k];
                    sumval+=mypoint.wtvalues[k];
                    sumsqval+=(mypoint.wtvalues[k]*mypoint.wtvalues[k]);
                }//end for k
                mypoint.wtpresence/=myreplicate.wtreplicates;
                mypoint.wtval=(float)(sumval/myreplicate.wtreplicates);
                mypoint.wtstdev=(float)java.lang.Math.sqrt((sumsqval/myreplicate.wtreplicates)-((mypoint.wtval)*(mypoint.wtval)));
                tmparr[j]=mypoint;
            }//end for j
            retvec.addElement(tmparr);
        }//end for i
        //now I have all the data, but ordered in the wrong manner (odered by file)
        //System.out.println("conditions:"+retvec.size()+" elements:"+java.lang.reflect.Array.getLength((datapoint[])retvec.elementAt(0)));
        //I need to reorder them by gene/spot-identifier
        datapoint[] tmparr2;
        nameshash.clear();
        for(int i=0;i<retvec.size();i++){
            tmparr=(datapoint[])retvec.elementAt(i);
            for(int j=java.lang.reflect.Array.getLength(tmparr);--j>=0;){
                if(nameshash.containsKey(tmparr[j].name)){
                    tmparr2=(datapoint[])nameshash.get(tmparr[j].name);
                    tmparr2[i]=tmparr[j];
                }else{
                    tmparr2=new datapoint[retvec.size()];
                    for(int k=retvec.size();--k>=0;){
                        tmparr2[k]=null;
                    }
                    tmparr2[i]=tmparr[j];
                    nameshash.put(tmparr[j].name,tmparr2);
                }
            }//end for j
        }//end for i
        //and now reconvert the data in nameshash to an array of datapoints
        tmpvec=new Vector();
        String[] names=(String[])(nameshash.keySet().toArray(new String[0]));
        java.util.Arrays.sort(names);
        for(int i=java.lang.reflect.Array.getLength(names);--i>=0;){
            tmpvec.addElement((datapoint[])nameshash.get(names[i]));
        }//end for i
        //now post-filter the data and assign a dummy entry to all datapoints assigned as "null" 
        //(can happen if certain files do not contain certain entries)
        //for(int i=tmpvec.size();--i>=0;){
        //    tmparr2=(datapoint[])tmpvec.get(i);
        //    for(int j=java.lang.reflect.Array.getLength(tmparr2);--j>=0;){
        //        if(tmparr2[j]==null){
        //            tmparr2[j]=new datapoint(true);
        //        }
        //    }//end for j
        //}//end for i
        //System.out.println("conditions:"+java.lang.reflect.Array.getLength((datapoint[])tmpvec.elementAt(0))+" elements:"+tmpvec.size());
        return tmpvec;
    }//end readdata
    
    //--------------------------------------------------------------------------
    
    static Vector readfile(File infile,float absentval,float marginalval,float presentval,float unknownval){
        //just read the data from the file; no filtering whatsoever
        //return the data as a Vector of indata elements
        //System.out.println("in readfile for "+infile);
        Vector retvec=new Vector();
        try{
            BufferedReader inread=new BufferedReader(new FileReader(infile));
            String inline;
            String[] tmparr;
            indata currpoint;
            boolean skipme=true;
            //skip all until line starting with probe set name
            while(((inline=inread.readLine())!=null)&&skipme){
                if(inline.matches("(?i)Probe Set Name.*")||inline.matches("(?i)start")){
                    //if this is a affymetrix data file
                    skipme=false;
                }
                //second possibility; this is a converted GCRMA file; then all entries are: name;value.
                //read GCRMA data here
                if(java.lang.reflect.Array.getLength(inline.trim().split("\\s+;\\s+"))==2){//then I have a converted GCRMA file
                    inread.close();
                    //System.out.println("GCRMA file");
                    System.out.print(".");
                    return readgcrmafile(infile,presentval);
                }
            }//end while reading file
            if(skipme==false){
                while((inline=inread.readLine())!=null){
                    inline=inline.trim();
                    if(inline.length()<1){
                        continue;
                    }
                    tmparr=inline.split("\\s+",0);
                    int arrsize=java.lang.reflect.Array.getLength(tmparr);
                    if((arrsize!=6)&&(arrsize!=3)){
                        if(arrsize==2){
                            System.out.println("missing name for '"+inline+"'");
                            continue;
                        }else{
                            javax.swing.JOptionPane.showMessageDialog(new javax.swing.JFrame(),"ERROR reading from file "+infile.getAbsolutePath()+" unreadable line '"+inline+"' elements="+arrsize);
                            break;
                        }
                    }else if(arrsize==3){//if preformatted
                        //System.out.println("read line"+inline);
                        currpoint=new indata();
                        currpoint.name=tmparr[0];
                        if((currpoint.name.length()<1)|(currpoint.name.equalsIgnoreCase("blank"))){
                            continue;
                        }
                        try{
                            currpoint.val=Float.parseFloat(tmparr[2]);
                        }catch (NumberFormatException ne){
                            javax.swing.JOptionPane.showMessageDialog(new javax.swing.JFrame(),"ERROR unable to parse float from '"+tmparr[2]+"'");
                            break;
                        }
                        if(tmparr[1].equals("0")){
                            currpoint.presentval=presentval;
                        }else if(tmparr[1].equals("3")){
                            currpoint.presentval=marginalval;
                        }else if(tmparr[1].equals("2")){
                            currpoint.presentval=absentval;
                        }else{
                            System.out.println("unknown elements for '"+tmparr[1]+"'");
                        }
                        retvec.addElement(currpoint);
                    }else{//standard ? affymetrix text output?
                        currpoint=new indata();
                        currpoint.name=tmparr[0];
                        try{
                            currpoint.val=Float.parseFloat(tmparr[3]);
                        }catch (NumberFormatException ne){
                            javax.swing.JOptionPane.showMessageDialog(new javax.swing.JFrame(),"ERROR unable to parse float from '"+tmparr[3]+"'");
                            break;
                        }
                        if(tmparr[4].equalsIgnoreCase("A")){
                            currpoint.presentval=absentval;
                        }else if(tmparr[4].equalsIgnoreCase("M")){
                            currpoint.presentval=marginalval;
                        }else if(tmparr[4].equalsIgnoreCase("P")){
                            currpoint.presentval=presentval;
                        }else{
                            currpoint.presentval=unknownval;
                        }
                        retvec.addElement(currpoint);
                    }
                }//end while reading
            }//end if skipme==false;
            inread.close();
        }catch (IOException ioe){
            javax.swing.JOptionPane.showMessageDialog(new javax.swing.JFrame(),"IOERROR reading from '"+infile.getAbsolutePath()+"'");
            return null;
        }
        return retvec;
    }//end readfile
    
    //--------------------------------------------------------------------------
    
    static Vector readgcrmafile(File infile,float presentval){
        Vector retvec=new Vector();
        try{
            BufferedReader inread=new BufferedReader(new FileReader(infile));
            String inline,tmp;
            String[] tmparr;
            indata currpoint;
            while((inline=inread.readLine())!=null){
                inline=inline.trim();
                if(inline.length()<1){
                    continue;
                }
                tmparr=inline.split(";",0);
                int arrsize=java.lang.reflect.Array.getLength(tmparr);
                if(arrsize!=2){
                    javax.swing.JOptionPane.showMessageDialog(new javax.swing.JFrame(),"ERROR reading from file "+infile.getAbsolutePath()+" unreadable line '"+inline+"' elements="+arrsize+" not in name;value format?");
                    break;
                }else{
                    currpoint=new indata();
                    currpoint.name=tmparr[0].trim();
                    if((currpoint.name.length()<1)|(currpoint.name.equalsIgnoreCase("blank"))){
                        continue;
                    }
                    try{
                        currpoint.val=Float.parseFloat(tmparr[1].trim());
                    }catch (NumberFormatException ne){
                        javax.swing.JOptionPane.showMessageDialog(new javax.swing.JFrame(),"ERROR unable to parse float from '"+tmparr[1]+"'");
                        break;
                    }
                    currpoint.presentval=presentval;
                    retvec.addElement(currpoint);
                }
            }//end while reading
            inread.close();
        }catch (IOException ioe){
            javax.swing.JOptionPane.showMessageDialog(new javax.swing.JFrame(),"IOERROR reading from '"+infile.getAbsolutePath()+"'");
            return null;
        }
        return retvec;
    }//end readgcrmafile
    
    //--------------------------------
    
    static class indata{
        
        public indata(){}
        
        String name=null;
        double val=-1;
        float presentval=0;
        boolean reference=false;
        
    }//end class indata
    
    //--------------------------------------------------------------------------
    
}
