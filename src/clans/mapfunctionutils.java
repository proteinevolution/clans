/*
 * mapfunctionutils.java
 *
 * Created on June 20, 2006, 3:13 PM
 */
package clans;
import java.io.*;
import java.util.*;
/**
 *
 * @author  tancred
 */
public class mapfunctionutils {
    
    /** Creates a new instance of mapfunctionutils */
    public mapfunctionutils() {
    }
    
    static HashMap loadlookup(File loadfile){
        //read data from a file containing affy_id classification_id lookup
        HashMap lookup=new HashMap();
        System.out.println("reading "+loadfile.getName());
        try{
            BufferedReader inread=new BufferedReader(new InputStreamReader(new javax.swing.ProgressMonitorInputStream(new javax.swing.JFrame(),"Reading "+loadfile.getName(),new FileInputStream(loadfile))));
            //BufferedReader inread=new BufferedReader(new FileReader(loadfile));
            String inline;
            String affyid,lookupid;
            String[] tmp;
            int count=0;
            while((inline=inread.readLine())!=null){
                count++;
                if(count%100==0){
                    System.out.print(".");
                }
                tmp=inline.split("\\s+");
                affyid=tmp[0].toLowerCase();
                lookupid=tmp[1].toLowerCase();
                lookup.put(lookupid,affyid);
                lookup.put(affyid,lookupid);
            }
            inread.close();
        }catch(IOException ioe){
            System.err.println("IOError reading from '"+loadfile.getAbsolutePath()+"'");
            return null;
        }
        System.out.println();
        return lookup;
    }//end loadlookup
    
    //--------------------------------------------------------------------------
    
    static void loadGOlookup(File loadfile, mapnode rootnode){
        System.out.println("reading "+loadfile.getName());
        //read data from a file containing affy_id classification_id lookup
        //in the GO-case I then have to add the classification names to the mapnode bins as "leaves"
        HashMap nodehash=new HashMap();
        //add each of the nodes to a hash to make a faster name-lookup
        addnodestohash(nodehash,rootnode);
        try{
            BufferedReader inread=new BufferedReader(new InputStreamReader(new javax.swing.ProgressMonitorInputStream(new javax.swing.JFrame(),"Reading "+loadfile.getName(),new FileInputStream(loadfile))));
            //BufferedReader inread=new BufferedReader(new FileReader(loadfile));
            String inline;
            String affyid,lookupid;
            String[] tmp;
            int count=0;
            mapnode curr;
            mapnode parent;
            while((inline=inread.readLine())!=null){
                count++;
                if(count%100==0){
                    System.out.print(".");
                }
                tmp=inline.split("\\s+");
                affyid=tmp[0].toLowerCase();
                lookupid=tmp[1].toLowerCase();
                if(nodehash.containsKey(affyid)){//add the other identifier as a mapnode
                    parent=((mapnode)nodehash.get(affyid));
                    curr=new mapnode();
                    curr.level=lookupid;
                    curr.child=null;
                    curr.info=parent.level+";"+parent.info;
                    curr.leaf=lookupid;
                    parent.add(curr);
                }else if(nodehash.containsKey(lookupid)){//add the other identifier as a mapnode
                    parent=((mapnode)nodehash.get(lookupid));
                    curr=new mapnode();
                    curr.level=affyid;
                    curr.child=null;
                    curr.info=parent.level+";"+parent.info;
                    curr.leaf=affyid;
                    parent.add(curr);
                }else{
                    System.err.println("WARNING: Neither of the id's '"+lookupid+"' or '"+affyid+"' could be found");
                }
            }
            inread.close();
        }catch(IOException ioe){
            System.err.println("IOError reading from '"+loadfile.getAbsolutePath()+"'");
            return;
        }
        System.out.println();
    }//end loadlookup
    
    //--------------------------------------------------------------------------
    
    static void addnodestohash(HashMap nodehash, mapnode curr){
        if(curr.child!=null){
            for(int i=java.lang.reflect.Array.getLength(curr.child);--i>=0;){
                addnodestohash(nodehash,curr.child[i]);
            }//end for i
        }
        if(curr.level.indexOf("|")>-1){
            String[] tmp=curr.level.split("\\|");
            for(int i=java.lang.reflect.Array.getLength(tmp);--i>=0;){
                nodehash.put(tmp[i],curr);
            }//end for i
        }else{
            nodehash.put(curr.level,curr);
        }
    }//end addnodestohash
    
    //--------------------------------------------------------------------------
    
    static mapnode loadunknownformat(File loadfile){
        //see if the file is in genebins/mapman format or GO.obo format
        //then use the corresponding load function to load the data
        try{
            BufferedReader inread=new BufferedReader(new FileReader(loadfile));
            String inline="";
            while((inline=inread.readLine())!=null){
                //only read until the first non-empty line is hit, then check the format
                inline=inline.trim();
                if(inline.length()>0){
                    String[] tmp=inline.split("[\"\']*[\t;][\"\']*",4);
                    if(java.lang.reflect.Array.getLength(tmp)>=4){
                        //then I am reading a valid MapMan or GeneBins file
                        inread.close();
                        return loadbins(loadfile);
                    }else{
                        //I suppose I am reading a G0 *.obo file
                        inread.close();
                        return loadGOobo(loadfile);
                    }
                }
            }
        }catch (IOException ioe){
            System.err.println("IOERROR trying to read from '"+loadfile.getAbsolutePath()+"'");
            return null;
        }
        return null;
    }//end loadunknownformat
    
    //--------------------------------------------------------------------------
    
    static mapnode loadGOobo(File loadfile){
        //load from a gene-ontology *.obo file
        //level should identify a hierarchichal node in a tree
        //name should match zero or one sequence id's in the clans map (or the lookup file)
        //description is optional and is used to describe a tree node
        HashMap nodehash=new HashMap();
        ArrayList noparentnodes=new ArrayList();
        try{
            System.out.println("Attempting to read Gene-Ontology *.obo format '"+loadfile.getName()+"'");
            BufferedReader inread=new BufferedReader(new FileReader(loadfile));
            String inline;
            mapnode currnode;
            String id="";//GO-number
            String name="";//name associated to GO-number
            //String annot="";//the definition of the go term
            //String namespace="";//the category this belongs to
            boolean obsolete=true;
            ArrayList parents=new ArrayList();//the parents of this node
            ArrayList altids=new ArrayList();
            int size;
            int count=0;
            mapnode curr,addcurr;
            String pid;
            while((inline=inread.readLine())!=null){
                inline=inline.trim();
                if(inline.equalsIgnoreCase("[Term]")){
                    count++;
                    if(count%1000==0){
                        System.out.print(".");
                    }
                    if(obsolete==false){
                        //then remember the data
                        //to do this, see if the current node is already defined in the hash.
                        if(nodehash.containsKey(id)==false){
                            //I don't need to check the alternate id's as the GO-terms always use the main id's
                            //check the alternate id's
                            //boolean hadalternate=false;
                            //for(int i=altids.size();--i>=0;){
                            //    if(nodehash.containsKey((String)altids.get(i))){
                            //        hadalternate=true;
                            //        curr=(mapnode)nodehash.get((String)altids.get(i));
                            //        nodehash.put(id,curr);
                            //    }
                            //}//end for i
                            //if it is not defined, create a new tmpnode
                            //if(hadalternate==false){
                            curr=new mapnode();
                            curr.child=new mapnode[0];
                            nodehash.put(id,curr);
                            //}else{//if I DID have an alternate name
                            //    curr=(mapnode)nodehash.get(id);
                            //}
                            //and now add the alternate names
                            //for(int i=altids.size();--i>=0;){
                            //    nodehash.put((String)altids.get(i),curr);
                            //}
                        }
                        //now add the data for this node
                        curr=(mapnode)nodehash.get(id);
                        curr.level=id;
                        for(int i=altids.size();--i>=0;){
                            curr.level+="|"+(String)altids.get(i);
                        }//end for i
                        curr.leaf=null;
                        curr.info=name;
                        //now see if this node has a parent assigned
                        if(parents.size()>0){
                            //see if each of the parents is defined
                            for(int i=parents.size();--i>=0;){
                                pid=(String)parents.get(i);
                                if(nodehash.containsKey(pid)==false){
                                    //if it is NOT defined, then define it and add the child
                                    nodehash.put(pid,new mapnode());
                                }
                                //and now add the child
                                ((mapnode)nodehash.get(pid)).add(curr);
                            }//end for i
                        }else{//if it has no parents, add it to the noparentshash
                            noparentnodes.add(curr);
                        }
                    }//end if obsolete==false
                    obsolete=false;
                    id="";
                    name="";
                    //annot="";
                    //namespace="";
                    parents.clear();
                    altids.clear();
                }else if(inline.startsWith("id: ")){
                    id=inline.substring(4).trim().toLowerCase();
                }else if(inline.startsWith("alt_id: ")){
                    altids.add(inline.substring(8).trim().toLowerCase());
                }else if(inline.startsWith("name: ")){
                    name=inline.substring(6).trim();
                }else if(inline.startsWith("is_obsolete: ")){
                    if(inline.indexOf("true")>-1){
                        obsolete=true;
                    }
                }else if(inline.startsWith("is_a: ")){
                    inline=inline.substring(6).trim().toLowerCase();
                    if(inline.indexOf("!")>-1){//get rid of anything except the GO number
                        inline=inline.substring(0,inline.indexOf("!")).trim();
                    }
                    parents.add(inline);
                }else if(inline.startsWith("relationship: part_of ")){
                    inline=inline.substring(22).trim().toLowerCase();
                    if(inline.indexOf("!")>-1){//get rid of anything except the GO number
                        inline=inline.substring(0,inline.indexOf("!")).trim();
                    }
                    parents.add(inline);
                }//else if(inline.startsWith("namespace: ")){
                //    namespace=inline.substring(11).trim();
                //}else if(inline.startsWith("def: ")){
                //    annot=inline.substring(5).trim();
                //}
            }
            inread.close();
        }catch(IOException ioe){
            System.err.println("IOERROR reading from "+loadfile.getName());
            return null;
        }
        System.out.println("finalizing bin structure; may take a while");
        //now I should have all the nodes in a hash more or less in a tree hierarchy.
        //Next, add the nodes without parents to rootnode
        mapnode rootnode=new mapnode();
        rootnode.level="GO-root";
        rootnode.leaf=null;
        for(int i=noparentnodes.size();--i>=0;){
            rootnode.add((mapnode)noparentnodes.get(i));
        }//end for i
        System.out.println("DONE");
        return rootnode;
    }//end loadobo
    
    //--------------------------------------------------------------------------
    
    static mapnode loadbins(File loadfile){
        //load from a file containing level;[description];name;[annot_info];
        //level should identify a hierarchichal node in a tree
        //name should match zero or one sequence id's in the clans map
        //description is optional and is used to describe a tree node
        HashMap nodehash=new HashMap();
        try{
            System.out.println("Attempting to read GeneBins/MapMan file '"+loadfile.getName()+"'");
            BufferedReader inread=new BufferedReader(new FileReader(loadfile));
            String inline;
            mapnode currnode;
            String[] tmp;
            String annot="";
            String idstr="";
            int size;
            int count=0;
            while((inline=inread.readLine())!=null){
                count++;
                if(count%1000==0){
                    System.out.print(".");
                }
                //System.err.println("reading "+inline);
                tmp=inline.split("[\"\']*[\t;][\"\']*",4);
                size=java.lang.reflect.Array.getLength(tmp);
                if(size<2){
                    System.err.println("ERROR parsing on line '"+inline.trim()+"'");
                    return null;
                }else if(size==2){
                    //i.e. I only have bincode and name
                    String[] tmptmp=new String[4];
                    tmptmp[0]=tmp[0];
                    tmptmp[1]=tmp[1];
                    tmptmp[2]="";
                    tmptmp[3]="";
                    tmp=tmptmp;
                }else if(size==3){
                    //i.e. I have bincode and name and identifiers, but no descriptions
                    String[] tmptmp=new String[4];
                    tmptmp[0]=tmp[0];
                    tmptmp[1]=tmp[1];
                    tmptmp[2]=tmp[2];
                    tmptmp[3]="";
                    tmp=tmptmp;
                }//els I have 4 or more elements; all that I need
                //tmp[0] is the level/node assignment (i.e. 3.2.14)
                tmp[0]=tmp[0].trim();
                //I might still have a string delimiter at the first position, check for that
                if(tmp[0].startsWith("'") || tmp[0].startsWith("\"")){
                    tmp[0]=tmp[0].substring(1);
                }
                tmp[1]=tmp[1].trim();//tmp[1] is the name
                if(tmp[1].length()==0){
                    tmp[1]="undefined";
                }
                if(nodehash.containsKey(tmp[0])==false){
                    //System.err.println("adding basenode "+tmp[0]);
                    //then create that node and add the annotation info for it
                    mapnode basenode=new mapnode();
                    basenode.level=tmp[0];
                    //tmp[1] is the level description (i.e. sugar metabolism)
                    basenode.info=tmp[1];
                    basenode.leaf=null;
                    nodehash.put(tmp[0],basenode);
                }
                //now add the current leaf-node
                //tmp[2] is the sequence identifier (i.e. at-number)
                tmp[2]=tmp[2].trim().toLowerCase();
                if(tmp[2].length()>0){//if this is a leaf node
                    idstr=tmp[0]+"."+tmp[2];
                    //System.err.println("adding newnode "+idstr);
                    mapnode newnode=new mapnode();
                    if(nodehash.containsKey(idstr)){
                        System.err.println(idstr+" already defined!");
                        continue;
                    }
                    newnode.level=idstr;
                    newnode.leaf=tmp[2];
                    //tmp[3] is annotation information for the sequence identifier
                    tmp[3]=tmp[3].trim();
                    //here too, I might still have a string delimiter at the end
                    if(tmp[3].endsWith("'") || tmp[3].endsWith("\"")){
                        tmp[3]=tmp[3].substring(0,tmp[3].length()-1);
                    }
                    if(tmp[3].length()>0){
                        newnode.info=tmp[3];
                    }
                    nodehash.put(idstr,newnode);
                }
                
            }
            inread.close();
        }catch(IOException ioe){
            System.err.println("IOERROR reading from "+loadfile.getName());
            return null;
        }
        System.out.println("finalizing, may take a while");
        //now I should have all the nodes in a hash; next, create the tree hierarchy
        String[] keys=(String[])nodehash.keySet().toArray(new String[0]);
        int num=java.lang.reflect.Array.getLength(keys);
        String[] tmp;
        String currlevel;
        mapnode lastnode;
        mapnode newnode;
        for(int i=num;--i>=0;){
            tmp=keys[i].split("\\.");
            //now create the level (string) array
            for(int j=1;j<java.lang.reflect.Array.getLength(tmp);j++){
                tmp[j]=tmp[j-1]+"."+tmp[j];
            }//end for j
            lastnode=(mapnode)nodehash.get(keys[i]);
            //now go through the levels and add nodes where necessary
            for(int j=java.lang.reflect.Array.getLength(tmp)-1;--j>=0;){//not the last element, as that is the node I am looking at!
                if(nodehash.containsKey(tmp[j])){
                    newnode=(mapnode)nodehash.get(tmp[j]);
                }else{
                    newnode=new mapnode();
                    newnode.level=tmp[j];
                }
                if(newnode.haschild(lastnode)==-1){
                    newnode.add(lastnode);
                    newnode.leaf=null;
                }
                nodehash.put(tmp[j],newnode);
                lastnode=newnode;
            }//end for j
        }//end for i
        //now I should have all the nodes (and the childnode assignments) in a hash
        //now go through all of the nodes and add the info of the parent node to the leaf-node if it's info is "undefined"
        keys=(String[])nodehash.keySet().toArray(new String[0]);
        mapnode currnode;
        String parentname;
        int lastdotindex;
        for(int i=java.lang.reflect.Array.getLength(keys);--i>=0;){
            currnode=(mapnode)nodehash.get(keys[i]);
            if(currnode.info.equalsIgnoreCase("undefined")&&currnode.leaf!=null){
                //then assign it the parent info
                if((lastdotindex=currnode.level.lastIndexOf("."))>-1){
                    parentname=currnode.level.substring(0,lastdotindex);
                    if(nodehash.containsKey(parentname)){
                        currnode.info=((mapnode)nodehash.get(parentname)).info;
                    }
                }
            }
        }//end for i
        //next, sort the childnodes, add the most basal nodes to "root", add root to the hash and return the hash
        keys=(String[])nodehash.keySet().toArray(new String[0]);
        mapnode rootnode=new mapnode();
        rootnode.level="root";
        rootnode.leaf=null;
        for(int i=java.lang.reflect.Array.getLength(keys);--i>=0;){
            if(keys[i].indexOf(".")==-1){//i.e. if this key has no dot in it
                rootnode.add((mapnode)nodehash.get(keys[i]));
            }
        }//end for i
        //and now sort the mapnodes
        sortnode(rootnode);
        //printnode(rootnode,"");
        System.out.println("DONE");
        return rootnode;
    }//end loadbins
    
    //--------------------------------------------------------------------------
    
    static void sortnode(mapnode innode){
        //sort the childnodes of this node
        if(innode.child==null){
            return;
        }
        for(int i=java.lang.reflect.Array.getLength(innode.child);--i>=0;){
            sortnode(innode.child[i]);
        }//end for i
        java.util.Arrays.sort(innode.child,new childnodecomparator());
    }//end sortnode
    
    //---------------------
    
    static class childnodecomparator implements java.util.Comparator{
        
        public int compare(Object o1, Object o2){
            String[] s1=((mapnode)o1).level.split("\\.");
            String[] s2=((mapnode)o2).level.split("\\.");
            int inum=java.lang.reflect.Array.getLength(s1);
            int jnum=java.lang.reflect.Array.getLength(s2);
            int lim=inum;
            if(lim>jnum){
                lim=jnum;
            }
            int i1,i2;
            try{
                for(int i=0;i<lim;i++){
                    i1=Integer.parseInt(s1[i]);
                    i2=Integer.parseInt(s2[i]);
                    if(i1>i2){
                        return -1;
                    }else if(i1<i2){
                        return 1;
                    }//else, check the next element
                }//end for i
            }catch (NumberFormatException ne){
                //i.e. they are not both only integers
                //then simply compare the last element of the array as string
                return -1*(s1[inum-1].compareTo(s2[jnum-1]));
            }
            //if they were the same up to here
            if(inum>jnum){//i.e. o1 is child of o2
                return 1;
            }else if(jnum>inum){
                return -1;
            }else{
                return 0;
            }
        }
        
    }//end comparator
    
    //--------------------------------------------------------------------------
    
    static void printnode(mapnode innode, String spacer){
        System.out.println(spacer+"'"+innode.level+"' "+innode.info);
        if(innode.leaf==null){
            for(int i=java.lang.reflect.Array.getLength(innode.child);--i>=0;){
                printnode(innode.child[i],spacer+"  ");
            }//end for i
        }
    }//end printnode
    
}
