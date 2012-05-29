/*
 * clusternode.java
 *
 * Created on April 29, 2003, 10:15 AM
 */
package clans;
/**
 *
 * @author  tancred
 */
public class cluster{
    
    /** Creates a new instance of clusternode */
    public cluster() {
    }
    
    String name="no name";
    int[] member=new int[0];
    float clusterconfidence=-1;
    float[] seqconfidence=null;
    
    void add(int newmember){
        int[] oldmember=member;
        int oldsize=java.lang.reflect.Array.getLength(oldmember);
        member=new int[oldsize+1];
        for(int i=0;i<oldsize;i++){
            member[i]=oldmember[i];
            if(newmember==oldmember[i]){//if what I want to add is already present as an element
                member=oldmember;
                return;//don't add anything
            }
        }//end for i
        member[oldsize]=newmember;
    }//end addmember
    
    void add(int[] newmembers){
        int[] oldmember=member;
        int oldsize=java.lang.reflect.Array.getLength(oldmember);
        int newsize=java.lang.reflect.Array.getLength(newmembers);
        member=new int[oldsize+newsize];
        int poscount=0;
        for(int i=0;i<oldsize;i++){
            member[poscount]=oldmember[i];
            poscount++;
        }
        boolean hasmember;
        for(int i=0;i<newsize;i++){
            hasmember=false;
            for(int j=0;j<oldsize;j++){
                if(oldmember[j]==newmembers[i]){
                    hasmember=true;
                    break;//exit for j
                }
            }//end for j
            if(hasmember==false){
                member[poscount]=newmembers[i];
                poscount++;
            }
        }//end for i
        oldmember=member;
        member=new int[poscount];
        for(int i=0;i<poscount;i++){
            member[i]=oldmember[i];
        }
    }//end addmember
    
    //-------------------------------------
    
    void remove(int xmember){
        int[] oldmember=member;
        int oldsize=java.lang.reflect.Array.getLength(oldmember);
        if(oldsize<1){
            return;
        }
        member=new int[oldsize-1];
        int offset=0;
        for(int i=0;i<oldsize;i++){
            if(oldmember[i]==xmember){
                offset=1;
                continue;
            }
            if((i==(oldsize-1))&&(offset==0)){
                //if I didn't find the element to remove
                member=oldmember;
                return;
            }
            member[i-offset]=oldmember[i];
        }//end for i
    }//end remove
    
    //--------------------------------------
    
    int members(){
        return java.lang.reflect.Array.getLength(member);
    }//end members()
    
    //---------------------------------------
    
}//end class
