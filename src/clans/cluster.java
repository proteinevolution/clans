package clans;

public class cluster {
    
    String name="no name";
    int[] member=new int[0];
    float clusterconfidence=-1;
    float[] seqconfidence=null;
    
    /**
     * 
     * @param newmember
     */
    void add(int newmember){
        int[] oldmember=member;
        int oldsize=oldmember.length;
        member=new int[oldsize+1];
        
        for(int i=0;i<oldsize;i++){
            member[i]=oldmember[i];
            
            if(newmember==oldmember[i]){//if what I want to add is already present as an element
                member=oldmember;
                return;
            }
        }
        
        member[oldsize]=newmember;
    }
    
    /**
     * 
     * @param newmembers
     */
    void add(int[] newmembers){
        int[] oldmember=member;
        int oldsize=oldmember.length;
        int newsize=newmembers.length;
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
    }
    
    /**
     * 
     * @param xmember
     */
    void remove(int xmember){
        int[] oldmember=member;
        int oldsize=oldmember.length;
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
        }
    }
}