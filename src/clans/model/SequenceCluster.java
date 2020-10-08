package clans.model;

public class SequenceCluster {
    
	public String name="no name";
	public int[] member=new int[0];
	public float clusterconfidence=-1;
	public float[] seqconfidence=null;

    /**
     * 
     * @param newMember
     */
    public void add(int newMember) {
        int[] oldMember = member;
        int oldSize = oldMember.length;
        member = new int[oldSize + 1];
        
        for(int i = 0; i < oldSize; i++) {
            member[i] = oldMember[i];
            
            if(newMember == oldMember[i]) {//if what I want to add is already present as an element
                member = oldMember;
                return;
            }
        }
        
        member[oldSize] = newMember;
    }
    
    /**
     * 
     * @param newmembers
     */
    public void add(int[] newMembers) {
        int[] oldMember = member;
        int oldSize = oldMember.length;
        int newSize = newMembers.length;
        member = new int[oldSize + newSize];

        int posCount=0;
        for(int i = 0;i < oldSize; i++){
            member[posCount] = oldMember[i];
            posCount++;
        }

        boolean hasMember;
        for(int i = 0; i < newSize; i++) {
            hasMember = false;
            for(int j = 0; j < oldSize; j++) {
                if(oldMember[j] == newMembers[i]) {
                    hasMember = true;
                    break;//exit for j
                }
            }//end for j

            if(!hasMember){
                member[posCount] = newMembers[i];
                posCount++;
            }
        }//end for i

        oldMember = member;
        member = new int[posCount];
        for(int i = 0; i < posCount; i++){
            member[i] = oldMember[i];
        }
    }
    
    /**
     * 
     * @param xMember
     */
    void remove(int xMember) {
        int[] oldMember = member;
        int oldSize = oldMember.length;
        if(oldSize < 1) {
            return;
        }

        member = new int[oldSize - 1];
        int offset = 0;
        for(int i = 0; i < oldSize; i++) {
            if(oldMember[i] == xMember) {
                offset=1;
                continue;
            }

            if((i == (oldSize - 1)) && (offset == 0)) {
                //if I didn't find the element to remove
                member = oldMember;
                return;
            }
            member[i - offset] = oldMember[i];
        }
    }
}
