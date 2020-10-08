package clans.model;

import java.util.*;

public class SequenceCluster {
    
	public String name = "no name";
	public int[] members = new int[0];
	public float clusterConfidence = -1;
	public float[] seqConfidence = null;

	public SequenceCluster()
	{
	}

	/**
	 * 
	 * @param newMember
	 */
	public SequenceCluster(Integer newMember)
	{
		this.members = new int[1];
		this.members[0] = newMember.intValue();
	}

	/**
	 * 
	 * @param newMembers
	 */
	public SequenceCluster(ArrayList<Integer> newMembers)
	{
		this.members = new int[newMembers.size()];
		for (int i = 0; i < newMembers.size(); i++) {
			this.members[i] = newMembers.get(i).intValue();
		}
	}

	/**
	 * 
	 * @param newMembers
	 */
	public SequenceCluster(int[] newMembers)
	{
		this.members = newMembers;
	}

    /**
     * 
     * @param newMember
     */
    public void add(int newMember) {
        int[] oldMembers = members;
        int oldSize = oldMembers.length;
        members = new int[oldSize + 1];
        
        for(int i = 0; i < oldSize; i++) {
            members[i] = oldMembers[i];
            
            if(newMember == oldMembers[i]) {//if what I want to add is already present as an element
                members = oldMembers;
                return;
            }
        }
        
        members[oldSize] = newMember;
    }
    
    /**
     * 
     * @param newMembers
     */
    public void add(int[] newMembers) {
        int[] oldMembers = members;
        int oldSize = oldMembers.length;
        int newSize = newMembers.length;
        members = new int[oldSize + newSize];

        int posCount=0;
        for(int i = 0;i < oldSize; i++){
            members[posCount] = oldMembers[i];
            posCount++;
        }

        boolean hasMember;
        for(int i = 0; i < newSize; i++) {
            hasMember = false;
            for(int j = 0; j < oldSize; j++) {
                if(oldMembers[j] == newMembers[i]) {
                    hasMember = true;
                    break;//exit for j
                }
            }//end for j

            if(!hasMember){
                members[posCount] = newMembers[i];
                posCount++;
            }
        }//end for i

        oldMembers = members;
        members = new int[posCount];
        for(int i = 0; i < posCount; i++){
            members[i] = oldMembers[i];
        }
    }
    
    /**
     * 
     * @param xMember
     */
    void remove(int xMember) {
        int[] oldMembers = members;
        int oldSize = oldMembers.length;
        if(oldSize < 1) {
            return;
        }

        members = new int[oldSize - 1];
        int offset = 0;
        for(int i = 0; i < oldSize; i++) {
            if(oldMembers[i] == xMember) {
                offset=1;
                continue;
            }

            if((i == (oldSize - 1)) && (offset == 0)) {
                //if I didn't find the element to remove
                members = oldMembers;
                return;
            }
            members[i - offset] = oldMembers[i];
        }
    }
}
