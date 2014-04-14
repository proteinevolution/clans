package clans;

import java.util.HashMap;

public class MovementComputerThread extends java.lang.Thread {

	/**
	 * Crates a MovementComputerThread with name {@code MovementComputerThread#<thread_id>} that only optimizes the
	 * selected subset of sequences.
	 * 
	 * @param data
	 *            The data
	 * @param thread_id
	 *            The numerical identifier for this thread.
	 * @param informParentOfCompletenessLock
	 *            The lock on which parents get informed about completed computations in a synchronized manner.
	 * @param selectedNamesHash
	 * @param selectedNamesIndices
	 */
	public MovementComputerThread(ClusterData data, int thread_id, final Object informParentOfCompletenessLock,
			HashMap<String, Integer> selectedNamesHash) {

		this(data, thread_id, informParentOfCompletenessLock);

		this.selectedNamesHash = selectedNamesHash;
		this.selectedNamesIndices = data.selectedSequencesIndicesStableCopy;
	}

	/**
	 * Crates a MovementComputerThread with name {@code MovementComputerThread#<thread_id>}.
	 * 
	 * @param data
	 *            The data
	 * @param thread_id
	 *            The numerical identifier for this thread.
	 * @param informParentOfCompletenessLock
	 *            The lock on which parents get informed about completed computations in a synchronized manner.
	 */
	public MovementComputerThread(ClusterData data, int thread_id, final Object informParentOfCompletenessLock) {

		this.setName("MovementComputerThread#" + thread_id);
		this.threadId = thread_id;
		this.totalThreads = data.cpu;

		this.data = data;
		this.informParentOfCompletenessLock = informParentOfCompletenessLock;

		this.positions = data.positions;
		this.attractions = data.attractionValues;
		this.movements = data.movements;
	}

	private int threadId;
	private int totalThreads;

	private ClusterData data;
	private final Object informParentOfCompletenessLock;

	private float[][] positions;
	private MinimalAttractionValue[] attractions;
	private float[][] movements;

	private HashMap<String, Integer> selectedNamesHash = null;
	private int[] selectedNamesIndices = null; // a local copy of parent.selectednames, as that may change during a calculation

	/**
	 * Notifies the parent of completed computations.
	 * <p>
	 * Note: Attempts to remove this and instead check .isAlive() of this thread in the parent have failed for unknown
	 * reasons. isAlive() seems to often return true if checked immediately after executing the last line of the run()
	 * method, which obviously can cause problems.. Unless you have too much time or know your Java threads very well,
	 * don't try again ;)
	 */
	protected boolean done = false;

	private boolean doSelectedOnly() {
		return selectedNamesHash != null && selectedNamesIndices != null;
	}
    
    /**
     * use the positions of all elements and their attraction/repulsion values to calculate a movement vector for each
     * (take into account the last movement). repulsion doesn't have a specific value as all evalues below a certain
     * point are simply regarded as insignificant. therefore use a one formula for all to compute the repulsive forces
     */
    @Override
    public void run() {

        int elements = positions.length;
        double[] currmoverep = new double[3];
        double[] currmoveatt = new double[3];
        double totaldist = 0;

        int attnum = attractions.length;
        double minattract = data.minattract;
        int repvalpow = data.repvalpow;
        int attvalpow = data.attvalpow;
        float repfactor = data.repfactor;
        float attfactor = data.attfactor;
        float maxmove = data.maxmove;

        if (doSelectedOnly()) {

            //now get from where to where I should do my calculations
        	int selected_sequences_count = selectedNamesIndices.length;
            int start = threadId * (int) (selected_sequences_count / totalThreads);
            int end;
            if (threadId == (totalThreads - 1)) {
                end = selected_sequences_count;
            } else {
                end = (threadId + 1) * (int) (selected_sequences_count / totalThreads);
            }
            
            if (data.cluster2d == true) {

                //cluster only the selected sequences in 2D
                for (int i = start; i < end; i++) {
                    ClusterMethods.addAttractionTowardsOrigin(positions[selectedNamesIndices[i]], currmoveatt, minattract);
                    movements[selectedNamesIndices[i]][0] += currmoveatt[0];
                    movements[selectedNamesIndices[i]][1] += currmoveatt[1];

                    for (int j = elements; --j >= 0;) {
                        ClusterMethods.addRepulsion(positions[selectedNamesIndices[i]], positions[j], currmoverep, repvalpow,
                                repfactor, ClusterMethods.rand, data.cluster2d);
                        movements[selectedNamesIndices[i]][0] += currmoverep[0];
                        movements[selectedNamesIndices[i]][1] += currmoverep[1];
                    }
                }
                
                // now add the attraction values, but only for the query or hit sequences for which I am responsible
                for (int i = attnum; --i >= 0;) {
                    if ((selectedNamesHash.containsKey(String.valueOf(attractions[i].query))) && (selectedNamesHash.get(String.valueOf(attractions[i].query)).intValue() == threadId)) {
                        ClusterMethods.addAttraction(positions[attractions[i].query], positions[attractions[i].hit], attractions[i].att,
                                currmoveatt, attvalpow, attfactor, data.cluster2d);

                        movements[attractions[i].query][0] += currmoveatt[0];
                        movements[attractions[i].query][1] += currmoveatt[1];
                        //movement[attvals[i].query][2]+=currmoveatt[2];
                        if ((selectedNamesHash.containsKey(String.valueOf(attractions[i].hit))) && (selectedNamesHash.get(String.valueOf(attractions[i].hit)).intValue() == threadId)) {
                            movements[attractions[i].hit][0] -= currmoveatt[0];
                            movements[attractions[i].hit][1] -= currmoveatt[1];
                            //movement[attvals[i].hit][2]-=currmoveatt[2];
                        }
                    } else if ((selectedNamesHash.containsKey(String.valueOf(attractions[i].hit))) && (selectedNamesHash.get(String.valueOf(attractions[i].hit)).intValue() == threadId)) {
                        ClusterMethods.addAttraction(positions[attractions[i].query], positions[attractions[i].hit], attractions[i].att,
                                currmoveatt, attvalpow, attfactor, data.cluster2d);

                        movements[attractions[i].hit][0] -= currmoveatt[0];
                        movements[attractions[i].hit][1] -= currmoveatt[1];
                        //movement[attvals[i].hit][2]-=currmoveatt[2];
                    }
                }//end for i
                //double totaldistsq;
                for (int i = start; i < end; i++) {
                    movements[selectedNamesIndices[i]][0] /= elements;
                    movements[selectedNamesIndices[i]][1] /= elements;
                    movements[selectedNamesIndices[i]][2] = 0;
                    totaldist = java.lang.Math.sqrt((movements[selectedNamesIndices[i]][0] * movements[selectedNamesIndices[i]][0]) + (movements[selectedNamesIndices[i]][1] * movements[selectedNamesIndices[i]][1]));
                    if (totaldist > maxmove) {
                        movements[selectedNamesIndices[i]][0] *= maxmove / totaldist;
                        movements[selectedNamesIndices[i]][1] *= maxmove / totaldist;
                    }
                }//end for i
            } else {
                
                //cluster only the selected sequences in 3D
                for (int i = start; i < end; i++) {

                    ClusterMethods.addAttractionTowardsOrigin(positions[selectedNamesIndices[i]], currmoveatt, minattract);
                    movements[selectedNamesIndices[i]][0] += currmoveatt[0];
                    movements[selectedNamesIndices[i]][1] += currmoveatt[1];
                    movements[selectedNamesIndices[i]][2] += currmoveatt[2];
                    for (int j = elements; --j >= 0;) {
                        ClusterMethods.addRepulsion(positions[selectedNamesIndices[i]], positions[j], currmoverep, repvalpow, repfactor, ClusterMethods.rand, data.cluster2d);

                        movements[selectedNamesIndices[i]][0] += currmoverep[0];
                        movements[selectedNamesIndices[i]][1] += currmoverep[1];
                        movements[selectedNamesIndices[i]][2] += currmoverep[2];
                    }//end for j
                }//end for i
                
                // now add the attraction values, but only for the query or hit sequences for which I am responsible
                for (int i = attnum; --i >= 0;) {
                    if ((selectedNamesHash.containsKey(String.valueOf(attractions[i].query))) && (selectedNamesHash.get(String.valueOf(attractions[i].query)).intValue() == threadId)) {
                        ClusterMethods.addAttraction(positions[attractions[i].query], positions[attractions[i].hit], attractions[i].att,
                                currmoveatt, attvalpow, attfactor, data.cluster2d);

                        movements[attractions[i].query][0] += currmoveatt[0];
                        movements[attractions[i].query][1] += currmoveatt[1];
                        movements[attractions[i].query][2] += currmoveatt[2];
                        if ((selectedNamesHash.containsKey(String.valueOf(attractions[i].hit))) && (selectedNamesHash.get(String.valueOf(attractions[i].hit)).intValue() == threadId)) {
                            movements[attractions[i].hit][0] -= currmoveatt[0];
                            movements[attractions[i].hit][1] -= currmoveatt[1];
                            movements[attractions[i].hit][2] -= currmoveatt[2];
                        }
                    } else if ((selectedNamesHash.containsKey(String.valueOf(attractions[i].hit))) && (selectedNamesHash.get(String.valueOf(attractions[i].hit)).intValue() == threadId)) {
                        ClusterMethods.addAttraction(positions[attractions[i].query], positions[attractions[i].hit], attractions[i].att,
                                currmoveatt, attvalpow, attfactor, data.cluster2d);
                        
                        movements[attractions[i].hit][0] -= currmoveatt[0];
                        movements[attractions[i].hit][1] -= currmoveatt[1];
                        movements[attractions[i].hit][2] -= currmoveatt[2];
                    }
                }//end for i
                //double totaldistsq;
                for (int i = start; i < end; i++) {
                    movements[selectedNamesIndices[i]][0] /= elements;
                    movements[selectedNamesIndices[i]][1] /= elements;
                    movements[selectedNamesIndices[i]][2] /= elements;
                    totaldist = java.lang.Math.sqrt((movements[selectedNamesIndices[i]][0] * movements[selectedNamesIndices[i]][0]) + (movements[selectedNamesIndices[i]][1] * movements[selectedNamesIndices[i]][1]) + (movements[selectedNamesIndices[i]][2] * movements[selectedNamesIndices[i]][2]));
                    if (totaldist > maxmove) {
                        movements[selectedNamesIndices[i]][0] *= maxmove / totaldist;
                        movements[selectedNamesIndices[i]][1] *= maxmove / totaldist;
                        movements[selectedNamesIndices[i]][2] *= maxmove / totaldist;
                    }
                }//end for i
            }
        } else { // if no sequences were selected or all should be used

            // get start and end index of data concerning this thread
            int start, end;
            if (threadId == (totalThreads - 1)) {
                // do everything from here to the last element to avoid rounding errors
                start = threadId * (int) (elements / totalThreads);
                end = elements;
            } else {
                start = threadId * (int) (elements / totalThreads);
                end = (threadId + 1) * (int) (elements / totalThreads);
            }
            
            // add repulsive movements
            for (int i = start; i < end; i++) {
                ClusterMethods.addAttractionTowardsOrigin(positions[i], currmoveatt, minattract);
                
                for (int j = 0; j < movements[i].length; j ++) {
                    movements[i][j] += currmoveatt[j];
                }

                for (int j = 0; j < elements; j++) {
                    if (i == j) {
                        continue;
                    } 

                    ClusterMethods.addRepulsion(positions[i], positions[j], currmoverep, repvalpow, repfactor, ClusterMethods.rand, data.cluster2d);
                    
                    for (int k = 0; k < movements[i].length; k ++) {
                        movements[i][k] += currmoverep[k];
                    }
                }
            }
            
            // add attractive movements
            for (int i = attnum; --i >= 0;) {
                if (attractions[i].query >= start && attractions[i].query < end) {
                    ClusterMethods.addAttraction(positions[attractions[i].query], positions[attractions[i].hit], attractions[i].att,
                            currmoveatt, attvalpow, attfactor, data.cluster2d);
                    
                    for (int j = 0; j < movements[attractions[i].query].length; j ++) { 
                        movements[attractions[i].query][j] += currmoveatt[j];
                    }
                }
                
                if (attractions[i].hit >= start && attractions[i].hit < end) {
                    ClusterMethods.addAttraction(positions[attractions[i].hit], positions[attractions[i].query], attractions[i].att,
                            currmoveatt, attvalpow, attfactor, data.cluster2d);

                    for (int j = 0; j < movements[attractions[i].hit].length; j ++) {
                        movements[attractions[i].hit][j] += currmoveatt[j];
                    }
                }
            }

            // restrict to the maximal allowed move distance
            for (int i = start; i < end; i++) {
                for (int j = 0; j < movements[i].length; j ++) {
                    movements[i][j] /= elements;
                }
                
                totaldist = java.lang.Math.sqrt((movements[i][0] * movements[i][0]) + (movements[i][1] * movements[i][1])
                        + (movements[i][2] * movements[i][2]));
                
                if (totaldist > maxmove) {
                    for (int j = 0; j < movements[i].length; j ++) {
                        movements[i][j] *= maxmove / totaldist;
                    }
                }
            }
        }

        synchronized (informParentOfCompletenessLock) {
            this.done = true;
            informParentOfCompletenessLock.notify();
        }
    }
}

