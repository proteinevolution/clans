package clans.algorithms;

import java.util.*;

import clans.model.SequenceCluster;
import clans.model.proteins.MinimalAttractionValue;

public class ClusterDetection {

    /**
     * linkage clustering - form clusters with at least each seq connected by minlinks to the cluster
     * 
     * @param attvals
     * @param minlinks
     * @param minseqnum
     * @param seqnum
     * @return
     */
    public static Vector<SequenceCluster> getlinkage(MinimalAttractionValue[] attvals, int minlinks, int minseqnum, int seqnum) {
		Vector<SequenceCluster> retvec = new Vector<SequenceCluster>();
		// simple cases, 0 and 1
		if (minlinks < 1) {// return all sequences
			int[] newClusterMembers = new int[seqnum];
			for (int i = 0; i < seqnum; i++) {
				newClusterMembers[i] = i;
			}// end for i
			SequenceCluster newcluster = new SequenceCluster(newClusterMembers);
			retvec.add(newcluster);
		} else if (minlinks == 1) {// single linkage
			singlelinkage(attvals, retvec, minseqnum, seqnum);
			// now sort the return vector
			// cluster is sorted in singlelinkage method
		} else {// otherwise it gets a bit more difficult
			multilinkage(attvals, retvec, minseqnum, minlinks, seqnum);
			// now sort the vector
			retvec.sort(new ClusterSizeComparator());
		}
		return retvec;
	}// end getlinkage

	// --------------------------------------------------------------------------
	static void multilinkage(MinimalAttractionValue[] attvals, Vector<SequenceCluster> retvec,
			int minseqnum, int minlinks, int seqnum) {
		// get clusters of linkage N
		// quite a bit more complex than single linkage
		/*
		 * take the sequence with the highest number of connections.from those
		 * connected to sequence0 select the N-1 with highest numbers of
		 * connections.see which sequences they ALL find. If none, change in
		 * turn all seqs selected in step 2.once a cluster "seed" of N-sequences
		 * is present, just add sequences connected with more than Ncontinue
		 * until no 1st step sequences have more than "N" connections
		 */

		// No HashMap needed here turns out that clusterhash is used as a
		// two-dimensional array. In this case also have such an array.
		// This also saves memory, since we can fill it with primitives.
		int[][] clusterhash = new int[seqnum][3];
		for (int i = 0; i < seqnum; i++) {
			clusterhash[i][0] = -1;	// groupvals; is a sequence assigned to a cluster and
									// if so, which cluster
			// groupval: -3=used as seed2 in last round, -2=don't use as seed,
			// -1=unassigned, 0=added to this cluster; >0=number of last cluster
			// this seq was assigned to
			clusterhash[i][1] = 0;	// used as seqshits (flag whether or not this sequence
									// is relevant in this calculation
			clusterhash[i][2] = 0;	// how many connections each sequence has to all
									// others
		}// end for i
		getmaxnum(attvals, clusterhash, 2);// get the maximum number
														boolean done = false;
		int[] tmpseqs = new int[minlinks];
		int[] tmpseqs2;
		int counter = 1;// because 0 and -1 are used for sth. else
		SequenceCluster newcluster;
		int maxconn = 0;// counter for the max number of connections
		int maxnum = -1;// the sequence number with highest maxconn
		GETCLUSTER: while (done == false) {
			// first reset all the unassigned cluster sequences
			maxconn = 0;
			for (int i = 0; i < seqnum; i++) {
				if (clusterhash[i][0] == 0
						|| clusterhash[i][0] == -3) {
					clusterhash[i][0] = -1;
				}
                if (clusterhash[i][0] == -1
                        && clusterhash[i][2] > maxconn) { // if
                                                              // this
                                                              // sequence
                                                              // is
                                                              // unassigned
                                                              // and
                                                              // has
                                                              // max
                                                              // connections
					maxconn = clusterhash[i][2];
					maxnum = i;
					// System.out.println("maxconn="+maxconn+" new seed="+maxnum);

				}
			}// end for i
				// now see which sequences were connected to maxnum
				// System.out.println("seed="+maxnum+" maxconn="+maxconn);
			System.out.print(".");
			int[] conarr = getconnecteds(maxnum, seqnum, attvals);
			if (conarr.length < minlinks) {	// if I
											// found
											// no
											// new
											// sequence
											// with
											// at
											// least
											// minlinks
											// connections
											// (no
											// new
											// clusters
											// possible)
				done = true;
				break;
			}
			for (int i = 0; i < seqnum; i++) {
				clusterhash[i][1] = 0;	// seqshits[i]=0;//the
											// number of
											// links
											// connecting
											// sequences to
											// this group
			}// end for i
			for (int i = conarr.length - 1; i >= 0; i--) {
				// remember which sequences had hits with current seed
				clusterhash[conarr[i]][1] = 1;// seqshits[maxnumarr[i]]=1;
			}
			clusterhash[maxnum][0] = 1;// counter;//groupvals[maxnum]=counter;
			newcluster = new SequenceCluster();
			while (newcluster.members.length == 0) {
				// get the minlinks other sequences that share the most hits
				// with this sequence
				tmpseqs = getmaxshared(attvals, clusterhash, seqnum, maxnum, minlinks);
				// now see that I found at least minlinks valid connected
				// sequences
				for (int i = 0; i <= minlinks; i++) {
					if (tmpseqs[i] == -1) {
						// remove the current seed from analysis until end
						// (where it might come back in)
						clusterhash[maxnum][0] = -2;	// skip
														// this
														// seq
														// for
														// the
														// rest
														// of
														// the
														// analysis
						continue GETCLUSTER;
					}
				}// end for i
				for (int i = tmpseqs.length - 1; i >= 0; i--) {
					clusterhash[tmpseqs[i]][0] = -3;	// groupvals[tmpseqs[i]]=-3;//am
														// checking,
														// remember
														// it
														// is
														// used
														// as
														// seed2
				}// end for i
					// now I have the N most connected sequences sharing the
					// most connections and directly connected
					// next, see if I can start a cluster (i.e. if there are any
					// sequence found by all of them)

				tmpseqs2 = getfound(attvals, clusterhash, seqnum, tmpseqs,
						minlinks);// get all sequences found by all of the seqs
									// in tmpseqs
				// System.out.println("hitseqs="+tmpseqs2.length);
				for (int i = tmpseqs2.length - 1; i >= 0; i--) {
					newcluster.add(tmpseqs2[i]);
					// System.out.println("\tseq:"+tmpseqs2[i]);
					if (clusterhash[tmpseqs2[i]][0] != -2) {
						clusterhash[tmpseqs2[i]][0] = 0;// groupvals[i]=0;//added
																				// to
																				// cluster
					}
				}// end for i
					// now see if I have members>minlinks
				int foundnew = 0;
				int newmembers;
                while (((newmembers = newcluster.members.length) < minlinks)
                        && (newmembers != foundnew)) {
					// System.out.println("inwhile");
					foundnew = newmembers;
					// System.out.println("in while");
					// if I don't, I need to re-check, including the new members
					// until I do have more than minlinks
					tmpseqs2 = getfound(attvals, clusterhash, seqnum,
							newcluster.members, minlinks);
					for (int i = tmpseqs2.length - 1; i >= 0; i--) {
						if (clusterhash[tmpseqs2[i]][0] != 0) {
							newcluster.add(tmpseqs2[i]);
							if (clusterhash[tmpseqs2[i]][0] != -2) {
								clusterhash[tmpseqs2[i]][0] = 0;	// groupvals[i]=0;//added
																	// to
																	// cluster
							}
						}
					}// end for i
				}// end while newcluster.memberss<minlinks
					// System.out.println("newcluster.memberss="+newcluster.members.length);
				if (newcluster.members.length < minlinks) {
					// System.out.println("setting "+maxnum+" to -2");
					// if my currnet cluster has less than minlinks elements
					clusterhash[maxnum][0] = -2;	// skip
													// this
													// seq
													// for
													// the
													// rest
													// of
													// the
													// analysis
					continue GETCLUSTER;
				}
			}// end while newcluster.memberss==0
				// now I have a new cluster with at least minlinks members
				// next, see which sequences they find with at least minlinks
				// connections
			int lastmembers = 0;
			int newmembers;

            while ((newmembers = newcluster.members.length) != lastmembers) {
				lastmembers = newmembers;

				// while I find sequences to add
				tmpseqs2 = getfound(attvals, clusterhash, seqnum,
						newcluster.members, minlinks);
				for (int i = tmpseqs2.length - 1; i >= 0; i--) {
					if (clusterhash[tmpseqs2[i]][0] != 0) {
						newcluster.add(tmpseqs2[i]);
						clusterhash[tmpseqs2[i]][0] = 0;	// groupvals[i]=0;//added
																// to
																// cluster
					}
				}// end for i
			}// end while
				// now I have all sequences I can add in this way
				// set the groups info
			for (int i = newcluster.members.length - 1; i >= 0; i--) {
				// if(clusterhash[newcluster.members[i]][0]==0){
				clusterhash[newcluster.members[i]][0] = counter;
				// }
			}// end for i
			counter++;
			if (newcluster.members.length > minseqnum) {
				// System.out.println("adding cluster "+counter);
				retvec.addElement(newcluster);
			}
		}// end while done==false
		return;
	}// end multilinkage

	// --------------------------------------------------------------------------
	static int[] getfound(MinimalAttractionValue[] attvals, int[][] clusterhash,
			 int seqnum, int[] tmpseqs, int minlinks) {
		// get all the sequences found by all of the tmpseqs
		int attnum = attvals.length;
		HashSet<Integer> tmphash = new HashSet<Integer>();
		for (int i = tmpseqs.length - 1; i >= 0; i--) {
			tmphash.add(tmpseqs[i]);
		}// end for i
			// initialize to zero
		for (int i = 0; i < seqnum; i++) {
			clusterhash[i][1] = 0;// seqshits[i]=0;
		}// end for i
			// now see which were found by the above
		for (int i = 0; i < attnum; i++) {
			if (tmphash.contains(attvals[i].query)) {
				clusterhash[attvals[i].hit][1]++;
			}
			if (tmphash.contains(attvals[i].hit)) {
				clusterhash[attvals[i].query][1]++;
			}
		}// end for i
		Vector<Integer> tmpvec = new Vector<Integer>();
		for (int i = 0; i < seqnum; i++) {
			if (clusterhash[i][1] >= minlinks) {
				tmpvec.addElement(new Integer(i));
			}
		}// end for i
		int[] retarr = new int[tmpvec.size()];
		for (int i = tmpvec.size() - 1; i >= 0; i--) {
			retarr[i] = tmpvec.elementAt(i).intValue();
		}// end for i
		return retarr;
	}

    /**
     * get the minlinks-1 other sequences that are connected to maxnum and have the most connections
     * 
     * @param attvals
     * @param clusterhash
     * @param seqnum
     * @param maxnum
     * @param minlinks
     * @return
     */
    static int[] getmaxshared(MinimalAttractionValue[] attvals, int[][] clusterhash, int seqnum, int maxnum,
            int minlinks) {
		int[] retarr = new int[minlinks + 1];
		retarr[0] = maxnum;
		int[] maxconn = new int[minlinks + 1];
		for (int i = 1; i <= minlinks; i++) {
			retarr[i] = -1;
			maxconn[i] = minlinks - 1;// make sure I only get sequences with at
										// least minlinks connections in the
										// output
		}// end for i
		for (int i = 0; i < seqnum; i++) {
			if ((clusterhash[i][1] == 1)
					&& (clusterhash[i][0] == -1)) {	// if
														// this
														// sequence
														// had
														// a
														// connection
														// to
														// maxnum
														// and
														// is
														// currently
														// unassigned
				for (int j = 1; j <= minlinks; j++) {
					if (clusterhash[i][2] > maxconn[j]) {
						// shift all subsequent values
						for (int k = minlinks; k > j; k--) {
							maxconn[k] = maxconn[k - 1];
							retarr[k] = retarr[k - 1];
						}// end for k
						maxconn[j] = clusterhash[i][2];
						retarr[j] = i;
						break;
					}
				}// end for j
			}
		}// end for i
		return retarr;
	}

	/**
	 * get the sequences with a connection to maxnum
	 * @param maxnum
	 * @param seqnum
	 * @param attvals
	 * @return
	 */
    static int[] getconnecteds(int maxnum,
			int seqnum, MinimalAttractionValue[] attvals) {
		int attnum = attvals.length;
		HashSet<Integer> tmphash = new HashSet<Integer>();
		for (int i = 0; i < attnum; i++) {
			if (attvals[i].query == maxnum) {
				// System.out.println("found hit:"+attvals[i].query+":"+attvals[i].hit);
				if (tmphash.contains(attvals[i].hit) == false) {
					// System.out.println("adding:"+attvals[i].hit+" as "+attvals[i].hit);
					tmphash.add(attvals[i].hit);
				}
			} else if (attvals[i].hit == maxnum) {
				// System.out.println("found hit:"+attvals[i].hit+":"+attvals[i].query);
				if (tmphash.contains(attvals[i].query) == false) {
					// System.out.println("adding:"+attvals[i].query+" as "+attvals[i].query);
					tmphash.add(attvals[i].query);
				}
			}
		}// end for i
			// now go through tmphash and get the integers of the sequences with
			// connection to maxnum
		int[] retarr = new int[tmphash.size()];
		int count = 0;
		for (int i = seqnum - 1; i >= 0; i--) {
			if (tmphash.contains(i)) {
				retarr[count] = i;
				count++;
			}
		}// end for i
		return retarr;
	}// end for getconnecteds

	// --------------------------------------------------------------------------
	static void getmaxnum(MinimalAttractionValue[] attvals, int[][] clusterhash,
			int hashpos) {
		// get how many connections to others each sequence has
		int attnum = attvals.length;
		for (int i = 0; i < attnum; i++) {
			clusterhash[attvals[i].query][hashpos]++;
			clusterhash[attvals[i].hit][hashpos]++;
		}// end for i
		return;
	}// end getmaxnum

	// --------------------------------------------------------------------------
	// --------------------------------------------------------------------------
	static void singlelinkage(MinimalAttractionValue[] attvals, Vector<SequenceCluster> retvec,
			int minseqnum, int seqnum) {
		// get all sequences connected via single linkage
		int[] groupvals = new int[seqnum];
		for (int i = 0; i < seqnum; i++) {
			groupvals[i] = -1;
		}// end for i
		boolean foundnew = true;
		int counter = 0;
		int attnum = attvals.length;
		while (foundnew) {
			// counter++;
			foundnew = false;
			int newnum = -1;
			for (int i = 0; i < seqnum; i++) {
				if (groupvals[i] == -1) {
					newnum = i;
					foundnew = true;
					break;
				}
			}// end for i
			if (newnum > -1) {
				System.out.print(".");
				tracesingle(attvals, seqnum, groupvals, counter, newnum, attnum);	// trace
																					// the
																					// single
																					// linkages
				counter++;
			}// end if newnum>-1
		}// end while foundnew
		SequenceCluster[] clusterarr = new SequenceCluster[counter];
		for (int i = 0; i < counter; i++) {
			clusterarr[i] = new SequenceCluster();
		}// end for i
		for (int i = 0; i < seqnum; i++) {
			// System.out.println("num="+i+" group:"+groupvals[i]);
			clusterarr[groupvals[i]].add(i);
		}// end for i
			// now sort the clusterarr by cluster size
		java.util.Arrays.sort(clusterarr, new ClusterSizeComparator());
		for (int i = 0; i < counter; i++) {
            if (clusterarr[i].members.length >= minseqnum) {
				retvec.addElement(clusterarr[i]);
			} else {
				break;
			}
		}// end for i
		return;
	}// end singlelinkage

	// --------------------------------------------------------------------------
	static void tracesingle(MinimalAttractionValue[] attvals, int seqnum, int[] groupvals,
			int counter, int newnum, int attnum) {
		// assign all single linked sequences the same counter value in
		// groupvals.
		groupvals[newnum] = counter;
		for (int i = 0; i < attnum; i++) {
			if ((attvals[i].query == newnum)
					&& (groupvals[attvals[i].hit] != counter)) {
				tracesingle(attvals, seqnum, groupvals, counter,
						attvals[i].hit, attnum);
			} else if ((attvals[i].hit == newnum)
					&& (groupvals[attvals[i].query] != counter)) {
				tracesingle(attvals, seqnum, groupvals, counter,
						attvals[i].query, attnum);
			}
		}// end for i
		return;
	}// end tracesingle

	// --------------------------------------------------------------------------
	// ---------------------convex
	// clustering------------------------------------
	// --------------------------------------------------------------------------
	public static Vector<SequenceCluster> getConvex(MinimalAttractionValue[] attractionValues,
			float sigmaFactor, int minSeqNum, int seqNum) {
		ConvexClustering convex = new ConvexClustering(attractionValues, sigmaFactor, minSeqNum, seqNum);
		return convex.getConvex();
	}

	private static class ConvexClustering {

		// Input variables
		private MinimalAttractionValue[] attractionValues;
		private float sigmaFactor;
		private int minSeqNum;
		private int seqNum;

		// Internal variables
		private ArrayList<Integer> remainingSeqIDs;
		private ArrayList<Integer> newClusterSeqIDs;
		private ArrayList<ArrayList<MinimalAttractionValue>> queryAttractionValues;
		private ArrayList<ArrayList<MinimalAttractionValue>> hitAttractionValues;
		private float avgAttraction;
		private float attractionVar;

		private ConvexClustering(
		                         MinimalAttractionValue[] attractionValues,
		                         float sigmaFactor,
		                         int minSeqNum,
		                         int seqNum
		                        )
		{
			this.attractionValues = attractionValues;
			this.sigmaFactor      = sigmaFactor;
			this.minSeqNum        = minSeqNum;
			this.seqNum           = seqNum;
		}

		private void initialize()
		{
			// Assign all nodes to one cluster to start out
			HashSet<Integer> initHash = new HashSet<Integer>(this.seqNum);
			this.remainingSeqIDs = new ArrayList<Integer>(this.seqNum);
			for (int i = 0; i < this.seqNum; i++) {
				Integer seqID = new Integer(i);
				this.remainingSeqIDs.add(seqID);
				initHash.add(seqID);
			}

			// Then take a node from this base cluster and add all those with
			// higher affinity to it than to the base
			this.newClusterSeqIDs = new ArrayList<Integer>(this.seqNum);

			// Get the average attraction for all nodes
			this.computeAverageAttraction(initHash);
			// Get the variance of attraction over all nodes
			this.computeAttractionVariance(initHash);

			// What a stuff cannot create a generic array
			this.queryAttractionValues = new ArrayList<ArrayList<MinimalAttractionValue>>(this.seqNum);
			this.hitAttractionValues = new ArrayList<ArrayList<MinimalAttractionValue>>(this.seqNum);

			for (int i = 0; i < this.seqNum; i++)
			{
				this.queryAttractionValues.add(new ArrayList<MinimalAttractionValue>());
				this.hitAttractionValues.add(new ArrayList<MinimalAttractionValue>());
			}

			for(int i = 0; i < attractionValues.length; i++)
			{
				this.queryAttractionValues.get(attractionValues[i].query).add(attractionValues[i]);
				this.hitAttractionValues.get(attractionValues[i].hit).add(attractionValues[i]);
			}
		}

		private Vector<SequenceCluster> getConvex() {
			// return a vector containing all clusters found by a polythenic
			// additive method
			// assign all sequences to one cluster. then take always the sequence
			// with highest overall attraction and seed it to a new cluster
			// next add all sequences with attractions higher to the new cluster
			// than average+factor*variance

			this.initialize();

			// It is recommended to use  ArrayList instead, however that would
			// change the interface so I leave it as ToDo.
			Vector<SequenceCluster> returnClusters = new Vector<SequenceCluster>();

			while (this.remainingSeqIDs.size() > 0) {
				// Get the remaining node with maximum attraction value
				int seed = this.getMaxAttraction();

				if (seed == -1 && this.minSeqNum <= 1) {
					// If there are no more attraction values left
					// Add all leftover sequences as separate clusters
					while (this.remainingSeqIDs.size() > 0) {
						SequenceCluster newCluster = new SequenceCluster(this.remainingSeqIDs.remove(0));
						returnClusters.addElement(newCluster);
					}
					break;
				}

				this.getOneCluster(seed);

				if(this.newClusterSeqIDs.size() >= this.minSeqNum) {
					SequenceCluster newCluster = new SequenceCluster(this.newClusterSeqIDs);
					returnClusters.add(newCluster);
				}
				this.newClusterSeqIDs.clear();
			}

			// Now sort the vector
			returnClusters.sort(new ClusterSizeComparator());

			return returnClusters;
		}// end getConvex

		// --------------------------------------------------------------------------
		/*
		* static float getMaxAttraction(cluster c1, cluster c2, minattvals[] attvals){
		* //get the maximum attraction between these two clusters int
		* seqnum1=c1.members.length; int
		* seqnum2=c2.members.length; float retval=0;
		* float skipped=0;//used in the bootstrapping procedure where some values
		* are changed to -1 (removed) for(int i=0;i<seqnum1;i++){ for(int
		* j=0;j<seqnum2;j++){ if(attvals[c1.members[i]][c2.members[j]]>=0){
		* retval+=attvals[c1.members[i]][c2.members[j]]; }else{ skipped++; } }//end
		* for j }//end for i return retval/((seqnum1*seqnum2)-skipped); }//end
		* getMaxAttraction
		*/// --------------------------------------------------------------------------
		private int getMaxAttraction() {
			// Get the sequence with overall highest attraction values from remainingSeqIDs
			int remainingSeqNum = this.remainingSeqIDs.size();

			HashMap<Integer, Integer> remaining2Index = new HashMap<Integer, Integer>(remainingSeqNum);
			float[] sumVals = new float[remainingSeqNum];
			for (int i = 0; i < remainingSeqNum; i++) {
				remaining2Index.put(this.remainingSeqIDs.get(i), new Integer(i));
				sumVals[i] = 0;
			}

			float maxVal = 0;
			int maxNum = -1;
			int attNum = this.attractionValues.length;

			for(int i = 0; i < this.queryAttractionValues.size(); i++) {
				if (remaining2Index.containsKey(i)) {
					ArrayList<MinimalAttractionValue> queryPosArray = this.queryAttractionValues.get(i);
					int index = remaining2Index.get(i).intValue();
					for(int j = 0; j < queryPosArray.size(); j++) {
						sumVals[index] += queryPosArray.get(j).att;
					}
				}
			}

			for(int i = 0; i < this.hitAttractionValues.size(); i++) {
				if (remaining2Index.containsKey(i)) {
					ArrayList<MinimalAttractionValue> hitPosArray = this.hitAttractionValues.get(i);
					int index = remaining2Index.get(i).intValue();
					for(int j = 0; j < hitPosArray.size(); j++) {
						sumVals[index] += hitPosArray.get(j).att;
					}
				}
			}

			for (int i = 0; i < remainingSeqNum; i++) {
				if (sumVals[i] > maxVal) {
					maxVal = sumVals[i];
					maxNum = i;
				}
			}

			return maxNum;
		}// end getMaxAttraction

		// --------------------------------------------------------------------------
		private void getOneCluster(int seed) {
			// Split one cluster off remainingSeqIDs and put the representatives in newClusterSeqIDs

			this.newClusterSeqIDs.add(this.remainingSeqIDs.remove(seed));
			Integer seedNum = this.newClusterSeqIDs.get(0);
			HashSet<Integer> newSeqHash = new HashSet<Integer>();
			newSeqHash.add(seedNum);
			// Now add all those values with attraction to newClusterSeqIDs greater than to
			// remainingSeqIDs.
			int remainingSeqs = this.remainingSeqIDs.size();

			boolean foundNew = true;
			float limit = (this.avgAttraction + (this.sigmaFactor * this.attractionVar));
			while (foundNew) {

				foundNew = false;

				remainingSeqs = this.remainingSeqIDs.size();
				if (remainingSeqs % 100 == 0) {
					System.out.print(remainingSeqs);
				}

				float maxAtt = -1;
				int maxNum = -1;

				// Now get the element with highest attraction to the new vector of elements
				for (int i = 0; i < remainingSeqs; i++) {

					int newPos = this.remainingSeqIDs.get(i).intValue();
					float currAtt = this.getAverageLocalAttraction(newPos, newSeqHash);

					if (currAtt > maxAtt) {
						maxAtt = currAtt;
						maxNum = i;
					}
				}

				System.out.print(".");

				if (maxNum > -1) {
					if (limit < maxAtt) {
						newSeqHash.add(this.remainingSeqIDs.get(maxNum));
						this.newClusterSeqIDs.add(this.remainingSeqIDs.remove(maxNum));
						remainingSeqs--;
						foundNew = true;
					}
					// else foundnew==false
				}
			}// end while foundnew
		}// end getcluster

		// --------------------------------------------------------------------------
		private float getAverageLocalAttraction(int newPos, HashSet<Integer> newSeqHash) {

			int newClusterSeqs = this.newClusterSeqIDs.size();
			int attNum = this.attractionValues.length;
			float retVal = 0;
			int skipped = 0;

			// Now get the average attraction of newPos to the current cluster
			ArrayList<MinimalAttractionValue> queryPosArray = this.queryAttractionValues.get(newPos);
			for(int i = 0; i < queryPosArray.size(); i++) {
				MinimalAttractionValue queryAttVal = queryPosArray.get(i);
				if(newSeqHash.contains(queryAttVal.hit)) {
					if (queryAttVal.att >= 0) {
						retVal += queryAttVal.att;
					} else {
						skipped++;
					}
				}
			}

			ArrayList<MinimalAttractionValue> hitPosArray = this.hitAttractionValues.get(newPos);
			for(int i = 0; i < hitPosArray.size(); i++) {
				MinimalAttractionValue hitAttVal = hitPosArray.get(i);
				if(newSeqHash.contains(hitAttVal.query)) {
					if (hitAttVal.att >= 0) {
						retVal += hitAttVal.att;
					} else {
						skipped++;
					}
				}
			}

			return retVal / (float) (newClusterSeqs - skipped);
		}// end getAverageLocalAttraction

		// --------------------------------------------------------------------------
		private void computeAverageAttraction(HashSet<Integer> initHash) {
			// Get the average attraction value for all sequences.
			// Note, the attraction values should be symmetrical and only those >=0
			// should be considered.

			int attNum = this.attractionValues.length;
			float sumVal = 0;
			int skipped = 0;
			for (int i = 0; i < attNum; i++) {
				if (initHash.contains(this.attractionValues[i].query)
				||  initHash.contains(this.attractionValues[i].hit)) {
					if (this.attractionValues[i].att >= 0) {
						sumVal += this.attractionValues[i].att;
					} else {
						skipped++;
					}
				}
			}

			this.avgAttraction = sumVal / (float) (((this.seqNum * (this.seqNum - 1)) / 2) - skipped);
		}// end getAverageAttraction

		// --------------------------------------------------------------------------
		private void computeAttractionVariance(HashSet<Integer> initHash) {
			// Get the variance of the attraction values for this cluster.

			int attNum = this.attractionValues.length;
			float sumVal = 0;
			int skipped = 0;
			for (int i = 0; i < attNum; i++) {
				if (initHash.contains(this.attractionValues[i].query)
				||  initHash.contains(this.attractionValues[i].hit)) {
					if (this.attractionValues[i].att >= 0) {
						float tmpVal = this.attractionValues[i].att - this.avgAttraction;
						sumVal += java.lang.Math.sqrt(tmpVal * tmpVal);
					} else {
						skipped++;
					}
				}
			}

			this.attractionVar = sumVal / (float) (((this.seqNum * (this.seqNum - 1)) / 2) - skipped);
		}// end getAttractionVariance
	} // End ConvexClustering

	// ------------------------------------------------------------------------------
	// --------------------------network
	// clustering----------------------------------
	// ------------------------------------------------------------------------------
	public static Vector<SequenceCluster> getnetwork(MinimalAttractionValue[] attvals, int minseqnum, boolean dooffset, boolean globalaverage, int seqnum, int maxrounds) {
		// the "network" based clustering approach
		// all of the links are fully described in the minattvals array (i.e.
		// hit, query and value)
		// I want to return a vector containing cluster objects
		// first convert the minattvals into a list of element I can work with
		element[] elements = new element[seqnum];
		for (int i = seqnum; --i >= 0;) {
			elements[i] = new element();
			elements[i].id = i;
			elements[i].index = i;
		}// end for i
		MinimalAttractionValue tmpdat;
		float avgval = 0;
		for (int i = attvals.length; --i >= 0;) {
			tmpdat = attvals[i];
			if (elements[tmpdat.hit].links == null) {
				elements[tmpdat.hit].links = new ArrayList<element>();
				elements[tmpdat.hit].weights = new ArrayList<Float>();
			}
			elements[tmpdat.hit].links.add(elements[tmpdat.query]);
			elements[tmpdat.hit].weights.add(new Float(tmpdat.att));
			avgval += tmpdat.att;
		}// end for i
			// now I have the thing in an arraylist of element objects with
			// directional links
		if (dooffset) {
			// if I want to offset the values by their average value or the
			// global average value
			if (globalaverage) {
				avgval /= (seqnum * seqnum);
			} else {
				// only divide by the actual number of entries
				avgval /= attvals.length;
			}
			// now offset all of the values accordingly
			element myelem;
			System.out.println("offsetting data by " + avgval);
			for (int i = seqnum; --i >= 0;) {
				myelem = elements[i];
				if (myelem.weights != null) {
					for (int j = myelem.weights.size(); --j >= 0;) {
						myelem.weights.set(j, new Float(myelem.weights.get(j)
								.floatValue() - avgval));
					}// end for i
				}
			}// end for i
		}// end if dooffset
			// now cluster the data
		cluster(elements, maxrounds);
		// now assign the entries to clusters and return
		HashMap<Integer, ArrayList<element>> tmphash = new HashMap<Integer, ArrayList<element>>();
		Integer currid;
		for (int i = elements.length; --i >= 0;) {
			currid = new Integer(elements[i].id);
			if (tmphash.containsKey(currid)) {
				tmphash.get(currid).add(elements[i]);
			} else {
				ArrayList<element> tmplist = new ArrayList<element>();
				tmplist.add(elements[i]);
				tmphash.put(currid, tmplist);
			}
		}// end for i
		ArrayList<element>[] clusters = tmphash.values().toArray(new ArrayList[0]);
        
		System.out.println("sorting clusters by size (" + clusters.length + ")");
		java.util.Arrays.sort(clusters, new sizecomparator());
		
		Vector<SequenceCluster> clustervec = new Vector<SequenceCluster>();
		SequenceCluster currcluster;
		ArrayList<element> tmp;
		for (int i = clusters.length; --i >= 0;) {
			tmp = clusters[i];
			currcluster = new SequenceCluster();
			currcluster.name = "cluster " + i;
			currcluster.members = new int[tmp.size()];
			// System.out.println("finalizing cluster "+i+" elements="+tmp.size());
			if (tmp.size() >= minseqnum) {
				for (int j = tmp.size(); --j >= 0;) {
					currcluster.members[j] = tmp.get(j).index;
				}// end for j
				clustervec.add(currcluster);
			}
		}// end for i
		System.out.println("Done network based clustering");
		return clustervec;
	}// end getnetwork

	static class sizecomparator implements Comparator<ArrayList<element>> {

		public int compare(ArrayList<element> a1, ArrayList<element> a2) {
            int a1s = a1.size();
            int a2s = a2.size();

			if (a1s > a2s) {
				return 1;
			} else if (a1s < a2s) {
				return -1;
			} else {
				return 0;
			}
		}
	}

	static public class element {

		public element() {
		}

		int id = -1;
		float weight = -1;
		int newid = -1;
		float newweight = -1;
		int index = -1;
		java.util.ArrayList<element> links = null;
		java.util.ArrayList<Float> weights = null;
		//java.util.ArrayList<fuzzydat> idlist = null;
	}// end class element

	public static void cluster(element[] elements, int maxrounds) {
		// go through the elements one by one and assign the id's!
		boolean changed = true;
		int currround = 0;
		float checkval;
		int myid, maxid = -1;
		float maxval = 0;
		HashMap<Integer, Float> assignhash = new HashMap<Integer,Float>();
		// now loop while change happens
		Integer currid;
		Float currweight;
		element myelement;
		while (changed && currround <= maxrounds) {
			changed = false;
			currround++;
			System.out.println("Iteration " + currround);
			for (int i = elements.length; --i >= 0;) {
				myelement = elements[i];
				// System.out.println(i + " elementname=" + myelement.name);
				// now loop through the emission values and get the lowest
				// highest emission id
				assignhash.clear();
				if (myelement.links != null) {
					// System.out.println("\tlinks=" + myelement.links.size());
					for (int j = myelement.links.size(); --j >= 0;) {
						currid = new Integer(myelement.links.get(j).id);
						currweight = new Float(myelement.weights.get(j)
								.floatValue());// make a copy of this as I might
												// otherwise be actually
												// changing the value of this
												// weight!
						// System.out.println("\tlink " + j + " weight=" +
						// currweight + " -->id:" + currid);
						if (assignhash.containsKey(currid)) {
							currweight += assignhash.get(currid);
						}
						assignhash.put(currid, currweight);
					}// end for j
						// now all I have to do is find the assigned id with the
						// highest sum value
					Integer[] keys = assignhash.keySet()
							.toArray(new Integer[0]);
					maxid = -1;
					maxval = 0;
					for (int j = keys.length; --j >= 0;) {
						checkval = assignhash.get(keys[j]).floatValue();
						if (checkval > maxval) {
							maxval = checkval;
							maxid = keys[j].intValue();
						} else if (checkval == maxval) {
							// then see which id is smaller
							myid = keys[j].intValue();
							if (maxid > myid) {
								maxid = myid;
							}
						}
					}// end for j
						// System.out.println("\tmaxid=" + maxid + " maxweight="
						// + maxval);
					myelement.newid = maxid;
					myelement.newweight = maxval;
				} else {// i.e. I have no input data whatsoever,
					myelement.newid = -1;
					myelement.newweight = -1;
				}
			}// end for i
				// now I know which id I would assign to the relative elements.
				// now for the difficult part. if I just assign things this way,
				// elements just continually swap assignments.
				// I need to take that into account with some sort of reverse
				// lookup...
			for (int i = elements.length; --i >= 0;) {
				myelement = elements[i];
				// System.err.println(i + "\tid:" + myelement.id + "(" +
				// myelement.weight + ")\tnewid:" + myelement.newid + "(" +
				// myelement.newweight + ")");
				myid = myelement.id;
				// now loop through the linked objects and see whether any of
				// them would change to this id
				// if they do, subtract their weight from the current assignment
				// weight.
				// if the new assignment weight drops below the old one, don't
				// change the assignment!
				// NOTE: only do this if the new assignment id is larger than
				// the old one! (for directionality and controling forever
				// loops)
				// System.out.println(i + " " + myelement.name + " currid=" +
				// myid + "/" + myelement.weight + " newid=" + myelement.newid +
				// "/" + myelement.newweight);
				if (myelement.newid > myid) {
					// System.out.println("\t newid > old id");
					for (int j = myelement.links.size(); --j >= 0;) {
						if (myelement.links.get(j).newid == myid) {
							// System.out.println("\t\tID swap between " + j +
							// " and " + i + " weight=" +
							// myelement.weights.get(j).floatValue());
							// then subtract this weight from the newid weight
							if ((myelement.weights.get(j)).floatValue() > 0) {// if
																				// this
																				// was
																				// a
																				// positive
																				// interaction
																				// leading
																				// to
																				// the
																				// clustering
								myelement.newweight -= (myelement.weights
										.get(j)).floatValue();
							}// else don't do anything
						}
					}// end for j
				}
				// System.out.println(i + " " + myelement.name + " id:" +
				// myelement.id + "(" + myelement.weight + ")\tnewid:" +
				// myelement.newid + "(" + myelement.newweight + ")");
				if (myelement.newweight > 0
						&& myelement.newweight > myelement.weight) {
					// I need the check for newweight>0 as otherwise unconnected
					// dots get added to a !random cluster.
					// System.out.println("\tNEW ASSIGNMENT: id=" +
					// myelement.newid + " weight=" + myelement.newweight);
					myelement.id = myelement.newid;
					myelement.weight = myelement.newweight;
					changed = true;
				}
			}// end for i
		}// end while changed or smaller maxrounds
		if (changed == false) {
			System.out.println("NO more chages performed after iteration :"
					+ currround);
		} else {
			System.out
					.println("Stopped clustering because the limit ot iterations was reached: "
							+ currround);
		}
		// that's it. the elements should have their new id assigned!
	}

	static int getmax(int pos, netnode[] nodes) {
		// get the cluster number with maximum attraction for this node
		netnode main = nodes[pos];
		Vector<float[]> conn = main.connections;
		float maxnum = -1, maxatt = 0;
		HashMap<String, float[]> tmphash = new HashMap<String, float[]>();
		String key;
		float[] tmparr, newarr;
		for (int i = conn.size() - 1; i >= 0; i--) {
			tmparr = conn.elementAt(i);
			key = String.valueOf(nodes[(int) tmparr[0]].number);
			if (tmphash.containsKey(key)) {
				newarr = tmphash.get(key);
				newarr[1] += tmparr[1];
				tmphash.put(key, newarr);
			} else {
				newarr = new float[2];
				newarr[0] = nodes[(int) tmparr[0]].number;
				newarr[1] = tmparr[1];
				tmphash.put(key, newarr);
			}
		}// end for i
		float[][] retarr = tmphash.values()
				.toArray(new float[0][2]);
		for (int i = retarr.length - 1; i >= 0; i--) {
			if (retarr[i][1] > maxatt) {
				maxatt = retarr[i][1];
				maxnum = retarr[i][0];
			}
		}// end for i
		return (int) maxnum;
	}// end getmax

	// --------------------------------------------------------------------------
	static class netnode {

		public netnode(int number) {
			this.number = number;
			connections = new Vector<float[]>();
		}

		int number;
		int oldnumber;
		Vector<float[]> connections;

		public void addconn(int hit, float conn) {
			float[] tmparr = new float[2];
			tmparr[0] = hit;
			tmparr[1] = conn;
			connections.addElement(tmparr);
		}
	}// end class netnode
}

class ClusterSizeComparator implements Comparator<SequenceCluster> {

    public int compare(SequenceCluster o1, SequenceCluster o2) {
         // sort largest first
        int num1 = o1.members.length;
        int num2 = o2.members.length;
		return (num1 < num2 ? 1 : (num1 == num2 ? 0 : -1));
	}
}
