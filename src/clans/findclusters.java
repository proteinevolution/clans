/*
 * findclusters.java
 *
 * Created on March 9, 2004, 4:05 PM
 */
package clans;

import java.util.*;

/**
 * 
 * @author tancred
 */
public class findclusters {

	/** Creates a new instance of findclusters */
	public findclusters() {
	}

	// --------------------------------------------------------------------------
	// ---------------------linkage
	// clustering-----------------------------------
	// --------------------------------------------------------------------------
	public static Vector<cluster> getlinkage(minattvals[] attvals,
			int minlinks, int minseqnum, int seqnum) {
		// form clusters with at least each seq connected by minlinks to the
		// cluster
		Vector<cluster> retvec = new Vector<cluster>();
		// simple cases, 0 and 1
		if (minlinks < 1) {// return all sequences
			cluster newcluster = new cluster();
			int[] tmpint = new int[seqnum];
			for (int i = 0; i < seqnum; i++) {
				tmpint[i] = i;
			}// end for i
			newcluster.add(tmpint);
			retvec.add(newcluster);
		} else if (minlinks == 1) {// single linkage
			singlelinkage(attvals, retvec, minseqnum, seqnum);
			// now sort the return vector
			// cluster is sorted in singlelinkage method
		} else {// otherwise it gets a bit more difficult
			multilinkage(attvals, retvec, minseqnum, minlinks, seqnum);
			// now sort the vector
			int clusternum = retvec.size();
			cluster[] tmparr = new cluster[clusternum];
			retvec.copyInto(tmparr);
			java.util.Arrays.sort(tmparr, new clustersizecomparator());
			for (int i = 0; i < clusternum; i++) {
				retvec.setElementAt(tmparr[i], i);
			}// end for i
		}
		return retvec;
	}// end getlinkage

	// --------------------------------------------------------------------------
	static void multilinkage(minattvals[] attvals, Vector<cluster> retvec,
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

		// make a Hash and an array of strings.
		// make the corresponding hash containing as value an array of int
		HashMap clusterhash = new HashMap((int) (seqnum / 0.8) + 1, 0.8f);
		String[] hashkeys = new String[seqnum];
		int[] tmparr;
		for (int i = 0; i < seqnum; i++) {
			hashkeys[i] = String.valueOf(i);
			tmparr = new int[3];
			tmparr[0] = -1;// groupvals; is a sequence assigned to a cluster and
							// if so, which cluster
			// groupval: -3=used as seed2 in last round, -2=don't use as seed,
			// -1=unassigned, 0=added to this cluster; >0=number of last cluster
			// this seq was assigned to
			tmparr[1] = 0;// used as seqshits (flag whether or not this sequence
							// is relevant in this calculation
			tmparr[2] = 0;// how many connections each sequence has to all
							// others
			clusterhash.put(hashkeys[i], tmparr);// groupvals,
		}// end for i
		getmaxnum(attvals, clusterhash, hashkeys, 2);// get the maximum number
														// of connections for
														// each sequence
		// for(int i=0;i<seqnum;i++){
		// System.out.println(i+" "+((int[])clusterhash.get(hashkeys[i]))[2]);
		// }//end for i
		int[] maxnumarr;
		boolean done = false;
		int attnum = java.lang.reflect.Array.getLength(attvals);
		int[] tmpseqs = new int[minlinks];
		int[] tmpseqs2;
		int counter = 1;// because 0 and -1 are used for sth. else
		cluster newcluster;
		int[] seqshits = new int[seqnum];
		int maxconn = 0;// counter for the max number of connections
		int maxnum = -1;// the sequence number with highest maxconn
		GETCLUSTER: while (done == false) {
			// first reset all the unassigned cluster sequences
			maxconn = 0;
			for (int i = 0; i < seqnum; i++) {
				if (((int[]) clusterhash.get(hashkeys[i]))[0] == 0
						|| ((int[]) clusterhash.get(hashkeys[i]))[0] == -3) {
					((int[]) clusterhash.get(hashkeys[i]))[0] = -1;
				}
				if (((int[]) clusterhash.get(hashkeys[i]))[0] == -1
						&& ((int[]) clusterhash.get(hashkeys[i]))[2] > maxconn) {// if
																					// this
																					// sequence
																					// is
																					// unassigned
																					// and
																					// has
																					// max
																					// connections
					maxconn = ((int[]) clusterhash.get(hashkeys[i]))[2];
					maxnum = i;
					// System.out.println("maxconn="+maxconn+" new seed="+maxnum);

				}
			}// end for i
				// now see which sequences were connected to maxnum
				// System.out.println("seed="+maxnum+" maxconn="+maxconn);
			System.out.print(".");
			int[] conarr = getconnecteds(maxnum, clusterhash, hashkeys, attvals);
			if (java.lang.reflect.Array.getLength(conarr) < minlinks) {// if I
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
				((int[]) clusterhash.get(hashkeys[i]))[1] = 0;// seqshits[i]=0;//the
																// number of
																// links
																// connecting
																// sequences to
																// this group
			}// end for i
			for (int i = java.lang.reflect.Array.getLength(conarr) - 1; i >= 0; i--) {
				// remember which sequences had hits with current seed
				((int[]) clusterhash.get(hashkeys[conarr[i]]))[1] = 1;// seqshits[maxnumarr[i]]=1;
			}
			((int[]) clusterhash.get(hashkeys[maxnum]))[0] = 1;// counter;//groupvals[maxnum]=counter;
			newcluster = new cluster();
			while (java.lang.reflect.Array.getLength(newcluster.member) == 0) {
				// get the minlinks other sequences that share the most hits
				// with this sequence
				tmpseqs = getmaxshared(attvals, clusterhash, hashkeys, maxnum,
						minlinks);
				// now see that I found at least minlinks valid connected
				// sequences
				for (int i = 0; i <= minlinks; i++) {
					if (tmpseqs[i] == -1) {
						// remove the current seed from analysis until end
						// (where it might come back in)
						((int[]) clusterhash.get(hashkeys[maxnum]))[0] = -2;// skip
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
				for (int i = java.lang.reflect.Array.getLength(tmpseqs) - 1; i >= 0; i--) {
					((int[]) clusterhash.get(hashkeys[tmpseqs[i]]))[0] = -3;// groupvals[tmpseqs[i]]=-3;//am
																			// checking,
																			// remeber
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

				tmpseqs2 = getfound(attvals, clusterhash, hashkeys, tmpseqs,
						minlinks);// get all sequences found by all of the seqs
									// in tmpseqs
				// System.out.println("hitseqs="+java.lang.reflect.Array.getLength(tmpseqs2));
				for (int i = java.lang.reflect.Array.getLength(tmpseqs2) - 1; i >= 0; i--) {
					newcluster.add(tmpseqs2[i]);
					// System.out.println("\tseq:"+tmpseqs2[i]);
					if (((int[]) clusterhash.get(hashkeys[tmpseqs2[i]]))[0] != -2) {
						((int[]) clusterhash.get(hashkeys[tmpseqs2[i]]))[0] = 0;// groupvals[i]=0;//added
																				// to
																				// cluster
					}
				}// end for i
					// now see if I have members>minlinks
				int foundnew = 0;
				int newmembers;
				while (((newmembers = java.lang.reflect.Array
						.getLength(newcluster.member)) < minlinks)
						&& (newmembers != foundnew)) {
					// System.out.println("inwhile");
					foundnew = newmembers;
					// System.out.println("in while");
					// if I don't, I need to re-check, including the new members
					// until I do have more than minlinks
					tmpseqs2 = getfound(attvals, clusterhash, hashkeys,
							newcluster.member, minlinks);
					for (int i = java.lang.reflect.Array.getLength(tmpseqs2) - 1; i >= 0; i--) {
						if (((int[]) clusterhash.get(hashkeys[tmpseqs2[i]]))[0] != 0) {
							newcluster.add(tmpseqs2[i]);
							if (((int[]) clusterhash.get(hashkeys[tmpseqs2[i]]))[0] != -2) {
								((int[]) clusterhash.get(hashkeys[tmpseqs2[i]]))[0] = 0;// groupvals[i]=0;//added
																						// to
																						// cluster
							}
						}
					}// end for i
				}// end while newcluster.members<minlinks
					// System.out.println("newcluster.members="+java.lang.reflect.Array.getLength(newcluster.member));
				if (java.lang.reflect.Array.getLength(newcluster.member) < minlinks) {
					// System.out.println("setting "+maxnum+" to -2");
					// if my currnet cluster has less than minlinks elements
					((int[]) clusterhash.get(hashkeys[maxnum]))[0] = -2;// skip
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
			}// end while newcluster.members==0
				// now I have a new cluster with at least minlinks members
				// next, see which sequences they find with at least minlinks
				// connections
			int lastmembers = 0;
			int newmembers;
			// System.out.println("startmembers="+java.lang.reflect.Array.getLength(newcluster.member));
			while ((newmembers = java.lang.reflect.Array
					.getLength(newcluster.member)) != lastmembers) {
				lastmembers = newmembers;
				// System.out.println("in while; members="+newmembers);
				// while I find sequences to add
				tmpseqs2 = getfound(attvals, clusterhash, hashkeys,
						newcluster.member, minlinks);
				for (int i = java.lang.reflect.Array.getLength(tmpseqs2) - 1; i >= 0; i--) {
					if (((int[]) clusterhash.get(hashkeys[tmpseqs2[i]]))[0] != 0) {
						newcluster.add(tmpseqs2[i]);
						((int[]) clusterhash.get(hashkeys[tmpseqs2[i]]))[0] = 0;// groupvals[i]=0;//added
																				// to
																				// cluster
					}
				}// end for i
			}// end while
				// now I have all sequences I can add in this way
				// set the groups info
			for (int i = java.lang.reflect.Array.getLength(newcluster.member) - 1; i >= 0; i--) {
				// if(((int[])clusterhash.get(hashkeys[newcluster.member[i]]))[0]==0){
				((int[]) clusterhash.get(hashkeys[newcluster.member[i]]))[0] = counter;
				// }
			}// end for i
			counter++;
			if (java.lang.reflect.Array.getLength(newcluster.member) > minseqnum) {
				// System.out.println("adding cluster "+counter);
				retvec.addElement(newcluster);
			}
		}// end while done==false
		return;
	}// end multilinkage

	// --------------------------------------------------------------------------
	static int[] getfound(minattvals[] attvals, HashMap clusterhash,
			String[] hashkeys, int[] tmpseqs, int minlinks) {
		// get all the sequences found by all of the tmpseqs
		int attnum = java.lang.reflect.Array.getLength(attvals);
		HashMap tmphash = new HashMap();
		int seqnum = java.lang.reflect.Array.getLength(hashkeys);
		for (int i = java.lang.reflect.Array.getLength(tmpseqs) - 1; i >= 0; i--) {
			tmphash.put(hashkeys[tmpseqs[i]], null);
		}// end for i
			// initialize to zero
		for (int i = 0; i < seqnum; i++) {
			((int[]) clusterhash.get(hashkeys[i]))[1] = 0;// seqshits[i]=0;
		}// end for i
			// now see which were found by the above
		for (int i = 0; i < attnum; i++) {
			if (tmphash.containsKey(hashkeys[attvals[i].query])) {
				((int[]) clusterhash.get(hashkeys[attvals[i].hit]))[1]++;
			}
			if (tmphash.containsKey(hashkeys[attvals[i].hit])) {
				((int[]) clusterhash.get(hashkeys[attvals[i].query]))[1]++;
			}
		}// end for i
		Vector<Integer> tmpvec = new Vector<Integer>();
		for (int i = 0; i < seqnum; i++) {
			if (((int[]) clusterhash.get(hashkeys[i]))[1] >= minlinks) {
				tmpvec.addElement(new Integer(i));
			}
		}// end for i
		int[] retarr = new int[tmpvec.size()];
		for (int i = tmpvec.size() - 1; i >= 0; i--) {
			retarr[i] = ((Integer) tmpvec.elementAt(i)).intValue();
		}// end for i
		return retarr;
	}// end getfound

	// --------------------------------------------------------------------------
	static int[] getmaxshared(minattvals[] attvals, HashMap clusterhash,
			String[] hashkeys, int maxnum, int minlinks) {
		// get the minlinks-1 other sequences that are connected to maxnum and
		// have the most connections
		int[] retarr = new int[minlinks + 1];
		int seqnum = java.lang.reflect.Array.getLength(hashkeys);
		retarr[0] = maxnum;
		int[] maxconn = new int[minlinks + 1];
		for (int i = 1; i <= minlinks; i++) {
			retarr[i] = -1;
			maxconn[i] = minlinks - 1;// make sure I only get sequences with at
										// least minlinks connections in the
										// output
		}// end for i
		for (int i = 0; i < seqnum; i++) {
			if ((((int[]) clusterhash.get(hashkeys[i]))[1] == 1)
					&& (((int[]) clusterhash.get(hashkeys[i]))[0] == -1)) {// if
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
					if (((int[]) clusterhash.get(hashkeys[i]))[2] > maxconn[j]) {
						// shift all subsequent values
						for (int k = minlinks; k > j; k--) {
							maxconn[k] = maxconn[k - 1];
							retarr[k] = retarr[k - 1];
						}// end for k
						maxconn[j] = ((int[]) clusterhash.get(hashkeys[i]))[2];
						retarr[j] = i;
						break;
					}
				}// end for j
			}
		}// end for i
		return retarr;
	}// end getmaxshared

	// --------------------------------------------------------------------------
	static int[] getconnecteds(int maxnum, HashMap clusterhash,
			String[] hashkeys, minattvals[] attvals) {
		// get the sequences with a connection to maxnum
		int attnum = java.lang.reflect.Array.getLength(attvals);
		HashMap tmphash = new HashMap();
		for (int i = 0; i < attnum; i++) {
			if (attvals[i].query == maxnum) {
				// System.out.println("found hit:"+attvals[i].query+":"+attvals[i].hit);
				if (tmphash.containsKey(hashkeys[attvals[i].hit]) == false) {
					// System.out.println("adding:"+attvals[i].hit+" as "+hashkeys[attvals[i].hit]);
					tmphash.put(hashkeys[attvals[i].hit], null);
				}
			} else if (attvals[i].hit == maxnum) {
				// System.out.println("found hit:"+attvals[i].hit+":"+attvals[i].query);
				if (tmphash.containsKey(hashkeys[attvals[i].query]) == false) {
					// System.out.println("adding:"+attvals[i].query+" as "+hashkeys[attvals[i].query]);
					tmphash.put(hashkeys[attvals[i].query], null);
				}
			}
		}// end for i
			// now go through tmphash and get the integers of the sequences with
			// connection to maxnum
		int[] retarr = new int[tmphash.size()];
		int count = 0;
		for (int i = java.lang.reflect.Array.getLength(hashkeys) - 1; i >= 0; i--) {
			if (tmphash.containsKey(hashkeys[i])) {
				retarr[count] = i;
				count++;
			}
		}// end for i
		return retarr;
	}// end for getconnecteds

	// --------------------------------------------------------------------------
	static void getmaxnum(minattvals[] attvals, HashMap clusterhash,
			String[] hashkeys, int hashpos) {
		// get how many connections to others each sequence has
		int attnum = java.lang.reflect.Array.getLength(attvals);
		for (int i = 0; i < attnum; i++) {
			((int[]) clusterhash.get(hashkeys[attvals[i].query]))[hashpos]++;
			((int[]) clusterhash.get(hashkeys[attvals[i].hit]))[hashpos]++;
		}// end for i
		return;
	}// end getmaxnum

	// --------------------------------------------------------------------------
	// --------------------------------------------------------------------------
	static void singlelinkage(minattvals[] attvals, Vector<cluster> retvec,
			int minseqnum, int seqnum) {
		// get all sequences connected via single linkage
		int[] groupvals = new int[seqnum];
		for (int i = 0; i < seqnum; i++) {
			groupvals[i] = -1;
		}// end for i
		boolean foundnew = true;
		int counter = 0;
		int attnum = java.lang.reflect.Array.getLength(attvals);
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
				tracesingle(attvals, seqnum, groupvals, counter, newnum, attnum);// trace
																					// the
																					// single
																					// linkages
				counter++;
			}// end if newnum>-1
		}// end while foundnew
		cluster[] clusterarr = new cluster[counter];
		for (int i = 0; i < counter; i++) {
			clusterarr[i] = new cluster();
		}// end for i
		for (int i = 0; i < seqnum; i++) {
			// System.out.println("num="+i+" group:"+groupvals[i]);
			clusterarr[groupvals[i]].add(i);
		}// end for i
			// now sort the clusterarr by cluster size
		java.util.Arrays.sort(clusterarr, new clustersizecomparator());
		for (int i = 0; i < counter; i++) {
			if (java.lang.reflect.Array.getLength(clusterarr[i].member) >= minseqnum) {
				retvec.addElement(clusterarr[i]);
			} else {
				break;
			}
		}// end for i
		return;
	}// end singlelinkage

	// --------------------------------------------------------------------------
	static void tracesingle(minattvals[] attvals, int seqnum, int[] groupvals,
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
	public static Vector<cluster> getconvex(minattvals[] attvals,
			float sigmafac, int minseqnum, int seqnum) {
		// return a vector containing all clusters found by a polythenic
		// additive method
		// assign all sequencs to one cluster. then take always the sequence
		// with highest overall attraction and seed it to a new cluster
		// next add all sequences with attractions higher to the new cluster
		// than average+factor*variance
		// System.out.println("in getconvex getclusters");
		Vector<cluster> retvec = new Vector<cluster>();
		// assign all nodes to one cluster to start out
		Vector basevec = new Vector(seqnum);
		String[] hashkeys = new String[seqnum];
		for (int i = 0; i < seqnum; i++) {
			hashkeys[i] = String.valueOf(i);
			basevec.add(new Integer(i));
		}// end for i
			// then take a node from this base cluster and add all those with
			// higher affinity to it than to the base
		Vector currvec = new Vector(seqnum);
		cluster currcluster;
		float avgatt = getavgatt(basevec, attvals, hashkeys);// get the average
																// attraction
																// for the
																// values in
																// that vector
		float varatt = getvaratt(basevec, attvals, avgatt, hashkeys);// get the
																		// variance
																		// of
																		// attraction
																		// values
		int seed;
		while (basevec.size() > 0) {
			// System.out.println("doing getmaxatt");
			seed = getmaxatt(basevec, attvals, hashkeys);// get the vector
															// element with the
															// highest overall
															// attraction
			// System.out.println("done getmaxatt");
			// System.out.println("seed="+seed+" basevec.size="+basevec.size());
			if (seed == -1) {// if I have no further attraction values in the
								// matrix
				// System.out.println("adding all leftover sequences");
				// add all leftover sequences as separate clusters
				while (basevec.size() > 0) {
					currcluster = new cluster();
					currcluster.add(((Integer) basevec.remove(0)).intValue());
					retvec.addElement(currcluster);
				}// end while basevec.size 2nd
				break;
			}// end if seed==-1
				// System.out.println("doing getcluster");
			getcluster(seed, basevec, attvals, currvec, avgatt, varatt,
					sigmafac, hashkeys);// get one more cluster
			// System.out.println("done getcluster");
			currcluster = new cluster();
			currcluster.member = new int[currvec.size()];
			for (int i = currvec.size() - 1; i >= 0; i--) {
				currcluster.member[i] = ((Integer) currvec.elementAt(i))
						.intValue();
			}// end for i
			retvec.add(currcluster);
			currvec.clear();
		}// end while basecluster>0
			// now sort the vector
		int clusterelements = retvec.size();
		cluster[] sortarr = new cluster[clusterelements];
		retvec.copyInto(sortarr);
		java.util.Arrays.sort(sortarr, new clustersizecomparator());
		retvec.clear();
		// now add them back in sorted order and filter for minimum sequence
		// number
		for (int i = 0; i < clusterelements; i++) {
			if (java.lang.reflect.Array.getLength(sortarr[i].member) >= minseqnum) {
				retvec.addElement(sortarr[i]);
			} else {
				break;
			}
		}// end for i
		return retvec;
	}// end get

	// --------------------------------------------------------------------------
	/*
	 * static float getavgatt(cluster c1, cluster c2, minattvals[] attvals){
	 * //get the average attraction between these two clusters int
	 * seqnum1=java.lang.reflect.Array.getLength(c1.member); int
	 * seqnum2=java.lang.reflect.Array.getLength(c2.member); float retval=0;
	 * float skipped=0;//used in the bootstrapping procedure where some values
	 * are changed to -1 (removed) for(int i=0;i<seqnum1;i++){ for(int
	 * j=0;j<seqnum2;j++){ if(attvals[c1.member[i]][c2.member[j]]>=0){
	 * retval+=attvals[c1.member[i]][c2.member[j]]; }else{ skipped++; } }//end
	 * for j }//end for i return retval/((seqnum1*seqnum2)-skipped); }//end
	 * getavgatt
	 */// --------------------------------------------------------------------------
	static int getmaxatt(Vector invec, minattvals[] attvals, String[] hashkeys) {
		// get the sequence with overall highest attvals
		int elements = invec.size();
		HashMap tmphash = new HashMap(elements);
		float[] sumvals = new float[elements];
		for (int i = 0; i < elements; i++) {
			tmphash.put(hashkeys[((Integer) invec.elementAt(i)).intValue()],
					new Integer(i));
			sumvals[i] = 0;
		}// end for i
		float maxval = 0;
		int maxnum = -1;
		int attnum = java.lang.reflect.Array.getLength(attvals);
		for (int i = 0; i < attnum; i++) {
			if (tmphash.containsKey(hashkeys[attvals[i].query])) {
				sumvals[((Integer) tmphash.get(hashkeys[attvals[i].query]))
						.intValue()] += attvals[i].att;
			}
			if (tmphash.containsKey(hashkeys[attvals[i].hit])) {
				sumvals[((Integer) tmphash.get(hashkeys[attvals[i].hit]))
						.intValue()] += attvals[i].att;
			}
		}// end for i
		for (int i = 0; i < elements; i++) {
			if (sumvals[i] > maxval) {
				maxval = sumvals[i];
				maxnum = i;
			}
		}// end for i
		return maxnum;
	}// end getmaxatt

	// --------------------------------------------------------------------------
	static void getcluster(int seed, Vector basevec, minattvals[] attvals,
			Vector currvec, float avgatt, float varatt, float sigmafac,
			String[] hashkeys) {
		// split one cluster off basevec and put the representatives in currvec
		// System.out.println("seed "+seed);//("in getcluster for seed "+seed);//"seed="+seed);
		currvec.add(basevec.remove(seed));
		int seednum = ((Integer) currvec.elementAt(0)).intValue();
		HashMap tmphash = new HashMap();
		tmphash.put(hashkeys[seednum], null);
		// now add all those values with attraction to currvec greater than to
		// basevec.
		int elements = basevec.size();
		int attnum = java.lang.reflect.Array.getLength(attvals);
		boolean foundnew = true;
		float maxatt, curratt;
		int maxnum = 0;
		float limit = (avgatt + (sigmafac * varatt));
		while (foundnew == true) {
			elements = basevec.size();
			if (elements % 100 == 0) {
				System.out.print(elements);
			}
			foundnew = false;
			maxatt = -1;
			maxnum = -1;
			// now get the element with highest attraction to the new vector of
			// elements
			for (int i = 0; i < elements; i++) {
				curratt = getavgatt(
						((Integer) basevec.elementAt(i)).intValue(), currvec,
						attvals, hashkeys, tmphash);
				// System.out.println("done getavgatt "+i);
				if (curratt > maxatt) {
					maxatt = curratt;
					maxnum = i;
				}
			}// end for i
			System.out.print(".");
			// System.out.println("elements="+elements+" maxnum="+maxnum+" maxatt="+maxatt);
			if (maxnum > -1) {
				if (limit < maxatt) {
					tmphash.put(hashkeys[((Integer) basevec.elementAt(maxnum))
							.intValue()], null);
					currvec.addElement(basevec.remove(maxnum));
					elements--;
					foundnew = true;
				}
				// else foundnew==false
			}
		}// end while foundnew
	}// end getcluster

	// --------------------------------------------------------------------------
	static float getavgatt(int newpos, Vector newgroup, minattvals[] attvals,
			String[] hashkeys, HashMap tmphash) {
		int elements = newgroup.size();
		int attnum = java.lang.reflect.Array.getLength(attvals);
		float retval = 0;
		float skipped = 0;
		// now get the average attraction of newpos to the current cluster
		for (int i = 0; i < attnum; i++) {
			if (attvals[i].hit == newpos
					&& tmphash.containsKey(hashkeys[attvals[i].query])) {
				if (attvals[i].att >= 0) {
					retval += attvals[i].att;
				} else {
					skipped++;
				}
			} else if (attvals[i].query == newpos
					&& tmphash.containsKey(hashkeys[attvals[i].hit])) {
				if (attvals[i].att >= 0) {
					retval += attvals[i].att;
				} else {
					skipped++;
				}
			}
		}// end for i
		return retval / (float) (elements - skipped);
	}// end getavgatt

	// --------------------------------------------------------------------------
	static float getavgatt(Vector invec, minattvals[] attvals, String[] hashkeys) {
		// get the average attraction value for all sequences in this vector
		// note, the attraction values should be symmetrical and only those >=0
		// should be considered
		int elements = invec.size();
		HashMap tmphash = new HashMap(elements);
		for (int i = 0; i < elements; i++) {
			tmphash.put(hashkeys[((Integer) invec.elementAt(i)).intValue()],
					null);
		}// end for i
		int attnum = java.lang.reflect.Array.getLength(attvals);
		float sumval = 0;
		int skipped = 0;
		for (int i = 0; i < attnum; i++) {
			if (tmphash.containsKey(hashkeys[attvals[i].query])
					|| tmphash.containsKey(hashkeys[attvals[i].hit])) {
				if (attvals[i].att >= 0) {
					sumval += attvals[i].att;
				} else {
					skipped++;
				}
			}
		}// end for i
		return sumval / (float) (((elements * (elements - 1)) / 2) - skipped);
	}// end getavgatt

	// --------------------------------------------------------------------------
	static float getvaratt(Vector invec, minattvals[] attvals, float avgval,
			String[] hashkeys) {
		// get the variance of the attraction values for this cluster
		int elements = invec.size();
		HashMap tmphash = new HashMap(elements);
		for (int i = 0; i < elements; i++) {
			tmphash.put(hashkeys[((Integer) invec.elementAt(i)).intValue()],
					null);
		}// end for i
		int attnum = java.lang.reflect.Array.getLength(attvals);
		float sumval = 0, tmpval;
		int skipped = 0;
		for (int i = 0; i < attnum; i++) {
			if (tmphash.containsKey(hashkeys[attvals[i].query])
					|| tmphash.containsKey(hashkeys[attvals[i].hit])) {
				if (attvals[i].att >= 0) {
					tmpval = attvals[i].att - avgval;
					sumval += java.lang.Math.sqrt(tmpval * tmpval);
				} else {
					skipped++;
				}
			}
		}// end for i
		return sumval / (float) (((elements * (elements - 1)) / 2) - skipped);
	}// end getvaratt

	// ------------------------------------------------------------------------------
	// --------------------------network
	// clustering----------------------------------
	// ------------------------------------------------------------------------------
	static Vector<cluster> getnetwork(minattvals[] attvals, int minseqnum, boolean dooffset, boolean globalaverage, int seqnum, int maxrounds) {
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
		minattvals tmpdat;
		float avgval = 0;
		for (int i = java.lang.reflect.Array.getLength(attvals); --i >= 0;) {
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
				avgval /= java.lang.reflect.Array.getLength(attvals);
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
		for (int i = java.lang.reflect.Array.getLength(elements); --i >= 0;) {
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
		System.out.println("sorting clusters by size ("
				+ java.lang.reflect.Array.getLength(clusters) + ")");
		java.util.Arrays.sort(clusters, new sizecomparator());
		element[] tmparr;
		int clustercount = 0;
		Vector<cluster> clustervec = new Vector<cluster>();
		cluster currcluster;
		ArrayList<element> tmp;
		for (int i = java.lang.reflect.Array.getLength(clusters); --i >= 0;) {
			tmp = clusters[i];
			currcluster = new cluster();
			currcluster.name = "cluster " + i;
			currcluster.member = new int[tmp.size()];
			// System.out.println("finalizing cluster "+i+" elements="+tmp.size());
			if (tmp.size() >= minseqnum) {
				for (int j = tmp.size(); --j >= 0;) {
					currcluster.member[j] = tmp.get(j).index;
				}// end for j
				clustervec.add(currcluster);
			}
		}// end for i
		System.out.println("Done network based clustering");
		return clustervec;
	}// end getnetwork

	static class sizecomparator implements Comparator {

		public int compare(Object a1, Object a2) {
			int a1s = ((ArrayList) a1).size();
			int a2s = ((ArrayList) a2).size();
			// parameter are of type Object, so we have to downcast it to
			// Employee objects
			if (a1s > a2s) {
				return 1;
			} else if (a1s < a2s) {
				return -1;
			} else {
				return 0;
			}
		}// end compare
	}// end sizecomparator

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
			for (int i = java.lang.reflect.Array.getLength(elements); --i >= 0;) {
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
					for (int j = java.lang.reflect.Array.getLength(keys); --j >= 0;) {
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
			for (int i = java.lang.reflect.Array.getLength(elements); --i >= 0;) {
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
	}// end cluster

/*	public static void fuzzycluster(element[] elements, int maxrounds) {
		// go through the elements one by one and assign the id's!
		// NOTE: this does not fully assign a node to either one group or
		// another but uses the difference between assignments as its weight (doesn't work very well.)
		//trying other stuff for the moment
		//try to get the avg + stdev assignment to each cluster and then take the one with the best overlap?
		boolean changed = true;
		int currround = 0;
		float checkval;
		int myid, bestlistsize, bestid = -1;
		float myavg,mystdev,bestavg,beststdev;
		HashMap<Integer, ArrayList<Float>> assignhash = new HashMap<Integer,ArrayList<Float>>();
		// now loop while change happens
		Integer currid;
		Float currweight;
		element myelement;
		float diff;
		ArrayList<Float> tmplist=null;
		while (changed && currround <= maxrounds) {
			changed = false;
			currround++;
			System.out.println("Fuzzy iteration " + currround);
			for (int i = java.lang.reflect.Array.getLength(elements); --i >= 0;) {
				myelement = elements[i];
				// System.out.println(i + " elementname=" + myelement.name);
				// now loop through the emission values and get the lowest
				// highest emission id
				assignhash.clear();
				if (myelement.links != null) {
					for (int j = myelement.links.size(); --j >= 0;) {
						currid = new Integer(myelement.links.get(j).id);
						currweight = new Float(myelement.weights.get(j).floatValue());// make a copy of this as I might otherwise be actually changing the
						if (assignhash.containsKey(currid)) {
							assignhash.get(currid).add(currweight);
						}else{
							ArrayList<Float> addlist=new ArrayList<Float>();
							addlist.add(currweight);
							assignhash.put(currid, addlist);
						}
					}// end for j
					Integer[] keys = assignhash.keySet().toArray(new Integer[0]);
					//for each of the keys now get the avg + var of the values
					bestid = -1;
					bestavg = -1;
					beststdev = -1;
					bestlistsize=-1;
					for (int j = java.lang.reflect.Array.getLength(keys); --j >= 0;) {
						tmplist=assignhash.get(keys[j]);
						myid = keys[j].intValue();
						myavg=0;
						mystdev=0;
						for(int k=tmplist.size();--k>=0;){
							myavg+=tmplist.get(k).floatValue();
						}
						myavg/=tmplist.size();
						for(int k=tmplist.size();--k>=0;){
							diff=myavg-tmplist.get(k).floatValue();
							mystdev+=diff*diff;
						}
						//special case for tmplist size==1
						if(tmplist.size()==1){
							mystdev=-1;
						}else{
							mystdev=(float) java.lang.Math.sqrt(mystdev/tmplist.size());
						}
						//now see whether this avg and stdev fit my datapoint better that what I already have
						if(bestid==-1){
							bestid=myid;
							bestavg=myavg;
							beststdev=mystdev;
							bestlistsize=tmplist.size();
						}else{
							//here then see whether the current cluster assignment fits this datapoint better than bestid, bestavg, beststdev
							//simple approach: divide average by stdev and take that value (kind of a Z-score relative to 0)
							//this works somewhat ok, but I am sure there is a better way
							if(beststdev==-1 || mystdev==-1){
								if(myavg>bestavg){
									bestid=myid;
									bestavg=myavg;
									beststdev=mystdev;
									bestlistsize=tmplist.size();
								}
							}else{
								if((myavg/mystdev)>(bestavg/beststdev)){
									bestid=myid;
									bestavg=myavg;
									beststdev=mystdev;		
									bestlistsize=tmplist.size();
								}								
							}
						}
					}
					myelement.newid = bestid;
					myelement.newweight = bestavg*bestlistsize;
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
			for (int i = java.lang.reflect.Array.getLength(elements); --i >= 0;) {
				myelement = elements[i];
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
				if (myelement.newid > myid) {
					for (int j = myelement.links.size(); --j >= 0;) {
						if (myelement.links.get(j).newid == myid) {
							// if these two would swap assignments then subtract this weight from the newid weight
							if ((myelement.weights.get(j)).floatValue() > 0) {
								myelement.newweight -= (myelement.weights.get(j)).floatValue();
							}// else don't do anything
						}
					}// end for j
				}
				if (myelement.newweight > 0	&& myelement.newweight > myelement.weight) {
					// I need the check for newweight>0 as otherwise unconnected dots get added to a !random cluster.
					myelement.id = myelement.newid;
					myelement.weight = myelement.newweight;
					changed = true;
				}
			}// end for i
		}// end while changed or smaller maxrounds
		if (changed == false) {
			System.out.println("NO more chages performed after fuzzy iteration :"+ currround);
		} else {
			System.out.println("Stopped fuzzy network clustering because the limit of iterations was reached: "	+ currround);
		}
		// that's it. the elements should have their new id assigned!
	}// end fuzzycluster
*/
	/*
	 * static Vector getnetwork(minattvals[] attvals, int minseqnum, boolean
	 * dooffset, boolean globalaverage, int seqnum, int maxrounds) { //use a
	 * "network approach //each element is assigned an input value //it emits
	 * that value to the next layer which sums the emissions of the same value
	 * //in the next round each node emits the value of the highest sum it
	 * encountered //repeat until no more changes int attnum =
	 * java.lang.reflect.Array.getLength(attvals); netnode[] nodes = new
	 * netnode[seqnum]; for (int i = 0; i < seqnum; i++) { nodes[i] = new
	 * netnode(i); }//end for i if (dooffset) { float average = 0; if
	 * (globalaverage) { //offset the values by the global average attraction
	 * value for (int i = 0; i < attnum; i++) { average += attvals[i].att;
	 * }//end for i average /= (seqnum * (seqnum - 1)) / 2; } else { //only get
	 * the average for the non-null attractions for (int i = 0; i < attnum; i++)
	 * { average += attvals[i].att; }//end for i average /= attnum; } //now
	 * assign the individual node objects their respective attractions for (int
	 * i = 0; i < attnum; i++) { nodes[attvals[i].query].addconn(attvals[i].hit,
	 * attvals[i].att - average);
	 * nodes[attvals[i].hit].addconn(attvals[i].query, attvals[i].att -
	 * average); }//end for i } else { //don't offset the values at all; just
	 * use as is. for (int i = 0; i < attnum; i++) {
	 * nodes[attvals[i].query].addconn(attvals[i].hit, attvals[i].att);
	 * nodes[attvals[i].hit].addconn(attvals[i].query, attvals[i].att); }//end
	 * for i } int round = 0; boolean changedval = true; int newnum; while
	 * ((changedval == true) && (round < maxrounds)) {//I don't want to do more
	 * than 100 rounds System.out.print("r" + round + ";"); round++; changedval
	 * = false; for (int i = 0; i < seqnum; i++) { newnum = getmax(i,
	 * nodes);//get the maximum weight cluster for this node if (newnum !=
	 * nodes[i].number) { changedval = true; if (newnum == nodes[i].oldnumber) {
	 * if (nodes[i].number < nodes[i].oldnumber) { //don't do anything; needed
	 * to keep nodes from swapping the same numbers over and over } else {
	 * nodes[i].oldnumber = nodes[i].number; nodes[i].number = newnum; } } else
	 * { nodes[i].oldnumber = nodes[i].number; nodes[i].number = newnum; } }
	 * }//end for i }//end while changedval //now I should have all the groups
	 * assigned; now convert all the data into a vector of clusters and return
	 * Vector retvec = new Vector(); changedval = true; int counter = 0; int
	 * currnum = -1; Vector tmpvec = new Vector(); cluster currcluster; while
	 * (changedval == true) { changedval = false; for (int i = 0; i < seqnum;
	 * i++) { if (changedval == false) { if (nodes[i].number != -1) { changedval
	 * = true; currnum = nodes[i].number; nodes[i].number = -1; if
	 * (tmpvec.size() > minseqnum) { currcluster = new cluster();
	 * currcluster.name = "cluster_" + counter; currcluster.member = new
	 * int[tmpvec.size()]; for (int j = tmpvec.size() - 1; j >= 0; j--) {
	 * currcluster.member[j] = ((Integer) tmpvec.elementAt(j)).intValue();
	 * }//end for i retvec.addElement(currcluster); counter++; } tmpvec.clear();
	 * tmpvec.addElement(new Integer(i)); } } else { if (nodes[i].number ==
	 * currnum) { nodes[i].number = -1; tmpvec.addElement(new Integer(i)); } }
	 * }//end for i }//end while if (tmpvec.size() > minseqnum) { currcluster =
	 * new cluster(); currcluster.member = new int[tmpvec.size()]; for (int j =
	 * tmpvec.size() - 1; j >= 0; j--) { currcluster.member[j] = ((Integer)
	 * tmpvec.elementAt(j)).intValue(); }//end for i
	 * retvec.addElement(currcluster); } int clusterelements = retvec.size();
	 * cluster[] sortarr = new cluster[clusterelements];
	 * retvec.copyInto(sortarr); java.util.Arrays.sort(sortarr, new
	 * clustersizecomparator()); retvec.clear(); //now add them back in sorted
	 * order and filter for minimum sequence number for (int i = 0; i <
	 * clusterelements; i++) { if
	 * (java.lang.reflect.Array.getLength(sortarr[i].member) >= minseqnum) {
	 * retvec.addElement(sortarr[i]); } else { break; } }//end for i return
	 * retvec; }//end getnetwork
	 */// --------------------------------------------------------------------------
	static int getmax(int pos, netnode[] nodes) {
		// get the cluster number with maximum attraction for this node
		netnode main = nodes[pos];
		Vector conn = main.connections;
		float maxnum = -1, maxatt = 0;
		HashMap tmphash = new HashMap();
		String key;
		float[] tmparr, newarr;
		for (int i = conn.size() - 1; i >= 0; i--) {
			tmparr = (float[]) conn.elementAt(i);
			key = String.valueOf(nodes[(int) tmparr[0]].number);
			if (tmphash.containsKey(key)) {
				newarr = ((float[]) tmphash.get(key));
				newarr[1] += tmparr[1];
				tmphash.put(key, newarr);
			} else {
				newarr = new float[2];
				newarr[0] = nodes[(int) tmparr[0]].number;
				newarr[1] = tmparr[1];
				tmphash.put(key, newarr);
			}
		}// end for i
		float[][] retarr = (float[][]) tmphash.values()
				.toArray(new float[0][2]);
		for (int i = java.lang.reflect.Array.getLength(retarr) - 1; i >= 0; i--) {
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
			connections = new Vector();
		}

		int number;
		int oldnumber;
		Vector connections;

		public void addconn(int hit, float conn) {
			float[] tmparr = new float[2];
			tmparr[0] = hit;
			tmparr[1] = conn;
			connections.addElement(tmparr);
		}
	}// end class netnode
}// end class

// --------------------------------------------------------------------------
// --------------------------comparator--------------------------------------
// --------------------------------------------------------------------------
class clustersizecomparator implements Comparator {
	// sort largest first
	public int compare(Object o1, Object o2) {
		int num1 = java.lang.reflect.Array.getLength(((cluster) o1).member);
		int num2 = java.lang.reflect.Array.getLength(((cluster) o2).member);
		return (num1 < num2 ? 1 : (num1 == num2 ? 0 : -1));
	}// end compare
}// end comparator class
