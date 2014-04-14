/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package clans;
import java.io.*;
/**
 *
 * @author tancred
 */
public class ThreadedStreamReader  extends Thread{

    BufferedReader inread;
    StringBuffer outread;
    boolean done=false;
    /** Creates new threadstreamreader */
    public ThreadedStreamReader(BufferedReader inbuf, StringBuffer outbuf) {
        this.inread=inbuf;
        this.outread=outbuf;
        done=false;
    }
    
    public void run(){
        int readchar;
        try{
            synchronized(outread){//synchr to be sure that all data is written before trying to read from
                while((readchar=inread.read())!=-1){
                    outread.append((char)readchar);
                }
            }// end synchronized
            inread.close();
            synchronized (outread){
                this.done=true;
                outread.notify();
            }
       }catch (Exception e){
            System.err.println("read error, not end of stream!");
            e.printStackTrace();
        }// end catch
    }// end run

    //public void finalize(){
    //    System.out.println("closing threadreader");
    //}
    
    
}
