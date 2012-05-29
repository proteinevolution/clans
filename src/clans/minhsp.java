/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package clans;

/**
 *
 * @author tancred
 */
public class minhsp {
 
    /** Creates a new instance of minhsp */
    public minhsp() {
    }
    
    public minhsp(int i1, int i2, double[] pv){
        this.query=i1;
        this.hit=i2;
        this.val=pv;
    }
    
    public minhsp(int i1, int i2, double pv){
        this.query=i1;
        this.hit=i2;
        this.val=new double[1];
        val[0]=pv;
    }
    
    void addpval(double pv){
        int length=java.lang.reflect.Array.getLength(val);
        double[] tmp=new double[length+1];
        for(int i=0;i<length;i++){
            tmp[i]=val[i];
        }//end for i
        tmp[length]=pv;
        val=tmp;
    }
    
    public int query=-1;
    public int hit=-1;
    public double[] val=new double[0];
    
}
