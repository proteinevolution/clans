/*
 * hashkey.java
 *
 * Created on September 9, 2005, 7:15 PM
 */
package clans;
/**
 *
 * @author  tancred
 */
public class hashkey_bak {
    
    /** Creates a new instance of hashkey */
    public hashkey_bak(int i) {
        this.i=i;
    }
    
    int i;
    
    public boolean equals(Object other){
        if(other==null){
            return false;
        }else if(other instanceof hashkey_bak){
            if(((hashkey_bak)other).i==this.i){
                return true;
            }else{
                return false;
            }
        }else{
            return false;
        }
    }
    
    public int hashCode(){
        return i;
    }
    
}
