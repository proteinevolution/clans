/*
 * mapnode.java
 *
 * Created on June 20, 2006, 3:15 PM
 */
package clans;
/**
 *
 * @author  tancred
 */
public class mapnode {
    
    /** Creates a new instance of mapnode */
    public mapnode() {
    }
    
    String level="undefined";
    String info="undefined";
    mapnode[] child=null;
    String leaf="undefined";
    
    //--------------------------------------------------------------------------
    
    public int haschild(mapnode check){
        //see if this node contains a specific child
        if(child==null){
            return -1;
        }
        for(int i=java.lang.reflect.Array.getLength(child);--i>=0;){
            if(check==child[i]){
                return i;
            }
        }//end for i
        return -1;
    }
    
    //--------------------------------------------------------------------------
    
    public void add(mapnode addchild){
        if(child==null){
            child=new mapnode[1];
            child[0]=addchild;
        }else{
            int size=java.lang.reflect.Array.getLength(child);
            mapnode[] tmp=child;
            child=new mapnode[size+1];
            System.arraycopy(tmp,0,child,0,size);
            child[size]=addchild;
        }
    }
    
    //--------------------------------------------------------------------------
    
}
